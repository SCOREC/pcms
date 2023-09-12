#include <numeric> // std::exclusive_scan
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_scalar.hpp> // divide_no_remainder
#include "pcms.h"
#include "test_support.h"

namespace ts = test_support;


//TODO - use attributes on the geometric model to
//       define which model entities are in the
//       buffer/blended/overlap regions.
//       This is currently hardcoded for the D3D
//       case in the coupling data repo.
/**
 * return 1 if the specificed model entity is part of the overlap region, 0
 * otherwise
 */
OMEGA_H_DEVICE Omega_h::I8 isModelEntInOverlap(const int dim, const int id) {
  //the TOMMS generated geometric model has
  //entity IDs that increase with the distance
  //from the magnetic axis
  if (dim == 2 && (id >= 22 && id <= 34) ) {
      return 1;
  } else if (dim == 1 && (id >= 21 && id <= 34) ) {
      return 1;
  } else if (dim == 0 && (id >= 21 && id <= 34) ) {
      return 1;
  }
  return 0;
}

/**
 * Create the tag 'isOverlap' for each mesh vertex whose value is 1 if the
 * vertex is classified on a model entity in the closure of the geometric model
 * faces forming the overlap region; the value is 0 otherwise.
 */
Omega_h::Read<Omega_h::I8> markOverlapMeshEntities(Omega_h::Mesh& mesh) {
  //transfer vtx classification to host
  auto classIds = mesh.get_array<Omega_h::ClassId>(0, "class_id");
  auto classDims = mesh.get_array<Omega_h::I8>(0, "class_dim");
  auto isOverlap = Omega_h::Write<Omega_h::I8>(classIds.size(), "isOverlap");
  auto markOverlap = OMEGA_H_LAMBDA(int i) {
    isOverlap[i] = isModelEntInOverlap(classDims[i], classIds[i]);
  };
  Omega_h::parallel_for(classIds.size(), markOverlap);
  auto isOverlap_r = Omega_h::read(isOverlap);
  mesh.add_tag(0, "isOverlap", 1, isOverlap_r);
  return isOverlap_r;
}

redev::ClassPtn setupClientPartition(Omega_h::Mesh& mesh) {
  ts::writeVtk(mesh,"appPartition",0);
  auto ptn = ts::ClassificationPartition();
  return redev::ClassPtn(MPI_COMM_WORLD, ptn.ranks, ptn.modelEnts);
}

redev::ClassPtn setupServerPartition(Omega_h::Mesh& mesh, std::string_view cpnFileName) {
  auto ohComm = mesh.comm();
  const auto facePartition = !ohComm->rank() ? ts::readClassPartitionFile(cpnFileName) :
                                               ts::ClassificationPartition();
  ts::migrateMeshElms(mesh, facePartition);
  auto ptn = ts::CreateClassificationPartition(mesh);
  return redev::ClassPtn(MPI_COMM_WORLD, ptn.ranks, ptn.modelEnts);
}

auto setupComms(redev::Redev& rdv, std::string_view name, const int clientId) {
  const bool isSST = false;
  adios2::Params params{ {"Streaming", "On"}, {"OpenTimeoutSecs", "12"}};
  REDEV_ALWAYS_ASSERT(clientId == 0 || clientId ==1);
  std::stringstream clientName;
  clientName << name << "Client" << clientId;
  auto channel = rdv.CreateAdiosChannel(clientName.str(), params,static_cast<redev::TransportType>(isSST) );
  return channel.CreateComm<redev::GO>(clientName.str(), rdv.GetMPIComm());
}

Omega_h::HostRead<Omega_h::I8> markMeshOverlapRegion(Omega_h::Mesh& mesh) {
  auto isOverlap = markOverlapMeshEntities(mesh);
  return Omega_h::HostRead(isOverlap);
}

void clientCheckIncomingMessages(Omega_h::Mesh& mesh,
    Omega_h::HostRead<Omega_h::I8> isOverlap_h,
    const std::vector<redev::GO>& msgsIn,
    const ts::OutMsg& appOut) {
  //check incoming messages are in the correct order
  auto gids = mesh.globals(0);
  auto gids_h = Omega_h::HostRead(gids);
  int j=0;
  for(size_t i=0; i<gids_h.size(); i++) {
    if( isOverlap_h[i] ) {
      REDEV_ALWAYS_ASSERT(msgsIn[appOut.permute[j++]] == gids_h[i]);
    }
  }
}

void client(Omega_h::Mesh& mesh, std::string fieldName, const int clientId) {
  const bool isRdv = false;
  auto ohComm = mesh.comm();
  const auto rank = ohComm->rank();
  if(!rank) fprintf(stderr, "clientId %d\n", clientId);

  auto partition = setupClientPartition(mesh);
  auto rdv = redev::Redev(MPI_COMM_WORLD,redev::Partition{std::move(partition)},static_cast<redev::ProcessType>(isRdv));
  auto comm = setupComms(rdv,fieldName,clientId);

  auto isOverlap_h = markMeshOverlapRegion(mesh);
  //////////////////////////////////////////////////////
  //the non-rendezvous app sends global vtx ids to rendezvous
  //////////////////////////////////////////////////////
  ts::OutMsg appOut = ts::prepareAppOutMessage(mesh, std::get<decltype(partition)>(rdv.GetPartition()));
  //Build the dest, offsets, and permutation arrays for the forward
  //send from client to rendezvous/server.
  comm.SetOutMessageLayout(appOut.dest, appOut.offset); //TODO - can this be moved to the AdiosComm ctor
  //fill message array
  auto gids = mesh.globals(0);
  auto gids_h = Omega_h::HostRead(gids);
  redev::GOs msgs(isOverlap_h.size(),0);
  int j=0;
  for(size_t i=0; i<gids_h.size(); i++) {
    if( isOverlap_h[i] ) {
      msgs[appOut.permute[j++]] = gids_h[i];
    }
  }
  auto start = std::chrono::steady_clock::now();
  comm.Send(msgs.data(),redev::Mode::Synchronous);
  ts::getAndPrintTime(start,fieldName + " appWrite",rank);
  //////////////////////////////////////////////////////
  //the rendezvous app sends global vtx ids to non-rendezvous
  //////////////////////////////////////////////////////
  start = std::chrono::steady_clock::now();
  const auto msgsIn = comm.Recv(redev::Mode::Synchronous);
  ts::getAndPrintTime(start,fieldName + " appRead",rank);
  clientCheckIncomingMessages(mesh,isOverlap_h,msgsIn,appOut);
  //////////////////////////////////////////////////////
  //communication loop
  //////////////////////////////////////////////////////
  for(int iter=0; iter<3; iter++) {
    if(!rank) fprintf(stderr, "isRdv %d iter %d\n", isRdv, iter);
    //send to server
    int j=0;
    for(size_t i=0; i<gids_h.size(); i++) {
      if( isOverlap_h[i] ) {
        msgs[appOut.permute[j++]] = gids_h[i];
      }
    }
    auto start = std::chrono::steady_clock::now();
    comm.Send(msgs.data(),redev::Mode::Synchronous);
    ts::getAndPrintTime(start,fieldName + " appWrite",rank);
    //receive from server
    start = std::chrono::steady_clock::now();
    const auto msgsIn = comm.Recv(redev::Mode::Synchronous);
    ts::getAndPrintTime(start,fieldName + " appRead",rank);
    clientCheckIncomingMessages(mesh,isOverlap_h,msgsIn,appOut);
  } //end iter loop
}

struct ClientMetaData {
  redev::GOs inPermute;
  ts::CSR outPermute;
  ts::OutMsg outMsg;
};

void serverReceiveFromClient(ClientMetaData& clientMeta,
    redev::BidirectionalComm<redev::GO>& comm, Omega_h::Mesh& mesh,
    std::string_view fieldName, const int rank, const int clientId) {
  std::stringstream ss;
  ss << fieldName << " rdvRead clientId " << clientId;
  auto start = std::chrono::steady_clock::now();
  const auto msgsIn = comm.Recv(redev::Mode::Synchronous);
  ts::getAndPrintTime(start,ss.str(),rank);
  auto rdvIn = comm.GetInMessageLayout();
  //setup outbound meta data
  // FIXME: this assumes we know the correct order of the data already!
  // If we transfer data/msg that's not GID we can't do this comparison on the
  // message directly, so the global ids may need to be sent during an initialization
  // phase
  clientMeta.inPermute = ts::getRdvPermutation(mesh, msgsIn);
  clientMeta.outMsg = ts::prepareRdvOutMessage(mesh,rdvIn);
  // FIXME is this necessary to set the outgoing message layout here? The
  // reciever shouldn't need to know the layout of the sender
  comm.SetOutMessageLayout(clientMeta.outMsg.dest,clientMeta.outMsg.offset);
  clientMeta.outPermute = ts::getRdvOutPermutation(mesh, msgsIn);
  //attach ids to the mesh
  ss.str("");
  ss << "inVtxGidsClientId" << clientId;
  ts::checkAndAttachIds(mesh, ss.str(), msgsIn, clientMeta.inPermute);
}

void serverSendToClient(ClientMetaData& clientMeta,
    redev::BidirectionalComm<redev::GO>& comm, Omega_h::Mesh& mesh,
    Omega_h::HostRead<Omega_h::I8>& isOverlap_h,
    std::string_view fieldName, const int rank, const int clientId) {
  //fill message array
  auto gids = mesh.globals(0);
  auto gids_h = Omega_h::HostRead(gids);
  redev::GOs msgs(clientMeta.outPermute.off.back());
  for(int i=0; i<gids_h.size(); i++) {
    for(int j=clientMeta.outPermute.off[i]; j<clientMeta.outPermute.off[i+1]; j++) {
      REDEV_ALWAYS_ASSERT(isOverlap_h[i]);
      msgs[clientMeta.outPermute.val[j]] = gids_h[i];
    }
  }
  std::stringstream ss;
  ss << fieldName << " rdvWrite clientId " << clientId;
  auto start = std::chrono::steady_clock::now();
  comm.Send(msgs.data());
  ts::getAndPrintTime(start,ss.str(),rank);
}

void server(Omega_h::Mesh& mesh, std::string fieldName, std::string_view cpnFileName) {
  const bool isRdv = true;
  auto ohComm = mesh.comm();
  const auto rank = ohComm->rank();
  auto partition = setupServerPartition(mesh,cpnFileName);
  auto rdv = redev::Redev(MPI_COMM_WORLD,std::move(partition),static_cast<redev::ProcessType>(isRdv));
  auto commClient0 = setupComms(rdv,fieldName,0);
  auto commClient1 = setupComms(rdv,fieldName,1);

  auto isOverlap_h = markMeshOverlapRegion(mesh);
  //TODO - Document why rendezvous needs two permutations but the app does not
  ClientMetaData client0;
  ClientMetaData client1;

  serverReceiveFromClient(client0,commClient0,mesh,fieldName,rank,0);
  serverReceiveFromClient(client1,commClient1,mesh,fieldName,rank,1);
  ts::writeVtk(mesh,"rdvInGids",0);

  serverSendToClient(client0,commClient0,mesh,isOverlap_h,fieldName,rank,0);
  serverSendToClient(client1,commClient1,mesh,isOverlap_h,fieldName,rank,1);
  //////////////////////////////////////////////////////
  //communication loop
  //////////////////////////////////////////////////////
  for(int iter=0; iter<3; iter++) {
    if(!rank) fprintf(stderr, "isRdv %d iter %d\n", isRdv, iter);
    //receive from clients
    auto start = std::chrono::steady_clock::now();
    const auto msgsIn0 = commClient0.Recv(redev::Mode::Synchronous);
    ts::getAndPrintTime(start,fieldName + " rdvRead clientId 0",rank);
    ts::checkAndAttachIds(mesh, "inVtxGidsClient0", msgsIn0, client0.inPermute);

    start = std::chrono::steady_clock::now();
    const auto msgsIn1 = commClient1.Recv(redev::Mode::Synchronous);
    ts::getAndPrintTime(start,fieldName + " rdvRead clientId 1",rank);
    ts::checkAndAttachIds(mesh, "inVtxGidsClient1", msgsIn1, client1.inPermute);
    //send to clients
    serverSendToClient(client0,commClient0,mesh,isOverlap_h,fieldName,rank,0);
    serverSendToClient(client1,commClient1,mesh,isOverlap_h,fieldName,rank,1);
  } //end iter loop
}


int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  if(argc != 4) {
    if(!rank) {
      std::cerr << "Usage: " << argv[0] << " <clientId=-1|0|1> /path/to/omega_h/mesh /path/to/partitionFile.cpn\n";
    }
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 4);
  const auto clientId = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <=1);
  const auto meshFile = argv[2];
  const auto classPartitionFile = argv[3];
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(meshFile, lib.world(), &mesh);
  const std::string name = "meshVtxIds";
  if(clientId == -1) { //rendezvous
    server(mesh,name,classPartitionFile);
  } else {
    client(mesh,name,clientId);
  }
  return 0;
}
