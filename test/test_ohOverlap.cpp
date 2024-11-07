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


int main(int argc, char** argv) {
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  if(argc != 4) {
    if(!rank) {
      std::cerr << "Usage: " << argv[0] << " <1=isRendezvousApp,0=isParticipant> /path/to/omega_h/mesh /path/to/partitionFile.cpn\n";
    }
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 4);
  const auto isRdv = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(isRdv==1 || isRdv==0);
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[2], lib.world(), &mesh);
  const std::string name = "meshVtxIds";
  if(isRdv) {
    ///////////////////  SERVER /////////////////////////
    std::string_view cpnFileName(argv[3]);
    auto partition = setupServerPartition(mesh,cpnFileName);
    auto rdv = redev::Redev(MPI_COMM_WORLD,redev::Partition{partition},static_cast<redev::ProcessType>(isRdv));
    adios2::Params params{ {"Streaming", "On"}, {"OpenTimeoutSecs", "60"}};
    auto channel = rdv.CreateAdiosChannel(std::string(name),params,redev::TransportType::BP4);
    auto comm = channel.CreateComm<redev::GO>(std::string(name),rdv.GetMPIComm());
    auto isOverlap_h = markMeshOverlapRegion(mesh);
    ts::OutMsg appOut = ts::OutMsg(); //FIXME is this needed?
    //TODO - Document why rendezvous needs two permutations but the app does not
    redev::GOs rdvInPermute;
    ts::CSR rdvOutPermute;
    ts::OutMsg rdvOut;
    //////////////////////////////////////////////////////
    //the non-rendezvous app sends global vtx ids to rendezvous
    //////////////////////////////////////////////////////
    auto start = std::chrono::steady_clock::now();
    channel.BeginReceiveCommunicationPhase();
    const auto msgsIn = comm.Recv(redev::Mode::Synchronous);
    channel.EndReceiveCommunicationPhase();
    ts::getAndPrintTime(start,name + " rdvRead",rank);
    //We have received the first input message in the rendezvous
    //processes.  Using the meta data of the incoming message we will:
    //- compute the permutation from the incoming vertex global ids to the
    //  on-process global ids
    //- set the message layout for the reverse (rendezvous->non-rendezvous) send by
    //  building the dest and offsets array.
    //- compute the reverse send's permutation array using the layout of
    //  global vertex ids in 'msgsIn'.
    //These operations only need to be done once per coupling as long as
    //the topology and partition of the rendezvous and non-rendezvous meshes
    //remains the same.
    auto rdvIn = comm.GetInMessageLayout();
    rdvInPermute = ts::getRdvPermutation(mesh, msgsIn);
    rdvOut = ts::prepareRdvOutMessage(mesh,rdvIn);
    comm.SetOutMessageLayout(rdvOut.dest,rdvOut.offset);
    rdvOutPermute = ts::getRdvOutPermutation(mesh, msgsIn);
    //attach ids to the mesh
    ts::checkAndAttachIds(mesh, "inVtxGids", msgsIn, rdvInPermute);
    ts::writeVtk(mesh,"rdvInGids",0);
    //////////////////////////////////////////////////////
    //the rendezvous app sends global vtx ids to the client
    //////////////////////////////////////////////////////
    //fill message array
    auto gids = mesh.globals(0);
    auto gids_h = Omega_h::HostRead(gids);
    redev::GOs msgs(rdvOutPermute.off.back());
    for(int i=0; i<gids_h.size(); i++) {
      for(int j=rdvOutPermute.off[i]; j<rdvOutPermute.off[i+1]; j++) {
        REDEV_ALWAYS_ASSERT(isOverlap_h[i]);
        msgs[rdvOutPermute.val[j]] = gids_h[i];
      }
    }
    start = std::chrono::steady_clock::now();
    channel.BeginSendCommunicationPhase();
    comm.Send(msgs.data(),redev::Mode::Synchronous);
    channel.EndSendCommunicationPhase();
    ts::getAndPrintTime(start,name + " rdvWrite",rank);
    //////////////////////////////////////////////////////
    //communication loop
    //////////////////////////////////////////////////////
    for(int iter=0; iter<3; iter++) {
      if(!rank) fprintf(stderr, "isRdv %d iter %d\n", isRdv, iter);
      //receive from client
      auto start = std::chrono::steady_clock::now();
      channel.BeginReceiveCommunicationPhase();
      const auto msgsIn = comm.Recv(redev::Mode::Synchronous);
      channel.EndReceiveCommunicationPhase();
      ts::getAndPrintTime(start,name + " rdvRead",rank);
      ts::checkAndAttachIds(mesh, "inVtxGids", msgsIn, rdvInPermute);
      //send to client
      for(int i=0; i<gids_h.size(); i++) {
        for(int j=rdvOutPermute.off[i]; j<rdvOutPermute.off[i+1]; j++) {
          REDEV_ALWAYS_ASSERT(isOverlap_h[i]);
          msgs[rdvOutPermute.val[j]] = gids_h[i];
        }
      }
      start = std::chrono::steady_clock::now();
      channel.BeginSendCommunicationPhase();
      comm.Send(msgs.data(),redev::Mode::Synchronous);
      channel.EndSendCommunicationPhase();
      ts::getAndPrintTime(start,name + " rdvWrite",rank);
    } //end iter loop
  } else {
    ///////////////////  CLIENT /////////////////////////
    auto partition = setupClientPartition(mesh);
    auto rdv = redev::Redev(MPI_COMM_WORLD,redev::Partition{std::move(partition)},static_cast<redev::ProcessType>(isRdv));
    adios2::Params params{ {"Streaming", "On"}, {"OpenTimeoutSecs", "60"}};
    auto channel = rdv.CreateAdiosChannel(std::string(name),params,redev::TransportType::BP4);
    auto comm = channel.CreateComm<redev::GO>(std::string(name),rdv.GetMPIComm());
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
    for(pcms::LO i=0; i<gids_h.size(); i++) {
      if( isOverlap_h[i] ) {
        msgs[appOut.permute[j++]] = gids_h[i];
      }
    }
    auto start = std::chrono::steady_clock::now();
    channel.BeginSendCommunicationPhase();
    comm.Send(msgs.data(),redev::Mode::Synchronous);
    channel.EndSendCommunicationPhase();
    ts::getAndPrintTime(start,name + " appWrite",rank);
    //////////////////////////////////////////////////////
    //the rendezvous app sends global vtx ids to non-rendezvous
    //////////////////////////////////////////////////////
    start = std::chrono::steady_clock::now();
    channel.BeginReceiveCommunicationPhase();
    const auto msgsIn = comm.Recv(redev::Mode::Synchronous);
    channel.EndReceiveCommunicationPhase();
    ts::getAndPrintTime(start,name + " appRead",rank);
    clientCheckIncomingMessages(mesh,isOverlap_h,msgsIn,appOut);
    //////////////////////////////////////////////////////
    //communication loop
    //////////////////////////////////////////////////////
    for(int iter=0; iter<3; iter++) {
      if(!rank) fprintf(stderr, "isRdv %d iter %d\n", isRdv, iter);
      //send to server
      int j=0;
      for(pcms::LO i=0; i<gids_h.size(); i++) {
        if( isOverlap_h[i] ) {
          msgs[appOut.permute[j++]] = gids_h[i];
        }
      }
      auto start = std::chrono::steady_clock::now();
      channel.BeginSendCommunicationPhase();
      comm.Send(msgs.data(),redev::Mode::Synchronous);
      channel.EndSendCommunicationPhase();
      ts::getAndPrintTime(start,name + " appWrite",rank);
      //receive from server
      start = std::chrono::steady_clock::now();
      channel.BeginReceiveCommunicationPhase();
      const auto msgsIn = comm.Recv(redev::Mode::Synchronous);
      channel.EndReceiveCommunicationPhase();
      ts::getAndPrintTime(start,name + " appRead",rank);
      clientCheckIncomingMessages(mesh,isOverlap_h,msgsIn,appOut);
    } //end iter loop
  }
  return 0;
}
