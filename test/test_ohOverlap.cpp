#include <numeric> // std::exclusive_scan
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_scalar.hpp> // divide_no_remainder
#include <redev_comm.h>
#include "wdmcpl.h"
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
  std::string_view cpnFileName(argv[3]);
  ts::ClassificationPartition classPartition;
  if(isRdv) {
    const auto facePartition = !rank ? ts::readClassPartitionFile(cpnFileName) : ts::ClassificationPartition();
    ts::migrateMeshElms(mesh, facePartition);
    classPartition = ts::CreateClassificationPartition(mesh);
    ts::writeVtk(mesh,"rdvClassPtn",0);
  } else {
    ts::writeVtk(mesh,"appPartition",0);
  }
  auto partition = redev::ClassPtn(classPartition.ranks,classPartition.modelEnts);
  partition.Gather(MPI_COMM_WORLD); //FIXME - move to redev::ClassPtn ctor
  redev::Redev rdv(MPI_COMM_WORLD,partition,isRdv);
  rdv.Setup(); //FIXME - move to redev ctor

  const std::string name = "meshVtxIds";
  const int rdvRanks = 4; //TODO - add the exchange of rank count to the redev::Setup call
  const int appRanks = 16;
  //TODO - name the endpoints in the rdv.get*Engine() APIs
  redev::AdiosComm<redev::GO> commA2R(MPI_COMM_WORLD, rdvRanks, rdv.getToEngine(), rdv.getToIO(), name+"_A2R");
  redev::AdiosComm<redev::GO> commR2A(MPI_COMM_WORLD, appRanks, rdv.getFromEngine(), rdv.getFromIO(), name+"_R2A");

  auto isOverlap = markOverlapMeshEntities(mesh);
  auto isOverlap_h = Omega_h::HostRead(isOverlap);
  if(isRdv) {
    ts::writeVtk(mesh,"rdvOverlap",0);
  } else {
    ts::writeVtk(mesh,"appOverlap",0);
  }

  //Build the dest, offsets, and permutation arrays for the forward
  //send from non-rendezvous to rendezvous.
  ts::OutMsg appOut = !isRdv ? ts::prepareAppOutMessage(mesh, partition) : ts::OutMsg();
  if(!isRdv) {
    commA2R.SetOutMessageLayout(appOut.dest, appOut.offset); //TODO - can this be moved to the AdiosComm ctor 
  }

  //TODO - Document why rendezvous needs two permutations but the app does not
  redev::GOs rdvInPermute;
  ts::CSR rdvOutPermute;
  ts::OutMsg rdvOut;

  for(int iter=0; iter<3; iter++) {
    if(!rank) fprintf(stderr, "isRdv %d iter %d\n", isRdv, iter);
    //////////////////////////////////////////////////////
    //the non-rendezvous app sends global vtx ids to rendezvous
    //////////////////////////////////////////////////////
    if(!isRdv) {
      //fill message array
      auto gids = mesh.globals(0);
      auto gids_h = Omega_h::HostRead(gids);
      redev::GOs msgs(isOverlap.size(),0);
      int j=0;
      for(size_t i=0; i<gids_h.size(); i++) {
        if( isOverlap_h[i] ) {
          msgs[appOut.permute[j++]] = gids_h[i];
        }
      }
      auto start = std::chrono::steady_clock::now();
      commA2R.Send(msgs.data());
      ts::getAndPrintTime(start,name + " appWrite",rank);
    } else {
      auto start = std::chrono::steady_clock::now();
      const auto msgs = commA2R.Recv();
      ts::getAndPrintTime(start,name + " rdvRead",rank);
      //attach the ids to the mesh
      if(iter==0) {
        //We have received the first input message in the rendezvous
        //processes.  Using the meta data of the incoming message we will:
        //- compute the permutation from the incoming vertex global ids to the
        //  on-process global ids
        //- set the message layout for the reverse (rendezvous->non-rendezvous) send by
        //  building the dest and offsets array.
        //- compute the reverse send's permutation array using the layout of
        //  global vertex ids in 'msgs'.
        //These operations only need to be done once per coupling as long as
        //the topology and partition of the rendezvous and non-rendezvous meshes
        //remains the same.
        auto rdvIn = commA2R.GetInMessageLayout();
        rdvInPermute = ts::getRdvPermutation(mesh, msgs);
        rdvOut = ts::prepareRdvOutMessage(mesh,rdvIn);
        commR2A.SetOutMessageLayout(rdvOut.dest,rdvOut.offset);
        rdvOutPermute = ts::getRdvOutPermutation(mesh, msgs);
      }
      ts::checkAndAttachIds(mesh, "inVtxGids", msgs, rdvInPermute);
      ts::writeVtk(mesh,"rdvInGids",iter);
    } //end non-rdv -> rdv
    //////////////////////////////////////////////////////
    //the rendezvous app sends global vtx ids to non-rendezvous
    //////////////////////////////////////////////////////
    if(isRdv) {
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
      auto start = std::chrono::steady_clock::now();
      commR2A.Send(msgs.data());
      ts::getAndPrintTime(start,name + " rdvWrite",rank);
    } else {
      auto start = std::chrono::steady_clock::now();
      const auto msgs = commR2A.Recv();
      ts::getAndPrintTime(start,name + " appRead",rank);
      { //check incoming messages are in the correct order
        auto gids = mesh.globals(0);
        auto gids_h = Omega_h::HostRead(gids);
        int j=0;
        for(size_t i=0; i<gids_h.size(); i++) {
          if( isOverlap_h[i] ) {
            REDEV_ALWAYS_ASSERT(msgs[appOut.permute[j++]] == gids_h[i]);
          }
        }
      }
    } //end rdv -> non-rdv
  } //end iter loop
  return 0;
}
