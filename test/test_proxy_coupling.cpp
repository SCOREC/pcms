#include <Omega_h_mesh.hpp>
#include <iostream>
#include <wdmcpl.h>
#include <wdmcpl/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <redev_variant_tools.h>
#include "test_support.h"
#include <wdmcpl/omega_h_field.h>

/**
 * return 1 if the specificed model entity is part of the overlap region, 0
 * otherwise
 */
OMEGA_H_DEVICE Omega_h::I8 isModelEntInOverlap(const int dim, const int id)
{
  // the TOMMS generated geometric model has
  // entity IDs that increase with the distance
  // from the magnetic axis
  if (dim == 2 && (id >= 22 && id <= 34)) {
    return 1;
  } else if (dim == 1 && (id >= 21 && id <= 34)) {
    return 1;
  } else if (dim == 0 && (id >= 21 && id <= 34)) {
    return 1;
  }
  return 0;
}

/**
 * Create the tag 'isOverlap' for each mesh vertex whose value is 1 if the
 * vertex is classified on a model entity in the closure of the geometric model
 * faces forming the overlap region; the value is 0 otherwise.
 * OnlyIncludesOverlapping and owned verts
 */
Omega_h::Read<Omega_h::I8> markOverlapMeshEntities(Omega_h::Mesh& mesh)
{
  // transfer vtx classification to host
  auto classIds = mesh.get_array<Omega_h::ClassId>(0, "class_id");
  auto classDims = mesh.get_array<Omega_h::I8>(0, "class_dim");
  auto isOverlap = Omega_h::Write<Omega_h::I8>(classIds.size(), "isOverlap");
  auto markOverlap = OMEGA_H_LAMBDA(int i)
  {
    isOverlap[i] = isModelEntInOverlap(classDims[i], classIds[i]);
  };
  Omega_h::parallel_for(classIds.size(), markOverlap);
  auto isOwned = mesh.owned(0);
  // try masking out to only owned entities
  Omega_h::parallel_for(
    isOverlap.size(),
    OMEGA_H_LAMBDA(int i) { isOverlap[i] = (isOwned[i] && isOverlap[i]); });

  auto isOverlap_r = Omega_h::read(isOverlap);
  mesh.add_tag(0, "isOverlap", 1, isOverlap_r);
  return isOverlap_r;
}

redev::ClassPtn setupServerPartition(Omega_h::Mesh& mesh,
                                     std::string_view cpnFileName)
{
  namespace ts = test_support;
  auto ohComm = mesh.comm();
  const auto facePartition = !ohComm->rank()
                               ? ts::readClassPartitionFile(cpnFileName)
                               : ts::ClassificationPartition();
  ts::migrateMeshElms(mesh, facePartition);
  auto ptn = ts::CreateClassificationPartition(mesh);
  return redev::ClassPtn(MPI_COMM_WORLD, ptn.ranks, ptn.modelEnts);
}

using PT = wdmcpl::ProcessType;
using wdmcpl::Copy;
using wdmcpl::GO;
using wdmcpl::Lagrange;
using wdmcpl::make_array_view;
using wdmcpl::OmegaHField;
using wdmcpl::OmegaHFieldShim;
using wdmcpl::Coupler;
using wdmcpl::FieldEvaluationMethod;
using wdmcpl::FieldTransferMethod;

static constexpr bool done = true;

void xgc_delta_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{

  Coupler<PT::Client> cpl("proxy_couple", comm, redev::ClassPtn{});
  auto is_overlap = markOverlapMeshEntities(mesh);
  auto* df_gid_field = cpl.AddField("delta_f_gids", OmegaHFieldShim<GO>("gids", mesh, is_overlap));

  do {
    cpl.SendField("_gids"); //(Alt) df_gid_field->Send();
    cpl.ReceiveField("delta_f_gids"); //(Alt) df_gid_field->Receive();
  } while(!done);
}
void xgc_total_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  wdmcpl::Coupler<PT::Client> cpl("proxy_couple", comm, redev::ClassPtn{});
  auto is_overlap = markOverlapMeshEntities(mesh);
  auto tf_gid_field = cpl.AddField("total_f_gids", OmegaHFieldShim<GO>("gids", mesh, is_overlap));
  do {
    cpl.SendField("total_f_gids"); //(Alt) tf_gid_field->Send();
    cpl.ReceiveField("total_f_gids"); //(Alt) tf_gid_field->Receive();
  } while(!done);
}
void coupler(MPI_Comm comm, Omega_h::Mesh& mesh, std::string_view cpn_file)
{
  wdmcpl::Coupler<PT::Server> cpl("proxy_couple", comm,
                                  setupServerPartition(mesh, cpn_file));
  auto is_overlap = markOverlapMeshEntities(mesh);
  // Note: coupler takes ownership of the field shim as well as
  cpl.AddField("total_f_gids",
               OmegaHFieldShim<GO>("total_f_gids", mesh, is_overlap),
               FieldTransferMethod::None, // to Omega_h
               FieldEvaluationMethod::None,
               FieldTransferMethod::None, // from Omega_h
               FieldEvaluationMethod::None);
  cpl.AddField(
    "delta_f_gids", OmegaHFieldShim<GO>("delta_f_gids", mesh, is_overlap),
    FieldTransferMethod::None, FieldEvaluationMethod::None,
    FieldTransferMethod::None, FieldEvaluationMethod::None);

  // Combiner is a functor that takes a vector of omega_h fields combines their values
  // and sets the combined values into the resultant field
  auto scatter = cpl.AddGatherFieldsOp("cpl1", {"total_f_gids", "delta_f_gids"},
                        "combined_gids", MeanCombiner{});
  auto gather = cpl.AddScatterFieldsOp("cpl1", "combined_gids",
                         {"total_f_gids", "delta_f_gids"});
  // for case with symmetric Gather/Scatter we have
  //auto[scatter, gather] = cpl.AddSymmetricGatherScatterOp("cpl1", {"total_f_gids", "delta_f_gids"},
  //                      "combined_gids", MeanCombiner{});

  do {
    //  Gather Field
    // 1. receives any member fields .Receive()
    // 2. field_transfer native to internal
    // 3. combine internal fields into combined internal field
    cpl.GatherFields("cpl1"); // (Alt) scatter->Run();
    // Scatter Field
    // 1. Field transfer internal to native
    // 2. Send data to members
    cpl.ScatterFields("cpl1"); // (Alt) gather->Run();
  } while(!done);
}

int main(int argc, char** argv)
{
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  if (argc != 4) {
    if (!rank) {
      std::cerr << "Usage: " << argv[0]
                << " <clientId=-1|0|1> /path/to/omega_h/mesh "
                   "/path/to/partitionFile.cpn\n";
    }
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 4);
  const auto clientId = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);
  const auto meshFile = argv[2];
  const auto classPartitionFile = argv[3];
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(meshFile, lib.world(), &mesh);
  MPI_Comm mpi_comm = lib.world()->get_impl();
  const std::string name = "meshVtxIds";
  switch (clientId) {
    case -1: coupler(mpi_comm, mesh, classPartitionFile); break;
    case 0: xgc_delta_f(mpi_comm, mesh); break;
    case 1: xgc_total_f(mpi_comm, mesh); break;
    default:
      std::cerr << "Unhandled client id (should be -1, 0,1)\n";
      exit(EXIT_FAILURE);
  }
  return 0;
}