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
 * On the server we mark the vertices on each process that are in the overlap
 * region and are owned by the process as defined by the Classification
 * Partition GetRank(modelEntity) function.  The ownership of mesh entities on
 * the mesh partition boundary is not guaranteed to match the ownership of the
 * geometric model entities defined by the Classification Partition so we must
 * use GetRank(modelEntity) to ensure consistency with data incoming from
 * the client processes.
 *
 * On the client side we only care that each mesh vertex is sent by exactly one
 * client process (to the server process returned by GetRank(modelEntity)) so
 * using the ownership of mesh entities following the mesh partition ownership
 * is OK. The function markMeshOverlapRegion(...) supports this.
 */
auto markServerOverlapRegion(Omega_h::Mesh& mesh,
                             const redev::ClassPtn& classPtn)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // transfer vtx classification to host
  auto classIds = mesh.get_array<Omega_h::ClassId>(0, "class_id");
  auto classIds_h = Omega_h::HostRead(classIds);
  auto classDims = mesh.get_array<Omega_h::I8>(0, "class_dim");
  auto classDims_h = Omega_h::HostRead(classDims);
  auto isOverlap = Omega_h::Write<Omega_h::I8>(classIds.size(), "isOverlap");
  auto markOverlap = OMEGA_H_LAMBDA(int i)
  {
    isOverlap[i] = isModelEntInOverlap(classDims[i], classIds[i]);
  };
  Omega_h::parallel_for(classIds.size(), markOverlap);
  auto owned_h = Omega_h::HostRead(mesh.owned(0));
  auto isOverlap_h = Omega_h::HostRead<Omega_h::I8>(isOverlap);
  // mask to only class partition owned entities
  auto isOverlapOwned = Omega_h::HostWrite<Omega_h::I8>(
    classIds.size(), "isOverlapAndOwnsModelEntInClassPartition");
  for (int i = 0; i < mesh.nverts(); i++) {
    redev::ClassPtn::ModelEnt ent(classDims_h[i], classIds_h[i]);
    auto destRank = classPtn.GetRank(ent);
    auto isModelEntOwned = (destRank == rank);
    isOverlapOwned[i] = isModelEntOwned && isOverlap_h[i];
    if (owned_h[i] && !isModelEntOwned) {
      fprintf(stderr, "%d owner conflict %d ent (%d,%d) owner %d owned %d\n",
              rank, i, classDims_h[i], classIds_h[i], destRank, owned_h[i]);
    }
  }
  // this is a crime: host -> device -> host
  auto isOverlapOwned_dr = Omega_h::read(Omega_h::Write(isOverlapOwned));
  // auto isOverlapOwned_hr = Omega_h::HostRead(isOverlapOwned_dr);
  mesh.add_tag(0, "isOverlap", 1, isOverlapOwned_dr);
  return isOverlapOwned_dr;
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
  REDEV_ALWAYS_ASSERT(mesh.nelems()); // all ranks should have elements
  auto ptn = ts::CreateClassificationPartition(mesh);
  return redev::ClassPtn(MPI_COMM_WORLD, ptn.ranks, ptn.modelEnts);
}

struct MeanCombiner
{
  // void operator()(const std::vector<std::reference_wrapper<OHField>>& fields,
  //                 OHField& combined) const
  void operator()(
    const nonstd::span<
      const std::reference_wrapper<wdmcpl::detail::InternalField>>& fields,
    wdmcpl::detail::InternalField& combined_variant) const
  {
    WDMCPL_ALWAYS_ASSERT(!combined_variant.valueless_by_exception());
    WDMCPL_ALWAYS_ASSERT(!fields[0].get().valueless_by_exception());
    /*
    std::visit(
      //[&fields](auto&& combined_field) {
        [](auto&& combined_field) {
        //using T = typename std::remove_reference_t<
        //  decltype(combined_field)>::value_type;
        //Omega_h::Write<T> combined_array(combined_field.Size());
        //for (auto& field_variant : fields) {
          //std::visit(
          //  [&combined_array](auto&& field) {
          //    WDMCPL_ALWAYS_ASSERT(field.Size() == combined_array.size());
          //    auto field_array = get_nodal_data(field);
          //    Omega_h::parallel_for(
          //      field_array.size(),
          //      OMEGA_H_LAMBDA(int i) { combined_array[i] += field_array[i]; });
          //  },
          //  field_variant.get());
        //}
        //auto num_fields = fields.size();
        //Omega_h::parallel_for(
        //  combined_array.size(),
        //  OMEGA_H_LAMBDA(int i) { combined_array[i] /= num_fields; });
        //set(combined_field, make_array_view(Omega_h::Read(combined_array)));
      },
      combined_variant);
      */
  }
};

using wdmcpl::Copy;
using wdmcpl::CouplerClient;
using wdmcpl::CouplerServer;
using wdmcpl::FieldEvaluationMethod;
using wdmcpl::FieldTransferMethod;
using wdmcpl::GO;
using wdmcpl::Lagrange;
using wdmcpl::make_array_view;
using wdmcpl::OmegaHField;
using wdmcpl::OmegaHFieldShim;

static constexpr bool done = true;

void xgc_delta_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  CouplerClient cpl("proxy_couple", comm, redev::ClassPtn{});
  auto is_overlap = markOverlapMeshEntities(mesh);
  cpl.AddField("delta_f_gids", OmegaHFieldShim<GO>("global", mesh, is_overlap));
  do {
    cpl.SendField("delta_f_gids");    //(Alt) df_gid_field->Send();
    cpl.ReceiveField("delta_f_gids"); //(Alt) df_gid_field->Receive();
  } while (!done);
}
void xgc_total_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  wdmcpl::CouplerClient cpl("proxy_couple", comm, redev::ClassPtn{});
  auto is_overlap = markOverlapMeshEntities(mesh);
  cpl.AddField("total_f_gids", OmegaHFieldShim<GO>("global", mesh, is_overlap));
  do {
    cpl.SendField("total_f_gids");    //(Alt) tf_gid_field->Send();
    cpl.ReceiveField("total_f_gids"); //(Alt) tf_gid_field->Receive();
  } while (!done);
}
void coupler(MPI_Comm comm, Omega_h::Mesh& mesh, std::string_view cpn_file)
{
  auto p = setupServerPartition(mesh, cpn_file);
  // coupling server using same mesh as application
  wdmcpl::CouplerServer cpl("proxy_couple", comm, std::move(p), mesh);
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  auto is_overlap = markServerOverlapRegion(mesh, partition);
  cpl.AddField("total_f_gids",
               OmegaHFieldShim<GO>("total_f_gids", mesh, is_overlap),
               FieldTransferMethod::None, // to Omega_h
               FieldEvaluationMethod::None,
               FieldTransferMethod::None, // from Omega_h
               FieldEvaluationMethod::None);
  cpl.AddField("delta_f_gids",
               OmegaHFieldShim<GO>("delta_f_gids", mesh, is_overlap),
               FieldTransferMethod::None, FieldEvaluationMethod::None,
               FieldTransferMethod::None, FieldEvaluationMethod::None);
  // CombinerFunction is a functor that takes a vector of omega_h fields
  // combines their values and sets the combined values into the resultant field
  auto* gather = cpl.AddGatherFieldsOp("cpl1", {"total_f_gids", "delta_f_gids"},
                                       "combined_gids", MeanCombiner{});
  auto* scatter = cpl.AddScatterFieldsOp("cpl1", "combined_gids",
                                         {"total_f_gids", "delta_f_gids"});
  // for case with symmetric Gather/Scatter we have
  // auto [gather, scatter] = cpl.AddSymmetricGatherScatterOp("cpl1",
  // {"total_f_gids", "delta_f_gids"},
  //                      "combined_gids", MeanCombiner{});
  //WDMCPL_ALWAYS_ASSERT(gather != nullptr);
  //WDMCPL_ALWAYS_ASSERT(scatter != nullptr);
  do {
    //  Gather OHField
    // 1. receives any member fields .Receive()
    // 2. field_transfer native to internal
    // 3. combine internal fields into combined internal field
     //cpl.GatherFields("cpl1"); // (Alt) gather->Run();
    gather->Run();
    //cpl.ReceiveField("total_f_gids");
    //cpl.ReceiveField("delta_f_gids");

    // Scatter OHField
    // 1. OHField transfer internal to native
    // 2. Send data to members
    //cpl.ScatterFields("cpl1"); // (Alt) scatter->Run();
    //cpl.SendField("total_f_gids");
    //cpl.SendField("delta_f_gids");
    scatter->Run();
  } while (!done);
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
