#include <wdmcpl.h>
#include <wdmcpl/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include "test_support.h"
#include <wdmcpl/omega_h_field.h>
#include <wdmcpl/xgc_field_adapter.h>
#include <chrono>

using wdmcpl::ConstructRCFromOmegaHMesh;
using wdmcpl::Copy;
using wdmcpl::CouplerClient;
using wdmcpl::CouplerServer;
using wdmcpl::FieldEvaluationMethod;
using wdmcpl::FieldTransferMethod;
using wdmcpl::GO;
using wdmcpl::LO;
using wdmcpl::Lagrange;
using wdmcpl::make_array_view;
using wdmcpl::OmegaHField;
using wdmcpl::OmegaHFieldAdapter;
using wdmcpl::ReadReverseClassificationVertex;
using wdmcpl::ReverseClassificationVertex;

namespace ts = test_support;

// TODO: we should communicate the geometric ids in the overlap regions.
// is there a way to use the isOverlap functor to do this. This allows for
// maximum flexibility moving forward

void omegah_coupler(MPI_Comm comm, Omega_h::Mesh& mesh,
                    std::string_view cpn_file, int nphi)
{
  std::chrono::duration<double> elapsed_seconds;
  double min, max, avg;
  int rank;
  MPI_Comm_rank(comm, &rank);
  auto time1 = std::chrono::steady_clock::now();
  // Should construct the global ids
  //Omega_h::reorder_by_hilbert(&mesh); 
  //auto ids = mesh.get_array<LO>(0,"simNumbering");
  //auto gids_array = Omega_h::Write<GO>(ids.size());
  //Omega_h::parallel_for(ids.size(),OMEGA_H_LAMBDA(int i){
  //  gids_array[i] = ids[i] - 1;
  //});
  //if(mesh.has_tag(0,"global")) {
  //  mesh.set_tag<GO>(0, "global", gids_array); 
  //} else {
  //  mesh.add_tag<GO>(0, "global", 1, gids_array); 
  //}


  wdmcpl::CouplerServer cpl("xgc_n0_coupling", comm,
                            redev::Partition{ts::setupServerPartition(mesh, cpn_file)}, mesh);
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  std::string numbering = "simNumbering";
  WDMCPL_ALWAYS_ASSERT(mesh.has_tag(0, numbering));

  auto is_overlap = ts::markServerOverlapRegion(
    mesh, partition, KOKKOS_LAMBDA(const int dim, const int id) {
      //if (id >= 1 && id <= 2) {
      //  return 1;
      //}
      if (id >= 100 && id <= 140) {
        return 1;
      }
      return 0;
    });
  auto time2 = std::chrono::steady_clock::now();
  elapsed_seconds = time2-time1;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if(!rank) ts::printTime("Initialize Coupler/Mesh", min, max, avg);
  std::vector<wdmcpl::ConvertibleCoupledField*> potential_fields(nphi * 2);
  std::vector<wdmcpl::ScatterOperation*> scatter_ops(nphi);
  std::vector<wdmcpl::GatherOperation*> gather_ops(nphi);
  std::cerr << "ADDING FIELDS\n";
  for (int i = 0; i < nphi; ++i) {
    std::stringstream field1_name;
    field1_name << "dpot_plane_" << i;
    std::cerr << field1_name.str() << "\n";
    potential_fields[2 * i] =
      cpl.AddField(field1_name.str(),
                   wdmcpl::OmegaHFieldAdapter<wdmcpl::Real>(
                     field1_name.str(), mesh, is_overlap, numbering),
                   FieldTransferMethod::Copy, // to Omega_h
                   FieldEvaluationMethod::None,
                   FieldTransferMethod::Copy, // from Omega_h
                   FieldEvaluationMethod::None, is_overlap, "deltaf1/");
    potential_fields[2 * i + 1] =
      cpl.AddField(field1_name.str(),
                   wdmcpl::OmegaHFieldAdapter<wdmcpl::Real>(
                     field1_name.str(), mesh, is_overlap, numbering),
                   FieldTransferMethod::Copy, // to Omega_h
                   FieldEvaluationMethod::None,
                   FieldTransferMethod::Copy, // from Omega_h
                   FieldEvaluationMethod::None, is_overlap, "deltaf2/");
    gather_ops[i] = cpl.AddGatherFieldsOp(
      field1_name.str(),
      {"deltaf1/" + field1_name.str(), "deltaf2/" + field1_name.str()},
      "combined" + field1_name.str(), ts::MeanCombiner{}, is_overlap);
    scatter_ops[i] = cpl.AddScatterFieldsOp(
      field1_name.str(), "combined" + field1_name.str(),
      {"deltaf1/" + field1_name.str(), "deltaf2/" + field1_name.str()},
      is_overlap);
    std::cerr << field1_name.str() << " DONE\n";
  }
  auto time3 = std::chrono::steady_clock::now();
  elapsed_seconds = time3-time2;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if(!rank) ts::printTime("Add Meshes", min, max, avg);
  while (true) {
    auto sr_time1 = std::chrono::steady_clock::now();
    for(auto* gather : gather_ops) {
      gather->Run();
    }
    auto sr_time2 = std::chrono::steady_clock::now();
    elapsed_seconds = sr_time2-sr_time1;
    ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
    if(!rank) ts::printTime("Gather", min, max, avg);
    for(auto* scatter : scatter_ops) {
      scatter->Run();
    }
    auto sr_time3 = std::chrono::steady_clock::now();
    elapsed_seconds = sr_time3-sr_time2;
    ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
    if(!rank) ts::printTime("Scatter", min, max, avg);
    // TODO: do this w/o blocking...
    //for (auto* field : potential_fields) {
    //  WDMCPL_ALWAYS_ASSERT(field != nullptr);
    //  field->Receive();
    //}
    //for (auto* field : potential_fields) {
    //  field->Send();
    //}
  }
}

int main(int argc, char** argv)
{
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  int size = world->size();
  if (argc != 4) {
    if (!rank) {
      std::cerr << "Usage: " << argv[0]
                << "</path/to/omega_h/mesh> "
                   "</path/to/partitionFile.cpn> "
                   "sml_nphi_total";
    }
    exit(EXIT_FAILURE);
  }

  const auto meshFile = argv[1];
  const auto classPartitionFile = argv[2];
  const int sml_nphi_total = std::atoi(argv[3]);

  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(meshFile, lib.world(), &mesh);
  MPI_Comm mpi_comm = lib.world()->get_impl();
  omegah_coupler(mpi_comm, mesh, classPartitionFile, sml_nphi_total);
  return 0;
}