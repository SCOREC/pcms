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
//

[[nodiscard]]
static wdmcpl::ConvertibleCoupledField* AddField(wdmcpl::Application *application, const std::string& name, const std::string& path, Omega_h::Read<Omega_h::I8> is_overlap, const std::string& numbering, Omega_h::Mesh& mesh, int plane) {
      WDMCPL_ALWAYS_ASSERT(application != nullptr);
      std::stringstream field_name;
      field_name << path << name << "_" << plane;
      return application->AddField(field_name.str(),
                   wdmcpl::OmegaHFieldAdapter<wdmcpl::Real>(
                     field_name.str(), mesh, is_overlap, numbering),
                   FieldTransferMethod::Copy, // to Omega_h
                   FieldEvaluationMethod::None,
                   FieldTransferMethod::Copy, // from Omega_h
                   FieldEvaluationMethod::None, is_overlap);
}

struct XGCAnalysis {
  using FieldVec = std::vector<wdmcpl::ConvertibleCoupledField*>;
  std::array<FieldVec,2> dpot;
  FieldVec pot0;
  std::array<FieldVec,2> edensity;
  std::array<FieldVec,2> idensity;
};

static void ReceiveFields(const std::vector<wdmcpl::ConvertibleCoupledField*> & fields) {
  for(auto* field : fields) {
    field->Receive();
  }
}
static void SendFields(const std::vector<wdmcpl::ConvertibleCoupledField*> & fields) {
  for(auto* field : fields) {
    field->Send();
  }
}
static void CopyFields(const std::vector<wdmcpl::ConvertibleCoupledField*> & from_fields,
                       const std::vector<wdmcpl::ConvertibleCoupledField*> & to_fields) {
  WDMCPL_ALWAYS_ASSERT(from_fields.size() == to_fields.size());
  for(size_t i=0; i<from_fields.size(); ++i) {
    const auto* from = from_fields[i]->GetFieldAdapter<wdmcpl::OmegaHFieldAdapter<wdmcpl::Real>>();
    auto* to = to_fields[i]->GetFieldAdapter<wdmcpl::OmegaHFieldAdapter<wdmcpl::Real>>();
    copy_field(from->GetField(), to->GetField());
  }
}

using GatherScatterFieldsVec = std::vector<std::reference_wrapper<wdmcpl::ConvertibleCoupledField>>;

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
  auto* core = cpl.AddApplication("core", "core/");
  auto* edge = cpl.AddApplication("edge", "core/");
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

  std::vector<wdmcpl::ScatterOperation*> scatter_density_ops;
  std::vector<wdmcpl::GatherOperation*> gather_density_ops(nphi);
  XGCAnalysis core_analysis;
  XGCAnalysis edge_analysis;
  std::cerr << "ADDING FIELDS\n";
  for (int i = 0; i < nphi; ++i) {
    core_analysis.dpot[0].push_back(AddField(core, "dpot_1_plane", "core/", 
                                             is_overlap, numbering, mesh, i));
    core_analysis.dpot[1].push_back(AddField(core, "dpot_2_plane", "core/", 
                                             is_overlap, numbering, mesh, i));
    core_analysis.pot0.push_back(AddField(core, "pot0_plane", "core/",
                                          is_overlap, numbering, mesh, i));
    core_analysis.edensity[0].push_back(AddField(core, "edensity_1_plane", "core/", 
                                             is_overlap, numbering, mesh, i));
    core_analysis.edensity[1].push_back(AddField(core, "edensity_2_plane", "core/", 
                                             is_overlap, numbering, mesh, i));
    core_analysis.idensity[0].push_back(AddField(core, "idensity_1_plane", "core/", 
                                             is_overlap, numbering, mesh, i));
    core_analysis.idensity[1].push_back(AddField(core, "idensity_2_plane", "core/", 
                                             is_overlap, numbering, mesh, i));

    edge_analysis.dpot[0].push_back(AddField(edge, "dpot_1_plane", "edge/", 
                                             is_overlap, numbering, mesh, i));
    edge_analysis.dpot[1].push_back(AddField(edge, "dpot_2_plane", "edge/", 
                                             is_overlap, numbering, mesh, i));
    edge_analysis.pot0.push_back(AddField(edge, "pot0_plane", "edge/",
                                          is_overlap, numbering, mesh, i));
    edge_analysis.edensity[0].push_back(AddField(edge, "edensity_1_plane", "edge/", 
                                             is_overlap, numbering, mesh, i));
    edge_analysis.edensity[1].push_back(AddField(edge, "edensity_2_plane", "edge/", 
                                             is_overlap, numbering, mesh, i));
    edge_analysis.idensity[0].push_back(AddField(edge, "idensity_1_plane", "edge/", 
                                             is_overlap, numbering, mesh, i));
    edge_analysis.idensity[1].push_back(AddField(edge, "idensity_2_plane", "edge/", 
                                             is_overlap, numbering, mesh, i));

    // gather density from core/edge to density
    //const std::string edensity_name = "edensity_1_plane_"+std::to_string(i);
    gather_density_ops.push_back(cpl.AddGatherFieldsOp(
      "edensity_1_plane_"+std::to_string(i),
      {*core_analysis.edensity[0].back(),
       *edge_analysis.edensity[0].back() },
      "combined_edensity_1_plane_"+std::to_string(i),
      ts::MeanCombiner{}, is_overlap));
    gather_density_ops.push_back(cpl.AddGatherFieldsOp(
      "edensity_2_plane_"+std::to_string(i),
      {*core_analysis.edensity[1].back(),
       *edge_analysis.edensity[1].back() },
      "combined_edensity_2_plane_"+std::to_string(i),
      ts::MeanCombiner{}, is_overlap));
    gather_density_ops.push_back(cpl.AddGatherFieldsOp(
      "idensity_1_plane_"+std::to_string(i),
      {*core_analysis.idensity[0].back(),
       *edge_analysis.idensity[0].back() },
      "combined_idensity_1_plane_"+std::to_string(i),
      ts::MeanCombiner{}, is_overlap));
    gather_density_ops.push_back(cpl.AddGatherFieldsOp(
      "idensity_2_plane_"+std::to_string(i),
      {*core_analysis.idensity[1].back(),
       *edge_analysis.idensity[1].back() },
      "combined_idensity_2_plane_"+std::to_string(i),
      ts::MeanCombiner{}, is_overlap));

    scatter_density_ops.push_back(cpl.AddScatterFieldsOp(
      "edensity_1_plane_"+std::to_string(i),
      "combined_edensity_1_plane_"+std::to_string(i),
      {*edge_analysis.edensity[0].back()},is_overlap));
    scatter_density_ops.push_back(cpl.AddScatterFieldsOp(
      "edensity_2_plane_"+std::to_string(i),
      "combined_edensity_2_plane_"+std::to_string(i),
      {*edge_analysis.edensity[1].back()},is_overlap));
    scatter_density_ops.push_back(cpl.AddScatterFieldsOp(
      "idensity_1_plane_"+std::to_string(i),
      "combined_idensity_1_plane_"+std::to_string(i),
      {*edge_analysis.idensity[0].back()},is_overlap));
    scatter_density_ops.push_back(cpl.AddScatterFieldsOp(
      "idensity_2_plane_"+std::to_string(i),
      "combined_idensity_2_plane_"+std::to_string(i),
      {*edge_analysis.idensity[1].back()},is_overlap));
  }
  auto time3 = std::chrono::steady_clock::now();
  elapsed_seconds = time3-time2;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if(!rank) ts::printTime("Add Meshes", min, max, avg);
  while (true) {
    auto sr_time1 = std::chrono::steady_clock::now();
    // gather density fields (Core+Edge)
    for(auto* gather : gather_density_ops) {
      gather->Run();
    }
    auto sr_time2 = std::chrono::steady_clock::now();
    elapsed_seconds = sr_time2-sr_time1;
    ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
    if(!rank) ts::printTime("Gather", min, max, avg);
    // Scatter density field (Edge)
    for(auto* scatter : scatter_density_ops) {
      scatter->Run();
    }
    auto sr_time3 = std::chrono::steady_clock::now();
    elapsed_seconds = sr_time3-sr_time2;
    ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
    if(!rank) ts::printTime("Scatter", min, max, avg);

    // deal with phi fields (pot0/dpot1/dpot2)
    // 1. reveive fields from Edge
    ReceiveFields(edge_analysis.dpot[0]);
    ReceiveFields(edge_analysis.dpot[1]);
    ReceiveFields(edge_analysis.pot0);
    // 2. Copy fields from Edge->Core
    // Since we know that both native fields are Omega_h,
    // we don't need to sync to the internal and combined
    // fields that would be created by using Gather/Scatter operation
    // This is a bit "hacky" and will use our intimate knowledge of how
    // OmegaHField adapters work. i.e., every field adapter with the same
    // name refers to the same OmegaH field
    CopyFields(edge_analysis.dpot[0], core_analysis.dpot[0]);
    CopyFields(edge_analysis.dpot[1], core_analysis.dpot[1]);
    CopyFields(edge_analysis.pot0, core_analysis.pot0);
    // 3. Send fields to Core
    SendFields(edge_analysis.dpot[0]);
    SendFields(edge_analysis.dpot[1]);
    SendFields(edge_analysis.pot0);

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
