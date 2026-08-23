#include <pcms.h>
#include <pcms/utility/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include "test_support.h"
#include "pcms/coupler/coupler.hpp"
#include "pcms/coupler/field_serializer.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/field/field.h"
#include "pcms/field/value_view.hpp"
#include "pcms/field/data/simple.h"
#include "pcms/transfer/copy.h"
#include <chrono>
#include <optional>
#include "pcms/field/coordinate_systems/cartesian.hpp"

using pcms::GO;
using pcms::LO;

namespace ts = test_support;

// TODO: we should communicate the geometric ids in the overlap regions.
// is there a way to use the isOverlap functor to do this. This allows for
// maximum flexibility moving forward
//

[[nodiscard]]
static std::string MakeFieldName(const std::string& name, int plane)
{
  std::stringstream field_name;
  field_name << name;
  if (plane >= 0) {
    field_name << "_" << plane;
  }
  return field_name.str();
}

struct RegisteredField
{
  pcms::FieldHandle<pcms::Real> handle;
};

[[nodiscard]]
static RegisteredField AddField(
  pcms::Application* application,
  const std::shared_ptr<pcms::LagrangeFunctionSpace>& function_space,
  const std::string& name, const std::string& path, int plane)
{
  PCMS_ALWAYS_ASSERT(application != nullptr);
  auto field_name = MakeFieldName(name, plane);
  auto field = function_space->CreateFunction<pcms::Real>(path + field_name);
  std::unique_ptr<pcms::FieldSerializer<pcms::Real>> serializer =
    std::make_unique<pcms::FieldSerializer<pcms::Real>>();
  auto handle = application->AddField(std::move(field), std::move(serializer));
  return {std::move(handle)};
}

struct XGCAnalysis
{
  using FieldVec = std::vector<RegisteredField>;
  std::array<FieldVec, 2> dpot;
  FieldVec pot0;
  std::array<FieldVec, 2> edensity;
  std::array<FieldVec, 2> idensity;
  std::optional<RegisteredField> psi;
  std::optional<RegisteredField> gids;
};

static void ReceiveFields(const std::vector<RegisteredField>& fields)
{
  for (const auto& field : fields) {
    field.handle.Receive();
  }
}
static void SendFields(const std::vector<RegisteredField>& fields)
{
  for (const auto& field : fields) {
    field.handle.Send();
  }
}
static void CopyFields(const std::vector<RegisteredField>& from_fields,
                       const std::vector<RegisteredField>& to_fields)
{
  PCMS_ALWAYS_ASSERT(from_fields.size() == to_fields.size());
  for (size_t i = 0; i < from_fields.size(); ++i) {
    auto& source = from_fields[i].handle.GetField();
    auto& target = to_fields[i].handle.GetField();
    target.SetDOFHolderDataHost(source.GetDOFHolderDataHost());
  }
}

void SendRecvDensity(pcms::Application* core, pcms::Application* edge,
                     XGCAnalysis& core_analysis, XGCAnalysis& edge_analysis,
                     int rank)
{

  std::chrono::duration<double> elapsed_seconds;
  double min, max, avg;
  if (!rank)
    std::cerr << "Send/Recv Density\n";
  auto sr_time1 = std::chrono::steady_clock::now();
  // gather density fields (Core+Edge)
  core->BeginReceivePhase();
  edge->BeginReceivePhase();
  // Gather
  ReceiveFields(core_analysis.edensity[0]);
  ReceiveFields(core_analysis.edensity[1]);
  ReceiveFields(edge_analysis.edensity[0]);
  ReceiveFields(edge_analysis.edensity[1]);
  ReceiveFields(core_analysis.idensity[0]);
  ReceiveFields(core_analysis.idensity[1]);
  ReceiveFields(edge_analysis.idensity[0]);
  ReceiveFields(edge_analysis.idensity[1]);

  core->EndReceivePhase();
  edge->EndReceivePhase();
  auto sr_time2 = std::chrono::steady_clock::now();
  elapsed_seconds = sr_time2 - sr_time1;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if (!rank)
    ts::printTime("Recv Density", min, max, avg);

  CopyFields(core_analysis.edensity[0], edge_analysis.edensity[0]);
  CopyFields(core_analysis.edensity[1], edge_analysis.edensity[1]);
  CopyFields(core_analysis.idensity[0], edge_analysis.idensity[0]);
  CopyFields(core_analysis.idensity[1], edge_analysis.idensity[1]);

  sr_time1 = std::chrono::steady_clock::now();
  elapsed_seconds = sr_time1 - sr_time2;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if (!rank)
    ts::printTime("Average Density", min, max, avg);
  edge->BeginSendPhase();
  SendFields(edge_analysis.edensity[0]);
  SendFields(edge_analysis.edensity[1]);
  SendFields(edge_analysis.idensity[0]);
  SendFields(edge_analysis.idensity[1]);
  edge->EndSendPhase();
  auto sr_time3 = std::chrono::steady_clock::now();
  elapsed_seconds = sr_time3 - sr_time1;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if (!rank)
    ts::printTime("Send Density", min, max, avg);
}
void SendRecvPotential(pcms::Application* core, pcms::Application* edge,
                       XGCAnalysis& core_analysis, XGCAnalysis& edge_analysis,
                       int rank)
{

  std::chrono::duration<double> elapsed_seconds;
  double min, max, avg;
  if (!rank)
    std::cerr << "Send/Recv Potential\n";
  auto sr_time3 = std::chrono::steady_clock::now();
  edge->BeginReceivePhase();
  // deal with phi fields (pot0/dpot1/dpot2)
  // 1. reveive fields from Edge
  for (auto& f : edge_analysis.dpot) {
    ReceiveFields(f);
  }
  ReceiveFields(edge_analysis.pot0);
  // core->EndReceivePhase();
  edge->EndReceivePhase();
  auto sr_time4 = std::chrono::steady_clock::now();
  elapsed_seconds = sr_time4 - sr_time3;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if (!rank)
    ts::printTime("Receive Potential", min, max, avg);
  // 2. Copy fields from Edge->Core
  for (int i = 0; i < edge_analysis.dpot.size(); ++i) {
    CopyFields(edge_analysis.dpot[i], core_analysis.dpot[i]);
    CopyFields(edge_analysis.dpot[i], core_analysis.dpot[i]);
  }
  CopyFields(edge_analysis.pot0, core_analysis.pot0);
  auto sr_time5 = std::chrono::steady_clock::now();
  elapsed_seconds = sr_time5 - sr_time4;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if (!rank)
    ts::printTime("Copy Potential", min, max, avg);
  core->BeginSendPhase();
  for (auto& f : core_analysis.dpot) {
    SendFields(f);
  }
  SendFields(core_analysis.pot0);
  core->EndSendPhase();
  auto sr_time6 = std::chrono::steady_clock::now();
  elapsed_seconds = sr_time6 - sr_time5;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if (!rank)
    ts::printTime("Send Potential", min, max, avg);
}

void omegah_coupler(MPI_Comm comm, Omega_h::Mesh& mesh,
                    std::string_view cpn_file, int nphi)
{
  std::chrono::duration<double> elapsed_seconds;
  double min, max, avg;
  int rank;
  MPI_Comm_rank(comm, &rank);
  auto time1 = std::chrono::steady_clock::now();

  pcms::Coupler cpl("xgc_n0_coupling", comm, true,
                    redev::Partition{ts::setupServerPartition(mesh, cpn_file)});
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  std::string numbering = "simNumbering";
  PCMS_ALWAYS_ASSERT(mesh.has_tag(0, numbering));
  auto* core = cpl.AddApplication("core", "core/");
  auto* edge = cpl.AddApplication("edge", "edge/");
  auto is_overlap = ts::markServerOverlapRegion(
    mesh, partition, KOKKOS_LAMBDA(const int dim, const int id) {
      // if (id >= 1 && id <= 2) {
      //   return 1;
      // }
      // if (id >= 100 && id <= 140) {
      //   return 1;
      // }
      // return 0;
      return 1;
    });
  auto function_space = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred(), is_overlap, numbering,
    pcms::LagrangeFunctionSpace::Backend::OmegaH, "n0_layout");
  auto time2 = std::chrono::steady_clock::now();
  elapsed_seconds = time2 - time1;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if (!rank)
    ts::printTime("Initialize Coupler/Mesh", min, max, avg);

  XGCAnalysis core_analysis;
  XGCAnalysis edge_analysis;
  std::cerr << "ADDING FIELDS\n";
  for (int i = 0; i < nphi; ++i) {
    // core_analysis.dpot[0].push_back(AddField(core, "dpot_m1_plane", "core/",
    //                                          is_overlap, numbering, mesh,
    //                                          i));
    core_analysis.dpot[0].push_back(
      AddField(core, function_space, "dpot_0_plane", "core/", i));
    core_analysis.dpot[1].push_back(
      AddField(core, function_space, "dpot_1_plane", "core/", i));
    // core_analysis.dpot[3].push_back(AddField(core, "dpot_2_plane", "core/",
    //                                          is_overlap, numbering, mesh,
    //                                          i));
    core_analysis.pot0.push_back(
      AddField(core, function_space, "pot0_plane", "core/", i));
    core_analysis.edensity[0].push_back(
      AddField(core, function_space, "edensity_1_plane", "core/", i));
    core_analysis.edensity[1].push_back(
      AddField(core, function_space, "edensity_2_plane", "core/", i));
    core_analysis.idensity[0].push_back(
      AddField(core, function_space, "idensity_1_plane", "core/", i));
    core_analysis.idensity[1].push_back(
      AddField(core, function_space, "idensity_2_plane", "core/", i));

    // edge_analysis.dpot[0].push_back(AddField(edge, "dpot_m1_plane", "edge/",
    //                                          is_overlap, numbering, mesh,
    //                                          i));
    edge_analysis.dpot[0].push_back(
      AddField(edge, function_space, "dpot_0_plane", "edge/", i));
    edge_analysis.dpot[1].push_back(
      AddField(edge, function_space, "dpot_1_plane", "edge/", i));
    // edge_analysis.dpot[3].push_back(AddField(edge, "dpot_2_plane", "edge/",
    //                                          is_overlap, numbering, mesh,
    //                                          i));
    edge_analysis.pot0.push_back(
      AddField(edge, function_space, "pot0_plane", "edge/", i));
    edge_analysis.edensity[0].push_back(
      AddField(edge, function_space, "edensity_1_plane", "edge/", i));
    edge_analysis.edensity[1].push_back(
      AddField(edge, function_space, "edensity_2_plane", "edge/", i));
    edge_analysis.idensity[0].push_back(
      AddField(edge, function_space, "idensity_1_plane", "edge/", i));
    edge_analysis.idensity[1].push_back(
      AddField(edge, function_space, "idensity_2_plane", "edge/", i));
  }
  core_analysis.psi = AddField(core, function_space, "psi", "core/", -1);
  edge_analysis.psi = AddField(edge, function_space, "psi", "edge/", -1);
  core_analysis.gids = AddField(core, function_space, "gid_debug", "core/", -1);
  edge_analysis.gids = AddField(edge, function_space, "gid_debug", "edge/", -1);
  auto time3 = std::chrono::steady_clock::now();
  elapsed_seconds = time3 - time2;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if (!rank)
    ts::printTime("Add Meshes", min, max, avg);

  Omega_h::vtk::write_parallel("initial.vtk", &mesh);
  edge->BeginReceivePhase();
  edge_analysis.psi->handle.Receive();
  edge_analysis.gids->handle.Receive();
  edge->EndReceivePhase();
  core->BeginReceivePhase();
  core_analysis.psi->handle.Receive();
  core_analysis.gids->handle.Receive();
  core->EndReceivePhase();
  Omega_h::vtk::write_parallel("psi-only.vtk", &mesh);
  auto time4 = std::chrono::steady_clock::now();
  elapsed_seconds = time4 - time3;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if (!rank)
    ts::printTime("Receive Psi", min, max, avg);
  int step = 0;
  while (true) {
    std::stringstream ss;
    SendRecvDensity(core, edge, core_analysis, edge_analysis, rank);
    SendRecvPotential(core, edge, core_analysis, edge_analysis, rank);
    ss << "step-" << step++ << ".vtk";
    Omega_h::vtk::write_parallel(ss.str(), &mesh);
  }
}

int main(int argc, char** argv)
{
  try {
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
  } catch (const std::exception& e) {
    std::cerr << "Exception caught in main: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception caught in main" << std::endl;
    return 1;
  }
}
