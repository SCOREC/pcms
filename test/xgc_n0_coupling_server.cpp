#include <pcms.h>
#include <pcms/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include "test_support.h"
#include <pcms/omega_h_field.h>
#include <pcms/xgc_field_adapter.h>
#include <chrono>

using pcms::Copy;
using pcms::CouplerClient;
using pcms::CouplerServer;
using pcms::GO;
using pcms::LO;
using pcms::OmegaHFieldAdapter;

namespace ts = test_support;

// TODO: we should communicate the geometric ids in the overlap regions.
// is there a way to use the isOverlap functor to do this. This allows for
// maximum flexibility moving forward
//

[[nodiscard]]
static pcms::ConvertibleCoupledField* AddField(pcms::Application *application, const std::string& name, const std::string& path, Omega_h::Read<Omega_h::I8> is_overlap, const std::string& numbering, Omega_h::Mesh& mesh, int plane) {
      PCMS_ALWAYS_ASSERT(application != nullptr);
      std::stringstream field_name;
      field_name << name;
      if(plane >=0) {
        field_name<< "_" << plane;
      }
      return application->AddField(field_name.str(),
                   pcms::OmegaHFieldAdapter<pcms::Real>(
                   path+field_name.str(), mesh, is_overlap, numbering),
                   is_overlap);
}

struct XGCAnalysis {
  using FieldVec = std::vector<pcms::ConvertibleCoupledField*>;
  std::array<FieldVec,2> dpot;
  FieldVec pot0;
  std::array<FieldVec,2> edensity;
  std::array<FieldVec,2> idensity;
  pcms::ConvertibleCoupledField* psi;
  pcms::ConvertibleCoupledField* gids;
};

static void ReceiveFields(const std::vector<pcms::ConvertibleCoupledField*> & fields) {
  for(auto* field : fields) {
    field->Receive();
  }
}
static void SendFields(const std::vector<pcms::ConvertibleCoupledField*> & fields) {
  for(auto* field : fields) {
    field->Send();
  }
}
static void CopyFields(const std::vector<pcms::ConvertibleCoupledField*> & from_fields,
                       const std::vector<pcms::ConvertibleCoupledField*> & to_fields) {
  PCMS_ALWAYS_ASSERT(from_fields.size() == to_fields.size());
  for(size_t i=0; i<from_fields.size(); ++i) {
    const auto* from = from_fields[i]->GetFieldAdapter<pcms::OmegaHFieldAdapter<pcms::Real>>();
    auto* to = to_fields[i]->GetFieldAdapter<pcms::OmegaHFieldAdapter<pcms::Real>>();
    copy_field(from->GetField(), to->GetField());
  }
}

template <typename T>
static void AverageAndSetField(const pcms::OmegaHField<T> & a, pcms::OmegaHField<T> & b) {
  const auto a_data = get_nodal_data(a);
  const auto b_data = get_nodal_data(b);
  Omega_h::Write<T> combined_data(a_data.size());
  Omega_h::parallel_for(combined_data.size(), OMEGA_H_LAMBDA(size_t i) {
    combined_data[i] = (a_data[i] + b_data[i]) / 2.0;
  });
  auto combined_view = pcms::make_array_view(Omega_h::Read<T>(combined_data));
  pcms::set_nodal_data(b, combined_view);
}

/*
 * Takes the average of each pair of fields and sets the results in the the second
 * argument
 */
static void AverageAndSetFields(const std::vector<pcms::ConvertibleCoupledField*> & from_fields,
                       const std::vector<pcms::ConvertibleCoupledField*> & to_fields) {
  PCMS_ALWAYS_ASSERT(from_fields.size() == to_fields.size());
  for(size_t i=0; i<from_fields.size(); ++i) {
    const auto* from = from_fields[i]->GetFieldAdapter<pcms::OmegaHFieldAdapter<pcms::Real>>();
    auto* to = to_fields[i]->GetFieldAdapter<pcms::OmegaHFieldAdapter<pcms::Real>>();
    AverageAndSetField(from->GetField(),to->GetField());
  }
}

void SendRecvDensity(pcms::Application* core, pcms::Application* edge, XGCAnalysis& core_analysis, XGCAnalysis& edge_analysis, int rank) {

    std::chrono::duration<double> elapsed_seconds;
    double min, max, avg;
    if(!rank) std::cerr<<"Send/Recv Density\n"; 
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
    elapsed_seconds = sr_time2-sr_time1;
    ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
    if(!rank) ts::printTime("Recv Density", min, max, avg);

    CopyFields(core_analysis.edensity[0], edge_analysis.edensity[0]);
    CopyFields(core_analysis.edensity[1], edge_analysis.edensity[1]);
    CopyFields(core_analysis.idensity[0], edge_analysis.idensity[0]);
    CopyFields(core_analysis.idensity[1], edge_analysis.idensity[1]);
    //AverageAndSetFields(core_analysis.edensity[0], edge_analysis.edensity[0]);
    //AverageAndSetFields(core_analysis.edensity[1], edge_analysis.edensity[1]);
    //AverageAndSetFields(core_analysis.idensity[0], edge_analysis.idensity[0]);
    //AverageAndSetFields(core_analysis.idensity[1], edge_analysis.idensity[1]);

    sr_time1 = std::chrono::steady_clock::now();
    elapsed_seconds = sr_time1-sr_time2;
    ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
    if(!rank) ts::printTime("Average Density", min, max, avg);
    edge->BeginSendPhase();
    SendFields(edge_analysis.edensity[0]);
    SendFields(edge_analysis.edensity[1]);
    SendFields(edge_analysis.idensity[0]);
    SendFields(edge_analysis.idensity[1]);
    edge->EndSendPhase();
    auto sr_time3 = std::chrono::steady_clock::now();
    elapsed_seconds = sr_time3-sr_time1;
    ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
    if(!rank) ts::printTime("Send Density", min, max, avg);
}
void SendRecvPotential(pcms::Application* core, pcms::Application* edge, XGCAnalysis& core_analysis, XGCAnalysis& edge_analysis, int rank) {

    std::chrono::duration<double> elapsed_seconds;
     double min, max, avg;
    if(!rank) std::cerr<<"Send/Recv Potential\n"; 
    auto sr_time3 = std::chrono::steady_clock::now();
    edge->BeginReceivePhase();
    // deal with phi fields (pot0/dpot1/dpot2)
    // 1. reveive fields from Edge
    for(auto& f: edge_analysis.dpot) {
      ReceiveFields(f);
    }
    ReceiveFields(edge_analysis.pot0);
    //core->EndReceivePhase();
    edge->EndReceivePhase();
    auto sr_time4 = std::chrono::steady_clock::now();
    elapsed_seconds = sr_time4-sr_time3;
    ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
    if(!rank) ts::printTime("Receive Potential", min, max, avg);
    // 2. Copy fields from Edge->Core
    for(int i=0; i<edge_analysis.dpot.size(); ++i){
      CopyFields(edge_analysis.dpot[i], core_analysis.dpot[i]);
      CopyFields(edge_analysis.dpot[i], core_analysis.dpot[i]);
    }
    CopyFields(edge_analysis.pot0, core_analysis.pot0);
    auto sr_time5 = std::chrono::steady_clock::now();
    elapsed_seconds = sr_time5-sr_time4;
    ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
    if(!rank) ts::printTime("Copy Potential", min, max, avg);
    core->BeginSendPhase();
    for(auto& f: core_analysis.dpot) {
      SendFields(f);
    }
    SendFields(core_analysis.pot0);
    core->EndSendPhase();
    auto sr_time6 = std::chrono::steady_clock::now();
    elapsed_seconds = sr_time6-sr_time5;
    ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
    if(!rank) ts::printTime("Send Potential", min, max, avg);
}

void omegah_coupler(MPI_Comm comm, Omega_h::Mesh& mesh,
                    std::string_view cpn_file, int nphi)
{
  std::chrono::duration<double> elapsed_seconds;
  double min, max, avg;
  int rank;
  MPI_Comm_rank(comm, &rank);
  auto time1 = std::chrono::steady_clock::now();


  pcms::CouplerServer cpl("xgc_n0_coupling", comm,
                            redev::Partition{ts::setupServerPartition(mesh, cpn_file)}, mesh);
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  std::string numbering = "simNumbering";
  PCMS_ALWAYS_ASSERT(mesh.has_tag(0, numbering));
  auto* core = cpl.AddApplication("core", "core/");
  auto* edge = cpl.AddApplication("edge", "edge/");
  auto is_overlap = ts::markServerOverlapRegion(
    mesh, partition, KOKKOS_LAMBDA(const int dim, const int id) {
      //if (id >= 1 && id <= 2) {
      //  return 1;
      //}
      //if (id >= 100 && id <= 140) {
      //  return 1;
      //}
      //return 0;
      return 1;
    });
  auto time2 = std::chrono::steady_clock::now();
  elapsed_seconds = time2-time1;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if(!rank) ts::printTime("Initialize Coupler/Mesh", min, max, avg);

  XGCAnalysis core_analysis;
  XGCAnalysis edge_analysis;
  std::cerr << "ADDING FIELDS\n";
  for (int i = 0; i < nphi; ++i) {
    //core_analysis.dpot[0].push_back(AddField(core, "dpot_m1_plane", "core/", 
    //                                         is_overlap, numbering, mesh, i));
    core_analysis.dpot[0].push_back(AddField(core, "dpot_0_plane", "core/", 
                                             is_overlap, numbering, mesh, i));
    core_analysis.dpot[1].push_back(AddField(core, "dpot_1_plane", "core/", 
                                             is_overlap, numbering, mesh, i));
    //core_analysis.dpot[3].push_back(AddField(core, "dpot_2_plane", "core/", 
    //                                         is_overlap, numbering, mesh, i));
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

    //edge_analysis.dpot[0].push_back(AddField(edge, "dpot_m1_plane", "edge/", 
    //                                         is_overlap, numbering, mesh, i));
    edge_analysis.dpot[0].push_back(AddField(edge, "dpot_0_plane", "edge/", 
                                             is_overlap, numbering, mesh, i));
    edge_analysis.dpot[1].push_back(AddField(edge, "dpot_1_plane", "edge/", 
                                             is_overlap, numbering, mesh, i));
    //edge_analysis.dpot[3].push_back(AddField(edge, "dpot_2_plane", "edge/", 
    //                                         is_overlap, numbering, mesh, i));
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

  }
  core_analysis.psi = AddField(core, "psi", "core/", 
                                           is_overlap, numbering, mesh, -1);
  edge_analysis.psi = AddField(edge, "psi", "edge/", 
                                           is_overlap, numbering, mesh, -1);
  core_analysis.gids = AddField(core, "gid_debug", "core/", 
                                           is_overlap, numbering, mesh, -1);
  edge_analysis.gids = AddField(edge, "gid_debug", "edge/", 
                                           is_overlap, numbering, mesh, -1);
  auto time3 = std::chrono::steady_clock::now();
  elapsed_seconds = time3-time2;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if(!rank) ts::printTime("Add Meshes", min, max, avg);

  Omega_h::vtk::write_parallel("initial.vtk", &mesh);
  edge->BeginReceivePhase();
  edge_analysis.psi->Receive();
  edge_analysis.gids->Receive();
  edge->EndReceivePhase();
  core->BeginReceivePhase();
  core_analysis.psi->Receive();
  core_analysis.gids->Receive();
  core->EndReceivePhase();
  Omega_h::vtk::write_parallel("psi-only.vtk", &mesh);
  auto time4 = std::chrono::steady_clock::now();
  elapsed_seconds = time4-time3;
  ts::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if(!rank) ts::printTime("Receive Psi", min, max, avg);
  int step = 0;
  while (true) {
    std::stringstream ss;
    SendRecvDensity(core, edge, core_analysis, edge_analysis, rank);
    SendRecvPotential(core, edge, core_analysis, edge_analysis, rank);
    ss <<"step-"<<step++ <<".vtk";
    Omega_h::vtk::write_parallel(ss.str(), &mesh);
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
