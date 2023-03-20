#include <iostream>
#include <wdmcpl.h>
#include <wdmcpl/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include "test_support.h"
#include <wdmcpl/omega_h_field.h>
#include <wdmcpl/xgc_field_adapter.h>

using wdmcpl::ConstructRCFromOmegaHMesh;
using wdmcpl::Copy;
using wdmcpl::CouplerClient;
using wdmcpl::CouplerServer;
using wdmcpl::FieldEvaluationMethod;
using wdmcpl::FieldTransferMethod;
using wdmcpl::GO;
using wdmcpl::Lagrange;
using wdmcpl::make_array_view;
using wdmcpl::OmegaHField;
using wdmcpl::OmegaHFieldAdapter;
using wdmcpl::ReadReverseClassificationVertex;
using wdmcpl::ReverseClassificationVertex;

static constexpr bool done = true;
namespace ts = test_support;

void xgc_coupler(MPI_Comm comm, Omega_h::Mesh& mesh, std::string_view cpn_file)
{
  // coupling server using same mesh as application
  // note the xgc_coupler stores a reference to the internal mesh and it is the
  // user responsibility to keep it alive!
  wdmcpl::CouplerServer cpl(
    "proxy_couple_server", comm,
    redev::Partition{ts::setupServerPartition(mesh, cpn_file)}, mesh);
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  ReverseClassificationVertex rc;
  if (mesh.has_tag(0, "simNumbering")) {
    rc = ConstructRCFromOmegaHMesh(mesh, "simNumbering");
  } else {
    rc = ConstructRCFromOmegaHMesh<GO>(mesh, "global", wdmcpl::IndexBase::Zero);
  }

  auto is_overlap =
    ts::markServerOverlapRegion(mesh, partition, ts::IsModelEntInOverlap{});
  std::vector<GO> data(mesh.nverts());
  auto* application = cpl.AddApplication("proxy_couple");

  auto field_adapter = wdmcpl::XGCFieldAdapter<GO>(
    "xgc_gids", comm, make_array_view(data), rc, ts::IsModelEntInOverlap{});
  application->AddField("xgc_gids", std::move(field_adapter),
                        FieldTransferMethod::Copy, // to Omega_h
                        FieldEvaluationMethod::None,
                        FieldTransferMethod::Copy, // from Omega_h
                        FieldEvaluationMethod::None, is_overlap);
  do {
    application->ReceivePhase([&]() { application->ReceiveField("xgc_gids"); });
    application->SendPhase([&]() { application->SendField("xgc_gids"); });
    application->ReceivePhase([&]() { application->ReceiveField("xgc_gids"); });
    application->SendPhase([&]() { application->SendField("xgc_gids"); });
  } while (!done);
  Omega_h::vtk::write_parallel("proxy_couple", &mesh, mesh.dim());
}
void omegah_coupler(MPI_Comm comm, Omega_h::Mesh& mesh,
                    std::string_view cpn_file)
{
  // coupling server using same mesh as application
  // note the xgc_coupler stores a reference to the internal mesh and it is the
  // user responsibility to keep it alive!
  wdmcpl::CouplerServer cpl(
    "proxy_couple_server", comm,
    redev::Partition{ts::setupServerPartition(mesh, cpn_file)}, mesh);
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  auto* application = cpl.AddApplication("proxy_couple");
  ReverseClassificationVertex rc;
  std::string numbering;
  if (mesh.has_tag(0, "simNumbering")) {
    rc = ConstructRCFromOmegaHMesh(mesh, "simNumbering");
    numbering = "simNumbering";
  } else {
    rc = ConstructRCFromOmegaHMesh<GO>(mesh, "global", wdmcpl::IndexBase::Zero);
    numbering = "global";
  }

  auto is_overlap =
    ts::markServerOverlapRegion(mesh, partition, ts::IsModelEntInOverlap{});
  std::vector<GO> data(mesh.nverts());
  constexpr int nplanes = 1;
  std::vector<wdmcpl::ConvertibleCoupledField*> fields;
  for(int i=0; i<nplanes; ++i) {
    std::stringstream ss;
    //ss << "xgc_gids_plane_"<<i;
    ss << "xgc_gids";
    auto field_adapter =
      wdmcpl::OmegaHFieldAdapter<GO>(ss.str(), mesh, is_overlap, numbering);
    fields.push_back(application->AddField(ss.str(), std::move(field_adapter),
                          FieldTransferMethod::Copy, // to Omega_h
                          FieldEvaluationMethod::None,
                          FieldTransferMethod::Copy, // from Omega_h
                          FieldEvaluationMethod::None, is_overlap));
  }
  do {
    application->ReceivePhase([&]() { std::for_each(fields.begin(), fields.end(), [](wdmcpl::ConvertibleCoupledField* f){f->Receive();});});
    application->SendPhase([&]() { std::for_each(fields.begin(), fields.end(), [](wdmcpl::ConvertibleCoupledField* f){f->Send();});});
    application->ReceivePhase([&]() { std::for_each(fields.begin(), fields.end(), [](wdmcpl::ConvertibleCoupledField* f){f->Receive();});});
    application->SendPhase([&]() { std::for_each(fields.begin(), fields.end(), [](wdmcpl::ConvertibleCoupledField* f){f->Send();});});
  } while (!done);
  Omega_h::vtk::write_parallel("proxy_couple", &mesh, mesh.dim());
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
                   "</path/to/partitionFile.cpn>"
                   "<coupler type (0 xgc, 1 omega-h)>";
    }
    exit(EXIT_FAILURE);
  }

  const auto meshFile = argv[1];
  const auto classPartitionFile = argv[2];
  int coupler_type = std::atoi(argv[3]);

  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(meshFile, lib.world(), &mesh);
  MPI_Comm mpi_comm = lib.world()->get_impl();
  if (coupler_type == 0) {
    if (size != 1) {
      if (!rank) {
        std::cerr << "XGC Adapter only works on 1 rank (not a distributed mesh "
                     "datastructure)"
                  << std::endl;
      }
      std::abort();
    }
    xgc_coupler(mpi_comm, mesh, classPartitionFile);
  } else if (coupler_type == 1) {
    omegah_coupler(mpi_comm, mesh, classPartitionFile);
  } else {
    std::cerr << "Invalid coupler type. Choose 1 for XGC, 2 for Omega-h\n";
    std::abort();
  }
  return 0;
}
