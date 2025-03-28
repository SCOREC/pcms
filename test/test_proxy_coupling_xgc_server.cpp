#include <iostream>
#include <pcms.h>
#include <pcms/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include "test_support.h"
#include <pcms/omega_h_field.h>
#include <pcms/xgc_field_adapter.h>

using pcms::ConstructRCFromOmegaHMesh;
using pcms::Copy;
using pcms::CouplerServer;
using pcms::GO;
using pcms::Lagrange;
using pcms::make_array_view;
using pcms::OmegaHField;
using pcms::OmegaHFieldAdapter;
using pcms::ReadReverseClassificationVertex;
using pcms::ReverseClassificationVertex;

static constexpr bool done = true;
namespace ts = test_support;

void xgc_coupler(MPI_Comm comm, Omega_h::Mesh& mesh, std::string_view cpn_file)
{
  // coupling server using same mesh as application
  // note the xgc_coupler stores a reference to the internal mesh and it is the
  // user responsibility to keep it alive!
  pcms::CouplerServer cpl(
    "proxy_couple_server", comm,
    redev::Partition{ts::setupServerPartition(mesh, cpn_file)});
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  ReverseClassificationVertex rc;
  if (mesh.has_tag(0, "simNumbering")) {
    rc = ConstructRCFromOmegaHMesh(mesh, "simNumbering");
  } else {
    rc = ConstructRCFromOmegaHMesh<GO>(mesh, "global", pcms::IndexBase::Zero);
  }

  auto is_overlap =
    ts::markServerOverlapRegion(mesh, partition, ts::IsModelEntInOverlap{});
  auto* application = cpl.AddApplication("proxy_couple");

  constexpr int nplanes = 2;
  std::array<std::vector<GO>, nplanes> data;
  std::vector<pcms::CoupledField*> fields;
  for (int i = 0; i < nplanes; ++i) {
    data[i].resize(mesh.nverts());
    std::stringstream ss;
    ss << "xgc_gids_plane_" << i;
    auto field_adapter = pcms::XGCFieldAdapter<GO>(
      ss.str(), comm, make_array_view(data[i]), rc, ts::IsModelEntInOverlap{});
    fields.push_back(
      application->AddField(ss.str(), std::move(field_adapter)));
  }

  do {
    application->ReceivePhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](pcms::CoupledField* f) { f->Receive(); });
    });
    application->SendPhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](pcms::CoupledField* f) { f->Send(); });
    });
    application->ReceivePhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](pcms::CoupledField* f) { f->Receive(); });
    });
    application->SendPhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](pcms::CoupledField* f) { f->Send(); });
    });
  } while (!done);

  Omega_h::vtk::write_parallel("proxy_couple", &mesh, mesh.dim());
}
void omegah_coupler(MPI_Comm comm, Omega_h::Mesh& mesh,
                    std::string_view cpn_file)
{
  // coupling server using same mesh as application
  // note the xgc_coupler stores a reference to the internal mesh and it is the
  // user responsibility to keep it alive!
  pcms::CouplerServer cpl(
    "proxy_couple_server", comm,
    redev::Partition{ts::setupServerPartition(mesh, cpn_file)});
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());
  auto* application = cpl.AddApplication("proxy_couple");
  std::string numbering;
  if (mesh.has_tag(0, "simNumbering")) {
    numbering = "simNumbering";
  } else {
    Omega_h::Write<GO> gids(mesh.nverts());
    auto globals = mesh.globals(0);
    Omega_h::parallel_for(
      mesh.nverts(), OMEGA_H_LAMBDA(int i) { gids[i] = globals[i] + 1; });
    mesh.add_tag<GO>(0, "simNumbering", 1, Omega_h::Read(gids));
    numbering = "simNumbering";
  }

  auto is_overlap =
    ts::markServerOverlapRegion(mesh, partition, ts::IsModelEntInOverlap{});
  constexpr int nplanes = 2;
  std::vector<pcms::CoupledField*> fields;
  for (int i = 0; i < nplanes; ++i) {
    std::stringstream ss;
    ss << "xgc_gids_plane_" << i;
    auto field_adapter =
      pcms::OmegaHFieldAdapter<GO>(ss.str(), mesh, is_overlap, numbering);
    fields.push_back(
      application->AddField(ss.str(), std::move(field_adapter)));
  }
  do {
    application->ReceivePhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](pcms::CoupledField* f) { f->Receive(); });
    });
    application->SendPhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](pcms::CoupledField* f) { f->Send(); });
    });
    application->ReceivePhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](pcms::CoupledField* f) { f->Receive(); });
    });
    application->SendPhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](pcms::CoupledField* f) { f->Send(); });
    });
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
  int coupler_type = std::stoi(argv[3]);

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
