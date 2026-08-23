#include <iostream>
#include <pcms.h>
#include <pcms/utility/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <data/simple.h>

#include "test_support.h"
#include "pcms/field/function_space/xgc.h"
#include "pcms/field/layout/xgc.h"
#include "pcms/coupler/serializer/xgc.h"
#include "pcms/coupler/coupler.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/coordinate_systems/cartesian.hpp"

using pcms::ConstructRCFromOmegaHMesh;
using pcms::GO;
using pcms::make_array_view;
using pcms::ReverseClassificationVertex;

static constexpr bool done = true;
namespace ts = test_support;

void xgc_coupler(MPI_Comm comm, Omega_h::Mesh& mesh, std::string_view cpn_file)
{
  // coupling server using same mesh as application
  // note the xgc_coupler stores a reference to the internal mesh and it is the
  // user responsibility to keep it alive!
  pcms::Coupler cpl("proxy_couple_server", comm, true,
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
  (void)is_overlap;
  auto* application = cpl.AddApplication("proxy_couple");

  constexpr int nplanes = 2;
  std::array<std::vector<GO>, nplanes> data;
  std::vector<pcms::FieldHandle<GO>> fields;
  for (int i = 0; i < nplanes; ++i) {
    data[i].resize(mesh.nverts());
    std::stringstream ss;
    ss << "xgc_gids_plane_" << i;
    // FIXME: The current C/Fortran proxy API couples layout registration to
    // field registration, so each XGC plane is registered as a separate layout
    // communicator even though the layouts are geometrically identical.
    auto function_space =
      pcms::XGCFieldFactory(rc, ts::IsModelEntInOverlap{},
                            static_cast<pcms::LO>(mesh.nverts()), ss.str());
    auto field = function_space.CreateField<pcms::GO>(
      ss.str(), std::make_unique<pcms::XGCFieldData<pcms::GO>>(
                  function_space.GetXGCLayout(), pcms::ValueBasis{},
                  make_array_view(data[i])));
    std::unique_ptr<pcms::FieldSerializer<GO>> serializer =
      std::make_unique<pcms::XGCFieldSerializer<GO>>(comm);
    fields.push_back(
      application->AddField(std::move(field), std::move(serializer)));
  }

  do {
    application->ReceivePhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](const pcms::FieldHandle<GO>& f) { f.Receive(); });
    });
    application->SendPhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](const pcms::FieldHandle<GO>& f) { f.Send(); });
    });
    application->ReceivePhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](const pcms::FieldHandle<GO>& f) { f.Receive(); });
    });
    application->SendPhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](const pcms::FieldHandle<GO>& f) { f.Send(); });
    });
  } while (!done);

  // Result verification
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    std::cout << "\n=== XGC Overlap Verification (OLD approach) ==="
              << std::endl;
    std::cout << "Total vertices: " << mesh.nverts() << std::endl;

    int overlap_count = 0;
    if (mesh.has_tag(0, "isOverlap")) {
      auto is_overlap = mesh.get_array<Omega_h::I8>(0, "isOverlap");
      auto is_overlap_h = Omega_h::HostRead(is_overlap);
      for (int i = 0; i < mesh.nverts(); ++i) {
        if (is_overlap_h[i]) {
          overlap_count++;
        }
      }
      std::cout << "Overlap vertices (from mesh tag): " << overlap_count
                << std::endl;
    } else {
      std::cout << "No 'isOverlap' tag found on mesh" << std::endl;
    }

    int overlap_count_func = 0;
    auto class_dims = mesh.get_array<Omega_h::I8>(0, "class_dim");
    auto class_ids = mesh.get_array<Omega_h::ClassId>(0, "class_id");
    auto class_dims_h = Omega_h::HostRead(class_dims);
    auto class_ids_h = Omega_h::HostRead(class_ids);

    for (int i = 0; i < mesh.nverts(); ++i) {
      if (ts::IsModelEntInOverlap{}(class_dims_h[i], class_ids_h[i])) {
        overlap_count_func++;
      }
    }
    std::cout << "Overlap vertices (from function): " << overlap_count_func
              << std::endl;
    std::cout << "Non-overlap vertices: "
              << (mesh.nverts() - overlap_count_func) << std::endl;
    std::cout << "Overlap ratio: "
              << (100.0 * overlap_count_func / mesh.nverts()) << "%"
              << std::endl;

    // Show field data statistics
    std::cout << "\nField data (plane 0) first 10 values: ";
    for (int i = 0; i < std::min(10, static_cast<int>(data[0].size())); ++i) {
      std::cout << data[0][i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Layout communicator count: "
              << application->GetLayoutCommunicatorCount() << std::endl;
    std::cout << "========================================\n" << std::endl;
  }

  Omega_h::vtk::write_parallel("proxy_couple", &mesh, mesh.dim());
}
void omegah_coupler(MPI_Comm comm, Omega_h::Mesh& mesh,
                    std::string_view cpn_file)
{
  // coupling server using same mesh as application
  // note the xgc_coupler stores a reference to the internal mesh and it is the
  // user responsibility to keep it alive!
  pcms::Coupler cpl("proxy_couple_server", comm, true,
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
  std::vector<pcms::FieldHandle<GO>> fields;
  for (int i = 0; i < nplanes; ++i) {
    std::stringstream ss;
    ss << "xgc_gids_plane_" << i;
    // FIXME: The current C/Fortran proxy API couples layout registration to
    // field registration, so each XGC plane is registered as a separate layout
    // communicator even though the layouts are geometrically identical.
    auto factory = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, 1, 1, pcms::csys::Cartesian::Deferred(), numbering,
      pcms::LagrangeFunctionSpace::Backend::OmegaH, ss.str());
    auto field = factory->CreateFunction<GO>(
      ss.str(), std::make_unique<pcms::SimpleFieldData<GO>>(
                  factory->GetLayout(), pcms::ValueBasis{}));
    std::unique_ptr<pcms::FieldSerializer<GO>> serializer =
      std::make_unique<pcms::FieldSerializer<GO>>();
    fields.push_back(
      application->AddField(std::move(field), std::move(serializer)));
  }
  do {
    application->ReceivePhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](const pcms::FieldHandle<GO>& f) { f.Receive(); });
    });
    application->SendPhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](const pcms::FieldHandle<GO>& f) { f.Send(); });
    });
    application->ReceivePhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](const pcms::FieldHandle<GO>& f) { f.Receive(); });
    });
    application->SendPhase([&]() {
      std::for_each(fields.begin(), fields.end(),
                    [](const pcms::FieldHandle<GO>& f) { f.Send(); });
    });
  } while (!done);
  Omega_h::vtk::write_parallel("proxy_couple", &mesh, mesh.dim());
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
          std::cerr
            << "XGC Adapter only works on 1 rank (not a distributed mesh "
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
  } catch (const std::exception& e) {
    std::cerr << "Exception caught in main: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Unknown exception caught in main" << std::endl;
    return 1;
  }
}
