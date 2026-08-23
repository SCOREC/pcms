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
#include "pcms/coupler/overlap_mask.h"
#include "pcms/field/value_view.hpp"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/coordinate_systems/cartesian.hpp"

using pcms::ConstructRCFromOmegaHMesh;
using pcms::GO;
using pcms::make_array_view;
using pcms::ReverseClassificationVertex;

static constexpr bool done = true;
namespace ts = test_support;

void xgc_coupler_with_overlap(MPI_Comm comm, Omega_h::Mesh& mesh,
                              std::string_view cpn_file)
{
  pcms::Coupler cpl("proxy_couple_server", comm, true,
                    redev::Partition{ts::setupServerPartition(mesh, cpn_file)});
  const auto partition = std::get<redev::ClassPtn>(cpl.GetPartition());

  ReverseClassificationVertex rc;
  if (mesh.has_tag(0, "simNumbering")) {
    rc = ConstructRCFromOmegaHMesh(mesh, "simNumbering");
  } else {
    rc = ConstructRCFromOmegaHMesh<GO>(mesh, "global", pcms::IndexBase::Zero);
  }

  auto* application = cpl.AddApplication("proxy_couple");

  constexpr int nplanes = 2;
  std::array<std::vector<GO>, nplanes> data;
  std::vector<pcms::FieldHandle<GO>> fields;

  for (int i = 0; i < nplanes; ++i) {
    data[i].resize(mesh.nverts());
    std::stringstream ss;
    ss << "xgc_gids_plane_" << i;

    auto overlap_mask = std::make_unique<pcms::OverlapMask>(
      mesh.nverts(), [](int dim, int id) -> int8_t {
        return ts::IsModelEntInOverlap{}(dim, id);
      });

    application->SetLayoutOverlapMask(ss.str(), std::move(overlap_mask));

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

  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    std::cout << "\n=== Field Communication Verification (NEW approach) ==="
              << std::endl;
    std::cout << "Total DOFs in field: " << data[0].size() << std::endl;

    int received_count = 0;
    int zero_count = 0;
    for (size_t i = 0; i < data[0].size(); ++i) {
      if (data[0][i] != 0) {
        received_count++;
      } else {
        zero_count++;
      }
    }

    std::cout << "Field values received (non-zero): " << received_count
              << std::endl;
    std::cout << "Field values not received (zero): " << zero_count
              << std::endl;

    std::cout << "\nSample field data (plane 0, first 20 values):" << std::endl;
    for (int i = 0; i < std::min(20, static_cast<int>(data[0].size())); ++i) {
      std::cout << "  DOF[" << i << "] = " << data[0][i] << std::endl;
    }

    int overlap_count = 0;
    auto class_dims = mesh.get_array<Omega_h::I8>(0, "class_dim");
    auto class_ids = mesh.get_array<Omega_h::ClassId>(0, "class_id");
    auto class_dims_h = Omega_h::HostRead(class_dims);
    auto class_ids_h = Omega_h::HostRead(class_ids);

    for (int i = 0; i < mesh.nverts(); ++i) {
      if (ts::IsModelEntInOverlap{}(class_dims_h[i], class_ids_h[i])) {
        overlap_count++;
      }
    }

    std::cout << "\nExpected overlap DOFs: " << overlap_count << std::endl;
    std::cout << "Layout communicator count: "
              << application->GetLayoutCommunicatorCount() << std::endl;
    std::cout << "======================================================\n"
              << std::endl;
  }

  Omega_h::vtk::write_parallel("proxy_couple_overlap", &mesh, mesh.dim());
}

void omegah_coupler_with_overlap(MPI_Comm comm, Omega_h::Mesh& mesh,
                                 std::string_view cpn_file)
{
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

  constexpr int nplanes = 2;
  std::vector<pcms::FieldHandle<GO>> fields;

  for (int i = 0; i < nplanes; ++i) {
    std::stringstream ss;
    ss << "xgc_gids_plane_" << i;

    auto overlap_mask = std::make_unique<pcms::OverlapMask>(
      mesh.nverts(), [](int dim, int id) -> int8_t {
        return ts::IsModelEntInOverlap{}(dim, id);
      });

    application->SetLayoutOverlapMask(ss.str(), std::move(overlap_mask));

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

  Omega_h::vtk::write_parallel("proxy_couple_overlap", &mesh, mesh.dim());
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
                  << " </path/to/omega_h/mesh> "
                     "</path/to/partitionFile.cpn> "
                     "<coupler type (0 xgc, 1 omega-h)>\n";
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
               "datastructure)\n";
        }
        std::abort();
      }
      xgc_coupler_with_overlap(mpi_comm, mesh, classPartitionFile);
    } else if (coupler_type == 1) {
      omegah_coupler_with_overlap(mpi_comm, mesh, classPartitionFile);
    } else {
      std::cerr << "Invalid coupler type. Choose 0 for XGC, 1 for Omega-h\n";
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
