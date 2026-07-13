// Two-application coupling test for the MFEM field adapter.
//
// A client application owns an MFEM order-1 vertex scalar field and sends it,
// restricted to a masked overlap domain, to a rendezvous "coupler" server.
// The overlap domain is selected from MFEM element attributes following the
// create_mask strategy in the mfem-pcms-example: a vertex participates if it is
// incident to an element with the target attribute. Routing uses a single-rank
// RCB (coordinate) partition.
//
// Usage: test_mfem_coupling <role: -1 server | 0 client>

#include <pcms/coupler/coupler.hpp>
#include <pcms/coupler/overlap_mask.h>
#include <pcms/field/function_space/mfem.h>
#include <pcms/field/layout/mfem.h>

#include <mfem.hpp>
#include <redev.h>

#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

using pcms::Real;

namespace
{

constexpr int TargetAttribute = 2;

// Seeded into receiver DOFs outside the overlap before receiving. Distinct from
// any sent value (gid + 1 >= 1), so a correct masked receive must leave it
// untouched.
constexpr Real OutsideOverlapSentinel = -7.0;

// Build a 4x4 quad mesh and tag the left half (centroid x < 0.5) with the
// target attribute; the rest keep attribute 1. The same construction runs on
// both apps so vertex global ids and coordinates match.
mfem::Mesh MakeAttributedMesh()
{
  auto mesh = mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::QUADRILATERAL);
  for (int e = 0; e < mesh.GetNE(); ++e) {
    mfem::Vector center;
    mesh.GetElementCenter(e, center);
    mesh.SetAttribute(e, center[0] < 0.5 ? TargetAttribute : 1);
  }
  mesh.SetAttributes();
  return mesh;
}

redev::Partition MakeRCBPartition(int dim)
{
  redev::LOs ranks(1);
  std::iota(ranks.begin(), ranks.end(), 0);
  redev::Reals cuts = {0};
  return redev::Partition{redev::RCBPtn{dim, ranks, cuts}};
}

// Expected sent value at a vertex: its global id + 1 (strictly positive so it
// is distinguishable from the receiver's initial zero state).
Real ExpectedValue(pcms::GO gid)
{
  return static_cast<Real>(gid + 1);
}

int RunClient(MPI_Comm comm)
{
  auto serial = MakeAttributedMesh();
  mfem::ParMesh pmesh(comm, serial);
  mfem::H1_FECollection fec(1, pmesh.Dimension());
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  mfem::ParGridFunction gf(&pfes);
  gf = 0.0;

  pcms::Coupler cpl("mfem_overlap_coupler", comm, false, redev::Partition{});
  auto* app = cpl.AddApplication("mfem_app");

  auto fs = pcms::MFEMFieldFactory(pmesh, pfes, gf,
                                       pcms::CoordinateSystem::Cartesian);
  auto layout = fs.GetLayout();

  auto overlap_view =
    pcms::MFEMLayout::OverlapMaskFromAttribute(pmesh, TargetAttribute);
  app->SetLayoutOverlapMask(
    "field", std::make_unique<pcms::OverlapMask>(
               static_cast<size_t>(layout->GetNumOwnedDofHolder()),
               overlap_view));
  app->AddLayout("field", layout);

  auto handle = app->AddField("field", fs.CreateField<Real>());

  // Seed the field so each vertex holds (gid + 1).
  auto gids = layout->GetGidsHost();
  const auto n = static_cast<size_t>(layout->GetNumOwnedDofHolder());
  Kokkos::View<Real*, pcms::HostMemorySpace> values("client_values", n);
  for (size_t v = 0; v < n; ++v) {
    values(v) = ExpectedValue(gids[v]);
  }
  handle.GetField().SetDOFHolderDataHost(pcms::make_const_array_view(values));

  app->SendPhase([&]() { handle.Send(); });
  return 0;
}

int RunServer(MPI_Comm comm)
{
  auto serial = MakeAttributedMesh();
  mfem::ParMesh pmesh(comm, serial);
  mfem::H1_FECollection fec(1, pmesh.Dimension());
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  mfem::ParGridFunction gf(&pfes);
  gf = 0.0;

  pcms::Coupler cpl("mfem_overlap_coupler", comm, true,
                    MakeRCBPartition(pmesh.SpaceDimension()));
  auto* app = cpl.AddApplication("mfem_app");

  auto fs = pcms::MFEMFieldFactory(pmesh, pfes, gf,
                                       pcms::CoordinateSystem::Cartesian);
  auto layout = fs.GetLayout();
  app->AddLayout("field", layout);
  auto handle = app->AddField("field", fs.CreateField<Real>());

  // Seed the whole receiver field with a sentinel. The masked send only carries
  // the overlap DOFs, so a correct receive must overwrite only those and leave
  // every DOF outside the overlap holding the sentinel.
  {
    const auto n = static_cast<size_t>(layout->GetNumOwnedDofHolder());
    Kokkos::View<Real*, pcms::HostMemorySpace> seed("server_seed", n);
    for (size_t v = 0; v < n; ++v) {
      seed(v) = OutsideOverlapSentinel;
    }
    handle.GetField().SetDOFHolderDataHost(pcms::make_const_array_view(seed));
  }

  app->ReceivePhase([&]() { handle.Receive(); });

  // Verify: every overlap vertex received the expected value. The overlap set
  // is recomputed from the identical mesh's attributes.
  auto overlap =
    pcms::MFEMLayout::OverlapMaskFromAttribute(pmesh, TargetAttribute);
  auto gids = layout->GetGidsHost();
  auto received = handle.GetField().GetDOFHolderDataHost();

  int overlap_count = 0;
  int mismatches = 0;
  int corrupted_outside = 0;
  for (size_t v = 0; v < received.size(); ++v) {
    if (overlap(v)) {
      ++overlap_count;
      const Real expected = ExpectedValue(gids[v]);
      if (received[v] != expected) {
        ++mismatches;
        std::cerr << "Mismatch at vertex " << v << ": expected " << expected
                  << " got " << received[v] << "\n";
      }
    } else if (received[v] != OutsideOverlapSentinel) {
      // A DOF outside the overlap was overwritten on receive.
      ++corrupted_outside;
    }
  }

  const int nv = pmesh.GetNV();
  std::cout << "MFEM overlap coupling: " << overlap_count << " / " << nv
            << " overlap vertices received; " << corrupted_outside
            << " DOFs outside the overlap corrupted\n";

  if (overlap_count == 0 || overlap_count == nv) {
    std::cerr << "Overlap mask is trivial (count=" << overlap_count
              << ", nv=" << nv << "); test is not meaningful\n";
    return 1;
  }
  if (mismatches != 0) {
    std::cerr << "MFEM overlap coupling FAILED with " << mismatches
              << " mismatches\n";
    return 1;
  }
  if (corrupted_outside != 0) {
    std::cerr << "MFEM overlap coupling FAILED: " << corrupted_outside
              << " DOFs outside the overlap were overwritten on receive "
                 "(expected to be preserved)\n";
    return 1;
  }
  std::cout << "MFEM overlap coupling PASSED\n";
  return 0;
}

} // namespace

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int rc = 0;
  {
    Kokkos::ScopeGuard kokkos(argc, argv);
    if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <role: -1 server | 0 client>\n";
      MPI_Finalize();
      return EXIT_FAILURE;
    }
    const int role = std::atoi(argv[1]);
    try {
      rc = (role == -1) ? RunServer(MPI_COMM_WORLD) : RunClient(MPI_COMM_WORLD);
    } catch (const std::exception& e) {
      std::cerr << "Exception: " << e.what() << "\n";
      rc = 1;
    }
  }
  MPI_Finalize();
  return rc;
}
