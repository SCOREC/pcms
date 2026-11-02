// Two-application coupling test for the MFEM field adapter.

#include <pcms/coupler/coupler.hpp>
#include <pcms/coupler/overlap_mask.h>
#include <pcms/field/function_space/mfem.hpp>
#include <pcms/field/layout/mfem.hpp>

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
constexpr const char* LayoutName = "mfem_vertex";
constexpr const char* FieldName = "field";

constexpr Real OutsideOverlapSentinel = -7.0;

// Build a 4x4 quad mesh and tag the left half (centroid x < 0.5) with the
// target attribute, everything else has tag 1.
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

void FillField(pcms::Field<Real>& field, const std::vector<Real>& values)
{
  Kokkos::View<Real**, pcms::HostMemorySpace> view("values", values.size(), 1);
  for (size_t v = 0; v < values.size(); ++v) {
    view(v, 0) = values[v];
  }
  field.SetDOFHolderDataHost(pcms::MakeConstRank2View(view));
}

int RunClient(MPI_Comm comm)
{
  auto serial = MakeAttributedMesh();
  mfem::ParMesh pmesh(comm, serial);
  mfem::H1_FECollection fec(1, pmesh.Dimension());
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  mfem::ParGridFunction gf(&pfes);
  gf = 0.0;

  pcms::CouplerComm cpl("mfem_overlap_coupler", comm, false,
                        redev::Partition{});
  auto* app = cpl.AddApplication("mfem_app");

  pcms::MFEMFieldFactory factory(pmesh, pfes, gf,
                                 pcms::CoordinateSystem::Cartesian, LayoutName);
  auto layout = factory.GetLayout();

  auto overlap_view =
    pcms::MFEMLayout::OverlapMaskFromAttribute(pmesh, TargetAttribute);
  app->SetLayoutOverlapMask(
    LayoutName,
    std::make_unique<pcms::OverlapMask>(
      static_cast<size_t>(layout->GetNumOwnedDofHolder()), overlap_view));

  auto handle = app->AddField(factory.CreateField<Real>(FieldName));

  // Seed the field so each vertex holds (gid + 1).
  auto gids = layout->GetGidsHost();
  std::vector<Real> values(static_cast<size_t>(layout->GetNumOwnedDofHolder()));
  for (size_t v = 0; v < values.size(); ++v) {
    values[v] = ExpectedValue(gids[v]);
  }
  FillField(handle.GetField(), values);

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

  pcms::CouplerComm cpl("mfem_overlap_coupler", comm, true,
                        MakeRCBPartition(pmesh.SpaceDimension()));
  auto* app = cpl.AddApplication("mfem_app");

  pcms::MFEMFieldFactory factory(pmesh, pfes, gf,
                                 pcms::CoordinateSystem::Cartesian, LayoutName);
  auto layout = factory.GetLayout();
  auto handle = app->AddField(factory.CreateField<Real>(FieldName));

  FillField(
    handle.GetField(),
    std::vector<Real>(static_cast<size_t>(layout->GetNumOwnedDofHolder()),
                      OutsideOverlapSentinel));

  app->ReceivePhase([&]() { handle.Receive(); });

  auto overlap =
    pcms::MFEMLayout::OverlapMaskFromAttribute(pmesh, TargetAttribute);
  auto gids = layout->GetGidsHost();
  auto received = handle.GetField().GetDOFHolderDataHost();

  int overlap_count = 0;
  int mismatches = 0;
  int corrupted_outside = 0;
  for (size_t v = 0; v < received.extent(0); ++v) {
    if (overlap(v)) {
      ++overlap_count;
      const Real expected = ExpectedValue(gids[v]);
      if (received(v, 0) != expected) {
        ++mismatches;
        std::cerr << "Mismatch at vertex " << v << ": expected " << expected
                  << " got " << received(v, 0) << "\n";
      }
    } else if (received(v, 0) != OutsideOverlapSentinel) {
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
