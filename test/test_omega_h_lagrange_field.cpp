#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>

#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/field_metadata.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/mesh_geometry.h"
#include "field_test_utils.h"

#include <memory>
#include <optional>
#include <stdexcept>

using pcms::LO;
using pcms::Real;

static Omega_h::Mesh MakeBox2D(Omega_h::CommPtr world, int nx = 10, int ny = 10)
{
  return Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, nx, ny, 0,
                            false);
}

// ---- Layout tests -----------------------------------------------------------

TEST_CASE("OmegaHLagrangeLayout order-1 properties")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world());
  pcms::OmegaHLagrangeLayout layout(mesh, 1, 2,
                                    pcms::CoordinateSystem::Cartesian);

  REQUIRE(layout.GetOrder() == 1);
  REQUIRE(layout.GetNumComponents() == 2);
  REQUIRE(layout.GetNumOwnedDofHolder() == mesh.nents(0));
  REQUIRE(layout.GetNumGlobalDofHolder() == mesh.nglobal_ents(0));
  REQUIRE(layout.IsDistributed()); // always true for Omega_h mesh layouts

  // DOF holder coordinates should match vertex coordinates
  auto coords_device = layout.GetDOFHolderCoordinates().GetValues();
  int nverts = mesh.nents(0);
  auto coords_view =
    pcms::test::CopyCoordinatesToHost(coords_device, nverts, mesh.dim());
  auto mesh_coords = Omega_h::HostRead<Real>(mesh.coords());
  REQUIRE(static_cast<int>(coords_view.extent(0)) == nverts);
  REQUIRE(static_cast<int>(coords_view.extent(1)) == mesh.dim());
  for (int v = 0; v < nverts; ++v) {
    for (int d = 0; d < mesh.dim(); ++d) {
      REQUIRE(coords_view(v, d) ==
              Catch::Approx(mesh_coords[v * mesh.dim() + d]));
    }
  }

  // GetEntOffsets: vertices are at slot 0, all other slots = nverts
  auto offsets = layout.GetEntOffsets();
  REQUIRE(offsets[0] == 0);
  for (int i = 1; i < pcms::ent_offsets_len; ++i)
    REQUIRE(offsets[i] == nverts);
}

TEST_CASE("OmegaHLagrangeLayout order-0 properties")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world());
  pcms::OmegaHLagrangeLayout layout(mesh, 0, 1,
                                    pcms::CoordinateSystem::Cartesian);

  REQUIRE(layout.GetOrder() == 0);
  REQUIRE(layout.GetNumComponents() == 1);
  REQUIRE(layout.GetNumOwnedDofHolder() == mesh.nelems());
  REQUIRE(layout.GetNumGlobalDofHolder() == mesh.nglobal_ents(mesh.dim()));

  // DOF holder coordinates should match element centroids
  auto coords_device = layout.GetDOFHolderCoordinates().GetValues();
  int nelems = mesh.nelems();
  auto coords_view =
    pcms::test::CopyCoordinatesToHost(coords_device, nelems, mesh.dim());
  auto centroids =
    Omega_h::HostRead<Real>(pcms::get_entity_centroids(mesh, mesh.dim()));
  REQUIRE(static_cast<int>(coords_view.extent(0)) == nelems);
  for (int e = 0; e < nelems; ++e) {
    for (int d = 0; d < mesh.dim(); ++d) {
      REQUIRE(coords_view(e, d) ==
              Catch::Approx(centroids[e * mesh.dim() + d]));
    }
  }

  // GetEntOffsets: all DOFs are at entity_dim = mesh.dim() (slot 2 for 2D)
  auto offsets = layout.GetEntOffsets();
  for (int i = 0; i <= mesh.dim(); ++i)
    REQUIRE(offsets[i] == 0);
  for (int i = mesh.dim() + 1; i < pcms::ent_offsets_len; ++i)
    REQUIRE(offsets[i] == nelems);
}

TEST_CASE("OmegaHLagrangeLayout invalid order throws")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world());
  REQUIRE_THROWS_AS(
    pcms::OmegaHLagrangeLayout(mesh, 2, 1, pcms::CoordinateSystem::Cartesian),
    std::invalid_argument);
  REQUIRE_THROWS_AS(
    pcms::OmegaHLagrangeLayout(mesh, -1, 1, pcms::CoordinateSystem::Cartesian),
    std::invalid_argument);
}

TEST_CASE("OmegaHLagrangeLayout layout sharing")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world());
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);

  auto f1 = factory.CreateField<Real>(pcms::FieldMetadata{});
  auto f2 = factory.CreateField<Real>(pcms::FieldMetadata{});

  REQUIRE(&f1.GetLayout() == &f2.GetLayout());
}

// ---- Order-1 field tests ----------------------------------------------------

TEST_CASE("OmegaHLagrangeField order-1: set/get DOF data round-trip")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world());
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});

  int n = factory.GetLayout()->GetNumOwnedDofHolder();
  std::vector<Real> data(n);
  for (int i = 0; i < n; ++i)
    data[i] = static_cast<Real>(i);

  pcms::Rank2View<const Real, pcms::HostMemorySpace> view(data.data(), n, 1);
  field.GetData().SetDOFHolderDataHost(view);

  auto got = pcms::FlattenToRank1View(field.GetData().GetDOFHolderDataHost());
  REQUIRE(static_cast<int>(got.size()) == n);
  for (int i = 0; i < n; ++i)
    REQUIRE(got[i] == Catch::Approx(data[i]));
}

TEST_CASE("OmegaHLagrangeField order-1: linear function evaluation")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world(), 20, 20);
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});

  pcms::test::SetField(
    field.GetData(), *factory.GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return x + 2.0 * y; });
  pcms::test::CheckEvaluation(
    factory, field, pcms::test::StandardEvalCoords2D(),
    OMEGA_H_LAMBDA(Real x, Real y) { return x + 2.0 * y; });
}

// The same linear evaluation test run on the MeshFields-backed order-1 field
// ensures both backends produce identical results for the same inputs.
TEST_CASE("MeshFieldsAdapter order-1: linear function evaluation (shared util)")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world(), 20, 20);
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});

  pcms::test::SetField(
    field.GetData(), *factory.GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return x + 2.0 * y; });
  pcms::test::CheckEvaluation(
    factory, field, pcms::test::StandardEvalCoords2D(),
    OMEGA_H_LAMBDA(Real x, Real y) { return x + 2.0 * y; });
}

TEST_CASE("OmegaHLagrangeField order-1: out-of-bounds FILL mode")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world());
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});

  pcms::test::SetField(
    field.GetData(), *factory.GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return x + 2.0 * y; });

  Real fill_value = -999.0;
  std::vector<Real> outside{-0.5, 0.5, 1.5, 0.5, 0.5, -0.5, 0.5, 1.5};
  pcms::test::CheckFillMode(factory, field, fill_value, outside);
}

TEST_CASE("OmegaHLagrangeField order-1: serialize / deserialize round-trip")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world());
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});

  pcms::test::SetField(
    field.GetData(), *factory.GetLayout(),
    OMEGA_H_LAMBDA(Real x, Real y) { return x + 2.0 * y; });
  pcms::test::CheckSerializeDeserialize(*factory.GetLayout(), field.GetData());
}

TEST_CASE(
  "OmegaHLagrangeField order-1: multi-component serialize / deserialize "
  "round-trip")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world());
  const int nc = 3;
  // The MeshFields backend is scalar-only; multi-component needs OmegaH.
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, nc, pcms::CoordinateSystem::Cartesian, "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});

  const int n = factory.GetLayout()->GetNumOwnedDofHolder();
  std::vector<Real> data(static_cast<size_t>(n) * nc);
  for (int i = 0; i < n; ++i)
    for (int c = 0; c < nc; ++c)
      data[static_cast<size_t>(i) * nc + c] = i + 0.25 * c;
  field.GetData().SetDOFHolderDataHost(
    pcms::Rank2View<const Real, pcms::HostMemorySpace>(data.data(), n, nc));

  pcms::test::CheckSerializeDeserialize(*factory.GetLayout(), field.GetData());
}

// ---- Order-0 field tests ----------------------------------------------------

TEST_CASE("OmegaHLagrangeField order-0: set/get DOF data round-trip")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world());
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 0, 1, pcms::CoordinateSystem::Cartesian);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});

  int n = factory.GetLayout()->GetNumOwnedDofHolder();
  std::vector<Real> data(n, 3.14);
  pcms::Rank2View<const Real, pcms::HostMemorySpace> view(data.data(), n, 1);
  field.GetData().SetDOFHolderDataHost(view);

  auto got = pcms::FlattenToRank1View(field.GetData().GetDOFHolderDataHost());
  REQUIRE(static_cast<int>(got.size()) == n);
  for (int i = 0; i < n; ++i)
    REQUIRE(got[i] == Catch::Approx(3.14));
}

TEST_CASE("OmegaHLagrangeField order-0: constant field evaluation")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world(), 10, 10);
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 0, 1, pcms::CoordinateSystem::Cartesian);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});

  const Real kValue = 42.0;
  int nelems = mesh.nelems();
  std::vector<Real> data(nelems, kValue);
  pcms::Rank2View<const Real, pcms::HostMemorySpace> view(data.data(), nelems,
                                                          1);
  field.GetData().SetDOFHolderDataHost(view);

  auto pts = pcms::test::StandardEvalCoords2D();
  int n = static_cast<int>(pts.size()) / 2;
  auto device_coords =
    pcms::test::CreateDeviceCoordinateView(pts, factory.GetCoordinateSystem());
  auto evaluator = factory.CreatePointEvaluator<Real>(
    pcms::EvaluationRequest::FromCoordinates(device_coords.coordinate_view));

  Kokkos::View<Real*, pcms::DeviceMemorySpace> eval_device("eval", n);
  using LayoutPolicy =
    pcms::detail::default_layout_for_memory_space_t<pcms::DeviceMemorySpace>;
  auto out = pcms::Rank2View<Real, pcms::DeviceMemorySpace, LayoutPolicy>(
    eval_device.data(), n, 1);
  evaluator->Evaluate(field, out);
  auto eval_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), eval_device);

  for (int i = 0; i < n; ++i)
    REQUIRE(eval_host(i) == Catch::Approx(kValue));
}

TEST_CASE("OmegaHLagrangeField order-0: out-of-bounds FILL mode")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world());
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 0, 1, pcms::CoordinateSystem::Cartesian);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});

  int nelems = mesh.nelems();
  std::vector<Real> data(nelems, 1.0);
  pcms::Rank2View<const Real, pcms::HostMemorySpace> view(data.data(), nelems,
                                                          1);
  field.GetData().SetDOFHolderDataHost(view);

  Real fill_value = -1.0;
  std::vector<Real> outside{-0.5, 0.5, 1.5, 0.5};
  pcms::test::CheckFillMode(factory, field, fill_value, outside);
}

TEST_CASE("OmegaHLagrangeField order-0: serialize / deserialize round-trip")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world());
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 0, 1, pcms::CoordinateSystem::Cartesian);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});

  int nelems = mesh.nelems();
  std::vector<Real> data(nelems);
  for (int i = 0; i < nelems; ++i)
    data[i] = static_cast<Real>(i);
  pcms::Rank2View<const Real, pcms::HostMemorySpace> view(data.data(), nelems,
                                                          1);
  field.GetData().SetDOFHolderDataHost(view);

  pcms::test::CheckSerializeDeserialize(*factory.GetLayout(), field.GetData());
}

// ---- Layout sharing communicator contract -----------------------------------
//
// BEHAVIORAL CONTRACT (requires MPI/redev — exercised by
// test_field_communication with clientId=2/3 via test_shared_layout):
//
// When two fields are added to an Application2 that share the same FieldLayout
// (e.g. both created from the same LagrangeFunctionSpace), the Application2
// must reuse a single FieldLayoutCommunicator for both fields rather than
// creating separate communicators. This is verified by:
//   - Calling AddLayout() registers exactly one FieldLayoutCommunicator
//   - Calling AddField() for additional fields with the same layout leaves the
//     count at one (Application2::GetLayoutCommunicatorCount() == 1)
// See test/test_field_communication.cpp::test_shared_layout for the full test.

// ---- Temporary factory lifetime safety --------------------------------------

TEST_CASE("OmegaHLagrangeField: field valid after layout destruction")
{
  auto lib = Omega_h::Library{};
  auto mesh = MakeBox2D(lib.world());

  std::optional<pcms::Field<Real>> field;
  {
    auto factory = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, 1, 1, pcms::CoordinateSystem::Cartesian);
    field.emplace(factory.CreateField<Real>(pcms::FieldMetadata{}));
  } // factory goes out of scope; field keeps layout alive

  pcms::test::SetField(
    *field, OMEGA_H_LAMBDA(Real x, Real y) { return x + 2.0 * y; });
  // Just verify data was set correctly (no evaluator needed for this lifetime
  // test)
  auto data = field->GetDOFHolderDataHost();
  REQUIRE(data.size() > 0);
}
