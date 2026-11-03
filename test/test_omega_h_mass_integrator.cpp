#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>
#include <pcms/field/function_space/lagrange.h>
#include <pcms/field/layout/omega_h_lagrange.h>
#include <pcms/transfer/mass_matrix_integrator.hpp>
#include <pcms/transfer/omega_h_mass_integrator.hpp>
#include "field_test_utils.h"
#include <KokkosController.hpp>
#include <MeshField.hpp>
#include <MeshField_Config.hpp>
#include <MeshField_Element.hpp>
#include <petscksp.h>
#include <map>
#include <utility>
#include <vector>
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

namespace
{

// dummy coordinate system that's public interfaces gives
// same result as Cartesian should fail equality
class OtherCartesian2 final : public pcms::CoordinateSystem
{
public:
  [[nodiscard]] std::string_view Kind() const noexcept override
  {
    return "cartesian";
  }
  [[nodiscard]] int Dimension() const noexcept override { return 2; }
  [[nodiscard]] bool HasOrthogonalBasis() const noexcept override
  {
    return true;
  }
  [[nodiscard]] bool HasUnitScaleFactors() const noexcept override
  {
    return true;
  }
  [[nodiscard]] bool operator==(
    const pcms::CoordinateSystem& other) const noexcept override
  {
    return dynamic_cast<const OtherCartesian2*>(&other) != nullptr;
  }
};

std::map<std::pair<pcms::GO, pcms::GO>, pcms::Real> BuildReferenceMassMap(
  Omega_h::Mesh& mesh, const pcms::FunctionSpace& space)
{
  using namespace pcms;

  MeshField::OmegahMeshField<DefaultExecutionSpace, 2,
                             MeshField::KokkosController>
    omf(mesh);
  auto coordField = omf.getCoordField();
  const auto [shp, map] = MeshField::Omegah::getTriangleElement<1>(mesh);
#if MeshFields_VERSION < 10000
  MeshField::FieldElement coordFe(mesh.nelems(), coordField, shp, map);
#else
  MeshField::FieldElement coordFe(mesh.nelems(), coordField.field, shp, map);
#endif
  auto elm_mass_dev = buildElementMassMatrix(mesh, coordFe);
  auto elm_mass_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, elm_mass_dev);

  const auto layout =
    std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(space.GetLayout());
  REQUIRE(layout != nullptr);
  const auto faces2nodes_dev = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  const auto faces2nodes = Omega_h::HostRead(faces2nodes_dev);
  const auto& gids = layout->GetGidsHost();

  std::map<std::pair<GO, GO>, Real> ref;
  for (int e = 0; e < mesh.nelems(); ++e) {
    const Omega_h::Few<Omega_h::LO, 3> verts = {
      faces2nodes[e * 3 + 0], faces2nodes[e * 3 + 1], faces2nodes[e * 3 + 2]};
    for (int i = 0; i < 3; ++i) {
      const GO row = static_cast<GO>(gids[verts[i]]);
      for (int j = 0; j < 3; ++j) {
        const GO col = static_cast<GO>(gids[verts[j]]);
        ref[{row, col}] += elm_mass_host(e * 9 + i * 3 + j);
      }
    }
  }
  return ref;
}

} // namespace

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE(
  "OmegaHMassIntegrator: assembled matrix matches buildElementMassMatrix",
  "[mass_integrator]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);

  auto space = pcms::test::MakeP1Space(mesh);

  const auto ref = BuildReferenceMassMap(mesh, *space);

  auto integrator = pcms::BuildOmegaHMassIntegrator(*space);
  Mat mat = integrator->GetMatrix();

  for (const auto& [key, ref_val] : ref) {
    PetscInt r = static_cast<PetscInt>(key.first);
    PetscInt c = static_cast<PetscInt>(key.second);
    PetscScalar got = 0.0;
    MatGetValues(mat, 1, &r, 1, &c, &got);
    CAPTURE(key.first, key.second, ref_val, got);
    CHECK(static_cast<pcms::Real>(got) ==
          Catch::Approx(ref_val).epsilon(1e-12));
  }
  // mat is owned by integrator; no MatDestroy here.
}

TEST_CASE("OmegaHMassIntegrator: row sums match lumped mass",
          "[mass_integrator]")
{
  // For a consistent mass matrix M, the row sums equal the lumped mass at
  // each node: row_sum(i) = integral(N_i dx). For a regular mesh of unit
  // area with uniform nodal distribution the sum over all nodes equals 1.
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);

  auto space = pcms::test::MakeP1Space(mesh);

  const auto layout =
    std::dynamic_pointer_cast<const pcms::OmegaHLagrangeLayout>(
      space->GetLayout());
  const auto& gids = layout->GetGidsHost();
  const int nverts = mesh.nverts();

  auto integrator = pcms::BuildOmegaHMassIntegrator(*space);
  Mat mat = integrator->GetMatrix();

  // Sum all row sums.
  pcms::Real grand_total = 0.0;
  for (int v = 0; v < nverts; ++v) {
    const pcms::GO row = static_cast<pcms::GO>(gids[v]);
    pcms::Real row_sum = 0.0;
    for (int u = 0; u < nverts; ++u) {
      const pcms::GO col = static_cast<pcms::GO>(gids[u]);
      PetscInt r = static_cast<PetscInt>(row);
      PetscInt c = static_cast<PetscInt>(col);
      PetscScalar val = 0.0;
      MatGetValues(mat, 1, &r, 1, &c, &val);
      row_sum += static_cast<pcms::Real>(val);
    }
    CAPTURE(v, row, row_sum);
    CHECK(row_sum > 0.0); // each row sum must be positive
    grand_total += row_sum;
  }

  // Total == area of domain == 1.
  CHECK(grand_total == Catch::Approx(1.0).epsilon(1e-10));
  // mat is owned by integrator; no MatDestroy here.
}

TEST_CASE("OmegaHMassIntegrator: lumped diagonal equals consistent row sums",
          "[mass_integrator]")
{
  // Row-sum lumping: M_L(i,i) = sum_j M_C(i,j) = integral(N_i dx), all
  // off-diagonal entries zero. The trace equals the domain area.
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);

  auto space = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);

  const auto layout =
    std::dynamic_pointer_cast<const pcms::OmegaHLagrangeLayout>(
      space->GetLayout());
  const auto& gids = layout->GetGidsHost();
  const int nverts = mesh.nverts();

  auto consistent = pcms::BuildOmegaHMassIntegrator(*space);
  auto lumped =
    pcms::BuildOmegaHMassIntegrator(*space, pcms::MassMatrixType::Lumped);
  CHECK_FALSE(consistent->IsDiagonal());
  CHECK(lumped->IsDiagonal());
  Mat consistent_mat = consistent->GetMatrix();
  Mat lumped_mat = lumped->GetMatrix();

  pcms::Real trace = 0.0;
  for (int v = 0; v < nverts; ++v) {
    const PetscInt row = static_cast<PetscInt>(gids[v]);
    pcms::Real row_sum = 0.0;
    for (int u = 0; u < nverts; ++u) {
      PetscInt r = row;
      PetscInt c = static_cast<PetscInt>(gids[u]);
      PetscScalar val = 0.0;
      MatGetValues(consistent_mat, 1, &r, 1, &c, &val);
      row_sum += static_cast<pcms::Real>(val);

      PetscScalar lumped_val = 0.0;
      MatGetValues(lumped_mat, 1, &r, 1, &c, &lumped_val);
      CAPTURE(v, u, row_sum, lumped_val);
      if (r == c) {
        trace += static_cast<pcms::Real>(lumped_val);
      } else {
        CHECK(static_cast<pcms::Real>(lumped_val) ==
              Catch::Approx(0.0).margin(1e-14));
      }
    }
    PetscInt r = row;
    PetscScalar diag = 0.0;
    MatGetValues(lumped_mat, 1, &r, 1, &r, &diag);
    CAPTURE(v, row_sum, diag);
    CHECK(static_cast<pcms::Real>(diag) ==
          Catch::Approx(row_sum).epsilon(1e-12));
  }

  // Trace == area of domain == 1.
  CHECK(trace == Catch::Approx(1.0).epsilon(1e-10));
  // matrices are owned by the integrators; no MatDestroy here.
}

TEST_CASE("OmegaHMassIntegrator: P0 lumped equals consistent",
          "[mass_integrator]")
{
  // P0 basis functions have disjoint support, so the consistent mass matrix
  // is already diagonal (M_ee = area(e)) and lumping is a no-op.
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);

  auto space = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 0, 1, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);

  auto consistent = pcms::BuildOmegaHMassIntegrator(*space);
  auto lumped =
    pcms::BuildOmegaHMassIntegrator(*space, pcms::MassMatrixType::Lumped);
  CHECK(consistent->IsDiagonal());
  CHECK(lumped->IsDiagonal());
  Mat consistent_mat = consistent->GetMatrix();
  Mat lumped_mat = lumped->GetMatrix();

  const auto elem_areas = Omega_h::measure_elements_real(&mesh);
  const auto elem_areas_h = Omega_h::HostRead<Omega_h::Real>(elem_areas);
  const auto layout =
    std::dynamic_pointer_cast<const pcms::OmegaHLagrangeLayout>(
      space->GetLayout());
  const auto& gids = layout->GetGidsHost();

  for (int e = 0; e < mesh.nelems(); ++e) {
    PetscInt r = static_cast<PetscInt>(gids[e]);
    PetscScalar consistent_val = 0.0;
    PetscScalar lumped_val = 0.0;
    MatGetValues(consistent_mat, 1, &r, 1, &r, &consistent_val);
    MatGetValues(lumped_mat, 1, &r, 1, &r, &lumped_val);
    CAPTURE(e, elem_areas_h[e]);
    CHECK(static_cast<pcms::Real>(consistent_val) ==
          Catch::Approx(elem_areas_h[e]).epsilon(1e-12));
    CHECK(static_cast<pcms::Real>(lumped_val) ==
          Catch::Approx(elem_areas_h[e]).epsilon(1e-12));
  }
  // matrices are owned by the integrators; no MatDestroy here.
}

TEST_CASE("OmegaHMassIntegrator: rejects invalid layouts", "[mass_integrator]")
{
  Omega_h::Library lib;
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);

  SECTION("multi-component space throws")
  {
    auto space = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, 1, 2, pcms::csys::Cartesian::Deferred(), "global",
      pcms::LagrangeFunctionSpace::Backend::OmegaH);
    REQUIRE_THROWS(pcms::BuildOmegaHMassIntegrator(*space));
  }

  SECTION("non-Cartesian coordinate system throws")
  {
    auto space = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, 1, 1, pcms::csys::CylindricalRZ::Create(), "global",
      pcms::LagrangeFunctionSpace::Backend::OmegaH);
    REQUIRE_THROWS(pcms::BuildOmegaHMassIntegrator(*space));
  }

  SECTION("distinct cartesian-kind system throws")
  {
    const std::shared_ptr<const pcms::CoordinateSystem> other =
      std::make_shared<OtherCartesian2>();
    auto space = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, 1, 1, other, "global",
      pcms::LagrangeFunctionSpace::Backend::OmegaH);
    REQUIRE_THROWS(pcms::BuildOmegaHMassIntegrator(*space));
  }
}
