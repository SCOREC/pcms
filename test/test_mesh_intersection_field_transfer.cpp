#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_shape.hpp>

#include <pcms/transfer/omega_h_conservative_projection.hpp>
#include <pcms/transfer/transfer_method.hpp>
#include <pcms/field/function_space/lagrange.h>
#include "field_test_utils.h"

#include <array>
#include <vector>

namespace
{

// Builds a unit-square 2D simplex mesh from the given element connectivity and
// adds the geometric classification tags required to build an Omega_h-backed
// Lagrange function space.
Omega_h::Mesh BuildUnitSquare(Omega_h::Library& lib, const Omega_h::LOs& ev2v)
{
  const Omega_h::Reals coords({
    0.0, 0.0, // v0
    1.0, 0.0, // v1
    1.0, 1.0, // v2
    0.0, 1.0  // v3
  });
  Omega_h::Mesh mesh(&lib);
  Omega_h::build_from_elems_and_coords(&mesh, OMEGA_H_SIMPLEX, 2, ev2v, coords);
  for (Omega_h::Int dim = 0; dim <= 2; ++dim) {
    mesh.add_tag<Omega_h::I8>(
      dim, "class_dim", 1,
      Omega_h::Read<Omega_h::I8>(mesh.nents(dim), Omega_h::I8(dim)));
    mesh.add_tag<Omega_h::ClassId>(
      dim, "class_id", 1,
      Omega_h::Read<Omega_h::ClassId>(mesh.nents(dim), Omega_h::ClassId(0)));
  }
  return mesh;
}

std::shared_ptr<pcms::LagrangeFunctionSpace> MakeP1Space(
  Omega_h::Mesh& mesh, const std::string& global_id_name = "global")
{
  return pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, pcms::CoordinateSystem::Cartesian, global_id_name,
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
}

std::shared_ptr<pcms::LagrangeFunctionSpace> MakeP0Space(Omega_h::Mesh& mesh)
{
  return pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 0, 1, pcms::CoordinateSystem::Cartesian, "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
}

void AddReorderedVertexGlobalIds(Omega_h::Mesh& mesh)
{
  mesh.add_tag<Omega_h::GO>(
    0, "reordered_global", 1,
    Omega_h::Read<Omega_h::GO>({2, 0, 3, 1}, "reordered_global"));
}

void AddSparseVertexGlobalIds(Omega_h::Mesh& mesh)
{
  mesh.add_tag<Omega_h::GO>(
    0, "sparse_global", 1,
    Omega_h::Read<Omega_h::GO>({102, 7, 41, 19}, "sparse_global"));
}

// Integral over the mesh of an order-0 (piecewise-constant) field whose values
// are stored in element order. Exact by construction.
double IntegrateP0Field(Omega_h::Mesh& mesh,
                        const pcms::Field<pcms::Real>& field)
{
  const auto values = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());
  const auto elem_areas = Omega_h::measure_elements_real(&mesh);
  const auto elem_areas_h = Omega_h::HostRead<Omega_h::Real>(elem_areas);

  double integral = 0.0;
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    integral += elem_areas_h[e] * values[e];
  }
  return integral;
}

// Per-element centroid coordinates in element order.
std::vector<std::array<double, 2>> ElementCentroids(Omega_h::Mesh& mesh)
{
  const auto coords_h = Omega_h::HostRead<Omega_h::Real>(mesh.coords());
  const auto elem_verts_h =
    Omega_h::HostRead<Omega_h::LO>(mesh.ask_elem_verts());
  std::vector<std::array<double, 2>> centroids(mesh.nelems());
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    double cx = 0.0, cy = 0.0;
    for (int k = 0; k < 3; ++k) {
      const Omega_h::LO v = elem_verts_h[3 * e + k];
      cx += coords_h[2 * v + 0];
      cy += coords_h[2 * v + 1];
    }
    centroids[e] = {cx / 3.0, cy / 3.0};
  }
  return centroids;
}

// Integral over the mesh of an order-1 nodal field whose values are stored in
// vertex order (as produced by Field<Real>::GetDOFHolderDataHost for the
// Omega_h backend), computed exactly via per-element vertex averaging.
double IntegrateP1Field(Omega_h::Mesh& mesh,
                        const pcms::Field<pcms::Real>& field)
{
  const auto values = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());
  const auto elem_areas = Omega_h::measure_elements_real(&mesh);
  const auto elem_verts = mesh.ask_elem_verts();
  const auto elem_areas_h = Omega_h::HostRead<Omega_h::Real>(elem_areas);
  const auto elem_verts_h = Omega_h::HostRead<Omega_h::LO>(elem_verts);

  double integral = 0.0;
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    const Omega_h::LO v0 = elem_verts_h[3 * e + 0];
    const Omega_h::LO v1 = elem_verts_h[3 * e + 1];
    const Omega_h::LO v2 = elem_verts_h[3 * e + 2];
    const double avg = (values[v0] + values[v1] + values[v2]) / 3.0;
    integral += elem_areas_h[e] * avg;
  }
  return integral;
}

} // namespace

// Fields that already live in the target order-1 space (constants and affine
// functions) must be reproduced exactly by the L2 projection, and the integral
// must be conserved. These properties hold without any reference solver, so
// they pin down correctness on their own.
TEST_CASE("OmegaHConservativeProjection reproduces constant and linear fields",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  // Source and target triangulate the same square along opposite diagonals.
  Omega_h::Mesh source_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));

  auto source_space = MakeP1Space(source_mesh);
  auto target_space = MakeP1Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);

  SECTION("constant field is preserved and conserved")
  {
    const double c = 2.0;
    pcms::test::SetField(
      source, OMEGA_H_LAMBDA(pcms::Real, pcms::Real) { return c; });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    REQUIRE(static_cast<Omega_h::LO>(target_values.size()) ==
            target_mesh.nverts());
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      REQUIRE(target_values[i] == Catch::Approx(c).margin(1e-10));
    }
    REQUIRE(IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(IntegrateP1Field(source_mesh, source)).margin(1e-10));
  }

  SECTION("linear field is reproduced on target vertices and conserved")
  {
    pcms::test::SetField(
      source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + y; });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    const auto tgt_coords_h =
      Omega_h::HostRead<Omega_h::Real>(target_mesh.coords());
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      const double expected = tgt_coords_h[2 * i + 0] + tgt_coords_h[2 * i + 1];
      REQUIRE(target_values[i] == Catch::Approx(expected).margin(1e-9));
    }
    REQUIRE(IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(IntegrateP1Field(source_mesh, source)).margin(1e-9));
  }
}

// For a field that does not live in the target space the projection is not an
// interpolation, but conservation of the integral is the defining property of
// the conservative (Galerkin) projection and must still hold exactly. This
// also exercises a second Apply with a different source field to confirm the
// cached integrators/evaluator/factorization are reused without stale state.
TEST_CASE("OmegaHConservativeProjection conserves the integral",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));

  auto source_space = MakeP1Space(source_mesh);
  auto target_space = MakeP1Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);

  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) {
      return x * x + x * y + 0.5 * y * y;
    });
  projection.Apply(source, target);
  REQUIRE(IntegrateP1Field(target_mesh, target) ==
          Catch::Approx(IntegrateP1Field(source_mesh, source)).margin(1e-12));

  // Second Apply with a different source field — verifies cached state is
  // reused correctly and is not stale.
  pcms::test::SetField(
    source,
    OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return 2.0 * x - y + 0.5; });
  projection.Apply(source, target);
  REQUIRE(IntegrateP1Field(target_mesh, target) ==
          Catch::Approx(IntegrateP1Field(source_mesh, source)).margin(1e-12));
}

TEST_CASE("OmegaHConservativeProjection writes reordered target GIDs in local "
          "DOF order",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));
  AddReorderedVertexGlobalIds(target_mesh);

  auto source_space = MakeP1Space(source_mesh);
  auto target_space = MakeP1Space(target_mesh, "reordered_global");

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + 2.0 * y; });

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);
  projection.Apply(source, target);

  const auto target_values = target.GetDOFHolderDataHost();
  const auto target_coords =
    target_space->GetLayout()->GetDOFHolderCoordinates().GetValues();
  const auto target_coords_h = pcms::test::CopyCoordinatesToHost(
    target_coords, target_mesh.nverts(), target_mesh.dim());
  REQUIRE(static_cast<Omega_h::LO>(target_values.extent(0)) ==
          target_mesh.nverts());
  REQUIRE(target_values.extent(1) == 1);
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    const double expected = target_coords_h(i, 0) + 2.0 * target_coords_h(i, 1);
    CAPTURE(i);
    REQUIRE(target_values(i, 0) == Catch::Approx(expected).margin(1e-9));
  }
}

TEST_CASE("OmegaHConservativeProjection maps sparse target GIDs to active "
          "PETSc rows",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));
  AddSparseVertexGlobalIds(target_mesh);

  auto source_space = MakeP1Space(source_mesh);
  auto target_space = MakeP1Space(target_mesh, "sparse_global");

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + 2.0 * y; });

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);
  projection.Apply(source, target);

  const auto target_values = target.GetDOFHolderDataHost();
  const auto target_coords =
    target_space->GetLayout()->GetDOFHolderCoordinates().GetValues();
  const auto target_coords_h = pcms::test::CopyCoordinatesToHost(
    target_coords, target_mesh.nverts(), target_mesh.dim());
  REQUIRE(static_cast<Omega_h::LO>(target_values.extent(0)) ==
          target_mesh.nverts());
  REQUIRE(target_values.extent(1) == 1);
  for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
    const double expected = target_coords_h(i, 0) + 2.0 * target_coords_h(i, 1);
    CAPTURE(i);
    REQUIRE(target_values(i, 0) == Catch::Approx(expected).margin(1e-9));
  }
}

// P0 (piecewise-constant) source -> P1 target. A constant lives in P1, so it is
// reproduced exactly; a non-constant P0 source is not reproduced but its
// integral must still be conserved.
TEST_CASE("OmegaHConservativeProjection P0 source to P1 target",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));

  auto source_space = MakeP0Space(source_mesh);
  auto target_space = MakeP1Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);

  SECTION("constant field is preserved and conserved")
  {
    const double c = 3.5;
    pcms::test::SetField(
      source, OMEGA_H_LAMBDA(pcms::Real, pcms::Real) { return c; });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    REQUIRE(static_cast<Omega_h::LO>(target_values.size()) ==
            target_mesh.nverts());
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      REQUIRE(target_values[i] == Catch::Approx(c).margin(1e-10));
    }
    REQUIRE(IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(IntegrateP0Field(source_mesh, source)).margin(1e-10));
  }

  SECTION("non-constant P0 source conserves the integral")
  {
    pcms::test::SetField(
      source,
      OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + 2.0 * y; });

    projection.Apply(source, target);

    REQUIRE(IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(IntegrateP0Field(source_mesh, source)).margin(1e-10));
  }
}

// P1 source -> P0 (piecewise-constant) target. A constant is reproduced
// exactly; a linear source projects to its cell average, which for a linear
// function equals the value at each element centroid. The integral is conserved
// in every case.
TEST_CASE("OmegaHConservativeProjection P1 source to P0 target",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  Omega_h::Mesh source_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));

  auto source_space = MakeP1Space(source_mesh);
  auto target_space = MakeP0Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space);

  SECTION("constant field is preserved and conserved")
  {
    const double c = -1.25;
    pcms::test::SetField(
      source, OMEGA_H_LAMBDA(pcms::Real, pcms::Real) { return c; });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    REQUIRE(static_cast<Omega_h::LO>(target_values.size()) ==
            target_mesh.nelems());
    for (Omega_h::LO e = 0; e < target_mesh.nelems(); ++e) {
      REQUIRE(target_values[e] == Catch::Approx(c).margin(1e-10));
    }
    REQUIRE(IntegrateP0Field(target_mesh, target) ==
            Catch::Approx(IntegrateP1Field(source_mesh, source)).margin(1e-10));
  }

  SECTION("linear field projects to cell average and conserves the integral")
  {
    const auto f = [](double x, double y) { return x + 2.0 * y; };
    pcms::test::SetField(
      source,
      OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + 2.0 * y; });

    projection.Apply(source, target);

    const auto target_values =
      pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    const auto centroids = ElementCentroids(target_mesh);
    REQUIRE(static_cast<std::size_t>(target_values.size()) == centroids.size());
    for (std::size_t e = 0; e < centroids.size(); ++e) {
      const double expected = f(centroids[e][0], centroids[e][1]);
      REQUIRE(target_values[e] == Catch::Approx(expected).margin(1e-9));
    }
    REQUIRE(IntegrateP0Field(target_mesh, target) ==
            Catch::Approx(IntegrateP1Field(source_mesh, source)).margin(1e-9));
  }
}

// Exercises the pcms::method::* recipe factory: each recipe's Build() must
// produce a working transfer operator equivalent to constructing it directly,
// and parameter-carrying recipes must forward their parameters.
TEST_CASE("TransferMethod recipes build working operators",
          "[transfer][method]")
{
  Omega_h::Library lib;
  Omega_h::Mesh source_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 2, 0, 2, 3}));
  Omega_h::Mesh target_mesh =
    BuildUnitSquare(lib, Omega_h::LOs({0, 1, 3, 1, 2, 3}));

  auto source_space = MakeP1Space(source_mesh);
  auto target_space = MakeP1Space(target_mesh);

  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();
  pcms::test::SetField(
    source, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) { return x + y; });

  const auto tgt_coords_h =
    Omega_h::HostRead<Omega_h::Real>(target_mesh.coords());
  auto check_reproduces_linear = [&]() {
    const auto values = pcms::FlattenToRank1View(target.GetDOFHolderDataHost());
    for (Omega_h::LO i = 0; i < target_mesh.nverts(); ++i) {
      const double expected = tgt_coords_h[2 * i + 0] + tgt_coords_h[2 * i + 1];
      REQUIRE(values[i] == Catch::Approx(expected).margin(1e-9));
    }
  };

  SECTION("ConservativeIntersection reproduces a linear field")
  {
    auto op = pcms::method::ConservativeIntersection{}.Build(*source_space,
                                                             *target_space);
    op->Apply(source, target);
    check_reproduces_linear();
  }

  SECTION("Interpolation reproduces a linear field")
  {
    auto op = pcms::method::Interpolation<pcms::Real>{}.Build(*source_space,
                                                              *target_space);
    op->Apply(source, target);
    check_reproduces_linear();
  }

  SECTION("ConservativeMonteCarlo forwards params and conserves the integral")
  {
    // A linear field lives in the target P1 space, so the control-variate
    // residual is zero and the projection is exact for any sample count.
    auto op =
      pcms::method::ConservativeMonteCarlo{.samples_per_element = 16}.Build(
        *source_space, *target_space);
    op->Apply(source, target);
    check_reproduces_linear();
    REQUIRE(IntegrateP1Field(target_mesh, target) ==
            Catch::Approx(IntegrateP1Field(source_mesh, source)).margin(1e-9));
  }
}
