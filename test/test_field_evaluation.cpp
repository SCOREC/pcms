#include <catch2/catch_test_macros.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/field_metadata.h"
#include "pcms/utility/assert.h"
#include "field_test_utils.h"
#include <cmath>

using pcms::Real;

// Non-linear test function used to exercise quadratic interpolation.
KOKKOS_INLINE_FUNCTION static Real sin_f(Real x, Real y)
{
  return std::sin(20 * x * y) / 2 + 0.5;
}

// Standard interior test points shared across all evaluation tests.
static const std::vector<Real> kEvalCoords = {
  0.7681, 0.886,  0.5337, 0.5205,   0.8088, 0.1513, 0.13,
  0.43,   0.5484, 0.8263, 0.006119, 0.8642, 0.5889, 0.5622,
  0.9268, 0.1749, 0.2615, 0.1468,   0.9793, 0.9612,
};

TEST_CASE("evaluate linear 2d omega_h_field")
{
  auto lib = Omega_h::Library{};
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 100,
                                 100, 0, false);
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

#ifdef PCMS_ENABLE_MESHFIELDS
TEST_CASE("evaluate quadratic 2d meshfields_field")
{
  auto lib = Omega_h::Library{};
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 100,
                                 100, 0, false);
  auto factory = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 2, 1, pcms::CoordinateSystem::Cartesian, "global",
    pcms::LagrangeFunctionSpace::Backend::MeshFields);

  // Quadratic DOF holders span vertices and edge midpoints; set them inline.
  const auto nverts = mesh.nents(0);
  const auto nedges = mesh.nents(1);
  auto mesh_coords = mesh.coords();
  auto edge_verts = mesh.ask_verts_of(1);

  Omega_h::Write<Real> test_f(nverts + nedges);
  Omega_h::parallel_for(
    nverts, OMEGA_H_LAMBDA(int i) {
      test_f[i] = sin_f(mesh_coords[2 * static_cast<size_t>(i)],
                        mesh_coords[2 * static_cast<size_t>(i) + 1]);
    });
  Omega_h::parallel_for(
    nedges, OMEGA_H_LAMBDA(int i) {
      auto ep = Omega_h::gather_verts<2>(edge_verts, i);
      Real cx = (mesh_coords[2 * static_cast<size_t>(ep[0])] +
                 mesh_coords[2 * static_cast<size_t>(ep[1])]) /
                2;
      Real cy = (mesh_coords[2 * static_cast<size_t>(ep[0]) + 1] +
                 mesh_coords[2 * static_cast<size_t>(ep[1]) + 1]) /
                2;
      test_f[nverts + i] = sin_f(cx, cy);
    });

  Omega_h::HostWrite<Real> test_f_host(test_f);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});
  field.GetData().SetDOFHolderDataHost(
    pcms::Rank2View<const Real, pcms::HostMemorySpace>(test_f_host.data(),
                                                       test_f_host.size(), 1));

  pcms::test::CheckEvaluation(
    factory, field, kEvalCoords,
    OMEGA_H_LAMBDA(Real x, Real y) { return std::sin(20 * x * y) / 2 + 0.5; },
    1.0e-2);
}
#endif

TEST_CASE("evaluate quadratic 2d omega_h_field throws")
{
  auto lib = Omega_h::Library{};
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1, 1, 0, 100,
                                 100, 0, false);

  REQUIRE_THROWS_AS(pcms::LagrangeFunctionSpace::FromMesh(
                      mesh, 2, 1, pcms::CoordinateSystem::Cartesian, "global",
                      pcms::LagrangeFunctionSpace::Backend::OmegaH),
                    pcms::pcms_error);
}
