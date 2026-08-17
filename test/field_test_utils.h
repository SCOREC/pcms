#ifndef PCMS_TEST_FIELD_TEST_UTILS_H
#define PCMS_TEST_FIELD_TEST_UTILS_H

#include <catch2/catch_approx.hpp>
#include <Kokkos_Core.hpp>
#include "pcms/configuration.h"
#include "pcms/field/field.h"
#include "pcms/field/field_data.h"
#include "pcms/field/field_evaluator_factory.h"
#include "pcms/field/evaluation_request.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/coordinate_system.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#ifdef PCMS_ENABLE_OMEGA_H
#include <Omega_h_build.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>
#include "pcms/field/function_space/lagrange.h"
#include "pcms/localization/adj_search.hpp"
#endif
#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
#include "pcms/transfer/linear_form_integrator.hpp"
#endif
#include <cmath>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

// Shared utilities for field tests that apply equally to MeshFields-backed
// and native OmegaH-backed field implementations.

namespace pcms::test
{

// Affine test function — exactly representable on linear elements.
KOKKOS_INLINE_FUNCTION Real linear_f(Real x, Real y)
{
  return x + 2.0 * y;
}

// Interior test points for a unit [0,1]^2 box mesh.
inline std::vector<Real> StandardEvalCoords2D()
{
  return {0.1, 0.2, 0.5, 0.5, 0.7, 0.3, 0.9, 0.1, 0.2, 0.8};
}

// Points strictly outside a unit [0,1]^2 box mesh.
inline std::vector<Real> StandardOutsideCoords2D()
{
  return {-0.5, 0.5, 1.5, 0.5, 0.5, -0.5, 0.5, 1.5};
}

#ifdef PCMS_ENABLE_OMEGA_H
// Builds a unit-square 2D simplex mesh from the given element connectivity and
// adds the geometric classification tags required to build an Omega_h-backed
// Lagrange function space.
inline Omega_h::Mesh BuildUnitSquare(Omega_h::Library& lib,
                                     const Omega_h::LOs& ev2v)
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

// Convenience overload: diagonal=0 splits the square along vertices (1,3),
// diagonal=1 along (0,2).
inline Omega_h::Mesh BuildUnitSquare(Omega_h::Library& lib, int diagonal)
{
  return BuildUnitSquare(lib, (diagonal == 0)
                                ? Omega_h::LOs({0, 1, 3, 1, 2, 3})
                                : Omega_h::LOs({0, 1, 2, 0, 2, 3}));
}

inline std::shared_ptr<LagrangeFunctionSpace> MakeP1Space(
  Omega_h::Mesh& mesh, const std::string& global_id_name = "global")
{
  return LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, CoordinateSystem::Cartesian, global_id_name,
    LagrangeFunctionSpace::Backend::OmegaH);
}

inline std::shared_ptr<LagrangeFunctionSpace> MakeP0Space(Omega_h::Mesh& mesh)
{
  return LagrangeFunctionSpace::FromMesh(
    mesh, 0, 1, CoordinateSystem::Cartesian, "global",
    LagrangeFunctionSpace::Backend::OmegaH);
}

inline bool AreArraysEqualUnordered(
  const Omega_h::HostRead<Omega_h::LO>& array1,
  const Omega_h::HostRead<Omega_h::LO>& array2, int start, int end)
{
  std::unordered_map<Omega_h::LO, int> freq1, freq2;
  for (int i = start; i < end; ++i) {
    freq1[array1[i]]++;
    freq2[array2[i]]++;
  }
  return freq1 == freq2;
}

inline void CheckSupportResultsEquivalent(const SupportResults& actual,
                                          const SupportResults& expected)
{
  auto actual_ptr = Omega_h::HostRead<Omega_h::LO>(actual.supports_ptr);
  auto actual_idx = Omega_h::HostRead<Omega_h::LO>(actual.supports_idx);
  auto expected_ptr = Omega_h::HostRead<Omega_h::LO>(expected.supports_ptr);
  auto expected_idx = Omega_h::HostRead<Omega_h::LO>(expected.supports_idx);

  REQUIRE(actual_ptr.size() == expected_ptr.size());
  REQUIRE(actual_idx.size() == expected_idx.size());

  for (int i = 0; i < actual_ptr.size(); ++i)
    REQUIRE(actual_ptr[i] == expected_ptr[i]);

  for (int i = 0; i < actual_ptr.size() - 1; ++i) {
    CAPTURE(i);
    REQUIRE(AreArraysEqualUnordered(actual_idx, expected_idx, actual_ptr[i],
                                    actual_ptr[i + 1]));
  }
}

inline std::vector<Real> CopyOmegaHRealsToVector(const Omega_h::Reals& coords)
{
  auto coords_read = Omega_h::HostRead<Omega_h::Real>(coords);
  return std::vector<Real>(coords_read.data(),
                           coords_read.data() + coords_read.size());
}

inline double IntegrateP0Field(Omega_h::Mesh& mesh, const Field<Real>& field)
{
  const auto values = FlattenToRank1View(field.GetDOFHolderDataHost());
  const auto measures = Omega_h::measure_elements_real(&mesh);
  const auto measures_h = Omega_h::HostRead<Omega_h::Real>(measures);

  double integral = 0.0;
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    integral += measures_h[e] * values[e];
  }
  return integral;
}

inline double IntegrateP1Field(Omega_h::Mesh& mesh, const Field<Real>& field)
{
  const auto values = FlattenToRank1View(field.GetDOFHolderDataHost());
  const auto measures = Omega_h::measure_elements_real(&mesh);
  const auto measures_h = Omega_h::HostRead<Omega_h::Real>(measures);
  const auto elem_verts_h =
    Omega_h::HostRead<Omega_h::LO>(mesh.ask_elem_verts());
  const int verts_per_elem = mesh.dim() + 1;

  double integral = 0.0;
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    double avg = 0.0;
    for (int k = 0; k < verts_per_elem; ++k) {
      avg += values[elem_verts_h[verts_per_elem * e + k]];
    }
    integral += measures_h[e] * (avg / verts_per_elem);
  }
  return integral;
}
#endif // PCMS_ENABLE_OMEGA_H

// Copy coordinates from device memory to a host view.
// This handles potential layout mismatches between host and device memory
// spaces.
template <typename CoordView>
inline Kokkos::View<const Real**, HostMemorySpace> CopyCoordinatesToHost(
  const CoordView& coords_device, int nents, int dim)
{
  auto coords_view =
    Kokkos::View<Real**, HostMemorySpace>("coords_view", nents, dim);
  auto coords_view_device =
    Kokkos::create_mirror_view(DeviceMemorySpace(), coords_view);
  ConvertMismatchLayoutView2D(coords_view_device, coords_device);
  Kokkos::deep_copy(coords_view, coords_view_device);
  return coords_view;
}

template <typename ExecutionSpace, typename Func>
inline std::vector<Real> EvaluateReferenceFunction(const std::vector<Real>& pts,
                                                   Func func)
{
  using MemorySpace = typename ExecutionSpace::memory_space;

  int n = static_cast<int>(pts.size()) / 2;
  auto pts_host = Kokkos::View<const Real**, HostMemorySpace>(pts.data(), n, 2);
  auto pts_device =
    Kokkos::create_mirror_view_and_copy(DeviceMemorySpace(), pts_host);

  Kokkos::View<Real*, MemorySpace> expected_device("expected_device", n);
  Kokkos::parallel_for(
    "field_test_utils_expected", Kokkos::RangePolicy<ExecutionSpace>(0, n),
    KOKKOS_LAMBDA(int i) {
      expected_device(i) = func(pts_device(i, 0), pts_device(i, 1));
    });

  auto expected_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), expected_device);
  return std::vector<Real>(expected_host.data(),
                           expected_host.data() + expected_host.extent(0));
}

template <typename DataView, typename CoordsView, typename Func>
struct SetFieldFunctor
{
  DataView data;
  CoordsView coords;
  Func f;

  SetFieldFunctor(DataView data_, CoordsView coords_, Func f_)
    : data(data_), coords(coords_), f(f_)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const
  {
    if constexpr (std::is_invocable_v<Func, Real, Real, Real>) {
      data(i) = f(coords(i, 0), coords(i, 1), coords(i, 2));
    } else {
      data(i) = f(coords(i, 0), coords(i, 1));
    }
  }
};

// Set scalar DOF data by sampling func at each DOF-holder coordinate. The
// arity of func selects the spatial dimension: func(x, y) for 2D layouts,
// func(x, y, z) for 3D.
template <typename ExecutionSpace = DefaultExecutionSpace, typename Func>
inline void SetField(const FieldLayout& layout, FieldData<Real>& field,
                     Func func)
{
  using MemorySpace = typename ExecutionSpace::memory_space;

  static_assert(std::is_invocable_v<Func, Real, Real> ||
                  std::is_invocable_v<Func, Real, Real, Real>,
                "SetField requires func(x, y) or func(x, y, z)");

  auto dof_coords = layout.GetDOFHolderCoordinates().GetValues();
  int n = static_cast<int>(dof_coords.extent(0));

  Kokkos::View<Real*, MemorySpace> data_device("data_device", n);
  Kokkos::parallel_for("field_test_utils_set_field",
                       Kokkos::RangePolicy<ExecutionSpace>(0, n),
                       SetFieldFunctor{data_device, dof_coords, func});

  auto data_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), data_device);
  field.SetDOFHolderDataHost(
    Rank2View<const Real, HostMemorySpace>(data_host.data(), n, 1));
}

template <typename ExecutionSpace = DefaultExecutionSpace, typename Func>
inline void SetField(Field<Real>& field, Func func)
{
  SetField<ExecutionSpace>(field.GetLayout(), field.GetData(), func);
}

template <typename ExecutionSpace = DefaultExecutionSpace, typename Func>
inline void SetField(FieldData<Real>& field, const FieldLayout& layout,
                     Func func)
{
  SetField<ExecutionSpace>(layout, field, func);
}

// Helper structure to hold device coordinates with proper lifetime management
struct DeviceCoordinates
{
  Kokkos::View<Real**, DeviceMemorySpace> view;
  CoordinateView<DeviceMemorySpace> coordinate_view;
};

// Helper function to create device CoordinateView from a vector of interleaved
// points Returns both the underlying View (to keep memory alive) and the
// CoordinateView pts contains interleaved coordinates: [x0, y0, x1, y1, ...]
// for 2D or [x0, y0, z0, x1, y1, z1, ...] for 3D
inline DeviceCoordinates CreateDeviceCoordinateView(
  const std::vector<Real>& pts, CoordinateSystem coord_system, int dim = 2)
{
  int n = static_cast<int>(pts.size()) / dim;
  // Create host view from input data
  Kokkos::View<Real**, HostMemorySpace> coords_host("coords_host", n, dim);
  for (int i = 0; i < n; ++i) {
    for (int d = 0; d < dim; ++d) {
      coords_host(i, d) = pts[dim * i + d];
    }
  }
  // Copy to device - using default layout for device
  auto coords_device =
    Kokkos::View<Real**, DeviceMemorySpace>("coords_device", n, dim);
  DeepCopyMismatchLayouts(coords_device, coords_host);
  auto coords_view = pcms::MakeRank2View(coords_device);
  return DeviceCoordinates{coords_device, CoordinateView<DeviceMemorySpace>{
                                            coord_system, coords_view}};
}

// Evaluate field at explicit test points using a PointEvaluator and check
// results.
template <typename ExecutionSpace = DefaultExecutionSpace, typename Func>
void CheckEvaluation(const PointEvaluator<Real>& evaluator,
                     const Field<Real>& field, const std::vector<Real>& pts,
                     Func func, double abs_tol = 1e-10)
{
  int n = static_cast<int>(pts.size()) / 2;

  Kokkos::View<Real**, DeviceMemorySpace> out_device("out_device", n, 1);
  evaluator.Evaluate(field, MakeRank2View(out_device));
  auto out_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out_device);

  auto expected = EvaluateReferenceFunction<ExecutionSpace>(pts, func);

  for (int i = 0; i < n; ++i) {
    INFO("Point " << i << " (" << pts[2 * static_cast<size_t>(i)] << ", "
                  << pts[2 * static_cast<size_t>(i) + 1] << ")"
                  << " got=" << out_host(i, 0) << " expected=" << expected[i]);
    REQUIRE(out_host(i, 0) == Catch::Approx(expected[i]).margin(abs_tol));
  }
}

// Overload that creates the evaluator from any factory with
// CreatePointEvaluator and GetCoordinateSystem (e.g. LagrangeFunctionSpace,
// FieldEvaluatorFactory<Real>).
template <typename Factory, typename ExecutionSpace = DefaultExecutionSpace,
          typename Func>
void CheckEvaluation(const Factory& factory, const Field<Real>& field,
                     const std::vector<Real>& pts, Func func,
                     double abs_tol = 1e-10)
{
  auto device_coords =
    CreateDeviceCoordinateView(pts, factory->GetCoordinateSystem());
  auto evaluator = factory->template CreatePointEvaluator<Real>(
    EvaluationRequest::FromCoordinates(device_coords.coordinate_view));
  CheckEvaluation<ExecutionSpace>(*evaluator, field, pts, func, abs_tol);
}

// Evaluate field at points known to be outside the mesh and verify fill value.
inline void CheckFillMode(const PointEvaluator<Real>& evaluator,
                          const Field<Real>& field, Real fill_value,
                          const std::vector<Real>& outside_pts)
{
  int n = static_cast<int>(outside_pts.size()) / 2;

  Kokkos::View<Real**, DeviceMemorySpace> out_device("out_device", n, 1);
  evaluator.Evaluate(field, MakeRank2View(out_device));
  auto out_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out_device);

  for (int i = 0; i < n; ++i) {
    REQUIRE(out_host(i, 0) == fill_value);
  }
}

// Overload that creates the evaluator from any factory with FILL policy.
template <typename Factory>
void CheckFillMode(const Factory& factory, const Field<Real>& field,
                   Real fill_value, const std::vector<Real>& outside_pts)
{
  auto device_coords =
    CreateDeviceCoordinateView(outside_pts, factory->GetCoordinateSystem());
  OutOfBoundsPolicy policy{OutOfBoundsMode::FILL, fill_value};
  auto evaluator = factory->template CreatePointEvaluator<Real>(
    EvaluationRequest::FromCoordinates(device_coords.coordinate_view, policy));
  CheckFillMode(*evaluator, field, fill_value, outside_pts);
}

// Check a mix of inside/outside points. inside points verified against func,
// outside points against fill_value.
template <typename Factory, typename ExecutionSpace = DefaultExecutionSpace,
          typename Func>
void CheckEvaluationWithFill(const Factory& factory, const Field<Real>& field,
                             const std::vector<Real>& pts,
                             const std::vector<bool>& is_inside, Func func,
                             Real fill_value, double abs_tol = 1e-10)
{
  int n = static_cast<int>(pts.size()) / 2;
  REQUIRE(static_cast<int>(is_inside.size()) == n);

  auto device_coords =
    CreateDeviceCoordinateView(pts, factory->GetCoordinateSystem());
  OutOfBoundsPolicy policy{OutOfBoundsMode::FILL, fill_value};
  auto evaluator = factory->template CreatePointEvaluator<Real>(
    EvaluationRequest::FromCoordinates(device_coords.coordinate_view, policy));

  Kokkos::View<Real**, DeviceMemorySpace> out_device("out_device", n, 1);
  evaluator->Evaluate(field, MakeRank2View(out_device));
  auto out_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out_device);

  auto expected = EvaluateReferenceFunction<ExecutionSpace>(pts, func);

  for (int i = 0; i < n; ++i) {
    INFO("Point " << i << " (" << pts[2 * static_cast<size_t>(i)] << ", "
                  << pts[2 * static_cast<size_t>(i) + 1] << ")"
                  << " got=" << out_host(i, 0));
    if (is_inside[i]) {
      INFO(" expected=" << expected[i]);
      REQUIRE(out_host(i, 0) == Catch::Approx(expected[i]).margin(abs_tol));
    } else {
      REQUIRE(out_host(i, 0) == fill_value);
    }
  }
}

#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
// Evaluates source_field at the integrator's sample points and assembles the
// load vector.
inline void EvaluateAndAssemble(
  LinearFormIntegrator& integrator,
  const std::shared_ptr<const FunctionSpace>& source_space,
  const Field<Real>& source_field)
{
  const auto& pts = integrator.GetIntegrationPoints();
  const std::size_t npts = pts.GetValues().extent(0);
  auto evaluator = source_space->CreatePointEvaluator<Real>(
    EvaluationRequest::FromCoordinates(pts));
  Kokkos::View<Real**, DeviceMemorySpace> sampled("sampled", npts, 1);
  evaluator->Evaluate(source_field, MakeRank2View(sampled));
  integrator.Assemble(MakeConstRank2View(sampled));
}
#endif // PCMS_ENABLE_PETSC && PCMS_ENABLE_MESHFIELDS

} // namespace pcms::test

#endif // PCMS_TEST_FIELD_TEST_UTILS_H
