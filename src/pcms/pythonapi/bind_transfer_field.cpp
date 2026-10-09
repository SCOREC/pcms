#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pcms/configuration.h"
#include "pcms/transfer/copy.h"
#include "pcms/field/field.h"
#include "pcms/field/function_space.h"
#include "pcms/field/point_evaluator.h"
#include "../transfer/interpolator.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "numpy_array_transform.h"
#include "pcms/utility/types.h"
#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
#include "pcms/transfer/mass_matrix_type.hpp"
#include "pcms/transfer/mass_smoother.hpp"
#include "pcms/transfer/mesh_intersection.hpp"
#include "pcms/transfer/omega_h_intersection_quadrature.hpp"
#include "pcms/transfer/omega_h_conservative_projection.hpp"
#include "pcms/transfer/omega_h_control_variate_projection.hpp"
#include "pcms/transfer/omega_h_mc_rhs_integrator.hpp"
#include <cstdint>
#endif

namespace py = pybind11;

namespace pcms
{

void bind_transfer_field_module(py::module& m)
{
  py::class_<PointEvaluator<Real>, std::unique_ptr<PointEvaluator<Real>>>(
    m, "PointEvaluator")
    .def(
      "evaluate",
      [](const PointEvaluator<Real>& self, const Field<Real>& field,
         py::array_t<Real> output) {
        auto output_view = numpy_to_kokkos_view_2d<Real>(output);
        auto output_device = Kokkos::View<Real**, DeviceMemorySpace>(
          "output_device", output_view.extent(0), output_view.extent(1));
        DeepCopyMismatchLayouts(output_device, output_view);
        auto output_rank2 = MakeRank2View(output_device);
        self.Evaluate(field, output_rank2);
        DeepCopyMismatchLayouts(output_view, output_device);
      },
      py::arg("field"), py::arg("output"),
      "Evaluate the given field at the cached query coordinates into a "
      "preallocated 2D numpy array of shape "
      "(num_query_points, num_components).");

  // Bind Interpolator<Real>: construct once per source×target function-space
  // pair (localization happens at construction), then call apply() repeatedly
  // for different field states at zero additional localization cost.
  py::class_<Interpolator<Real>>(m, "Interpolator")
    .def(
      py::init([](const FunctionSpace& source_space,
                  const FunctionSpace& target_space, OutOfBoundsPolicy policy) {
        return Interpolator<Real>(source_space, target_space, policy);
      }),
      py::arg("source_space"), py::arg("target_space"),
      py::arg("policy") = OutOfBoundsPolicy{},
      "Construct an interpolator. Localization is performed here and cached. "
      "Call apply() repeatedly without re-localizing.")
    .def(
      "apply",
      [](const Interpolator<Real>& self, const Field<Real>& source,
         Field<Real>& target) { self.Apply(source, target); },
      py::arg("source"), py::arg("target"),
      "Interpolate source field to target DOF locations (cheap; reuses cached "
      "localization).");

  py::class_<Copy<Real>>(m, "Copy")
    .def(py::init([](const FunctionSpace& source_space,
                     const FunctionSpace& target_space) {
           return Copy<Real>(source_space, target_space);
         }),
         py::arg("source_space"), py::arg("target_space"),
         "Construct a copy operator for compatible function spaces.")
    .def(
      "apply",
      [](const Copy<Real>& self, const Field<Real>& source,
         Field<Real>& target) { self.Apply(source, target); },
      py::arg("source"), py::arg("target"),
      "Copy source field data to target field (same layout required).");

#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
  // Mass-matrix formulation for the Galerkin projections. Bound before the
  // projection classes so it can be used as a default argument value.
  py::enum_<MassMatrixType>(m, "MassMatrixType")
    .value("Consistent", MassMatrixType::Consistent,
           "Full Galerkin mass matrix")
    .value("Lumped", MassMatrixType::Lumped,
           "Row-sum lumped diagonal mass matrix (exact diagonal solve; "
           "conserves the integral but adds lumping error)")
    .export_values();

  // Source-target overlap map, built once per mesh pair and shared by every
  // projection on that pair.
  py::class_<MeshIntersection, std::shared_ptr<MeshIntersection>>(
    m, "MeshIntersection",
    "Overlap map between a source and a target mesh (for each target element, "
    "the source elements it intersects). Build once per mesh pair and pass to "
    "every projection between spaces on those meshes.");
  py::class_<OmegaHMeshIntersection, MeshIntersection,
             std::shared_ptr<OmegaHMeshIntersection>>(m,
                                                      "OmegaHMeshIntersection")
    .def(py::init([](const FunctionSpace& source_space,
                     const FunctionSpace& target_space) {
           return std::make_shared<OmegaHMeshIntersection>(source_space,
                                                           target_space);
         }),
         py::arg("source_space"), py::arg("target_space"),
         "Intersect the two spaces' Omega_h meshes, reusing the source "
         "space's point search. Any spaces on the same two meshes may share "
         "the result, whatever their order.");
  m.def(
    "intersect_meshes",
    [](const FunctionSpace& source_space, const FunctionSpace& target_space) {
      return IntersectMeshes(source_space, target_space);
    },
    py::arg("source_space"), py::arg("target_space"),
    "MeshIntersection for the two spaces' meshes, dispatched on their type.");

  // Quadrature of the intersection for one pair of spaces: points, target
  // weights, source elements and a search-free source evaluator. Shared by
  // every projection between the same two spaces.
  py::class_<OmegaHIntersectionQuadrature,
             std::shared_ptr<OmegaHIntersectionQuadrature>>(
    m, "OmegaHIntersectionQuadrature")
    .def(py::init([](const FunctionSpace& source_space,
                     const FunctionSpace& target_space,
                     std::shared_ptr<MeshIntersection> intersection) {
           return std::make_shared<OmegaHIntersectionQuadrature>(
             source_space, target_space, std::move(intersection));
         }),
         py::arg("source_space"), py::arg("target_space"),
         py::arg("intersection") = nullptr,
         "Integration points of the source-target intersection on the target "
         "mesh, with the target basis weights, the source element of each "
         "point and a source evaluator that needs no point search. Depends on "
         "both spaces (orders included) but not on the mass matrix; pass to "
         "every OmegaHConservativeProjection between the same two spaces.")
    .def("num_points", &OmegaHIntersectionQuadrature::GetNumPoints);

  // Conservative L2 (Galerkin) projection between order-1 Lagrange spaces on
  // Omega_h 2D simplex meshes. The RHS is integrated exactly over the
  // intersection of the source and target meshes.
  py::class_<OmegaHConservativeProjection>(m, "OmegaHConservativeProjection")
    .def(py::init([](const FunctionSpace& source_space,
                     const FunctionSpace& target_space,
                     MassMatrixType mass_matrix_type,
                     std::shared_ptr<MeshIntersection> intersection,
                     std::shared_ptr<OmegaHIntersectionQuadrature> quadrature) {
           return std::make_unique<OmegaHConservativeProjection>(
             source_space, target_space, mass_matrix_type,
             std::move(intersection), std::move(quadrature));
         }),
         py::arg("source_space"), py::arg("target_space"),
         py::arg("mass_matrix_type") = MassMatrixType::Consistent,
         py::arg("intersection") = nullptr, py::arg("quadrature") = nullptr,
         "Construct a mesh-intersection conservative projection. Mesh "
         "intersection and quadrature (unless passed in) and the mass-matrix "
         "factorization happen here and are cached; call apply() repeatedly.")
    .def(
      "apply",
      [](const OmegaHConservativeProjection& self, const Field<Real>& source,
         Field<Real>& target) { self.Apply(source, target); },
      py::arg("source"), py::arg("target"),
      "Conservatively project the source field onto the target space using "
      "exact mesh-intersection quadrature.");

  // Same-space conservative smoother; one apply() is one sweep of inv(M_L) M.
  py::class_<MassSmoother>(m, "MassSmoother")
    .def(py::init(
           [](const FunctionSpace& space) { return BuildMassSmoother(space); }),
         py::arg("space"),
         "Conservative, bounds-preserving smoother on one function space: "
         "apply() performs one sweep u <- inv(M_L) M u (Farrell et al. 2009, "
         "Eqs. 37-38, applied to the whole field). The mass matrix is "
         "assembled once here; call apply() in a loop for more sweeps.")
    .def(
      "apply",
      [](const MassSmoother& self, Field<Real>& field) { self.Apply(field); },
      py::arg("field"), "Smooth the field in place (one sweep).");

  // Sampling strategy for the Monte Carlo RHS integrator.
  py::enum_<MonteCarloSampling>(m, "MonteCarloSampling")
    .value("UniformRandom", MonteCarloSampling::UniformRandom,
           "Independent pseudo-random uniform draws per element")
    .export_values();

  // Variance-reduced Monte Carlo Galerkin projection. The source field is
  // first interpolated onto the target space as a control variate; only the
  // residual is integrated stochastically, so sampling noise is small (and
  // exactly zero when the source already lives in the target space).
  py::class_<OmegaHControlVariateProjection>(m,
                                             "OmegaHControlVariateProjection")
    .def(py::init([](const FunctionSpace& source_space,
                     const FunctionSpace& target_space, int samples_per_element,
                     MonteCarloSampling sampling, std::uint64_t seed,
                     MassMatrixType mass_matrix_type) {
           return std::make_unique<OmegaHControlVariateProjection>(
             source_space, target_space, samples_per_element, sampling, seed,
             mass_matrix_type);
         }),
         py::arg("source_space"), py::arg("target_space"),
         py::arg("samples_per_element"),
         py::arg("sampling") = MonteCarloSampling::UniformRandom,
         py::arg("seed") = std::uint64_t(8675309),
         py::arg("mass_matrix_type") = MassMatrixType::Consistent,
         "Construct a Monte Carlo (control-variate) conservative projection. "
         "Sample-point generation, the control-variate interpolator, and the "
         "mass-matrix factorization are cached; call apply() repeatedly.")
    .def(
      "apply",
      [](const OmegaHControlVariateProjection& self, const Field<Real>& source,
         Field<Real>& target) { self.Apply(source, target); },
      py::arg("source"), py::arg("target"),
      "Project the source field onto the target space using Monte Carlo "
      "integration with control-variate variance reduction.");
#endif
}

} // namespace pcms
