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
        // The evaluator writes the source field's stored basis; numpy buffers
        // carry no tag of their own, so claim that basis on the caller's
        // behalf.
        self.Evaluate(field, ValueView<Real, DeviceMemorySpace>(
                               field.GetData().GetValueBasis(),
                               MakeRank2View(output_device)));
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

  // Conservative L2 (Galerkin) projection between order-1 Lagrange spaces on
  // Omega_h 2D simplex meshes. The RHS is integrated exactly over the
  // intersection of the source and target meshes.
  py::class_<OmegaHConservativeProjection>(m, "OmegaHConservativeProjection")
    .def(py::init([](const FunctionSpace& source_space,
                     const FunctionSpace& target_space,
                     MassMatrixType mass_matrix_type) {
           return std::make_unique<OmegaHConservativeProjection>(
             source_space, target_space, mass_matrix_type);
         }),
         py::arg("source_space"), py::arg("target_space"),
         py::arg("mass_matrix_type") = MassMatrixType::Consistent,
         "Construct a mesh-intersection conservative projection. Mesh "
         "intersection, quadrature setup, and mass-matrix factorization happen "
         "here and are cached; call apply() repeatedly.")
    .def(
      "apply",
      [](const OmegaHConservativeProjection& self, const Field<Real>& source,
         Field<Real>& target) { self.Apply(source, target); },
      py::arg("source"), py::arg("target"),
      "Conservatively project the source field onto the target space using "
      "exact mesh-intersection quadrature.");

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
