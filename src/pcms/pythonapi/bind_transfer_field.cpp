#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pcms/configuration.hpp"
#include "pcms/transfer/copy.hpp"
#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/transfer/interpolator.hpp"
#include "pcms/transfer/projection.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/pythonapi/numpy_array_transform.hpp"
#include "pcms/pythonapi/py_value_types.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/geometry/quadrature.hpp"
#include "pcms/transfer/mass_matrix_type.hpp"
#include <cstdint>
#include <optional>
#if defined(PCMS_ENABLE_OMEGA_H)
#include "pcms/transfer/omega_h_intersection_quadrature.hpp"
#endif

namespace py = pybind11;

namespace pcms
{

void bind_transfer_field_module(py::module& m)
{
  py::class_<PointEvaluator<Real>, std::unique_ptr<PointEvaluator<Real>>>(
    m, "PointEvaluator",
    "Evaluates fields at the query points it was created for")
    .def(
      "evaluate",
      [](const PointEvaluator<Real>& self, const Field<Real>& field,
         PyValueView output) {
        WithDeviceOutput(output.values, [&](auto out) {
          self.Evaluate(field,
                        ValueView<Real, DeviceMemorySpace>(output.basis, out));
        });
      },
      py::arg("field"), py::arg("output"),
      "Evaluate the field into `output`, a ValueView over a writeable "
      "C-contiguous float64 array (num_points × num_components) whose basis "
      "must be output_basis() of the field's stored basis")
    .def(
      "evaluate",
      [](const PointEvaluator<Real>& self, const Field<Real>& field,
         py::array output) {
        const ValueBasis basis =
          self.OutputBasis(field.GetData().GetValueBasis());
        RequireScalarForPlainArray("evaluate", basis);
        WithDeviceOutput(output, [&](auto out) {
          self.Evaluate(field, ValueView<Real, DeviceMemorySpace>(basis, out));
        });
      },
      py::arg("field"), py::arg("output").noconvert(),
      "Evaluate a scalar field into a writeable C-contiguous float64 array "
      "(num_points × num_components). A 1D array (num_points,) is accepted "
      "when the output has one component.")
    .def("num_points", &PointEvaluator<Real>::NumPoints)
    .def("output_basis", &PointEvaluator<Real>::OutputBasis, py::arg("stored"),
         "The basis this evaluator writes for a field stored in `stored`")
    .def(
      "filled_points",
      [](py::object self) {
        return IndicesToNumpy(
          self.cast<const PointEvaluator<Real>&>().FilledPoints(), self);
      },
      "Indices of the query points that get the out-of-bounds fill value, "
      "read-only");

  py::class_<TransferOperator<Real>>(m, "TransferOperator")
    .def("apply", &TransferOperator<Real>::Apply, py::arg("source"),
         py::arg("target"), "Transfer the source field onto the target");

  m.def(
    "create_interpolator",
    [](const Immutable<PointEvaluatorFactory>& source,
       const Immutable<FunctionSpace>& target_space, OutOfBoundsPolicy policy) {
      return CreateInterpolator<Real>(source.ptr, target_space.ptr, policy);
    },
    py::arg("source"), py::arg("target_space"),
    py::arg("policy") = OutOfBoundsPolicy{},
    "Create an interpolator from any PointEvaluatorFactory (a FunctionSpace or "
    "an operator source such as gradient()) to the target space's DOF "
    "holders. Localization is performed here and cached; call apply() "
    "repeatedly. To transfer across coordinate systems, bind with a map "
    "(source.get_geometry().bind(target_space, policy, map)) and pass the "
    "BoundPoints.");
  m.def(
    "create_interpolator",
    [](const Immutable<PointEvaluatorFactory>& source,
       const Immutable<FunctionSpace>& target_space,
       const Immutable<BoundPoints>& bound) {
      return CreateInterpolator<Real>(source.ptr, target_space.ptr, bound.ptr);
    },
    py::arg("source"), py::arg("target_space"), py::arg("bound"),
    "Create an interpolator that reuses `bound`, a binding of the target "
    "space's DOF holders by `source`'s geometry "
    "(source.get_geometry().bind(target_space)), so a field and its gradient "
    "can share one localization.");

  m.def(
    "create_copy",
    [](const Immutable<FunctionSpace>& source_space,
       const Immutable<FunctionSpace>& target_space) {
      return CreateCopy<Real>(source_space.ptr, target_space.ptr);
    },
    py::arg("source_space"), py::arg("target_space"),
    "Create a copy operator for function spaces sharing a layout.");

  // Bound before the projections, which take it as a default argument value.
  py::enum_<MassMatrixType>(m, "MassMatrixType")
    .value("Consistent", MassMatrixType::Consistent,
           "Full Galerkin mass matrix")
    .value("Lumped", MassMatrixType::Lumped,
           "Row-sum lumped diagonal mass matrix (exact diagonal solve; "
           "conserves the integral but adds lumping error)")
    .export_values();

  py::class_<Quadrature>(m, "Quadrature",
                         "Points bound on a geometry, with one weight each")
    .def("num_points",
         [](const Quadrature& self) { return self.points->NumPoints(); });

  m.def(
    "create_gauss_quadrature",
    [](const Immutable<Geometry>& geometry, int order) {
      return CreateGaussQuadrature(geometry.ptr, order);
    },
    py::arg("geometry"), py::arg("order"),
    "The Gauss rule of each cell of the geometry, exact for polynomials of "
    "degree `order`.");

  m.def(
    "create_monte_carlo_quadrature",
    [](const Immutable<Geometry>& geometry, int samples_per_cell,
       std::uint64_t seed) {
      return CreateMonteCarloQuadrature(geometry.ptr, samples_per_cell, seed);
    },
    py::arg("geometry"), py::arg("samples_per_cell"),
    py::arg("seed") = DefaultQuadratureSeed,
    "Uniform random samples in each cell of the geometry.");

#if defined(PCMS_ENABLE_OMEGA_H)
  m.def(
    "create_intersection_quadrature",
    [](const Immutable<Geometry>& source, const Immutable<Geometry>& target,
       int order) {
      return omega_h::CreateIntersectionQuadrature(source.ptr, target.ptr,
                                                   order);
    },
    py::arg("source"), py::arg("target"), py::arg("order"),
    "The quadrature of the intersection of two Omega_h meshes, bound on the "
    "target geometry, exact for polynomials of degree `order` on the pieces.");
#endif

  m.def(
    "create_projection",
    [](const Immutable<PointEvaluatorFactory>& source,
       const Immutable<FunctionSpace>& target_space,
       std::optional<Quadrature> quadrature, bool control_variate,
       MassMatrixType mass_matrix_type,
       const std::optional<Immutable<CoordinateMap>>& to_source) {
      return CreateProjection(
        source.ptr, target_space.ptr,
        ProjectionOptions{.quadrature = std::move(quadrature),
                          .control_variate = control_variate,
                          .mass_matrix_type = mass_matrix_type,
                          .to_source = to_source ? to_source->ptr : nullptr});
    },
    py::arg("source"), py::arg("target_space"),
    py::arg("quadrature") = py::none(), py::arg("control_variate") = false,
    py::arg("mass_matrix_type") = MassMatrixType::Consistent,
    py::arg("to_source") = py::none(),
    "Create a Galerkin projection of `source` onto the target space, with "
    "the right-hand side on `quadrature` (a quadrature of the target space's "
    "geometry; None is the Gauss rule of degree 2). control_variate projects "
    "only the residual of the source from its interpolation. to_source maps "
    "the target's coordinate system to the source's, or None if they are the "
    "same. Setup is cached; call apply() repeatedly.");

#if defined(PCMS_ENABLE_OMEGA_H_MESHFIELDS)

  m.def(
    "create_conservative_projection",
    [](const Immutable<FunctionSpace>& source_space,
       const Immutable<FunctionSpace>& target_space,
       MassMatrixType mass_matrix_type) {
      return CreateConservativeProjection(
        source_space.ptr, target_space.ptr,
        ConservativeProjectionOptions{.mass_matrix_type = mass_matrix_type});
    },
    py::arg("source_space"), py::arg("target_space"),
    py::arg("mass_matrix_type") = MassMatrixType::Consistent,
    "Create a mesh-intersection conservative projection between order-1 "
    "Lagrange spaces on Omega_h meshes. Mesh intersection, quadrature setup, "
    "and mass-matrix factorization happen here and are cached; call apply() "
    "repeatedly.");

  m.def(
    "create_control_variate_projection",
    [](const Immutable<PointEvaluatorFactory>& source,
       const Immutable<FunctionSpace>& target_space, int samples_per_element,
       std::uint64_t seed, MassMatrixType mass_matrix_type,
       const std::optional<Immutable<CoordinateMap>>& to_source) {
      return CreateControlVariateProjection(
        source.ptr, target_space.ptr,
        ControlVariateProjectionOptions{
          .samples_per_element = samples_per_element,
          .seed = seed,
          .mass_matrix_type = mass_matrix_type,
          .to_source = to_source ? to_source->ptr : nullptr});
    },
    py::arg("source"), py::arg("target_space"), py::arg("samples_per_element"),
    py::arg("seed") = DefaultQuadratureSeed,
    py::arg("mass_matrix_type") = MassMatrixType::Consistent,
    py::arg("to_source") = py::none(),
    "Create a Monte Carlo (control-variate) conservative projection onto an "
    "order-1 Omega_h Lagrange space. Sample-point generation, the "
    "control-variate interpolator, and the mass-matrix factorization are "
    "cached; call apply() repeatedly. to_source maps the target's coordinate "
    "system to the source's, or None if they are the same.");
#endif
}

} // namespace pcms
