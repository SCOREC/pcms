#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pcms/configuration.h"
#include "pcms/transfer/copy.h"
#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/field.h"
#include "pcms/field/function_space.h"
#include "pcms/field/point_evaluator.h"
#include "../transfer/interpolator.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "numpy_array_transform.h"
#include "py_value_types.h"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/transfer/transfer_operator.hpp"
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

  // Interpolator<Real>: construct once per source × target pair (localization
  // happens at construction), then apply() repeatedly for different field
  // states at zero additional localization cost.
  py::class_<Interpolator<Real>, TransferOperator<Real>>(m, "Interpolator")
    .def(
      py::init([](const Immutable<PointEvaluatorFactory>& source,
                  const Immutable<FunctionSpace>& target_space,
                  OutOfBoundsPolicy policy) {
        return std::make_unique<Interpolator<Real>>(*source.ptr,
                                                    *target_space.ptr, policy);
      }),
      py::arg("source"), py::arg("target_space"),
      py::arg("policy") = OutOfBoundsPolicy{}, py::keep_alive<1, 2>(),
      py::keep_alive<1, 3>(),
      "Construct an interpolator from any PointEvaluatorFactory (a "
      "FunctionSpace or an operator source such as gradient()) to the target "
      "space's DOF holders. Localization is performed here and cached. To "
      "transfer across coordinate systems, bind with a map "
      "(source.get_geometry().bind(target_space, policy, map)) and pass the "
      "bound source.")
    .def(py::init([](const Immutable<PointEvaluatorFactory>& source,
                     const Immutable<FunctionSpace>& target_space,
                     const Immutable<BoundSource>& bound) {
           return std::make_unique<Interpolator<Real>>(
             *source.ptr, *target_space.ptr, bound.ptr);
         }),
         py::arg("source"), py::arg("target_space"), py::arg("bound"),
         py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
         "Construct an interpolator that reuses `bound`, a binding of the "
         "target space's DOF holders by `source`'s geometry "
         "(source.get_geometry().bind(target_space)), so a field and its "
         "gradient can share one localization.");

  py::class_<Copy<Real>, TransferOperator<Real>>(m, "Copy").def(
    py::init([](const Immutable<FunctionSpace>& source_space,
                const Immutable<FunctionSpace>& target_space) {
      return std::make_unique<Copy<Real>>(*source_space.ptr, *target_space.ptr);
    }),
    py::arg("source_space"), py::arg("target_space"),
    "Construct a copy operator for function spaces sharing a layout.");

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
  py::class_<OmegaHConservativeProjection, TransferOperator<Real>>(
    m, "OmegaHConservativeProjection")
    .def(py::init([](const Immutable<FunctionSpace>& source_space,
                     const Immutable<FunctionSpace>& target_space,
                     MassMatrixType mass_matrix_type) {
           return std::make_unique<OmegaHConservativeProjection>(
             *source_space.ptr, *target_space.ptr, mass_matrix_type);
         }),
         py::arg("source_space"), py::arg("target_space"),
         py::arg("mass_matrix_type") = MassMatrixType::Consistent,
         py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
         "Construct a mesh-intersection conservative projection. Mesh "
         "intersection, quadrature setup, and mass-matrix factorization happen "
         "here and are cached; call apply() repeatedly.");

  // Variance-reduced Monte Carlo Galerkin projection. The source field is
  // first interpolated onto the target space as a control variate; only the
  // residual is integrated stochastically, so sampling noise is small (and
  // exactly zero when the source already lives in the target space).
  py::class_<OmegaHControlVariateProjection, TransferOperator<Real>>(
    m, "OmegaHControlVariateProjection")
    .def(py::init([](const Immutable<FunctionSpace>& source_space,
                     const Immutable<FunctionSpace>& target_space,
                     int samples_per_element, std::uint64_t seed,
                     MassMatrixType mass_matrix_type,
                     const std::optional<Immutable<CoordinateMap>>& to_source) {
           return std::make_unique<OmegaHControlVariateProjection>(
             *source_space.ptr, *target_space.ptr, samples_per_element, seed,
             mass_matrix_type, to_source ? to_source->ptr : nullptr);
         }),
         py::arg("source_space"), py::arg("target_space"),
         py::arg("samples_per_element"),
         py::arg("seed") = std::uint64_t(8675309),
         py::arg("mass_matrix_type") = MassMatrixType::Consistent,
         py::arg("to_source") = py::none(), py::keep_alive<1, 2>(),
         py::keep_alive<1, 3>(),
         "Construct a Monte Carlo (control-variate) conservative projection. "
         "Sample-point generation, the control-variate interpolator, and the "
         "mass-matrix factorization are cached; call apply() repeatedly. "
         "to_source maps the target's coordinate system to the source's, or "
         "None if they are the same.");
#endif
}

} // namespace pcms
