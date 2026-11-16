#include <pybind11/pybind11.h>
#include "pcms/utility/arrays.h"
#include "pcms/field/out_of_bounds_policy.h"

namespace py = pybind11;

namespace pcms
{

void bind_omega_h_field(py::module& m)
{
  // Bind OutOfBoundsMode enum
  py::enum_<OutOfBoundsMode>(m, "OutOfBoundsMode")
    .value("ERROR", OutOfBoundsMode::ERROR,
           "Raise error when points are out of bounds")
    .value("FILL", OutOfBoundsMode::FILL,
           "Fill with a specified value when points are out of bounds")
    .value("NEAREST_BOUNDARY", OutOfBoundsMode::NEAREST_BOUNDARY,
           "Clamp to nearest boundary cell (extrapolate)")
    .export_values();

  // Bind OutOfBoundsPolicy struct
  py::class_<OutOfBoundsPolicy>(m, "OutOfBoundsPolicy")
    .def(py::init<>(), "Default constructor (mode=ERROR, fill_value=0)")
    .def(py::init([](OutOfBoundsMode mode, Real fill_value) {
           OutOfBoundsPolicy p;
           p.mode = mode;
           p.fill_value = fill_value;
           return p;
         }),
         py::arg("mode"), py::arg("fill_value") = Real(0),
         "Construct with explicit mode and optional fill value")
    .def_readwrite("mode", &OutOfBoundsPolicy::mode,
                   "Out-of-bounds handling mode")
    .def_readwrite("fill_value", &OutOfBoundsPolicy::fill_value,
                   "Fill value (used only when mode == FILL)");

  // NOTE: MeshFieldsAdapter2<Real> and FieldT<Real> have been removed from
  // the C++ API. Python now works with Field<Real> objects created from
  // FunctionSpace-derived factories such as LagrangeFunctionSpace.
}

} // namespace pcms
