#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pcms/transfer_field2.h"
#include "pcms/field.h"

namespace py = pybind11;

namespace pcms
{

void bind_transfer_field2_module(py::module& m)
{
  // Only bind for double since that's the only FieldT type exposed to Python
  // (FieldT<float> and FieldT<Real> are not bound in bind_field_base.cpp)

  m.def(
    "copy_field",
    [](const FieldT<double>& source, FieldT<double>& target) {
      copy_field2(source, target);
    },
    py::arg("source"), py::arg("target"),
    "Copy field data from source to target. Source and target must be of the "
    "same type.");

  m.def(
    "interpolate_field",
    [](const FieldT<double>& source, FieldT<double>& target) {
      interpolate_field2(source, target);
    },
    py::arg("source"), py::arg("target"),
    "Interpolate field from source to target. Coordinate systems must match.");
}

} // namespace pcms