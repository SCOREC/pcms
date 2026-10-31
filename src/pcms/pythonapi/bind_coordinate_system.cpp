#include <functional>
#include <memory>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "numpy_view.h"
#include "py_value_types.h"

namespace py = pybind11;

namespace pcms
{

namespace
{

using SystemHolder = std::shared_ptr<CoordinateSystem>;

std::string Repr(const CoordinateSystem& system)
{
  return "<pcms.CoordinateSystem '" + std::string(system.Kind()) + "' dim " +
         std::to_string(system.Dimension()) + ">";
}

PyCoordinateView MakeCoordinateView(const SystemHolder& system,
                                    const NumpyInput<Real>& coordinates)
{
  const auto input = NumpyToHostView<2, Real>(coordinates);
  // pcms may keep the coordinates, so the view owns a copy
  auto values = ToNumpy(input);
  detail::SetReadOnly(values);
  return PyCoordinateView{
    ResolveCoordinateSystem(CartesianIfNone(system), input.extent(1)), values};
}

} // namespace

void bind_coordinate_system_module(py::module& m)
{
  // Concrete systems are registered as subclasses so Python sees the dynamic
  // type of any system pcms returns (isinstance(cs, pcms.csys.Cartesian)).
  py::class_<CoordinateSystem, SystemHolder>(
    m, "CoordinateSystem",
    "What coordinate numbers mean. Create one through the classes in "
    "pcms.csys.")
    .def(
      "kind",
      [](const CoordinateSystem& self) { return std::string(self.Kind()); },
      "Family label, e.g. 'cartesian'")
    .def(
      "dimension",
      [](const CoordinateSystem& self) { return self.Dimension(); },
      "Number of coordinate columns")
    .def("has_orthogonal_basis",
         [](const CoordinateSystem& self) { return self.HasOrthogonalBasis(); })
    .def(
      "has_unit_scale_factors",
      [](const CoordinateSystem& self) { return self.HasUnitScaleFactors(); })
    .def(
      "__eq__",
      [](const SystemHolder& a, const SystemHolder& b) {
        return SameCoordinateSystem(a, b);
      },
      py::is_operator())
    .def("__hash__",
         [](const CoordinateSystem& self) {
           return std::hash<std::string>{}(std::string(self.Kind()) + "/" +
                                           std::to_string(self.Dimension()));
         })
    .def("__repr__", [](const CoordinateSystem& self) { return Repr(self); });

  auto csys = m.def_submodule("csys", "Concrete coordinate systems");
  py::class_<csys::Cartesian, CoordinateSystem,
             std::shared_ptr<csys::Cartesian>>(csys, "Cartesian")
    .def_static(
      "create",
      [](int dim) { return ToPyHolder(csys::Cartesian::Create(dim)); },
      py::arg("dim"), "The Cartesian coordinate system of dimension 1..3");
  py::class_<csys::CylindricalRThetaZ, CoordinateSystem,
             std::shared_ptr<csys::CylindricalRThetaZ>>(csys,
                                                        "CylindricalRThetaZ")
    .def_static(
      "create", []() { return ToPyHolder(csys::CylindricalRThetaZ::Create()); },
      "Cylindrical (r, theta, z) with theta = atan2(y, x); dimension 3");
  py::class_<csys::CylindricalRZ, CoordinateSystem,
             std::shared_ptr<csys::CylindricalRZ>>(csys, "CylindricalRZ")
    .def_static(
      "create", []() { return ToPyHolder(csys::CylindricalRZ::Create()); },
      "(R, Z) poloidal section; dimension 2");

  m.def(
    "same_coordinate_system",
    [](const SystemHolder& a, const SystemHolder& b) {
      return SameCoordinateSystem(a, b);
    },
    py::arg("a"), py::arg("b"));
  m.def(
    "has_identity_metric",
    [](const CoordinateSystem& system) { return HasIdentityMetric(system); },
    py::arg("system"));
  py::class_<PyCoordinateView>(m, "CoordinateView",
                               "Coordinates tagged with their coordinate "
                               "system")
    .def(py::init(&MakeCoordinateView), py::arg("coordinate_system"),
         py::arg("coordinates"),
         "Coordinates from a 2D array (num_points × dim), which is copied. A "
         "coordinate_system of None is Cartesian in the array's dimension.")
    .def(py::init([](const NumpyInput<Real>& coordinates) {
           return MakeCoordinateView(nullptr, coordinates);
         }),
         py::arg("coordinates"),
         "Cartesian coordinates from a 2D array (num_points × dim), which is "
         "copied.")
    .def("get_coordinate_system",
         [](const PyCoordinateView& self) { return ToPyHolder(self.system); })
    .def(
      "get_values", [](const PyCoordinateView& self) { return self.values; },
      "The coordinates as a read-only 2D array (num_points × dim)")
    .def("__repr__", [](const PyCoordinateView& self) {
      return "<pcms.CoordinateView " + std::to_string(self.values.shape(0)) +
             " points in " + Repr(*self.system) + ">";
    });
}

} // namespace pcms
