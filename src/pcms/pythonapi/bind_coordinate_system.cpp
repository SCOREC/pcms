#include <functional>
#include <memory>
#include <optional>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "numpy_view.hpp"
#include "py_value_types.hpp"

namespace py = pybind11;

namespace pcms
{

namespace
{

std::string Repr(const CoordinateSystem& system)
{
  return "<pcms.CoordinateSystem '" + std::string(system.Kind()) + "' dim " +
         std::to_string(system.Dimension()) + ">";
}

PyCoordinateView MakeCoordinateView(
  const std::optional<Immutable<CoordinateSystem>>& system,
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
  py::class_<Immutable<CoordinateSystem>>(
    m, "CoordinateSystem",
    "What coordinate numbers mean. Create one through the classes in "
    "pcms.csys.")
    .def(
      "kind",
      [](const Immutable<CoordinateSystem>& self) {
        return std::string(self.ptr->Kind());
      },
      "Family label, e.g. 'cartesian'")
    .def(
      "dimension",
      [](const Immutable<CoordinateSystem>& self) {
        return self.ptr->Dimension();
      },
      "Number of coordinate columns")
    .def("has_orthogonal_basis",
         [](const Immutable<CoordinateSystem>& self) {
           return self.ptr->HasOrthogonalBasis();
         })
    .def("has_unit_scale_factors",
         [](const Immutable<CoordinateSystem>& self) {
           return self.ptr->HasUnitScaleFactors();
         })
    .def(
      "__eq__",
      [](const Immutable<CoordinateSystem>& a,
         const Immutable<CoordinateSystem>& b) {
        return SameCoordinateSystem(a.ptr, b.ptr);
      },
      py::is_operator())
    .def("__hash__",
         [](const Immutable<CoordinateSystem>& self) {
           return std::hash<std::string>{}(
             std::string(self.ptr->Kind()) + "/" +
             std::to_string(self.ptr->Dimension()));
         })
    .def("__repr__", [](const Immutable<CoordinateSystem>& self) {
      return Repr(*self.ptr);
    });

  auto csys = m.def_submodule("csys", "Concrete coordinate systems");
  py::class_<CartesianHandle, Immutable<CoordinateSystem>>(csys, "Cartesian")
    .def_static(
      "create",
      [](int dim) { return CartesianHandle{{csys::Cartesian::Create(dim)}}; },
      py::arg("dim"), "The Cartesian coordinate system of dimension 1..3");
  py::class_<CylindricalRThetaZHandle, Immutable<CoordinateSystem>>(
    csys, "CylindricalRThetaZ")
    .def_static(
      "create",
      []() {
        return CylindricalRThetaZHandle{{csys::CylindricalRThetaZ::Create()}};
      },
      "Cylindrical (r, theta, z) with theta = atan2(y, x); dimension 3");
  py::class_<CylindricalRZHandle, Immutable<CoordinateSystem>>(csys,
                                                               "CylindricalRZ")
    .def_static(
      "create",
      []() { return CylindricalRZHandle{{csys::CylindricalRZ::Create()}}; },
      "(R, Z) poloidal section; dimension 2");

  m.def(
    "same_coordinate_system",
    [](const Immutable<CoordinateSystem>& a,
       const Immutable<CoordinateSystem>& b) {
      return SameCoordinateSystem(a.ptr, b.ptr);
    },
    py::arg("a"), py::arg("b"));
  m.def(
    "has_identity_metric",
    [](const Immutable<CoordinateSystem>& system) {
      return HasIdentityMetric(*system.ptr);
    },
    py::arg("system"));
  py::class_<PyCoordinateView>(m, "CoordinateView",
                               "Coordinates tagged with their coordinate "
                               "system")
    .def(py::init(&MakeCoordinateView), py::arg("coordinate_system"),
         py::arg("coordinates"),
         "Coordinates from a 2D array (num_points × dim), which is copied. A "
         "coordinate_system of None is Cartesian in the array's dimension.")
    .def(py::init([](const NumpyInput<Real>& coordinates) {
           return MakeCoordinateView(std::nullopt, coordinates);
         }),
         py::arg("coordinates"),
         "Cartesian coordinates from a 2D array (num_points × dim), which is "
         "copied.")
    .def("get_coordinate_system",
         [](const PyCoordinateView& self) { return ToPython(self.system); })
    .def(
      "get_values", [](const PyCoordinateView& self) { return self.values; },
      "The coordinates as a read-only 2D array (num_points × dim)")
    .def("__repr__", [](const PyCoordinateView& self) {
      return "<pcms.CoordinateView " + std::to_string(self.values.shape(0)) +
             " points in " + Repr(*self.system) + ">";
    });
}

} // namespace pcms
