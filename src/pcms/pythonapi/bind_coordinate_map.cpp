#include <memory>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/pythonapi/py_value_types.hpp"

namespace py = pybind11;

namespace pcms
{

void bind_coordinate_map_module(py::module& m)
{
  py::class_<BoundCoordinateMap, std::unique_ptr<BoundCoordinateMap>>(
    m, "BoundCoordinateMap", "A CoordinateMap bound to one set of query points")
    .def(
      "mapped_points",
      [](py::object self) {
        return ToPyCoordinateView(
          self.cast<const BoundCoordinateMap&>().MappedPoints(), self);
      },
      "The query points in the map's target coordinate system, read-only")
    .def("num_points", &BoundCoordinateMap::NumPoints)
    .def(
      "output_basis",
      [](const BoundCoordinateMap& self, const ValueBasis& stored) {
        return self.ValueTransform().OutputBasis(stored);
      },
      py::arg("stored"),
      "The basis values stored in `stored` are re-expressed in")
    .def(
      "transform_values",
      [](const BoundCoordinateMap& self, PyValueView values) {
        const auto& transform = self.ValueTransform();
        const ValueBasis transformed = transform.OutputBasis(values.basis);
        WithDeviceOutput(values.values, [&](auto out) {
          transform.Apply(ValueView<Real, DeviceMemorySpace>(values.basis, out),
                          ValueView<Real, DeviceMemorySpace>(transformed, out));
        });
        return PyValueView{transformed, values.values};
      },
      py::arg("values"),
      "Re-express `values` (one row per query point, over a writeable "
      "C-contiguous float64 array) in place, and return a ValueView of the "
      "same array tagged with the new basis");

  ImmutableClass<CoordinateMap>(
    m, "CoordinateMap",
    "Maps points from a source to a target coordinate system; values are "
    "re-expressed the other way")
    .def("get_source_coordinate_system",
         [](const Immutable<CoordinateMap>& self) {
           return ToPython(self.ptr->GetSourceCoordinateSystem());
         })
    .def("get_target_coordinate_system",
         [](const Immutable<CoordinateMap>& self) {
           return ToPython(self.ptr->GetTargetCoordinateSystem());
         })
    .def(
      "bind",
      [](const Immutable<CoordinateMap>& self,
         const PyCoordinateView& query_points) {
        return self.ptr->Bind(CopyCoordinatesToDevice(query_points).view);
      },
      py::arg("query_points"),
      "Do every per-point computation for points in the source coordinate "
      "system");

  m.def(
    "create_cartesian_to_cylindrical",
    []() { return Immutable<CoordinateMap>{CreateCartesianToCylindrical()}; },
    "Map from Cartesian (x, y, z) to cylindrical (R, theta, Z) points");
  m.def(
    "create_cylindrical_to_cartesian",
    []() { return Immutable<CoordinateMap>{CreateCylindricalToCartesian()}; },
    "Map from cylindrical (R, theta, Z) to Cartesian (x, y, z) points");
}

} // namespace pcms
