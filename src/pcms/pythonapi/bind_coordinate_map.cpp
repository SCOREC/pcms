#include <memory>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/transformed_source.hpp"
#include "pcms/field/value_view.hpp"
#include "py_value_types.h"

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
    .def("output_basis", &BoundCoordinateMap::OutputBasis, py::arg("stored"),
         "The basis values stored in `stored` are re-expressed in")
    .def(
      "transform_values",
      [](const BoundCoordinateMap& self, PyValueView values) {
        ValueBasis transformed;
        WithDeviceOutput(values.values, [&](auto out) {
          transformed = self
                          .TransformValues(ValueView<Real, DeviceMemorySpace>(
                            values.basis, out))
                          .GetBasis();
        });
        return PyValueView{transformed, values.values};
      },
      py::arg("values"),
      "Re-express `values` (one row per query point, over a writeable "
      "C-contiguous float64 array) in place, and return a ValueView of the "
      "same array tagged with the new basis");

  py::class_<CoordinateMap, std::shared_ptr<CoordinateMap>>(
    m, "CoordinateMap",
    "Maps points from a source to a target coordinate system; values are "
    "re-expressed the other way")
    .def("get_source_coordinate_system",
         [](const CoordinateMap& self) {
           return ToPyHolder(self.GetSourceCoordinateSystem());
         })
    .def("get_target_coordinate_system",
         [](const CoordinateMap& self) {
           return ToPyHolder(self.GetTargetCoordinateSystem());
         })
    .def(
      "bind",
      [](const CoordinateMap& self, const PyCoordinateView& query_points) {
        return self.Bind(CopyCoordinatesToDevice(query_points).view);
      },
      py::arg("query_points"),
      "Do every per-point computation for points in the source coordinate "
      "system");

  py::class_<CartesianToCylindrical, CoordinateMap,
             std::shared_ptr<CartesianToCylindrical>>(m,
                                                      "CartesianToCylindrical")
    .def(py::init<>());
  py::class_<CylindricalToCartesian, CoordinateMap,
             std::shared_ptr<CylindricalToCartesian>>(m,
                                                      "CylindricalToCartesian")
    .def(py::init<>());

  py::class_<TransformedSource, PointEvaluatorFactory,
             std::shared_ptr<TransformedSource>>(
    m, "TransformedSource",
    "Presents a PointEvaluatorFactory as one in another coordinate system")
    .def(py::init([](std::shared_ptr<PointEvaluatorFactory> source,
                     std::shared_ptr<CoordinateMap> to_source) {
           return std::make_shared<TransformedSource>(std::move(source),
                                                      std::move(to_source));
         }),
         py::arg("source"), py::arg("to_source"),
         "`to_source` maps this object's coordinate system to the source's");
}

} // namespace pcms
