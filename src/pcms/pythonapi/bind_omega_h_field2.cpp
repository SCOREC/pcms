#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include "pcms/utility/arrays.h"
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include "pcms/adapter/omega_h/omega_h_field_layout.h"
#include "numpy_array_transform.h"

namespace py = pybind11;

namespace pcms
{

void bind_omega_h_field2(py::module& m)
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

  // Bind OmegaHField2 class
  py::class_<OmegaHField2, FieldT<Real>, std::shared_ptr<OmegaHField2>>(
    m, "OmegaHField2")
    .def(py::init<const OmegaHFieldLayout&>(), py::arg("layout"),
         "Constructor for OmegaHField2")

    .def(
      "get_localization_hint",
      [](const OmegaHField2& self, py::array_t<Real> coordinates,
         const CoordinateSystem& coord_system) {
        // Create CoordinateView from numpy array
        auto coords_view = numpy_to_view_2d<const Real>(coordinates);
        CoordinateView<HostMemorySpace> coord_view(coord_system, coords_view);
        return self.GetLocalizationHint(coord_view);
      },
      py::arg("coordinates"), py::arg("coordinate_system"),
      "Get localization hint for given coordinates")

    .def(
      "evaluate",
      [](const OmegaHField2& self, const LocalizationHint& location,
         py::array_t<Real> values, const CoordinateSystem& coord_system) {
        // Create FieldDataView
        auto values_view = numpy_to_view<Real>(values);
        FieldDataView<Real, HostMemorySpace> field_data_view(values_view,
                                                             coord_system);
        self.Evaluate(location, field_data_view);
      },
      py::arg("location"), py::arg("values"), py::arg("coordinate_system"),
      "Evaluate field at locations specified by localization hint")

    .def(
      "evaluate_gradient",
      [](OmegaHField2& self, py::array_t<Real> gradients,
         const CoordinateSystem& coord_system) {
        auto gradients_view = numpy_to_view<Real>(gradients);
        FieldDataView<Real, HostMemorySpace> field_data_view(gradients_view,
                                                             coord_system);
        self.EvaluateGradient(field_data_view);
      },
      py::arg("gradients"), py::arg("coordinate_system"),
      "Evaluate gradient of the field")

    .def("get_layout", &OmegaHField2::GetLayout,
         py::return_value_policy::reference, "Get the field layout")

    .def("can_evaluate_gradient", &OmegaHField2::CanEvaluateGradient,
         "Check if gradient evaluation is supported")

    .def(
      "serialize",
      [](const OmegaHField2& self, py::array_t<Real> buffer,
         py::array_t<const pcms::LO> permutation) {
        auto buffer_view = numpy_to_view<Real>(buffer);
        auto perm_view = numpy_to_view<const pcms::LO>(permutation);
        return self.Serialize(buffer_view, perm_view);
      },
      py::arg("buffer"), py::arg("permutation"),
      "Serialize field data into buffer")

    .def(
      "deserialize",
      [](OmegaHField2& self, py::array_t<const Real> buffer,
         py::array_t<const pcms::LO> permutation) {
        auto buffer_view = numpy_to_view<const Real>(buffer);
        auto perm_view = numpy_to_view<const pcms::LO>(permutation);
        self.Deserialize(buffer_view, perm_view);
      },
      py::arg("buffer"), py::arg("permutation"),
      "Deserialize field data from buffer")

    .def(
      "get_dof_holder_data",
      [](const OmegaHField2& self) {
        auto const_data = self.GetDOFHolderData();
        // Create a numpy array that owns its own data
        return view_to_numpy<const Real>(const_data);
      },
      "Get the DOF holder data")

    .def(
      "set_dof_holder_data",
      [](OmegaHField2& self, py::array_t<Real> data) {
        // Ensure array is contiguous
        auto contiguous_data = py::array_t<Real>(data);
        auto data_view = numpy_to_view<Real>(contiguous_data);
        // Create const view wrapper
        Rank1View<const Real, HostMemorySpace> const_view(
          data_view.data_handle(), data_view.size());
        self.SetDOFHolderData(const_view);
      },
      py::arg("data"), "Set the DOF holder data")

    .def("set_out_of_bounds_mode", &OmegaHField2::SetOutOfBoundsMode,
         py::arg("mode"), py::arg("fill_value") = 0.0,
         "Set the out-of-bounds mode and fill value")

    .def("get_out_of_bounds_mode", &OmegaHField2::GetOutOfBoundsMode,
         "Get the current out-of-bounds mode")

    .def("get_fill_value", &OmegaHField2::GetFillValue,
         "Get the current fill value for out-of-bounds points");

  // Helper functions for creating views (if needed for testing)
  m.def(
    "create_coordinate_view",
    [](py::array_t<Real> coordinates, const CoordinateSystem& coord_system) {
      auto coords_view = numpy_to_view_2d<const Real>(coordinates);
      return CoordinateView<HostMemorySpace>(coord_system, coords_view);
    },
    py::arg("coordinates"), py::arg("coordinate_system"),
    "Create a CoordinateView from numpy array");

  m.def(
    "create_field_data_view",
    [](py::array_t<Real> values, const CoordinateSystem& coord_system) {
      auto values_view = numpy_to_view<Real>(values);
      return FieldDataView<Real, HostMemorySpace>(values_view, coord_system);
    },
    py::arg("values"), py::arg("coordinate_system"),
    "Create a FieldDataView from numpy array");
}

} // namespace pcms