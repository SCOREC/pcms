#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include "pcms/utility/arrays.h"
#include "pcms/adapter/uniform_grid/uniform_grid_field.h"
#include "pcms/adapter/uniform_grid/uniform_grid_field_layout.h"
#include "numpy_array_transform.h"

namespace py = pybind11;

namespace pcms
{

void bind_uniform_grid_field_module(py::module& m)
{

  // Bind UniformGridField class for 2D
  py::class_<UniformGridField<2>, FieldT<Real>,
             std::shared_ptr<UniformGridField<2>>>(m, "UniformGridField2D")
    .def(py::init<const UniformGridFieldLayout<2>&>(), py::arg("layout"),
         "Constructor for UniformGridField2D")

    .def(
      "get_localization_hint",
      [](const UniformGridField<2>& self, py::array_t<Real> coordinates,
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
      [](const UniformGridField<2>& self, const LocalizationHint& location,
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
      [](UniformGridField<2>& self, py::array_t<Real> gradients,
         const CoordinateSystem& coord_system) {
        auto gradients_view = numpy_to_view<Real>(gradients);
        FieldDataView<Real, HostMemorySpace> field_data_view(gradients_view,
                                                             coord_system);
        self.EvaluateGradient(field_data_view);
      },
      py::arg("gradients"), py::arg("coordinate_system"),
      "Evaluate gradient of the field")

    .def("get_layout", &UniformGridField<2>::GetLayout,
         py::return_value_policy::reference, "Get the field layout")

    .def("can_evaluate_gradient", &UniformGridField<2>::CanEvaluateGradient,
         "Check if gradient evaluation is supported")

    .def(
      "serialize",
      [](const UniformGridField<2>& self, py::array_t<Real> buffer,
         py::array_t<const pcms::LO> permutation) {
        auto buffer_view = numpy_to_view<Real>(buffer);
        auto perm_view = numpy_to_view<const pcms::LO>(permutation);
        return self.Serialize(buffer_view, perm_view);
      },
      py::arg("buffer"), py::arg("permutation"),
      "Serialize field data into buffer")

    .def(
      "deserialize",
      [](UniformGridField<2>& self, py::array_t<const Real> buffer,
         py::array_t<const pcms::LO> permutation) {
        auto buffer_view = numpy_to_view<const Real>(buffer);
        auto perm_view = numpy_to_view<const pcms::LO>(permutation);
        self.Deserialize(buffer_view, perm_view);
      },
      py::arg("buffer"), py::arg("permutation"),
      "Deserialize field data from buffer")

    .def(
      "get_dof_holder_data",
      [](const UniformGridField<2>& self) {
        auto const_data = self.GetDOFHolderData();
        // Create a numpy array that owns its own data
        return view_to_numpy<const Real>(const_data);
      },
      "Get the DOF holder data")

    .def(
      "set_dof_holder_data",
      [](UniformGridField<2>& self, py::array_t<Real> data) {
        // Ensure array is contiguous
        auto contiguous_data = py::array_t<Real>(data);
        auto data_view = numpy_to_view<Real>(contiguous_data);
        // Create const view wrapper
        Rank1View<const Real, HostMemorySpace> const_view(
          data_view.data_handle(), data_view.size());
        self.SetDOFHolderData(const_view);
      },
      py::arg("data"), "Set the DOF holder data")

    .def(
      "to_numpy",
      [](const UniformGridField<2>& self) {
        return mdspan_view_to_numpy(self.to_mdspan());
      },
      "Get the DOF holder data as a 2D numpy array (copy)")

    .def("set_out_of_bounds_mode", &UniformGridField<2>::SetOutOfBoundsMode,
         py::arg("mode"), py::arg("fill_value") = 0.0,
         "Set the out-of-bounds mode and fill value")

    .def("get_out_of_bounds_mode", &UniformGridField<2>::GetOutOfBoundsMode,
         "Get the current out-of-bounds mode")

    .def("get_fill_value", &UniformGridField<2>::GetFillValue,
         "Get the current fill value for out-of-bounds points");

  // Bind UniformGridField class for 3D
  py::class_<UniformGridField<3>, FieldT<Real>,
             std::shared_ptr<UniformGridField<3>>>(m, "UniformGridField3D")
    .def(py::init<const UniformGridFieldLayout<3>&>(), py::arg("layout"),
         "Constructor for UniformGridField3D")

    .def(
      "get_localization_hint",
      [](const UniformGridField<3>& self, py::array_t<Real> coordinates,
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
      [](const UniformGridField<3>& self, const LocalizationHint& location,
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
      [](UniformGridField<3>& self, py::array_t<Real> gradients,
         const CoordinateSystem& coord_system) {
        auto gradients_view = numpy_to_view<Real>(gradients);
        FieldDataView<Real, HostMemorySpace> field_data_view(gradients_view,
                                                             coord_system);
        self.EvaluateGradient(field_data_view);
      },
      py::arg("gradients"), py::arg("coordinate_system"),
      "Evaluate gradient of the field")

    .def("get_layout", &UniformGridField<3>::GetLayout,
         py::return_value_policy::reference, "Get the field layout")

    .def("can_evaluate_gradient", &UniformGridField<3>::CanEvaluateGradient,
         "Check if gradient evaluation is supported")

    .def(
      "serialize",
      [](const UniformGridField<3>& self, py::array_t<Real> buffer,
         py::array_t<const pcms::LO> permutation) {
        auto buffer_view = numpy_to_view<Real>(buffer);
        auto perm_view = numpy_to_view<const pcms::LO>(permutation);
        return self.Serialize(buffer_view, perm_view);
      },
      py::arg("buffer"), py::arg("permutation"),
      "Serialize field data into buffer")

    .def(
      "deserialize",
      [](UniformGridField<3>& self, py::array_t<const Real> buffer,
         py::array_t<const pcms::LO> permutation) {
        auto buffer_view = numpy_to_view<const Real>(buffer);
        auto perm_view = numpy_to_view<const pcms::LO>(permutation);
        self.Deserialize(buffer_view, perm_view);
      },
      py::arg("buffer"), py::arg("permutation"),
      "Deserialize field data from buffer")

    .def(
      "get_dof_holder_data",
      [](const UniformGridField<3>& self) {
        auto const_data = self.GetDOFHolderData();
        // Create a numpy array that owns its own data
        return view_to_numpy<const Real>(const_data);
      },
      "Get the DOF holder data")

    .def(
      "set_dof_holder_data",
      [](UniformGridField<3>& self, py::array_t<Real> data) {
        // Ensure array is contiguous
        auto contiguous_data = py::array_t<Real>(data);
        auto data_view = numpy_to_view<Real>(contiguous_data);
        // Create const view wrapper
        Rank1View<const Real, HostMemorySpace> const_view(
          data_view.data_handle(), data_view.size());
        self.SetDOFHolderData(const_view);
      },
      py::arg("data"), "Set the DOF holder data")

    .def(
      "to_numpy",
      [](const UniformGridField<3>& self) {
        return mdspan_view_to_numpy(self.to_mdspan());
      },
      "Get the DOF holder data as a 3D numpy array (copy)")

    .def("set_out_of_bounds_mode", &UniformGridField<3>::SetOutOfBoundsMode,
         py::arg("mode"), py::arg("fill_value") = 0.0,
         "Set the out-of-bounds mode and fill value")

    .def("get_out_of_bounds_mode", &UniformGridField<3>::GetOutOfBoundsMode,
         "Get the current out-of-bounds mode")

    .def("get_fill_value", &UniformGridField<3>::GetFillValue,
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
