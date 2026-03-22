#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pcms/coordinate_system.h"
#include "pcms/coordinate.h"
#include "pcms/adapter/uniform_grid/uniform_grid_field.h"
#include "pcms/adapter/uniform_grid/uniform_grid_field_layout.h"
#include "pcms/adapter/uniform_grid/uniform_grid_binary_field.h"
#include "pcms/lagrange_field_factory.h"
#include "pcms/utility/uniform_grid.h"
#include "numpy_array_transform.h"

namespace py = pybind11;

namespace pcms
{

void bind_coordinate_system_module(py::module& m)
{
  // Bind CoordinateSystem enum
  py::enum_<CoordinateSystem>(m, "CoordinateSystem")
    .value("Cartesian", CoordinateSystem::Cartesian)
    .value("Cylindrical", CoordinateSystem::Cylindrical)
    .value("XGC", CoordinateSystem::XGC)
    .value("GNET", CoordinateSystem::GNET)
    .value("BEAMS3D", CoordinateSystem::BEAMS3D)
    .export_values();

  // Bind CoordinateView for HostMemorySpace
  py::class_<CoordinateView<HostMemorySpace>>(m, "CoordinateView")
    .def(py::init([](CoordinateSystem cs, py::array_t<Real> coords) {
           // Convert 2D numpy array to Rank2View
           auto buf = coords.request();
           if (buf.ndim != 2) {
             throw std::runtime_error("Coordinates must be a 2D array");
           }
           // Create a view from the numpy array
           Rank2View<const Real, HostMemorySpace> coords_view(
             static_cast<Real*>(buf.ptr), buf.shape[0], buf.shape[1]);
           return CoordinateView<HostMemorySpace>(cs, coords_view);
         }),
         py::arg("coordinate_system"), py::arg("coordinates"),
         "Constructor for CoordinateView")

    .def("get_coordinate_system",
         &CoordinateView<HostMemorySpace>::GetCoordinateSystem,
         "Get the coordinate system")

    .def("set_coordinate_system",
         &CoordinateView<HostMemorySpace>::SetCoordinateSystem, py::arg("cs"),
         "Set the coordinate system")

    .def(
      "get_coordinates",
      [](const CoordinateView<HostMemorySpace>& self) {
        auto coords = self.GetCoordinates();
        // Convert to numpy array
        py::array_t<Real> result({static_cast<py::ssize_t>(coords.extent(0)),
                                  static_cast<py::ssize_t>(coords.extent(1))});
        auto buf = result.request();
        Real* ptr = static_cast<Real*>(buf.ptr);
        for (size_t i = 0; i < coords.extent(0); ++i) {
          for (size_t j = 0; j < coords.extent(1); ++j) {
            ptr[i * coords.extent(1) + j] = coords(i, j);
          }
        }
        return result;
      },
      "Get the coordinates as numpy array");

  // Bind CoordinateTransformation abstract base class
  py::class_<CoordinateTransformation,
             std::shared_ptr<CoordinateTransformation>>(
    m, "CoordinateTransformation")
    .def("evaluate", &CoordinateTransformation::Evaluate, py::arg("from"),
         py::arg("to"), "Evaluate the coordinate transformation");
}

void bind_coordinate_module(py::module& m)
{
  // Bind Coordinate template for common cases
  // 3D Cartesian coordinates with Real type
  py::class_<Coordinate<CoordinateSystem, Real, 3>>(m, "Coordinate3D")
    .def(py::init<Real, Real, Real>(), py::arg("x"), py::arg("y"), py::arg("z"),
         "Constructor for 3D coordinate")

    .def(
      "values",
      [](const Coordinate<CoordinateSystem, Real, 3>& self) {
        auto vals = self.Values();
        py::array_t<Real> result(3);
        auto buf = result.request();
        Real* ptr = static_cast<Real*>(buf.ptr);
        ptr[0] = vals[0];
        ptr[1] = vals[1];
        ptr[2] = vals[2];
        return result;
      },
      "Get coordinate values as numpy array")

    .def("__getitem__", &Coordinate<CoordinateSystem, Real, 3>::operator[],
         py::arg("index"), "Get coordinate value by index");

  // 2D coordinates
  py::class_<Coordinate<CoordinateSystem, Real, 2>>(m, "Coordinate2D")
    .def(py::init<Real, Real>(), py::arg("x"), py::arg("y"),
         "Constructor for 2D coordinate")

    .def(
      "values",
      [](const Coordinate<CoordinateSystem, Real, 2>& self) {
        auto vals = self.Values();
        py::array_t<Real> result(2);
        auto buf = result.request();
        Real* ptr = static_cast<Real*>(buf.ptr);
        ptr[0] = vals[0];
        ptr[1] = vals[1];
        return result;
      },
      "Get coordinate values as numpy array")

    .def("__getitem__", &Coordinate<CoordinateSystem, Real, 2>::operator[],
         py::arg("index"), "Get coordinate value by index");

  // Bind CoordinateElement for common types
  py::class_<CoordinateElement<CoordinateSystem, Real>>(m, "CoordinateElement")
    .def(py::init<Real>(), py::arg("data"), "Constructor for CoordinateElement")

    .def("underlying",
         py::overload_cast<>(
           &CoordinateElement<CoordinateSystem, Real>::underlying, py::const_),
         "Get the underlying value");
}

void bind_create_field_module(py::module& m)
{
  // Bind LagrangeFieldFactory
  py::class_<LagrangeFieldFactory>(m, "LagrangeFieldFactory")
    .def_static(
      "from_mesh",
      [](Omega_h::Mesh& mesh, int order, int num_components,
         CoordinateSystem coordinate_system) {
        return LagrangeFieldFactory::FromMesh(mesh, order, num_components,
                                             coordinate_system);
      },
      py::arg("mesh"), py::arg("order"), py::arg("num_components") = 1,
      py::arg("coordinate_system") = CoordinateSystem::Cartesian,
      "Create a LagrangeFieldFactory from an Omega_h mesh")

    .def_static(
      "from_uniform_grid",
      [](const UniformGrid<2>& grid, int num_components, CoordinateSystem cs,
         int order) {
        auto el = grid.edge_length;
        auto bl = grid.bot_left;
        auto div = grid.divisions;
        Rank1View<Real, HostMemorySpace> el_view(el.data(), 2);
        Rank1View<Real, HostMemorySpace> bl_view(bl.data(), 2);
        Rank1View<LO, HostMemorySpace> div_view(div.data(), 2);
        return LagrangeFieldFactory::FromUniformGrid(el_view, bl_view, div_view,
                                                     num_components, cs, order);
      },
      py::arg("grid"), py::arg("num_components") = 1,
      py::arg("coordinate_system") = CoordinateSystem::Cartesian,
      py::arg("order") = 1,
      "Create a LagrangeFieldFactory from a 2D uniform grid")

    .def_static(
      "from_uniform_grid",
      [](const UniformGrid<3>& grid, int num_components, CoordinateSystem cs,
         int order) {
        auto el = grid.edge_length;
        auto bl = grid.bot_left;
        auto div = grid.divisions;
        Rank1View<Real, HostMemorySpace> el_view(el.data(), 3);
        Rank1View<Real, HostMemorySpace> bl_view(bl.data(), 3);
        Rank1View<LO, HostMemorySpace> div_view(div.data(), 3);
        return LagrangeFieldFactory::FromUniformGrid(el_view, bl_view, div_view,
                                                     num_components, cs, order);
      },
      py::arg("grid"), py::arg("num_components") = 1,
      py::arg("coordinate_system") = CoordinateSystem::Cartesian,
      py::arg("order") = 1,
      "Create a LagrangeFieldFactory from a 3D uniform grid")

    .def("get_layout", &LagrangeFieldFactory::GetLayout,
         "Get the field layout (shared_ptr)")

    .def(
      "create_field_real",
      [](const LagrangeFieldFactory& self) {
        return std::shared_ptr<FieldT<Real>>(self.CreateFieldReal());
      },
      "Create a real-valued field");

  // Bind CreateUniformGridFromMesh for 2D
  m.def(
    "create_uniform_grid_from_mesh",
    [](Omega_h::Mesh& mesh, const std::array<LO, 2>& divisions) {
      return CreateUniformGridFromMesh<2>(mesh, divisions);
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 2D uniform grid from an Omega_h mesh");

  // Bind CreateUniformGridFromMesh for 3D
  m.def(
    "create_uniform_grid_from_mesh",
    [](Omega_h::Mesh& mesh, const std::array<LO, 3>& divisions) {
      return CreateUniformGridFromMesh<3>(mesh, divisions);
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 3D uniform grid from an Omega_h mesh");

  m.def(
    "create_uniform_grid_binary_field",
    [](Omega_h::Mesh& mesh, const std::array<LO, 2>& divisions) {
      auto [layout, field] =
        CreateUniformGridBinaryField<2>(mesh, divisions);
      return py::make_tuple(layout,
                            std::shared_ptr<FieldT<Real>>(std::move(field)));
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 2D vertex mask field indicating inside/outside mesh");

}

template <typename T>
void bind_field_t(py::module& m, const std::string& type_suffix)
{
  std::string class_name = "FieldT_" + type_suffix;

  py::class_<FieldT<T>, std::shared_ptr<FieldT<T>>>(m, class_name.c_str())
    .def("get_coordinate_system", &FieldT<T>::GetCoordinateSystem,
         "Get the coordinate system of the field")

    .def("get_localization_hint", &FieldT<T>::GetLocalizationHint,
         py::arg("coordinates"),
         "Get a localization hint for a set of coordinates")

    .def(
      "evaluate",
      [](const FieldT<T>& self, LocalizationHint hint,
         py::array_t<T> results_array, CoordinateSystem coord_sys) {
        auto results_view = numpy_to_view<T>(results_array);
        FieldDataView<T, HostMemorySpace> results(results_view, coord_sys);
        self.Evaluate(hint, results);
      },
      py::arg("hint"), py::arg("results"), py::arg("coordinate_system"),
      "Evaluate the field at given locations")

    .def(
      "evaluate_gradient",
      [](FieldT<T>& self, py::array_t<T> results_array,
         CoordinateSystem coord_sys) {
        auto results_view = numpy_to_view<T>(results_array);
        FieldDataView<T, HostMemorySpace> results(results_view, coord_sys);
        self.EvaluateGradient(results);
      },
      py::arg("results"), py::arg("coordinate_system"),
      "Evaluate the gradient of the field")

    .def(
      "get_dof_holder_data",
      [](const FieldT<T>& self) {
        auto data = self.GetDOFHolderData();
        return view_to_numpy(data);
      },
      "Get the DOF holder data")

    .def(
      "set_dof_holder_data",
      [](FieldT<T>& self, py::array_t<const T> data) {
        auto data_view = numpy_to_view<const T>(data);
        self.SetDOFHolderData(data_view);
      },
      py::arg("data"), "Set the DOF holder data")

    .def("get_layout", &FieldT<T>::GetLayout,
         py::return_value_policy::reference, "Get the field layout")

    .def("can_evaluate_gradient", &FieldT<T>::CanEvaluateGradient,
         "Check if the field can evaluate gradients")

    .def(
      "serialize",
      [](const FieldT<T>& self, py::array_t<T> buffer,
         py::array_t<const LO> permutation) {
        auto buffer_view = numpy_to_view<T>(buffer);
        auto perm_view = numpy_to_view<const LO>(permutation);
        return self.Serialize(buffer_view, perm_view);
      },
      py::arg("buffer"), py::arg("permutation"), "Serialize the field data")

    .def(
      "deserialize",
      [](FieldT<T>& self, py::array_t<const T> buffer,
         py::array_t<const LO> permutation) {
        auto buffer_view = numpy_to_view<const T>(buffer);
        auto perm_view = numpy_to_view<const LO>(permutation);
        self.Deserialize(buffer_view, perm_view);
      },
      py::arg("buffer"), py::arg("permutation"),
      "Deserialize field data from buffer");
}

void bind_field_module(py::module& m)
{
  // Bind LocalizationHint (opaque type - data member is internal only)
  py::class_<LocalizationHint>(m, "LocalizationHint").def(py::init<>());

  // Bind FieldDataView for common types
  py::class_<FieldDataView<Real, HostMemorySpace>>(m, "FieldDataView_Real")
    .def(py::init<Rank1View<Real, HostMemorySpace>, CoordinateSystem>(),
         py::arg("values"), py::arg("coordinate_system"))
    .def("size", &FieldDataView<Real, HostMemorySpace>::Size)
    .def("get_coordinate_system",
         &FieldDataView<Real, HostMemorySpace>::GetCoordinateSystem)
    // TODO: const GetValues?
    .def("get_values", [](FieldDataView<Real, HostMemorySpace>& self) {
      return view_to_numpy(self.GetValues());
    });

  // Bind FieldT only for double (Real is defined as double in types.h)
  bind_field_t<double>(m, "Double");
}

} // namespace pcms
