#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pcms/field/coordinate_system.h"
#include "pcms/field/coordinate.h"
#include "pcms/field/field.h"
#include "pcms/field/field_evaluator_factory.h"
#include "pcms/field/field_layout.h"
#include "pcms/field/field_metadata.h"
#include "pcms/field/function_space.h"
#include "pcms/field/data/simple.h"
#include "pcms/field/layout/uniform_grid.h"
#include "pcms/field/uniform_grid_binary_field.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/function_space/polynomial_reconstruction.hpp"
#include "pcms/field/evaluator/mls_options.h"
#include "pcms/utility/uniform_grid.h"
#include "pcms/utility/arrays.h"
#include "numpy_array_transform.h"

namespace py = pybind11;

namespace pcms
{

namespace
{

struct PythonEvaluationRequest
{
  EvaluationRequest request;
  py::object owner;
  // Keep the device view alive to prevent dangling pointer in mdspan
  Kokkos::View<Real**, DeviceMemorySpace> device_coords;
};

} // namespace

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
           using LayoutPolicy =
             detail::default_layout_for_memory_space_t<HostMemorySpace>;
           Rank2View<const Real, HostMemorySpace, LayoutPolicy> coords_view(
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
        auto coords = self.GetValues();
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
  py::class_<PythonEvaluationRequest>(m, "EvaluationRequest")
    .def_static(
      "from_coordinates",
      [](py::array_t<Real> coords, CoordinateSystem coordinate_system,
         OutOfBoundsPolicy policy) {
        auto coords_view = numpy_to_view_2d<const Real>(coords);
        // Create a Kokkos::View from the host data and deep copy to device
        auto coords_host = Kokkos::View<Real**, HostMemorySpace>(
          "coords_host", coords_view.extent(0), coords_view.extent(1));
        for (size_t i = 0; i < coords_view.extent(0); ++i) {
          for (size_t j = 0; j < coords_view.extent(1); ++j) {
            coords_host(i, j) = coords_view(i, j);
          }
        }
        auto coords_device = Kokkos::View<Real**, DeviceMemorySpace>(
          "coords_device", coords_view.extent(0), coords_view.extent(1));
        DeepCopyMismatchLayouts(coords_device, coords_host);
        auto coords_device_view = MakeRank2View(coords_device);
        return PythonEvaluationRequest{
          EvaluationRequest::FromCoordinates(
            CoordinateView<DeviceMemorySpace>(coordinate_system,
                                              coords_device_view),
            policy),
          py::reinterpret_borrow<py::object>(coords),
          coords_device}; // Keep the View alive!
      },
      py::arg("coordinates"),
      py::arg("coordinate_system") = CoordinateSystem::Cartesian,
      py::arg("policy") = OutOfBoundsPolicy{},
      "Create an EvaluationRequest from an explicit coordinate array.")
    .def_static(
      "from_function_space",
      [](const FunctionSpace& function_space, OutOfBoundsPolicy policy) {
        return PythonEvaluationRequest{
          EvaluationRequest::FromFunctionSpace(function_space, policy),
          py::none(),
          Kokkos::View<Real**, DeviceMemorySpace>()}; // Empty view for
                                                      // function_space case
      },
      py::arg("function_space"), py::arg("policy") = OutOfBoundsPolicy{},
      "Create an EvaluationRequest from a FunctionSpace's DOF-holder sites.");

  // Bind Field<Real>: composed per-field object returned by
  // FunctionSpace-backed factories' create_field(). Move-only in C++; Python
  // holds it by value in a heap-allocated wrapper.
  py::class_<Field<Real>>(m, "Field")
    .def(
      "get_dof_holder_data",
      [](const Field<Real>& self) {
        auto data = self.GetDOFHolderDataHost();
        // Flatten the [dof][comp] data into node-major order for the 1D numpy
        // array (Python-facing DOF data stays 1D).
        const auto num_dof = data.extent(0);
        const auto num_comp = data.extent(1);
        py::array_t<Real> result(static_cast<py::ssize_t>(data.size()));
        auto buf = result.request();
        Real* ptr = static_cast<Real*>(buf.ptr);
        for (size_t i = 0; i < num_dof; ++i)
          for (size_t c = 0; c < num_comp; ++c)
            ptr[i * num_comp + c] = data(i, c);
        return result;
      },
      "Get the DOF holder data as a 1D numpy array")

    .def(
      "set_dof_holder_data",
      [](Field<Real>& self, py::array_t<const Real> arr) {
        auto buf = arr.request();
        if (buf.ndim != 1) {
          throw std::runtime_error("DOF holder data must be a 1D array");
        }
        // Reshape the flat 1D input into the field's [dof][comp] layout.
        const int nc = self.GetLayout().GetNumComponents();
        const auto total = static_cast<size_t>(buf.shape[0]);
        self.SetDOFHolderDataHost(Rank2View<const Real, HostMemorySpace>(
          static_cast<const Real*>(buf.ptr), static_cast<LO>(total / nc), nc));
      },
      py::arg("data"), "Set the DOF holder data from a 1D numpy array")

    .def(
      "get_num_dof_holders",
      [](const Field<Real>& self) {
        return self.GetLayout().GetNumOwnedDofHolder();
      },
      "Number of owned DOF holders (nodes/elements)")

    .def(
      "get_num_components",
      [](const Field<Real>& self) {
        return self.GetLayout().GetNumComponents();
      },
      "Number of field components per DOF holder")

    .def(
      "get_dof_holder_coordinates",
      [](const Field<Real>& self) {
        auto cv = self.GetLayout().GetDOFHolderCoordinates();
        auto coords = cv.GetValues();
        Kokkos::View<Real**, DeviceMemorySpace> coords_device(
          "coords_device", coords.extent(0), coords.extent(1));
        Kokkos::parallel_for(
          Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
            0, coords.extent(0)),
          KOKKOS_LAMBDA(size_t i) {
            for (size_t j = 0; j < coords.extent(1); ++j) {
              coords_device(i, j) = coords(i, j);
            }
          });
        Kokkos::View<Real**, HostMemorySpace> coords_host(
          "coords_host", coords.extent(0), coords.extent(1));
        DeepCopyMismatchLayouts(coords_host, coords_device);
        py::array_t<Real> result(
          {static_cast<py::ssize_t>(coords_host.extent(0)),
           static_cast<py::ssize_t>(coords_host.extent(1))});
        auto buf = result.request();
        Real* ptr = static_cast<Real*>(buf.ptr);
        for (size_t i = 0; i < coords_host.extent(0); ++i)
          for (size_t j = 0; j < coords_host.extent(1); ++j)
            ptr[i * coords_host.extent(1) + j] = coords_host(i, j);
        return result;
      },
      "DOF holder coordinates as a 2D numpy array (num_dof_holders × dim)");

  py::class_<FunctionSpace>(m, "FunctionSpace")
    .def(
      "create_field",
      [](const FunctionSpace& self) { return self.CreateField<Real>(); },
      "Create a Field<Real> for this function space.")
    .def(
      "create_point_evaluator",
      [](const FunctionSpace& self, const PythonEvaluationRequest& request) {
        return self.CreatePointEvaluator<Real>(request.request);
      },
      py::arg("request"),
      "Create a reusable point evaluator from an EvaluationRequest.")
    .def("get_coordinate_system", &FunctionSpace::GetCoordinateSystem,
         "Get the coordinate system for this function space");

  // Bind LagrangeFunctionSpace as a concrete FunctionSpace subtype.
  py::class_<LagrangeFunctionSpace, FunctionSpace>(m, "LagrangeFunctionSpace")
    .def_static(
      "from_mesh",
      [](Omega_h::Mesh& mesh, int order, int num_components,
         CoordinateSystem coordinate_system) {
        return LagrangeFunctionSpace::FromMesh(mesh, order, num_components,
                                               coordinate_system);
      },
      py::arg("mesh"), py::arg("order"), py::arg("num_components") = 1,
      py::arg("coordinate_system") = CoordinateSystem::Cartesian,
      "Create a LagrangeFunctionSpace from an Omega_h mesh")

    .def_static(
      "from_uniform_grid",
      [](const UniformGrid<2>& grid, int num_components, CoordinateSystem cs,
         int order) {
        return LagrangeFunctionSpace::FromUniformGrid(grid, num_components, cs,
                                                      order);
      },
      py::arg("grid"), py::arg("num_components") = 1,
      py::arg("coordinate_system") = CoordinateSystem::Cartesian,
      py::arg("order") = 1,
      "Create a LagrangeFunctionSpace from a 2D uniform grid")

    .def_static(
      "from_uniform_grid",
      [](const UniformGrid<3>& grid, int num_components, CoordinateSystem cs,
         int order) {
        return LagrangeFunctionSpace::FromUniformGrid(grid, num_components, cs,
                                                      order);
      },
      py::arg("grid"), py::arg("num_components") = 1,
      py::arg("coordinate_system") = CoordinateSystem::Cartesian,
      py::arg("order") = 1,
      "Create a LagrangeFunctionSpace from a 3D uniform grid");

  // Bind MLSOptions: configuration struct for
  // PolynomialReconstructionFunctionSpace MLS
  // evaluation.
  py::class_<MLSOptions>(m, "MLSOptions")
    .def(py::init<>(), "Default MLSOptions")
    .def_readwrite("radius", &MLSOptions::radius)
    .def_readwrite("min_req_supports", &MLSOptions::min_req_supports)
    .def_readwrite("degree", &MLSOptions::degree)
    .def_readwrite("adapt_radius", &MLSOptions::adapt_radius)
    .def_readwrite("lambda_reg", &MLSOptions::lambda)
    .def_readwrite("tol", &MLSOptions::tol)
    .def_readwrite("decay_factor", &MLSOptions::decay_factor)
    .def_readwrite("basis", &MLSOptions::basis);

  // Bind PolynomialReconstructionFunctionSpace as a concrete FunctionSpace
  // subtype.
  py::class_<PolynomialReconstructionFunctionSpace, FunctionSpace>(
    m, "PolynomialReconstructionFunctionSpace")
    .def_static(
      "from_coords",
      [](py::array_t<Real> coords, CoordinateSystem cs, MLSOptions opts) {
        auto view = numpy_to_view_2d<Real>(coords);
        return PolynomialReconstructionFunctionSpace::Create(view, cs, opts);
      },
      py::arg("coords"),
      py::arg("coordinate_system") = CoordinateSystem::Cartesian,
      py::arg("options") = MLSOptions{},
      "Create a PolynomialReconstructionFunctionSpace from a 2D array of "
      "source coordinates (shape: num_points × dim).")
    .def_static(
      "from_mesh",
      [](Omega_h::Mesh& mesh, int source_entity_dim, CoordinateSystem cs,
         MLSOptions opts) {
        return PolynomialReconstructionFunctionSpace::FromMesh(
          mesh, source_entity_dim, cs, opts);
      },
      py::arg("mesh"), py::arg("source_entity_dim"),
      py::arg("coordinate_system") = CoordinateSystem::Cartesian,
      py::arg("options") = MLSOptions{},
      "Create a PolynomialReconstructionFunctionSpace from mesh entity "
      "coordinates.");

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
      auto [layout, field] = CreateUniformGridBinaryField<2>(mesh, divisions);
      static_cast<void>(layout);
      return field;
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 2D vertex mask field indicating inside/outside mesh");
}

// bind_field_module is kept for compatibility but now registers nothing that
// refers to the deleted FieldT / LocalizationHint / FieldDataView types.
void bind_field_module(py::module& /*m*/) {}

} // namespace pcms
