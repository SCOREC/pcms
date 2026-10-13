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
#include "pcms/discretization/discretization/omega_h.hpp"
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
  py::class_<CoordinateView<HostMemorySpace>,
             std::shared_ptr<CoordinateView<HostMemorySpace>>>(m,
                                                               "CoordinateView")
    .def(py::init([](CoordinateSystem cs, const NumpyInput<Real>& coords) {
           // CoordinateView does not own its coordinates, so keep a copy
           // alive in the holder's deleter.
           auto storage = CopyNumpyToKokkosView<2, Real, HostMemorySpace>(
             coords, "CoordinateView");
           return std::shared_ptr<CoordinateView<HostMemorySpace>>(
             new CoordinateView<HostMemorySpace>(cs, MakeRank2View(storage)),
             [storage](CoordinateView<HostMemorySpace>* view) { delete view; });
         }),
         py::arg("coordinate_system"), py::arg("coordinates"),
         "Create a CoordinateView from a 2D array (num_points × dim). The "
         "coordinates are copied.")

    .def("get_coordinate_system",
         &CoordinateView<HostMemorySpace>::GetCoordinateSystem,
         "Get the coordinate system")

    .def("set_coordinate_system",
         &CoordinateView<HostMemorySpace>::SetCoordinateSystem, py::arg("cs"),
         "Set the coordinate system")

    .def(
      "get_coordinates",
      [](py::object self) {
        return ViewToNumpy(
          self.cast<const CoordinateView<HostMemorySpace>&>().GetValues(),
          self);
      },
      "The coordinates as a read-only 2D numpy array (num_points × dim) that "
      "shares memory with this CoordinateView.");

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
      [](const NumpyInput<Real>& coords, CoordinateSystem coordinate_system,
         OutOfBoundsPolicy policy) {
        auto coords_device =
          CopyNumpyToKokkosView<2, Real, DeviceMemorySpace>(coords, "coords");
        return PythonEvaluationRequest{
          EvaluationRequest::FromCoordinates(
            CoordinateView<DeviceMemorySpace>(coordinate_system,
                                              MakeRank2View(coords_device)),
            policy),
          coords_device};
      },
      py::arg("coordinates"),
      py::arg("coordinate_system") = CoordinateSystem::Cartesian,
      py::arg("policy") = OutOfBoundsPolicy{},
      "Create an EvaluationRequest from a 2D coordinate array (num_points × "
      "dim). The coordinates are copied.")
    .def_static(
      "from_function_space",
      [](const FunctionSpace& function_space, OutOfBoundsPolicy policy) {
        return PythonEvaluationRequest{
          EvaluationRequest::FromFunctionSpace(function_space, policy),
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
      [](py::object self) {
        return ViewToNumpy(self.cast<const Field<Real>&>().GetDOFHolderData(),
                           self);
      },
      "The DOF holder data as a read-only 2D numpy array (num_dof_holders × "
      "num_components). When the field's data is in host memory the array "
      "shares that memory, so it reflects later set_dof_holder_data calls; "
      "otherwise it is a copy.")

    .def(
      "set_dof_holder_data",
      [](Field<Real>& self, const NumpyInput<Real>& arr) {
        auto data = NumpyToHostView<2, Real>(arr);
        const LO num_dof_holders = self.GetLayout().GetNumOwnedDofHolder();
        const LO num_components = self.GetLayout().GetNumComponents();
        if (data.extent(0) != num_dof_holders ||
            data.extent(1) != num_components) {
          throw std::invalid_argument("DOF holder data must have shape (" +
                                      std::to_string(num_dof_holders) + ", " +
                                      std::to_string(num_components) +
                                      "), got (" +
                                      std::to_string(data.extent(0)) + ", " +
                                      std::to_string(data.extent(1)) + ")");
        }
        self.SetDOFHolderDataHost(data);
      },
      py::arg("data"),
      "Set the DOF holder data from a 2D numpy array (num_dof_holders × "
      "num_components)")

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
      [](py::object self) {
        return ViewToNumpy(self.cast<const Field<Real>&>()
                             .GetLayout()
                             .GetDOFHolderCoordinates()
                             .GetValues(),
                           self);
      },
      "Read-only DOFHolder coords (num_dof_holders × dim).");

  py::class_<FunctionSpace, std::shared_ptr<FunctionSpace>>(m, "FunctionSpace")
    .def(
      "create_field",
      [](const FunctionSpace& self) -> Field<Real> {
        return self.CreateFunction<Real>();
      },
      "Create a Field<Real> for this function space.")
    .def(
      "create_point_evaluator",
      [](const FunctionSpace& self, const PythonEvaluationRequest& request) {
        return self.CreatePointEvaluator<Real>(request.request);
      },
      py::arg("request"),
      "Create a reusable point evaluator from an EvaluationRequest.")
    .def("get_coordinate_system", &FunctionSpace::GetCoordinateSystem,
         "Get the coordinate system for this function space")
    .def(
      "mesh",
      [](const FunctionSpace& self) -> Omega_h::Mesh& {
        auto disc = std::dynamic_pointer_cast<const OmegaHDiscretization>(
          self.GetDiscretization());
        if (!disc) {
          throw std::runtime_error(
            "FunctionSpace::mesh: this function space is not backed by an "
            "Omega_h mesh");
        }
        return disc->GetMesh();
      },
      py::return_value_policy::reference,
      "The Omega_h mesh this function space is defined on (Omega_h backend "
      "only).");

  // Bind LagrangeFunctionSpace as a concrete FunctionSpace subtype.
  py::class_<LagrangeFunctionSpace, FunctionSpace,
             std::shared_ptr<LagrangeFunctionSpace>>
    lagrange_space(m, "LagrangeFunctionSpace");

  // Backend selecting the underlying field layout. The conservative-projection
  // transfer operators require the native Omega_h backend.
  py::enum_<LagrangeFunctionSpace::Backend>(lagrange_space, "Backend")
    .value("MeshFields", LagrangeFunctionSpace::Backend::MeshFields,
           "MeshFields-backed layout (default when MeshFields is enabled)")
    .value("OmegaH", LagrangeFunctionSpace::Backend::OmegaH,
           "Native Omega_h Lagrange layout (required by the conservative and "
           "Monte Carlo projection transfer operators)")
    .export_values();

  lagrange_space
    .def_static(
      "from_mesh",
      [](Omega_h::Mesh& mesh, int order, int num_components,
         CoordinateSystem coordinate_system, std::string global_id_name,
         LagrangeFunctionSpace::Backend backend) {
        return LagrangeFunctionSpace::FromMesh(
          mesh, order, num_components, coordinate_system,
          std::move(global_id_name), backend);
      },
      py::arg("mesh"), py::arg("order"), py::arg("num_components") = 1,
      py::arg("coordinate_system") = CoordinateSystem::Cartesian,
      py::arg("global_id_name") = "global",
      py::arg("backend") = LagrangeFunctionSpace::DefaultBackend,
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
  py::class_<PolynomialReconstructionFunctionSpace, FunctionSpace,
             std::shared_ptr<PolynomialReconstructionFunctionSpace>>(
    m, "PolynomialReconstructionFunctionSpace")
    .def_static(
      "from_coords",
      [](const NumpyInput<Real>& coords, CoordinateSystem cs, MLSOptions opts) {
        return PolynomialReconstructionFunctionSpace::Create(
          NumpyToHostView<2, Real>(coords), cs, opts);
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
