#include <optional>
#include <stdexcept>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pcms/discretization/discretization/omega_h.hpp"
#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/evaluation_request.h"
#include "pcms/field/evaluator/mls_options.h"
#include "pcms/field/field.h"
#include "pcms/field/field_data.h"
#include "pcms/field/field_layout.h"
#include "pcms/field/function_space.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/function_space/polynomial_reconstruction.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/uniform_grid_binary_field.h"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/uniform_grid.h"
#include "numpy_array_transform.h"
#include "py_value_types.h"

namespace py = pybind11;

namespace pcms
{

namespace
{

struct PyEvaluationRequest
{
  EvaluationRequest request;
  Kokkos::View<Real**, DeviceMemorySpace> device_coords;
};

using SystemHolder = std::shared_ptr<CoordinateSystem>;

void SetDOFHolderData(Field<Real>& self, const PyValueView& v)
{
  auto array = NumpyInput<Real>::ensure(v.values);
  if (!array) {
    throw py::type_error("set_dof_holder_data: values must be float64");
  }
  self.SetDOFHolderDataHost(ValueView<const Real, HostMemorySpace>(
    v.basis, NumpyToHostView<2, Real>(array)));
}

} // namespace

void bind_create_field_module(py::module& m)
{
  py::class_<FieldLayout, std::shared_ptr<FieldLayout>>(
    m, "FieldLayout", "How a field's DOF holders are laid out")
    .def("get_num_components", &FieldLayout::GetNumComponents)
    .def("get_num_owned_dof_holder", &FieldLayout::GetNumOwnedDofHolder)
    .def("get_num_global_dof_holder", &FieldLayout::GetNumGlobalDofHolder)
    .def("get_dimension", &FieldLayout::GetDimension)
    .def("is_distributed", &FieldLayout::IsDistributed)
    .def(
      "get_coordinate_system",
      [](const FieldLayout& self) {
        return ToPyHolder(self.GetCoordinateSystem());
      },
      "What this layout's DOF-holder coordinates mean")
    .def(
      "get_dof_holder_coordinates",
      [](py::object self) {
        return ToPyCoordinateView(
          self.cast<const FieldLayout&>().GetDOFHolderCoordinates(), self);
      },
      "The DOF-holder coordinates (num_dof_holders × dim), read-only");

  py::class_<FieldData<Real>>(m, "FieldData", "A field's stored values")
    .def("get_value_type", &FieldData<Real>::GetValueType)
    .def("get_value_basis", &FieldData<Real>::GetValueBasis,
         "The basis the stored components are in");

  py::class_<Field<Real>>(
    m, "Field", "DOF-holder values together with the layout that orders them")
    .def("get_name", &Field<Real>::GetName)
    .def("get_layout", &Field<Real>::GetLayout,
         py::return_value_policy::reference_internal)
    .def(
      "get_data",
      [](Field<Real>& self) -> FieldData<Real>& { return self.GetData(); },
      py::return_value_policy::reference_internal)
    .def(
      "get_dof_holder_data",
      [](py::object self) {
        return ToPyValueView(self.cast<const Field<Real>&>().GetDOFHolderData(),
                             self);
      },
      "The DOF-holder values (num_dof_holders × num_components) as a "
      "ValueView whose array is read-only. When the field's data is in host "
      "memory the array shares it, so it reflects later writes; otherwise it "
      "is a copy.")
    .def("set_dof_holder_data", &SetDOFHolderData, py::arg("data"),
         "Set the DOF-holder values from a ValueView (num_dof_holders × "
         "num_components) whose basis must match the field's")
    .def(
      "set_dof_holder_data",
      [](Field<Real>& self, const NumpyInput<Real>& data) {
        RequireScalarForPlainArray("set_dof_holder_data",
                                   self.GetData().GetValueBasis());
        self.SetDOFHolderDataHost(ValueView<const Real, HostMemorySpace>(
          ValueBasis{},
          NumpyToHostComponentView<Real>(data, PerEntityComponents(data))));
      },
      py::arg("data"),
      "Set a scalar field's DOF-holder values from a plain 2D array "
      "(num_dof_holders × num_components). A 1D array (num_dof_holders,) is "
      "accepted when the field has one component.")
    .def(
      "set_dof_holder_data_unchecked",
      [](Field<Real>& self, const NumpyInput<Real>& data) {
        const auto& layout = self.GetLayout();
        const LO num_dof_holders = layout.GetNumOwnedDofHolder();
        auto view =
          NumpyToHostComponentView<Real>(data, layout.GetNumComponents());
        if (view.extent(0) != static_cast<size_t>(num_dof_holders)) {
          throw std::invalid_argument(
            "set_dof_holder_data_unchecked: data must have " +
            std::to_string(num_dof_holders) + " rows, got " +
            std::to_string(view.extent(0)));
        }
        self.SetDOFHolderDataUncheckedHost(view);
      },
      py::arg("data"),
      "Set the DOF-holder values without checking their basis (the shape is "
      "still checked). A 1D array (num_dof_holders,) is accepted when the "
      "field has one component.");

  py::class_<Function<Real>, Field<Real>>(
    m, "Function", "A Field that knows its FunctionSpace")
    .def("get_space", [](const Function<Real>& self) {
      return ToPyHolder(self.GetSpacePtr());
    });

  py::class_<PyEvaluationRequest>(m, "EvaluationRequest")
    .def_static(
      "from_coordinates",
      [](const PyCoordinateView& coords, OutOfBoundsPolicy policy) {
        auto device = CopyCoordinatesToDevice(coords);
        return PyEvaluationRequest{
          EvaluationRequest::FromCoordinates(device.view, policy),
          device.storage};
      },
      py::arg("coordinates"), py::arg("policy") = OutOfBoundsPolicy{},
      "Request evaluation at the given CoordinateView's points")
    .def_static(
      "from_function_space",
      [](const FunctionSpace& function_space, OutOfBoundsPolicy policy) {
        return PyEvaluationRequest{
          EvaluationRequest::FromFunctionSpace(function_space, policy), {}};
      },
      py::arg("function_space"), py::arg("policy") = OutOfBoundsPolicy{},
      "Request evaluation at a FunctionSpace's DOF-holder coordinates")
    .def_property_readonly(
      "coords",
      [](py::object self) {
        return ToPyCoordinateView(
          self.cast<const PyEvaluationRequest&>().request.coords, self);
      })
    .def_property_readonly("policy", [](const PyEvaluationRequest& self) {
      return self.request.policy;
    });

  py::class_<PointEvaluatorFactory, std::shared_ptr<PointEvaluatorFactory>>(
    m, "PointEvaluatorFactory", "Creates evaluators bound to query points")
    .def(
      "get_coordinate_system",
      [](const PointEvaluatorFactory& self) {
        return ToPyHolder(self.GetCoordinateSystem());
      },
      "The coordinate system query points must be expressed in")
    .def(
      "create_point_evaluator",
      [](const PointEvaluatorFactory& self,
         const PyEvaluationRequest& request) {
        return self.CreatePointEvaluator<Real>(request.request);
      },
      py::arg("request"), py::keep_alive<0, 1>(),
      "Create a reusable point evaluator for the request's points");

  py::class_<FunctionSpace, PointEvaluatorFactory,
             std::shared_ptr<FunctionSpace>>(m, "FunctionSpace")
    .def(
      "create_function",
      [](const FunctionSpace& self, std::string name,
         VarianceSignature variance, const SystemHolder& system,
         std::optional<ComponentScaling> scaling) {
        return self.CreateFunction<Real>(std::move(name), variance, system,
                                         scaling);
      },
      py::arg("name") = "", py::arg("variance") = values::Scalar,
      py::arg("system") = py::none(), py::arg("scaling") = py::none(),
      "Create a Function in this space. Vector and tensor values are in "
      "`system`'s basis (default: this space's coordinate system), with "
      "`scaling` required when that system admits both scalings.")
    .def("get_layout",
         [](const FunctionSpace& self) { return ToPyHolder(self.GetLayout()); })
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

  // The mesh-backed spaces refer to the mesh, so it is kept alive with them.
  lagrange_space
    .def_static(
      "from_mesh",
      [](Omega_h::Mesh& mesh, int order, int num_components,
         const SystemHolder& mesh_coordinate_system, std::string global_id_name,
         LagrangeFunctionSpace::Backend backend) {
        return LagrangeFunctionSpace::FromMesh(
          mesh, order, num_components, CartesianIfNone(mesh_coordinate_system),
          std::move(global_id_name), backend);
      },
      py::arg("mesh"), py::arg("order"), py::arg("num_components"),
      py::arg("mesh_coordinate_system") = py::none(),
      py::arg("global_id_name") = "global",
      py::arg("backend") = LagrangeFunctionSpace::DefaultBackend,
      py::keep_alive<0, 1>(),
      "Create a LagrangeFunctionSpace from an Omega_h mesh. The coordinate "
      "system defaults to Cartesian.")
    .def_static(
      "from_uniform_grid",
      [](const UniformGrid<2>& grid, int num_components,
         const SystemHolder& mesh_coordinate_system, int order) {
        return LagrangeFunctionSpace::FromUniformGrid(
          grid, num_components, CartesianIfNone(mesh_coordinate_system), order);
      },
      py::arg("grid"), py::arg("num_components"),
      py::arg("mesh_coordinate_system") = py::none(), py::arg("order") = 1,
      "Create a LagrangeFunctionSpace from a 2D uniform grid. The coordinate "
      "system defaults to Cartesian.")
    .def_static(
      "from_uniform_grid",
      [](const UniformGrid<3>& grid, int num_components,
         const SystemHolder& mesh_coordinate_system, int order) {
        return LagrangeFunctionSpace::FromUniformGrid(
          grid, num_components, CartesianIfNone(mesh_coordinate_system), order);
      },
      py::arg("grid"), py::arg("num_components"),
      py::arg("mesh_coordinate_system") = py::none(), py::arg("order") = 1,
      "Create a LagrangeFunctionSpace from a 3D uniform grid. The coordinate "
      "system defaults to Cartesian.");

  // Configuration for PolynomialReconstructionFunctionSpace MLS evaluation.
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

  py::class_<PolynomialReconstructionFunctionSpace, FunctionSpace,
             std::shared_ptr<PolynomialReconstructionFunctionSpace>>(
    m, "PolynomialReconstructionFunctionSpace")
    .def_static(
      "create",
      [](const PyCoordinateView& coords, MLSOptions options) {
        auto values = NumpyInput<Real>::ensure(coords.values);
        return PolynomialReconstructionFunctionSpace::Create(
          CoordinateView<HostMemorySpace>(coords.system,
                                          NumpyToHostView<2, Real>(values)),
          options);
      },
      py::arg("coords"), py::arg("options") = MLSOptions{},
      "Create a PolynomialReconstructionFunctionSpace over the source points "
      "of a CoordinateView")
    .def_static(
      "from_mesh",
      [](Omega_h::Mesh& mesh, int source_entity_dim,
         const SystemHolder& coordinate_system, MLSOptions options) {
        return PolynomialReconstructionFunctionSpace::FromMesh(
          mesh, source_entity_dim, CartesianIfNone(coordinate_system), options);
      },
      py::arg("mesh"), py::arg("source_entity_dim"),
      py::arg("coordinate_system") = py::none(),
      py::arg("options") = MLSOptions{}, py::keep_alive<0, 1>(),
      "Create a PolynomialReconstructionFunctionSpace from mesh entity "
      "coordinates.");

  m.def(
    "create_uniform_grid_from_mesh",
    [](Omega_h::Mesh& mesh, const std::array<LO, 2>& divisions) {
      return CreateUniformGridFromMesh<2>(mesh, divisions);
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 2D uniform grid from an Omega_h mesh");

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
      return CreateUniformGridBinaryField<2>(mesh, divisions).second;
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 2D vertex mask field indicating inside/outside mesh");
}

} // namespace pcms
