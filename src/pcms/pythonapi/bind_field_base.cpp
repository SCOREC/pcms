#include <optional>
#include <stdexcept>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/field/evaluation_request.hpp"
#include "pcms/localization/mls_options.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/discretization/discretization.hpp"
#include "pcms/field/function_spaces.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/uniform_grid.hpp"
#include "pcms/pythonapi/numpy_array_transform.hpp"
#include "pcms/pythonapi/py_value_types.hpp"

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
      "get_discretization",
      [](const FunctionSpace& self) {
        return ToPyHolder(self.GetDiscretization());
      },
      "The discretization this function space is defined on");

  py::class_<Discretization, std::shared_ptr<Discretization>>(
    m, "Discretization", "The entities (mesh, grid, points) a space lives on")
    .def("get_dimension", &Discretization::GetDimension)
    .def("get_num_entities", &Discretization::GetNumEntities,
         py::arg("entity_dim"));

  m.def(
    "create_uniform_grid_discretization",
    [](const UniformGrid<2>& grid) {
      return ToPyHolder(CreateUniformGridDiscretization(grid));
    },
    py::arg("grid"), "Discretization of a 2D uniform grid");
  m.def(
    "create_uniform_grid_discretization",
    [](const UniformGrid<3>& grid) {
      return ToPyHolder(CreateUniformGridDiscretization(grid));
    },
    py::arg("grid"), "Discretization of a 3D uniform grid");
  m.def(
    "create_empty_discretization",
    []() { return ToPyHolder(CreateEmptyDiscretization()); },
    "Discretization without entities, for evaluation-only spaces");

  // Which implementation evaluates Lagrange fields on Omega_h meshes. The
  // conservative-projection transfer operators require OmegaH.
  py::enum_<LagrangeImplementation>(m, "LagrangeImplementation")
    .value("Default", LagrangeImplementation::Default,
           "MeshFields when PCMS is built with MeshFields, OmegaH otherwise")
    .value("MeshFields", LagrangeImplementation::MeshFields,
           "MeshFields-backed layout")
    .value("OmegaH", LagrangeImplementation::OmegaH,
           "Native Omega_h Lagrange layout (required by the conservative and "
           "Monte Carlo projection transfer operators)")
    .export_values();

  // A space refers to its discretization (and through it to the mesh), so
  // the discretization is kept alive with the space.
  m.def(
    "create_lagrange_space",
    [](std::shared_ptr<Discretization> discretization, int order,
       int num_components, const SystemHolder& coordinate_system,
       LagrangeImplementation implementation, std::string layout_name) {
      return CreateLagrangeSpace(
        std::move(discretization), order, num_components,
        CartesianIfNone(coordinate_system),
        LagrangeOptions{implementation, std::move(layout_name)});
    },
    py::arg("discretization"), py::arg("order"), py::arg("num_components"),
    py::arg("coordinate_system") = py::none(),
    py::arg("implementation") = LagrangeImplementation::Default,
    py::arg("layout_name") = "", py::keep_alive<0, 1>(),
    "Lagrange space on an Omega_h or uniform-grid discretization. The "
    "coordinate system defaults to Cartesian.");

  m.def(
    "create_spline_space",
    [](std::shared_ptr<Discretization> discretization,
       const SystemHolder& coordinate_system) {
      return CreateSplineSpace(std::move(discretization),
                               CartesianIfNone(coordinate_system));
    },
    py::arg("discretization"), py::arg("coordinate_system") = py::none(),
    py::keep_alive<0, 1>(),
    "Bicubic spline space on a 2D uniform-grid discretization");

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

  m.def(
    "create_polynomial_reconstruction_space",
    [](const PyCoordinateView& coords, MLSOptions options) {
      auto values = NumpyInput<Real>::ensure(coords.values);
      return CreatePolynomialReconstructionSpace(
        CoordinateView<HostMemorySpace>(coords.system,
                                        NumpyToHostView<2, Real>(values)),
        options);
    },
    py::arg("coords"), py::arg("options") = MLSOptions{},
    "Polynomial-reconstruction (MLS) space over the points of a "
    "CoordinateView");
  m.def(
    "create_polynomial_reconstruction_space",
    [](std::shared_ptr<Discretization> discretization, int source_entity_dim,
       const SystemHolder& coordinate_system, MLSOptions options) {
      return CreatePolynomialReconstructionSpace(
        std::move(discretization), source_entity_dim,
        CartesianIfNone(coordinate_system), options);
    },
    py::arg("discretization"), py::arg("source_entity_dim"),
    py::arg("coordinate_system") = py::none(),
    py::arg("options") = MLSOptions{}, py::keep_alive<0, 1>(),
    "Polynomial-reconstruction (MLS) space over the entities of dimension "
    "source_entity_dim of an Omega_h discretization");

  // Omega_h backend entry points
  auto omega_h_module =
    m.def_submodule("backend", "Backend-specific entry points")
      .def_submodule("omega_h", "Omega_h backend");
  omega_h_module.def(
    "create_discretization",
    [](Omega_h::Mesh& mesh, std::string global_id_name) {
      return ToPyHolder(pcms::backend::omega_h::CreateDiscretization(
        mesh, std::move(global_id_name)));
    },
    py::arg("mesh"), py::arg("global_id_name") = "global",
    py::keep_alive<0, 1>(),
    "Discretization on an Omega_h mesh (the mesh is kept alive with it)");
  omega_h_module.def(
    "get_mesh",
    [](const Discretization& discretization) -> Omega_h::Mesh& {
      return pcms::backend::omega_h::GetMesh(discretization);
    },
    py::arg("discretization"), py::return_value_policy::reference,
    py::keep_alive<0, 1>(), "The Omega_h mesh of an Omega_h discretization");
  omega_h_module.def(
    "create_uniform_grid_binary_field",
    [](Omega_h::Mesh& mesh, const std::array<LO, 2>& divisions) {
      return pcms::backend::omega_h::CreateUniformGridBinaryField(mesh,
                                                                  divisions);
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 2D vertex mask field indicating inside/outside mesh");

  omega_h_module.def(
    "create_uniform_grid_from_mesh",
    [](Omega_h::Mesh& mesh, const std::array<LO, 2>& divisions) {
      return pcms::backend::omega_h::CreateUniformGridFromMesh<2>(mesh,
                                                                  divisions);
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 2D uniform grid from an Omega_h mesh");

  omega_h_module.def(
    "create_uniform_grid_from_mesh",
    [](Omega_h::Mesh& mesh, const std::array<LO, 3>& divisions) {
      return pcms::backend::omega_h::CreateUniformGridFromMesh<3>(mesh,
                                                                  divisions);
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 3D uniform grid from an Omega_h mesh");
}

} // namespace pcms
