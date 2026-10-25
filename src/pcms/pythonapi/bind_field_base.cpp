#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/field.h"
#include "pcms/field/field_evaluator_factory.h"
#include "pcms/field/field_layout.h"
#include "pcms/field/value_view.hpp"
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
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"

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

// pybind11 holders cannot be shared_ptr<const T>, so the python-facing handle
// is a small value wrapper around the shared identity object.
struct PyCoordinateSystem
{
  std::shared_ptr<const CoordinateSystem> system;
};

} // namespace

void bind_coordinate_system_module(py::module& m)
{
  auto system_cls =
    py::class_<PyCoordinateSystem>(m, "CoordinateSystem")
      .def_static(
        "Cartesian",
        []() { return PyCoordinateSystem{csys::Cartesian::Deferred()}; },
        "The dimension-deferred Cartesian placeholder, resolved against the "
        "coordinate data by the constructor that receives it")
      .def_static(
        "Cartesian",
        [](int dim) {
          return PyCoordinateSystem{csys::Cartesian::Create(dim)};
        },
        py::arg("dim"),
        "The canonical Cartesian coordinate system of the given dimension "
        "(1..3)")
      .def_property_readonly(
        "kind",
        [](const PyCoordinateSystem& self) {
          return std::string(self.system->Kind());
        },
        "Family label")
      .def_property_readonly(
        "dimension",
        [](const PyCoordinateSystem& self) { return self.system->Dimension(); },
        "Number of stored coordinate columns")
      .def_property_readonly("has_orthogonal_basis",
                             [](const PyCoordinateSystem& self) {
                               return self.system->HasOrthogonalBasis();
                             })
      .def_property_readonly("has_unit_scale_factors",
                             [](const PyCoordinateSystem& self) {
                               return self.system->HasUnitScaleFactors();
                             })
      .def_static(
        "cylindrical_rthetaz",
        []() { return PyCoordinateSystem{csys::CylindricalRThetaZ::Create()}; },
        "Cylindrical (r, theta, z), theta = atan2(y, x); dimension 3")
      .def_static(
        "cylindrical_rz",
        []() { return PyCoordinateSystem{csys::CylindricalRZ::Create()}; },
        "(R, Z) poloidal section; dimension 2")
      .def(
        "__eq__",
        [](const PyCoordinateSystem& a, const PyCoordinateSystem& b) {
          return SameCoordinateSystem(a.system, b.system);
        },
        py::is_operator())
      .def("__repr__", [](const PyCoordinateSystem& self) {
        return "<pcms.CoordinateSystem '" + std::string(self.system->Kind()) +
               "' dim " + std::to_string(self.system->Dimension()) + ">";
      });
  // Bind CoordinateView for HostMemorySpace
  py::class_<CoordinateView<HostMemorySpace>,
             std::shared_ptr<CoordinateView<HostMemorySpace>>>(m,
                                                               "CoordinateView")
    .def(py::init([](const PyCoordinateSystem& coordinate_system,
                     const NumpyInput<Real>& coords) {
           // CoordinateView does not own its coordinates, so keep a copy
           // alive in the holder's deleter.
           auto storage = CopyNumpyToKokkosView<2, Real, HostMemorySpace>(
             coords, "CoordinateView");
           return std::shared_ptr<CoordinateView<HostMemorySpace>>(
             new CoordinateView<HostMemorySpace>(
               ResolveCoordinateSystem(coordinate_system.system,
                                       static_cast<int>(storage.extent(1))),
               MakeRank2View(storage)),
             [storage](CoordinateView<HostMemorySpace>* view) { delete view; });
         }),
         py::arg("coordinate_system"), py::arg("coordinates"),
         "Create a CoordinateView from a 2D array (num_points × dim). The "
         "coordinates are copied.")

    .def(
      "get_coordinate_system",
      [](const CoordinateView<HostMemorySpace>& self) {
        return PyCoordinateSystem{self.GetCoordinateSystem()};
      },
      "The coordinate system declaring what the coordinate numbers mean")

    .def(
      "get_coordinates",
      [](py::object self) {
        return ViewToNumpy(
          self.cast<const CoordinateView<HostMemorySpace>&>().GetValues(),
          self);
      },
      "The coordinates as a read-only 2D numpy array (num_points × dim) that "
      "shares memory with this CoordinateView.");
}

void bind_create_field_module(py::module& m)
{
  py::class_<PythonEvaluationRequest>(m, "EvaluationRequest")
    .def_static(
      "from_coordinates",
      [](const NumpyInput<Real>& coords,
         const PyCoordinateSystem& coordinate_system,
         OutOfBoundsPolicy policy) {
        auto coords_device =
          CopyNumpyToKokkosView<2, Real, DeviceMemorySpace>(coords, "coords");
        return PythonEvaluationRequest{
          EvaluationRequest::FromCoordinates(
            CoordinateView<DeviceMemorySpace>(
              ResolveCoordinateSystem(
                coordinate_system.system,
                static_cast<int>(coords_device.extent(1))),
              MakeRank2View(coords_device)),
            policy),
          coords_device};
      },
      py::arg("coordinates"),
      py::arg("coordinate_system") =
        PyCoordinateSystem{csys::Cartesian::Deferred()},
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
        return ViewToNumpy(
          self.cast<const Field<Real>&>().GetDOFHolderData().GetValues(), self);
      },
      "The DOF holder data as a read-only 2D numpy array (num_dof_holders × "
      "num_components). When the field's data is in host memory the array "
      "shares that memory, so it reflects later set_dof_holder_data calls; "
      "otherwise it is a copy.")

    .def(
      "set_dof_holder_data",
      [](Field<Real>& self, const NumpyInput<Real>& arr) {
        const LO num_dof_holders = self.GetLayout().GetNumOwnedDofHolder();
        auto data = NumpyToHostComponentView<Real>(
          arr, self.GetLayout().GetNumComponents());
        if (data.extent(0) != num_dof_holders) {
          throw std::invalid_argument(
            "DOF holder data must have " + std::to_string(num_dof_holders) +
            " rows, got " + std::to_string(data.extent(0)));
        }
        self.SetDOFHolderDataUncheckedHost(data);
      },
      py::arg("data"),
      "Set the DOF holder data from a 2D numpy array (num_dof_holders × "
      "num_components). A 1D array (num_dof_holders,) is accepted when the "
      "field has one component.")

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
    .def(
      "get_coordinate_system",
      [](const FunctionSpace& self) {
        return PyCoordinateSystem{self.GetCoordinateSystem()};
      },
      "The domain coordinate system: what this space's mesh coordinates mean")
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
         const PyCoordinateSystem& coordinate_system,
         std::string global_id_name, LagrangeFunctionSpace::Backend backend) {
        return LagrangeFunctionSpace::FromMesh(
          mesh, order, num_components, coordinate_system.system,
          std::move(global_id_name), backend);
      },
      py::arg("mesh"), py::arg("order"), py::arg("num_components") = 1,
      py::arg("mesh_coordinate_system") =
        PyCoordinateSystem{csys::Cartesian::Deferred()},
      py::arg("global_id_name") = "global",
      py::arg("backend") = LagrangeFunctionSpace::DefaultBackend,
      "Create a LagrangeFunctionSpace from an Omega_h mesh")

    .def_static(
      "from_uniform_grid",
      [](const UniformGrid<2>& grid, int num_components,
         const PyCoordinateSystem& coordinate_system, int order) {
        return LagrangeFunctionSpace::FromUniformGrid(
          grid, num_components, coordinate_system.system, order);
      },
      py::arg("grid"), py::arg("num_components") = 1,
      py::arg("mesh_coordinate_system") =
        PyCoordinateSystem{csys::Cartesian::Deferred()},
      py::arg("order") = 1,
      "Create a LagrangeFunctionSpace from a 2D uniform grid")

    .def_static(
      "from_uniform_grid",
      [](const UniformGrid<3>& grid, int num_components,
         const PyCoordinateSystem& coordinate_system, int order) {
        return LagrangeFunctionSpace::FromUniformGrid(
          grid, num_components, coordinate_system.system, order);
      },
      py::arg("grid"), py::arg("num_components") = 1,
      py::arg("mesh_coordinate_system") =
        PyCoordinateSystem{csys::Cartesian::Deferred()},
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
      [](const NumpyInput<Real>& coords,
         const PyCoordinateSystem& coordinate_system, MLSOptions opts) {
        auto view = NumpyToHostView<2, Real>(coords);
        return PolynomialReconstructionFunctionSpace::Create(
          CoordinateView<HostMemorySpace>(
            ResolveCoordinateSystem(coordinate_system.system,
                                    static_cast<int>(view.extent(1))),
            view),
          opts);
      },
      py::arg("coords"),
      py::arg("coordinate_system") =
        PyCoordinateSystem{csys::Cartesian::Deferred()},
      py::arg("options") = MLSOptions{},
      "Create a PolynomialReconstructionFunctionSpace from a 2D array of "
      "source coordinates (shape: num_points × dim).")
    .def_static(
      "from_mesh",
      [](Omega_h::Mesh& mesh, int source_entity_dim,
         const PyCoordinateSystem& coordinate_system, MLSOptions opts) {
        return PolynomialReconstructionFunctionSpace::FromMesh(
          mesh, source_entity_dim, coordinate_system.system, opts);
      },
      py::arg("mesh"), py::arg("source_entity_dim"),
      py::arg("coordinate_system") =
        PyCoordinateSystem{csys::Cartesian::Deferred()},
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
      return CreateUniformGridBinaryField<2>(mesh, divisions).second;
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 2D vertex mask field indicating inside/outside mesh");
}

// bind_field_module is kept for compatibility but now registers nothing that
// refers to the deleted FieldT / LocalizationHint / FieldDataView types.
void bind_field_module(py::module& /*m*/) {}

} // namespace pcms
