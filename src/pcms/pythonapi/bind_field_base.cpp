#include <algorithm>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include <Omega_h_library.hpp>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pcms/field/coordinate_map.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/geometry/bound_source.hpp"
#include "pcms/field/differential_operator.hpp"
#include "pcms/localization/mls_options.hpp"
#include "pcms/backends/omega_h/geometry.hpp"
#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/backends/uniform_grid/geometry.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/function_spaces.hpp"
#include "pcms/backends/point_cloud/geometry.hpp"
#include "pcms/field/gradient_source.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/query_points.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/uniform_grid.hpp"
#include "pcms/pythonapi/numpy_array_transform.hpp"
#include "pcms/pythonapi/py_value_types.hpp"

namespace pcms
{
using backend::omega_h::CreateGeometry;
using backend::omega_h::OmegaHGeometry;
} // namespace pcms

namespace py = pybind11;

namespace pcms
{

namespace
{

// The Python classes of the Omega_h geometries derive from Geometry, so Python
// sees the kind of any geometry pcms returns.
struct OmegaHGeometryHandle : Immutable<Geometry>
{};

void SetDOFHolderData(Field<Real>& self, const PyValueView& v)
{
  auto array = NumpyInput<Real>::ensure(v.values);
  if (!array) {
    throw py::type_error("set_dof_holder_data: values must be float64");
  }
  self.SetDOFHolderDataHost(ValueView<const Real, HostMemorySpace>(
    v.basis, NumpyToHostView<2, Real>(array)));
}

/// Return the Python handle matching the geometry kind.
py::object GeometryToPython(const std::shared_ptr<const Geometry>& geometry)
{
  if (geometry == nullptr) {
    return py::none();
  }
  if (dynamic_cast<const OmegaHGeometry*>(geometry.get()) != nullptr) {
    return py::cast(OmegaHGeometryHandle{{geometry}});
  }
  return py::cast(Immutable<Geometry>{geometry});
}

/// The methods a source of evaluators has, on the handle of a
/// PointEvaluatorFactory or of a FunctionSpace (which is one, and converts to
/// it where a source is taken).
template <typename Handle>
void BindEvaluatorFactoryMethods(py::class_<Handle>& cls)
{
  cls
    .def(
      "get_geometry",
      [](const Handle& self) {
        return GeometryToPython(self.ptr->GetGeometry());
      },
      "The geometry whose bound sources this source accepts; bind query points "
      "with its bind()")
    .def(
      "create_point_evaluator",
      [](const Handle& self, const Immutable<BoundSource>& bound,
         DifferentialOperator op) {
        return self.ptr->template CreatePointEvaluator<Real>(bound.ptr, op);
      },
      py::arg("bound"), py::arg("op") = DifferentialOperator::Value,
      py::keep_alive<0, 1>(),
      "Create an evaluator at the points of `bound` (a result of bind() on "
      "this source's geometry) that applies `op`: the values for "
      "DifferentialOperator.Value, the gradient for "
      "DifferentialOperator.Gradient")
    .def(
      "supports",
      [](const Handle& self, DifferentialOperator op) {
        return self.ptr->Supports(op);
      },
      py::arg("op"), "Whether create_point_evaluator accepts `op`");
}

} // namespace

void bind_create_field_module(py::module& m)
{
  ImmutableClass<FieldLayout>(m, "FieldLayout",
                              "How a field's DOF holders are laid out")
    .def("get_num_components",
         [](const Immutable<FieldLayout>& self) {
           return self.ptr->GetNumComponents();
         })
    .def(
      "get_num_dof_holder",
      [](const Immutable<FieldLayout>& self) {
        return self.ptr->GetNumDofHolder();
      },
      "Number of DOF holders this rank owns")
    .def("get_num_global_dof_holder",
         [](const Immutable<FieldLayout>& self) {
           return self.ptr->GetNumGlobalDofHolder();
         })
    .def("get_dimension",
         [](const Immutable<FieldLayout>& self) {
           return self.ptr->GetDimension();
         })
    .def("is_distributed",
         [](const Immutable<FieldLayout>& self) {
           return self.ptr->IsDistributed();
         })
    .def(
      "get_coordinate_system",
      [](const Immutable<FieldLayout>& self) {
        return ToPython(self.ptr->GetCoordinateSystem());
      },
      "What this layout's DOF-holder coordinates mean")
    .def(
      "get_dof_holder_coordinates",
      [](py::object self) {
        return ToPyCoordinateView(self.cast<const Immutable<FieldLayout>&>()
                                    .ptr->GetDOFHolderCoordinates(),
                                  self);
      },
      "The owned DOF-holder coordinates (num_dof_holders × dim), read-only");

  py::class_<FieldData<Real>>(m, "FieldData", "A field's stored values")
    .def("get_value_type", &FieldData<Real>::GetValueType)
    .def("get_value_basis", &FieldData<Real>::GetValueBasis,
         "The basis the stored components are in");

  py::class_<Field<Real>>(
    m, "Field", "DOF-holder values together with the layout that orders them")
    .def("get_name", &Field<Real>::GetName)
    .def("get_layout",
         [](const Field<Real>& self) {
           return Immutable<FieldLayout>{self.GetLayoutPtr()};
         })
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
        const LO num_dof_holders = layout.GetNumDofHolder();
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
      "field has one component.")
    .def(
      "synchronize", [](Field<Real>& self) { self.Synchronize(); },
      "Bring the copies of DOF holders other ranks own (e.g. Omega_h ghosts) "
      "up to date with their owners' values");

  py::class_<Function<Real>, Field<Real>>(
    m, "Function", "A Field that knows its FunctionSpace")
    .def("get_space", [](const Function<Real>& self) {
      return Immutable<FunctionSpace>{self.GetSpacePtr()};
    });

  ImmutableClass<Geometry>(
    m, "Geometry",
    "The coordinates of a set of entities (a mesh, a grid, a point cloud) and "
    "the point search; function spaces are created on it")
    .def("get_coordinate_system",
         [](const Immutable<Geometry>& g) {
           return ToPython(g.ptr->GetCoordinateSystem());
         })
    .def(
      "bind",
      [](const Immutable<Geometry>& g, const PyCoordinateView& coords,
         OutOfBoundsPolicy policy,
         const std::optional<Immutable<CoordinateMap>>& map) {
        // the bound source copies the points
        auto device = CopyCoordinatesToDevice(coords);
        return Immutable<BoundSource>{
          g.ptr->Bind(device.view, policy, map ? map->ptr : nullptr)};
      },
      py::arg("coordinates"), py::arg("policy") = OutOfBoundsPolicy{},
      py::arg("map") = py::none(),
      "Localize the points of a CoordinateView (in get_coordinate_system(), "
      "or in the source system of `map`, a CoordinateMap to "
      "get_coordinate_system(): the points are then located at their images "
      "and evaluated values re-expressed through it); the search is done "
      "once here")
    .def(
      "bind",
      [](const Immutable<Geometry>& g, const Immutable<FunctionSpace>& space,
         OutOfBoundsPolicy policy,
         const std::optional<Immutable<CoordinateMap>>& map) {
        auto layout = space.ptr->GetLayout();
        return Immutable<BoundSource>{g.ptr->Bind(
          layout->GetDOFHolderCoordinates(), DOFHolderOwners(*layout), policy,
          map ? map->ptr : nullptr)};
      },
      py::arg("space"), py::arg("policy") = OutOfBoundsPolicy{},
      py::arg("map") = py::none(),
      "Bind the DOF-holder coordinates of `space`, through `map` if given; "
      "without a map, a space whose DOF holders are entities of this geometry "
      "binds them with owner hints, which are verified (with a map the hints "
      "are ignored)");

  [[maybe_unused]] const auto omega_h_geometry =
    py::class_<OmegaHGeometryHandle, Immutable<Geometry>>(
      m, "OmegaHGeometry", "The geometry of an Omega_h simplex mesh");

  auto omega_h_module =
    m.def_submodule("backend", "Backend-specific entry points")
      .def_submodule("omega_h", "Omega_h backend");
  omega_h_module.def(
    "create_geometry",
    [](std::shared_ptr<Omega_h::Mesh> mesh,
       const std::optional<Immutable<CoordinateSystem>>& system, int divisions,
       std::string global_id_name,
       const std::optional<py::array_t<
         std::int8_t, py::array::c_style | py::array::forcecast>>& owned_mask) {
      auto library =
        py::cast(mesh->library(), py::return_value_policy::reference)
          .cast<std::shared_ptr<Omega_h::Library>>();
      backend::omega_h::OmegaHGeometryOptions options;
      options.search.divisions = divisions;
      options.global_id_name = std::move(global_id_name);
      if (owned_mask) {
        if (owned_mask->ndim() != 1) {
          throw py::value_error(
            "create_geometry: owned_mask must be a 1D array");
        }
        Omega_h::HostWrite<Omega_h::I8> host(
          static_cast<Omega_h::LO>(owned_mask->size()));
        std::copy(owned_mask->data(), owned_mask->data() + owned_mask->size(),
                  host.data());
        options.owned_mask = Omega_h::Read<Omega_h::I8>(host.write());
      }
      return OmegaHGeometryHandle{
        {CreateGeometry(std::move(mesh), std::move(library),
                        CartesianIfNone(system), std::move(options))}};
    },
    py::arg("mesh").none(false), py::arg("system") = py::none(),
    py::arg("divisions") = 0, py::arg("global_id_name") = "global",
    py::arg("owned_mask") = py::none(),
    "Create the geometry of an Omega_h mesh: its coordinate system (default "
    "Cartesian) and point search (`divisions` grid divisions per axis, 0 to "
    "derive them from the element count); `global_id_name` names the tag "
    "that holds the global entity ids; `owned_mask` is a per-vertex (int8) "
    "ownership override, by default the mesh's ownership");

  m.def(
    "create_uniform_grid_geometry",
    [](const UniformGrid<2>& grid,
       const std::optional<Immutable<CoordinateSystem>>& system) {
      return Immutable<Geometry>{
        backend::uniform_grid::CreateUniformGridGeometry(
          grid, CartesianIfNone(system))};
    },
    py::arg("grid"), py::arg("system") = py::none(),
    "Create the geometry of a 2D uniform grid; the coordinate system defaults "
    "to Cartesian");

  m.def(
    "create_uniform_grid_geometry",
    [](const UniformGrid<3>& grid,
       const std::optional<Immutable<CoordinateSystem>>& system) {
      return Immutable<Geometry>{
        backend::uniform_grid::CreateUniformGridGeometry(
          grid, CartesianIfNone(system))};
    },
    py::arg("grid"), py::arg("system") = py::none(),
    "Create the geometry of a 3D uniform grid; the coordinate system defaults "
    "to Cartesian");

  ImmutableClass<BoundSource>(
    m, "BoundSource",
    "A set of query points localized on a geometry; evaluators created from "
    "it share its search")
    .def(
      "get_geometry",
      [](const Immutable<BoundSource>& self) {
        return GeometryToPython(self.ptr->GetGeometry());
      },
      "The geometry the points were localized on")
    .def(
      "num_points",
      [](const Immutable<BoundSource>& self) { return self.ptr->NumPoints(); })
    .def_property_readonly(
      "policy",
      [](const Immutable<BoundSource>& self) { return self.ptr->Policy(); })
    .def(
      "found_points",
      [](py::object self) {
        return IndicesToNumpy(
          self.cast<const Immutable<BoundSource>&>().ptr->FoundPoints(), self);
      },
      "Indices of the query points the source evaluates, read-only")
    .def(
      "missing_points",
      [](py::object self) {
        return IndicesToNumpy(
          self.cast<const Immutable<BoundSource>&>().ptr->MissingPoints(),
          self);
      },
      "Indices of the query points the source does not evaluate, read-only");

  py::enum_<DifferentialOperator>(m, "DifferentialOperator",
                                  "The quantity a point evaluator computes")
    .value("Value", DifferentialOperator::Value)
    .value("Gradient", DifferentialOperator::Gradient)
    .value("Divergence", DifferentialOperator::Divergence)
    .value("Curl", DifferentialOperator::Curl)
    .value("Hessian", DifferentialOperator::Hessian);

  auto factory = ImmutableClass<PointEvaluatorFactory>(
    m, "PointEvaluatorFactory",
    "Creates evaluators at the points of bound sources of its geometry");
  BindEvaluatorFactoryMethods(factory);

  py::class_<QueryPoints>(
    m, "QueryPoints",
    "Evaluation at points that move between calls: locate() binds the points "
    "and rebinds them in place after; evaluate() creates the evaluator of a "
    "function on first use after each locate()")
    .def(
      py::init([](const Immutable<Geometry>& geometry, OutOfBoundsPolicy policy,
                  const std::optional<Immutable<CoordinateMap>>& map) {
        return std::make_unique<QueryPoints>(geometry.ptr, policy,
                                             map ? map->ptr : nullptr);
      }),
      py::arg("geometry"), py::arg("policy") = OutOfBoundsPolicy{},
      py::arg("map") = py::none())
    .def(
      "locate",
      [](QueryPoints& self, const PyCoordinateView& coords) {
        // the bound source copies the points
        auto device = CopyCoordinatesToDevice(coords);
        self.Locate(device.view);
      },
      py::arg("coordinates"),
      "Locate the points of a CoordinateView (as Geometry.bind takes them), "
      "reusing the previous points' storage when nothing else holds it")
    .def(
      "evaluate",
      [](QueryPoints& self, const Function<Real>& function, PyValueView output,
         DifferentialOperator op) {
        WithDeviceOutput(output.values, [&](auto out) {
          self.Evaluate(function,
                        ValueView<Real, DeviceMemorySpace>(output.basis, out),
                        op);
        });
      },
      py::arg("function"), py::arg("output"),
      py::arg("op") = DifferentialOperator::Value,
      "Evaluate `op` of the function at the located points into `output`, a "
      "ValueView as PointEvaluator.evaluate takes it")
    .def(
      "evaluate",
      [](QueryPoints& self, const Function<Real>& function, py::array output,
         DifferentialOperator op) {
        const ValueBasis basis = function.GetData().GetValueBasis();
        RequireScalarForPlainArray("evaluate", basis);
        WithDeviceOutput(output, [&](auto out) {
          self.Evaluate(function,
                        ValueView<Real, DeviceMemorySpace>(basis, out), op);
        });
      },
      py::arg("function"), py::arg("output").noconvert(),
      py::arg("op") = DifferentialOperator::Value,
      "Evaluate the values of a scalar function at the located points into a "
      "writeable C-contiguous float64 array (num_points × num_components)")
    .def_property_readonly(
      "bound",
      [](const QueryPoints& self) -> std::optional<Immutable<BoundSource>> {
        if (self.Bound() == nullptr) {
          return std::nullopt;
        }
        return Immutable<BoundSource>{self.Bound()};
      },
      "The BoundSource of the located points, or None before locate(); while "
      "it is held, the next locate() binds anew instead of reusing it");

  m.def(
    "evaluate",
    [](const Function<Real>& function, const PyCoordinateView& points,
       PyValueView output, DifferentialOperator op, OutOfBoundsPolicy policy,
       const std::optional<Immutable<CoordinateMap>>& map) {
      auto device = CopyCoordinatesToDevice(points);
      WithDeviceOutput(output.values, [&](auto out) {
        Evaluate(function, device.view,
                 ValueView<Real, DeviceMemorySpace>(output.basis, out), op,
                 policy, map ? map->ptr : nullptr);
      });
    },
    py::arg("function"), py::arg("points"), py::arg("output"),
    py::arg("op") = DifferentialOperator::Value,
    py::arg("policy") = OutOfBoundsPolicy{}, py::arg("map") = py::none(),
    "Evaluate `op` of the function at the points of a CoordinateView into "
    "`output`, a ValueView; binds the points on every call (use QueryPoints "
    "for repeated evaluation)");
  m.def(
    "evaluate",
    [](const Function<Real>& function, const PyCoordinateView& points,
       py::array output, DifferentialOperator op, OutOfBoundsPolicy policy,
       const std::optional<Immutable<CoordinateMap>>& map) {
      const ValueBasis basis = function.GetData().GetValueBasis();
      RequireScalarForPlainArray("evaluate", basis);
      auto device = CopyCoordinatesToDevice(points);
      WithDeviceOutput(output, [&](auto out) {
        Evaluate(function, device.view,
                 ValueView<Real, DeviceMemorySpace>(basis, out), op, policy,
                 map ? map->ptr : nullptr);
      });
    },
    py::arg("function"), py::arg("points"), py::arg("output").noconvert(),
    py::arg("op") = DifferentialOperator::Value,
    py::arg("policy") = OutOfBoundsPolicy{}, py::arg("map") = py::none(),
    "Evaluate the values of a scalar function at the points of a "
    "CoordinateView into a writeable C-contiguous float64 array; binds the "
    "points on every call (use QueryPoints for repeated evaluation)");

  m.def(
    "gradient",
    [](const Immutable<PointEvaluatorFactory>& source) {
      return Immutable<PointEvaluatorFactory>{
        pcms::detail::Gradient(source.ptr)};
    },
    py::arg("source"),
    "A source whose evaluators write the gradient of the field they are given; "
    "it accepts the bound sources of `source`'s geometry");

  m.def(
    "divergence",
    [](const Immutable<PointEvaluatorFactory>& source) {
      return Immutable<PointEvaluatorFactory>{
        pcms::detail::Divergence(source.ptr)};
    },
    py::arg("source"),
    "A source whose evaluators write the divergence of the field they are "
    "given; it accepts the bound sources of `source`'s geometry");

  m.def(
    "curl",
    [](const Immutable<PointEvaluatorFactory>& source) {
      return Immutable<PointEvaluatorFactory>{pcms::detail::Curl(source.ptr)};
    },
    py::arg("source"),
    "A source whose evaluators write the curl of the field they are given; it "
    "accepts the bound sources of `source`'s geometry");

  auto space = ImmutableClass<FunctionSpace>(m, "FunctionSpace");
  BindEvaluatorFactoryMethods(space);
  // a FunctionSpace is accepted where a PointEvaluatorFactory is
  factory.def(py::init([](const Immutable<FunctionSpace>& s) {
    return Immutable<PointEvaluatorFactory>{s.ptr};
  }));
  py::implicitly_convertible<Immutable<FunctionSpace>,
                             Immutable<PointEvaluatorFactory>>();
  space
    .def(
      "create_function",
      [](const Immutable<FunctionSpace>& self, std::string name,
         VarianceSignature variance,
         const std::optional<Immutable<CoordinateSystem>>& system,
         std::optional<ComponentScaling> scaling) {
        return self.ptr->CreateFunction<Real>(
          std::move(name), variance, system ? system->ptr : nullptr, scaling);
      },
      py::arg("name") = "", py::arg("variance") = values::Scalar,
      py::arg("system") = py::none(), py::arg("scaling") = py::none(),
      "Create a Function in this space. Vector and tensor values are in "
      "`system`'s basis (default: this space's coordinate system), with "
      "`scaling` required when that system admits both scalings.")
    .def("get_layout",
         [](const Immutable<FunctionSpace>& self) {
           return Immutable<FieldLayout>{self.ptr->GetLayout()};
         })
    .def(
      "get_coordinate_system",
      [](const Immutable<FunctionSpace>& self) {
        return ToPython(self.ptr->GetCoordinateSystem());
      },
      "The coordinate system of the space's geometry");

  // Backend selecting the underlying field layout. The conservative-projection
  // transfer operators require the native Omega_h backend.
  py::enum_<LagrangeBackend>(m, "LagrangeBackend")
    .value("MeshFields", LagrangeBackend::MeshFields,
           "MeshFields-backed layout (default when MeshFields is enabled)")
    .value("OmegaH", LagrangeBackend::OmegaH,
           "Native Omega_h Lagrange layout (required by the conservative and "
           "Monte Carlo projection transfer operators)");

  // The spaces keep their geometry alive, and with it the mesh of an Omega_h
  // geometry.
  m.def(
    "create_lagrange_space",
    [](const Immutable<Geometry>& geometry, int order, int num_components,
       LagrangeBackend backend, std::string layout_name) {
      return Immutable<FunctionSpace>{CreateLagrangeSpace(
        geometry.ptr, order, num_components,
        LagrangeOptions{.backend = backend,
                        .layout_name = std::move(layout_name)})};
    },
    py::arg("geometry"), py::arg("order"), py::arg("num_components"),
    py::arg("backend") = DefaultLagrangeBackend, py::arg("layout_name") = "",
    "Create a Lagrange FunctionSpace of order `order` on an Omega_h geometry "
    "(backend OmegaH: orders 0 and 1; MeshFields: orders 1 and 2) or a "
    "uniform-grid geometry (orders 0 and 1, `backend` ignored); spaces on one "
    "geometry accept each other's bound sources.");

  m.def(
    "create_spline_space",
    [](const Immutable<Geometry>& geometry) {
      return Immutable<FunctionSpace>{CreateSplineSpace(geometry.ptr)};
    },
    py::arg("geometry"),
    "Create a bicubic spline FunctionSpace on a 2D uniform-grid geometry");

  // Configuration for MLS evaluation of polynomial reconstruction spaces.
  py::class_<MLSOptions>(m, "MLSOptions")
    .def(py::init<>(), "Default MLSOptions")
    .def_readwrite("degree", &MLSOptions::degree)
    .def_readwrite("lambda_reg", &MLSOptions::lambda)
    .def_readwrite("tol", &MLSOptions::tol)
    .def_readwrite("decay_factor", &MLSOptions::decay_factor)
    .def_readwrite("basis", &MLSOptions::basis);

  py::class_<NeighbourhoodOptions>(m, "NeighbourhoodOptions",
                                   "Support search of a point-cloud geometry")
    .def(py::init<>(), "Default NeighbourhoodOptions")
    .def_readwrite("radius", &NeighbourhoodOptions::radius)
    .def_readwrite("min_supports", &NeighbourhoodOptions::min_supports)
    .def_readwrite("adapt_radius", &NeighbourhoodOptions::adapt_radius);

  m.def(
    "create_polynomial_reconstruction_space",
    [](const Immutable<Geometry>& geometry, MLSOptions options) {
      return Immutable<FunctionSpace>{
        CreatePolynomialReconstructionSpace(geometry.ptr, options)};
    },
    py::arg("geometry"), py::arg("options") = MLSOptions{},
    "Create an MLS reconstruction FunctionSpace on a point-cloud geometry "
    "(its sites are the DOF holders); `options` are the reconstruction "
    "parameters, the support search is the geometry's");

  m.def(
    "create_point_cloud_geometry",
    [](const PyCoordinateView& coords, NeighbourhoodOptions neighbourhood) {
      auto values = NumpyInput<Real>::ensure(coords.values);
      return Immutable<Geometry>{backend::point_cloud::CreatePointCloudGeometry(
        CoordinateView<HostMemorySpace>(coords.system,
                                        NumpyToHostView<2, Real>(values)),
        std::move(neighbourhood))};
    },
    py::arg("coords"), py::arg("neighbourhood") = NeighbourhoodOptions{},
    "Create a point-cloud geometry whose sites are the points of a "
    "CoordinateView; `neighbourhood` is its support search");

  m.def(
    "create_point_cloud_geometry",
    [](const Immutable<FieldLayout>& sites,
       NeighbourhoodOptions neighbourhood) {
      return Immutable<Geometry>{backend::point_cloud::CreatePointCloudGeometry(
        sites.ptr, std::move(neighbourhood))};
    },
    py::arg("sites"), py::arg("neighbourhood") = NeighbourhoodOptions{},
    "Create a point-cloud geometry whose sites are the DOF holders of a "
    "layout; `neighbourhood` is its support search");

  omega_h_module.def(
    "create_point_cloud_geometry",
    [](const Immutable<Geometry>& geometry, int entity_dim,
       NeighbourhoodOptions neighbourhood) {
      const auto omega_h =
        std::dynamic_pointer_cast<const backend::omega_h::OmegaHGeometry>(
          geometry.ptr);
      if (omega_h == nullptr) {
        throw py::value_error(
          "create_point_cloud_geometry: the geometry must be an Omega_h "
          "geometry");
      }
      try {
        return Immutable<Geometry>{backend::omega_h::CreatePointCloudGeometry(
          omega_h, entity_dim, std::move(neighbourhood))};
      } catch (const pcms_error& e) {
        throw py::value_error(e.what());
      }
    },
    py::arg("geometry"), py::arg("entity_dim"),
    py::arg("neighbourhood") = NeighbourhoodOptions{}, py::keep_alive<0, 1>(),
    "A point-cloud geometry whose sites are the vertices (entity_dim 0) or "
    "the element centroids (the mesh dimension) of the mesh of an Omega_h "
    "geometry; the support search uses the mesh adjacency. Spaces made on the "
    "same geometry object share its entities, so binding the cloud at their "
    "DOF holders (bind(space)) starts from those entities.");
  omega_h_module.def(
    "get_mesh",
    [](const Immutable<Geometry>& geometry) -> std::shared_ptr<Omega_h::Mesh> {
      const auto* omega_h =
        dynamic_cast<const backend::omega_h::OmegaHGeometry*>(
          geometry.ptr.get());
      if (omega_h == nullptr) {
        throw py::value_error(
          "get_mesh: the geometry must be an Omega_h geometry");
      }
      // owning handle aliased to the geometry: the mesh (and the Omega_h
      // library the geometry holds) stay alive even when the geometry is a
      // temporary and the original Python mesh wrapper is gone
      return std::shared_ptr<Omega_h::Mesh>(geometry.ptr, &omega_h->GetMesh());
    },
    py::arg("geometry"), "The Omega_h mesh of an Omega_h geometry");

  omega_h_module.def(
    "create_uniform_grid_binary_field",
    [](Omega_h::Mesh& mesh, const std::array<LO, 2>& divisions) {
      return backend::omega_h::CreateUniformGridBinaryField(mesh, divisions);
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
