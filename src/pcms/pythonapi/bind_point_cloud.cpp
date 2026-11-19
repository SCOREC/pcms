#include <memory>
#include <utility>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pcms/configuration.hpp"
#include "pcms/backends/point_cloud/geometry.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/function_spaces.hpp"
#include "pcms/field/mls_options.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/geometry/supports.hpp"
#include "pcms/pythonapi/numpy_array_transform.hpp"
#include "pcms/pythonapi/py_value_types.hpp"

namespace py = pybind11;

namespace pcms
{

// pcms::point_cloud -> pcms.point_cloud, and the polynomial-reconstruction
// space when it is built
void bind_point_cloud_module(py::module& m)
{
  auto point_cloud_module =
    m.def_submodule("point_cloud", "Point-cloud backend");

  py::class_<NeighbourhoodOptions>(m, "NeighbourhoodOptions",
                                   "Support search of a point-cloud geometry")
    .def(py::init<>(), "Default NeighbourhoodOptions")
    .def_readwrite("radius", &NeighbourhoodOptions::radius)
    .def_readwrite("min_supports", &NeighbourhoodOptions::min_supports)
    .def_readwrite("adapt_radius", &NeighbourhoodOptions::adapt_radius);

#ifdef PCMS_ENABLE_POLYNOMIAL_RECONSTRUCTION
  // RadialBasisFunction is part of the MLSOptions API used by
  // create_polynomial_reconstruction_space.
  py::enum_<pcms::RadialBasisFunction>(m, "RadialBasisFunction")
    .value("RBF_GAUSSIAN", pcms::RadialBasisFunction::RBF_GAUSSIAN)
    .value("RBF_C4", pcms::RadialBasisFunction::RBF_C4)
    .value("RBF_CONST", pcms::RadialBasisFunction::RBF_CONST)
    .value("NO_OP", pcms::RadialBasisFunction::NO_OP)
    .export_values();

  // Configuration for MLS evaluation of polynomial reconstruction spaces.
  py::class_<MLSOptions>(m, "MLSOptions")
    .def(py::init<>(), "Default MLSOptions")
    .def_readwrite("degree", &MLSOptions::degree)
    .def_readwrite("lambda_reg", &MLSOptions::lambda)
    .def_readwrite("tol", &MLSOptions::tol)
    .def_readwrite("decay_factor", &MLSOptions::decay_factor)
    .def_readwrite("basis", &MLSOptions::basis);

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
#endif

  point_cloud_module.def(
    "create_geometry",
    [](const PyCoordinateView& coords, NeighbourhoodOptions neighbourhood) {
      auto values = NumpyInput<Real>::ensure(coords.values);
      return Immutable<Geometry>{point_cloud::CreateGeometry(
        CoordinateView<HostMemorySpace>(coords.system,
                                        NumpyToHostView<2, Real>(values)),
        std::move(neighbourhood))};
    },
    py::arg("coords"), py::arg("neighbourhood") = NeighbourhoodOptions{},
    "Create a point-cloud geometry whose sites are the points of a "
    "CoordinateView; `neighbourhood` is its support search");

  point_cloud_module.def(
    "create_geometry",
    [](const Immutable<FieldLayout>& sites,
       NeighbourhoodOptions neighbourhood) {
      return Immutable<Geometry>{
        point_cloud::CreateGeometry(sites.ptr, std::move(neighbourhood))};
    },
    py::arg("sites"), py::arg("neighbourhood") = NeighbourhoodOptions{},
    "Create a point-cloud geometry whose sites are the DOF holders of a "
    "layout; `neighbourhood` is its support search");
}

} // namespace pcms
