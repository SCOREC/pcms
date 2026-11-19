#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>

#include <Omega_h_library.hpp>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pcms/backends/omega_h/geometry.hpp"
#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/configuration.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/pythonapi/py_value_types.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/uniform_grid.hpp"

namespace py = pybind11;

namespace pcms
{

using omega_h::CreateGeometry;

// pcms::omega_h -> pcms.omega_h. Bound after the point-cloud module, whose
// NeighbourhoodOptions the Omega_h point-cloud geometry takes.
void bind_omega_h_module(py::module& m)
{
  auto omega_h_module = m.def_submodule("omega_h", "Omega_h backend");

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
      omega_h::OmegaHGeometryOptions options;
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
      return Immutable<Geometry>{
        CreateGeometry(std::move(mesh), std::move(library),
                       CartesianIfNone(system), std::move(options))};
    },
    py::arg("mesh").none(false), py::arg("system") = py::none(),
    py::arg("divisions") = 0, py::arg("global_id_name") = "global",
    py::arg("owned_mask") = py::none(),
    "Create the geometry of an Omega_h mesh: its coordinate system (default "
    "Cartesian) and point search (`divisions` grid divisions per axis, 0 to "
    "derive them from the element count); `global_id_name` names the tag "
    "that holds the global entity ids; `owned_mask` is a per-vertex (int8) "
    "ownership override, by default the mesh's ownership");

  omega_h_module.def(
    "create_point_cloud_geometry",
    [](const Immutable<Geometry>& geometry, int entity_dim,
       NeighbourhoodOptions neighbourhood) {
      try {
        return Immutable<Geometry>{omega_h::CreatePointCloudGeometry(
          geometry.ptr, entity_dim, std::move(neighbourhood))};
      } catch (const pcms_error& e) {
        throw py::value_error(e.what());
      }
    },
    py::arg("geometry"), py::arg("entity_dim"),
    py::arg("neighbourhood") = NeighbourhoodOptions{},
    "A point-cloud geometry whose sites are the vertices (entity_dim 0) or "
    "the element centroids (the mesh dimension) of the mesh of an Omega_h "
    "geometry; the support search uses the mesh adjacency. Spaces made on the "
    "same geometry object share its entities, so binding the cloud at their "
    "DOF holders (bind(space)) starts from those entities.");

  omega_h_module.def(
    "create_uniform_grid_binary_field",
    [](Omega_h::Mesh& mesh, const std::array<LO, 2>& divisions) {
      return omega_h::CreateUniformGridBinaryField(mesh, divisions);
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 2D vertex mask field indicating inside/outside mesh");

  omega_h_module.def(
    "create_uniform_grid_from_mesh",
    [](Omega_h::Mesh& mesh, const std::array<LO, 2>& divisions) {
      return pcms::omega_h::CreateUniformGridFromMesh<2>(mesh, divisions);
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 2D uniform grid from an Omega_h mesh");

  omega_h_module.def(
    "create_uniform_grid_from_mesh",
    [](Omega_h::Mesh& mesh, const std::array<LO, 3>& divisions) {
      return pcms::omega_h::CreateUniformGridFromMesh<3>(mesh, divisions);
    },
    py::arg("mesh"), py::arg("divisions"),
    "Create a 3D uniform grid from an Omega_h mesh");
}

} // namespace pcms
