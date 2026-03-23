#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pcms/adapter/meshfields/mesh_fields_adapter_layout.h"
#include "pcms/field.h"
#include "pcms/field_layout.h"
#include "numpy_array_transform.h"

namespace py = pybind11;

namespace pcms
{

void bind_omega_h_field_layout_module(py::module& m)
{
  // Bind the MeshFieldsAdapterLayout class
  py::class_<MeshFieldsAdapterLayout, FieldLayout,
             std::shared_ptr<MeshFieldsAdapterLayout>>(
    m, "MeshFieldsAdapterLayout")
    .def(py::init<Omega_h::Mesh&, std::array<int, 4>, int, CoordinateSystem,
                  std::string>(),
         py::arg("mesh"), py::arg("nodes_per_dim"), py::arg("num_components"),
         py::arg("coordinate_system"), py::arg("global_id_name") = "global",
         "Constructor for MeshFieldsAdapterLayout")

    .def(
      "create_field",
      [](MeshFieldsAdapterLayout& self) {
        return std::shared_ptr<FieldT<Real>>(self.CreateFieldReal());
      },
      "Create a field with this layout")

    .def("get_num_components", &MeshFieldsAdapterLayout::GetNumComponents,
         "Get the number of components in the field")

    .def("get_num_owned_dof_holder",
         &MeshFieldsAdapterLayout::GetNumOwnedDofHolder,
         "Get the number of owned DOF holders")

    .def("get_num_global_dof_holder",
         &MeshFieldsAdapterLayout::GetNumGlobalDofHolder,
         "Get the number of global DOF holders")

    .def(
      "get_owned",
      [](const MeshFieldsAdapterLayout& self) {
        auto owned = self.GetOwned();
        // Convert to numpy array
        py::array_t<bool> result(owned.extent(0));
        auto buf = result.request();
        bool* ptr = static_cast<bool*>(buf.ptr);
        for (size_t i = 0; i < owned.extent(0); ++i) {
          ptr[i] = owned[i];
        }
        return result;
      },
      "Get the owned mask array")

    .def(
      "get_gids",
      [](const MeshFieldsAdapterLayout& self) {
        auto gids = self.GetGids();
        // Convert to numpy array
        py::array_t<GO> result(gids.extent(0));
        auto buf = result.request();
        GO* ptr = static_cast<GO*>(buf.ptr);
        for (size_t i = 0; i < gids.extent(0); ++i) {
          ptr[i] = gids[i];
        }
        return result;
      },
      "Get the global IDs array")

    .def(
      "get_dof_holder_coordinates",
      [](const MeshFieldsAdapterLayout& self) {
        auto coords = self.GetDOFHolderCoordinates().GetCoordinates();
        // Convert to numpy array (2D)
        py::array_t<Real> result({coords.extent(0), coords.extent(1)});
        auto buf = result.request();
        Real* ptr = static_cast<Real*>(buf.ptr);
        for (size_t i = 0; i < coords.extent(0); ++i) {
          for (size_t j = 0; j < coords.extent(1); ++j) {
            ptr[i * coords.extent(1) + j] = coords(i, j);
          }
        }
        return result;
      },
      "Get the DOF holder coordinates")

    .def("is_distributed", &MeshFieldsAdapterLayout::IsDistributed,
         "Check if the field layout is distributed")

    .def(
      "get_ent_offsets",
      [](const MeshFieldsAdapterLayout& self) {
        auto offsets = self.GetEntOffsets();
        // Convert std::array to list/tuple
        py::list result;
        for (size_t i = 0; i < offsets.size(); ++i) {
          result.append(offsets[i]);
        }
        return result;
      },
      "Get the entity offsets array")

    .def(
      "get_nodes_per_dim",
      [](const MeshFieldsAdapterLayout& self) {
        auto nodes = self.GetNodesPerDim();
        py::list result;
        for (size_t i = 0; i < nodes.size(); ++i) {
          result.append(nodes[i]);
        }
        return result;
      },
      "Get the nodes per dimension array")

    .def("get_num_ents", &MeshFieldsAdapterLayout::GetNumEnts,
         "Get the total number of entities")

    .def(
      "get_mesh",
      [](MeshFieldsAdapterLayout& self) -> Omega_h::Mesh& {
        return self.GetMesh();
      },
      py::return_value_policy::reference, "Get the underlying Omega_h mesh")

    .def("owned_size", &MeshFieldsAdapterLayout::OwnedSize,
         "Get the owned size (num_components * num_owned_dof_holder)")

    .def("global_size", &MeshFieldsAdapterLayout::GlobalSize,
         "Get the global size (num_components * num_global_dof_holder)");
}

} // namespace pcms