#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pcms/field_layout.h"
#include "numpy_array_transform.h"

namespace py = pybind11;

namespace pcms
{

void bind_field_layout_module(py::module& m)
{
  // Bind the base FieldLayout class as an abstract base
  py::class_<FieldLayout, std::shared_ptr<FieldLayout>>(m, "FieldLayout")
    .def("create_field", &FieldLayout::CreateField,
         "Create a field with this layout")

    .def("get_num_components", &FieldLayout::GetNumComponents,
         "Get the number of components in the field")

    .def("get_num_owned_dof_holder", &FieldLayout::GetNumOwnedDofHolder,
         "Get the number of owned DOF holders")

    .def("get_num_global_dof_holder", &FieldLayout::GetNumGlobalDofHolder,
         "Get the number of global DOF holders")

    .def("owned_size", &FieldLayout::OwnedSize,
         "Get the owned size (num_components * num_owned_dof_holder)")

    .def("global_size", &FieldLayout::GlobalSize,
         "Get the global size (num_components * num_global_dof_holder)")

    .def(
      "get_owned",
      [](const FieldLayout& self) {
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
      [](const FieldLayout& self) {
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

    .def("is_distributed", &FieldLayout::IsDistributed,
         "Check if the field layout is distributed")

    .def(
      "get_ent_offsets",
      [](const FieldLayout& self) {
        auto offsets = self.GetEntOffsets();
        // Convert std::array to list
        py::list result;
        for (size_t i = 0; i < offsets.size(); ++i) {
          result.append(offsets[i]);
        }
        return result;
      },
      "Get the entity offsets array")

    .def(
      "get_dof_holder_coordinates",
      [](const FieldLayout& self) {
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
      "Get the DOF holder coordinates");
}

} // namespace pcms