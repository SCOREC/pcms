#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pcms/adapter/uniform_grid/uniform_grid_field_layout.h"
#include "pcms/field.h"
#include "pcms/field_layout.h"
#include "pcms/uniform_grid.h"
#include "numpy_array_transform.h"

namespace py = pybind11;

namespace pcms
{

void bind_uniform_grid_field_layout_module(py::module& m)
{
  // Bind UniformGrid structure for 2D
  py::class_<UniformGrid<2>>(m, "UniformGrid2D")
    .def(py::init<>())
    .def_readwrite("edge_length", &UniformGrid<2>::edge_length,
                   "Edge length in each dimension")
    .def_readwrite("bot_left", &UniformGrid<2>::bot_left,
                   "Bottom-left corner coordinates")
    .def_readwrite("divisions", &UniformGrid<2>::divisions,
                   "Number of divisions in each dimension")
    .def("get_num_cells", &UniformGrid<2>::GetNumCells,
         "Get the total number of cells in the grid")
    .def(
      "closest_cell_id",
      [](const UniformGrid<2>& self, py::array_t<Real> point) {
        if (point.size() != 2) {
          throw std::runtime_error("Point must have 2 coordinates for 2D grid");
        }
        auto buf = point.request();
        Real* ptr = static_cast<Real*>(buf.ptr);
        Omega_h::Vector<2> omega_point({ptr[0], ptr[1]});
        return self.ClosestCellID(omega_point);
      },
      py::arg("point"),
      "Get the cell ID that contains or is closest to the given point")
    .def("get_cell_bbox", &UniformGrid<2>::GetCellBBOX, py::arg("cell_index"),
         "Get the bounding box of a cell");

  // Bind UniformGrid structure for 3D
  py::class_<UniformGrid<3>>(m, "UniformGrid3D")
    .def(py::init<>())
    .def_readwrite("edge_length", &UniformGrid<3>::edge_length,
                   "Edge length in each dimension")
    .def_readwrite("bot_left", &UniformGrid<3>::bot_left,
                   "Bottom-left corner coordinates")
    .def_readwrite("divisions", &UniformGrid<3>::divisions,
                   "Number of divisions in each dimension")
    .def("get_num_cells", &UniformGrid<3>::GetNumCells,
         "Get the total number of cells in the grid")
    .def(
      "closest_cell_id",
      [](const UniformGrid<3>& self, py::array_t<Real> point) {
        if (point.size() != 3) {
          throw std::runtime_error("Point must have 3 coordinates for 3D grid");
        }
        auto buf = point.request();
        Real* ptr = static_cast<Real*>(buf.ptr);
        Omega_h::Vector<3> omega_point({ptr[0], ptr[1], ptr[2]});
        return self.ClosestCellID(omega_point);
      },
      py::arg("point"),
      "Get the cell ID that contains or is closest to the given point")
    .def("get_cell_bbox", &UniformGrid<3>::GetCellBBOX, py::arg("cell_index"),
         "Get the bounding box of a cell");

  // Bind the UniformGridFieldLayout class for 2D
  py::class_<UniformGridFieldLayout<2>, FieldLayout,
             std::shared_ptr<UniformGridFieldLayout<2>>>(
    m, "UniformGridFieldLayout2D")
    .def(py::init<UniformGrid<2>&, int, CoordinateSystem>(), py::arg("grid"),
         py::arg("num_components"), py::arg("coordinate_system"),
         "Constructor for UniformGridFieldLayout2D")

    .def(
      "create_field",
      [](UniformGridFieldLayout<2>& self) {
        return std::shared_ptr<FieldT<Real>>(self.CreateFieldReal());
      },
      "Create a field with this layout")

    .def("get_num_components", &UniformGridFieldLayout<2>::GetNumComponents,
         "Get the number of components in the field")

    .def("get_num_owned_dof_holder",
         &UniformGridFieldLayout<2>::GetNumOwnedDofHolder,
         "Get the number of owned DOF holders")

    .def("get_num_global_dof_holder",
         &UniformGridFieldLayout<2>::GetNumGlobalDofHolder,
         "Get the number of global DOF holders")

    .def(
      "get_owned",
      [](const UniformGridFieldLayout<2>& self) {
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
      [](const UniformGridFieldLayout<2>& self) {
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
      [](const UniformGridFieldLayout<2>& self) {
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

    .def("is_distributed", &UniformGridFieldLayout<2>::IsDistributed,
         "Check if the field layout is distributed")

    .def("get_grid", &UniformGridFieldLayout<2>::GetGrid,
         py::return_value_policy::reference, "Get the underlying uniform grid")

    .def("get_num_cells", &UniformGridFieldLayout<2>::GetNumCells,
         "Get the number of cells in the grid")

    .def("get_num_vertices", &UniformGridFieldLayout<2>::GetNumVertices,
         "Get the number of vertices in the grid");

  // Bind the UniformGridFieldLayout class for 3D
  py::class_<UniformGridFieldLayout<3>, FieldLayout,
             std::shared_ptr<UniformGridFieldLayout<3>>>(
    m, "UniformGridFieldLayout3D")
    .def(py::init<UniformGrid<3>&, int, CoordinateSystem>(), py::arg("grid"),
         py::arg("num_components"), py::arg("coordinate_system"),
         "Constructor for UniformGridFieldLayout3D")

    .def(
      "create_field",
      [](UniformGridFieldLayout<3>& self) {
        return std::shared_ptr<FieldT<Real>>(self.CreateFieldReal());
      },
      "Create a field with this layout")

    .def("get_num_components", &UniformGridFieldLayout<3>::GetNumComponents,
         "Get the number of components in the field")

    .def("get_num_owned_dof_holder",
         &UniformGridFieldLayout<3>::GetNumOwnedDofHolder,
         "Get the number of owned DOF holders")

    .def("get_num_global_dof_holder",
         &UniformGridFieldLayout<3>::GetNumGlobalDofHolder,
         "Get the number of global DOF holders")

    .def(
      "get_owned",
      [](const UniformGridFieldLayout<3>& self) {
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
      [](const UniformGridFieldLayout<3>& self) {
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
      [](const UniformGridFieldLayout<3>& self) {
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

    .def("is_distributed", &UniformGridFieldLayout<3>::IsDistributed,
         "Check if the field layout is distributed")

    .def("get_grid", &UniformGridFieldLayout<3>::GetGrid,
         py::return_value_policy::reference, "Get the underlying uniform grid")

    .def("get_num_cells", &UniformGridFieldLayout<3>::GetNumCells,
         "Get the number of cells in the grid")

    .def("get_num_vertices", &UniformGridFieldLayout<3>::GetNumVertices,
         "Get the number of vertices in the grid");
}

} // namespace pcms
