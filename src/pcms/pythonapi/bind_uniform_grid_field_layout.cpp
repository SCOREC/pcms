#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "numpy_array_transform.h"
#include "pcms/field/layout/uniform_grid.h"
#include "pcms/utility/uniform_grid.h"

namespace py = pybind11;

namespace pcms
{

void bind_uniform_grid_field_layout_module(py::module& m)
{
  // Bind UniformGrid setup types. The corresponding layout types remain
  // internal to the Python API; users construct function spaces from grids
  // and work with Field objects.
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
      [](const UniformGrid<2>& self, const NumpyInput<Real>& point_arr) {
        auto point = NumpyToHostView<1, Real>(point_arr);
        if (point.extent(0) != 2) {
          throw std::invalid_argument(
            "Point must have 2 coordinates for 2D grid");
        }
        Omega_h::Vector<2> omega_point({point(0), point(1)});
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
      [](const UniformGrid<3>& self, const NumpyInput<Real>& point_arr) {
        auto point = NumpyToHostView<1, Real>(point_arr);
        if (point.extent(0) != 3) {
          throw std::invalid_argument(
            "Point must have 3 coordinates for 3D grid");
        }
        Omega_h::Vector<3> omega_point({point(0), point(1), point(2)});
        return self.ClosestCellID(omega_point);
      },
      py::arg("point"),
      "Get the cell ID that contains or is closest to the given point")
    .def("get_cell_bbox", &UniformGrid<3>::GetCellBBOX, py::arg("cell_index"),
         "Get the bounding box of a cell");
}

} // namespace pcms
