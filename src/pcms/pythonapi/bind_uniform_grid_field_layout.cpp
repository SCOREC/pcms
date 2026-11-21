#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pcms/pythonapi/numpy_array_transform.hpp"
#include "pcms/backends/uniform_grid/layout.hpp"
#include "pcms/utility/uniform_grid.hpp"
#include <Kokkos_Core.hpp>

namespace py = pybind11;

namespace pcms
{

namespace
{

template <unsigned Dim>
void BindUniformGrid(py::module& m)
{
  using Grid = UniformGrid<Dim>;
  py::class_<Grid>(m, ("UniformGrid" + std::to_string(Dim) + "D").c_str())
    .def(py::init<>())
    .def_readwrite("edge_length", &Grid::edge_length,
                   "Edge length in each dimension")
    .def_readwrite("bot_left", &Grid::bot_left,
                   "Bottom-left corner coordinates")
    .def_readwrite("divisions", &Grid::divisions,
                   "Number of divisions in each dimension")
    .def("get_num_cells", &Grid::GetNumCells,
         "Get the total number of cells in the grid")
    .def(
      "closest_cell_id",
      [](const Grid& self, const NumpyInput<Real>& point_arr) {
        auto point = NumpyToHostView<1, Real>(point_arr);
        if (point.extent(0) != Dim) {
          throw std::invalid_argument("Point must have " + std::to_string(Dim) +
                                      " coordinates for " +
                                      std::to_string(Dim) + "D grid");
        }
        Kokkos::Array<Real, Dim> cell_point;
        for (unsigned d = 0; d < Dim; ++d) {
          cell_point[d] = point(d);
        }
        return self.ClosestCellID(cell_point);
      },
      py::arg("point"),
      "Get the cell ID that contains or is closest to the given point")
    .def("get_cell_bbox", &Grid::GetCellBBOX, py::arg("cell_index"),
         "Get the bounding box of a cell");
}

} // namespace

void bind_uniform_grid_field_layout_module(py::module& m)
{
  // Bind UniformGrid setup types. The corresponding layout types remain
  // internal to the Python API; users construct function spaces from grids
  // and work with Field objects.
  BindUniformGrid<1>(m);
  BindUniformGrid<2>(m);
  BindUniformGrid<3>(m);
}

} // namespace pcms
