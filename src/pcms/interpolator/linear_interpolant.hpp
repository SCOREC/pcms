
#ifndef INTERPOLANT_HPP
#define INTERPOLANT_HPP

#include <iostream>
#include "multidimarray.hpp"
#define MAX_DIM 10

KOKKOS_INLINE_FUNCTION
void find_indices(const IntVecView& num_bins, const RealVecView& range,
                  const RealVecView& point, int* indices)
{
  int dim = point.extent(0);
  // IntVecView indices("parametric coordinate", dim);
  for (int i = 0; i < dim; ++i) {
    int id = i * 2;
    double length = range(id + 1) - range(id);
    double dlen = length / num_bins(i);
    indices[i] = ((point(i) - range(id)) / dlen);
  }
}

KOKKOS_INLINE_FUNCTION
void find_limits(const IntVecView& num_bins, const RealVecView& range,
                 const int* indices, double* limits)
{
  int dim = num_bins.extent(0);
  // RealVecView limits("limits", 2*dim);
  for (int i = 0; i < dim; ++i) {
    int ptr = i * 2;
    double dlen = (range(ptr + 1) - range(ptr)) / num_bins(i);
    limits[ptr] = indices[i] * dlen;
    limits[ptr + 1] = limits[ptr] + dlen;
  }
}

KOKKOS_INLINE_FUNCTION
void evaluateParametricCoord(const RealVecView& point, const double* limits,
                             double* parametric_coord)
{
  int dim = point.extent(0);
  //  RealVecView parametric_coord("parametric coordinate", dim);
  for (int i = 0; i < dim; ++i) {
    int index = i * 2;
    parametric_coord[i] =
      (point(i) - limits[index]) / (limits[index + 1] - limits[index]);
  }
}

KOKKOS_INLINE_FUNCTION
void basis_function(const RealVecView& parametric_coord,
                    double* linear_basis_each_dir)
{
  int dim = parametric_coord.extent(0);
  // RealVecView linear_basis_each_dir("basis function each direction", 2*dim);
  for (int i = 0; i < dim; ++i) {
    int index = i * 2;

    linear_basis_each_dir[index] = 1 - parametric_coord(i);

    linear_basis_each_dir[index + 1] = parametric_coord(i);
  }
}

KOKKOS_INLINE_FUNCTION
double linear_interpolant(const IntVecView& dimensions,
                          const double* linear_basis_each_dir,
                          const IntVecView& indices, const RealVecView& values)
{
  int dim = dimensions.extent(0);
  double sum = 0;
  int ids[MAX_DIM];

  int array_size = 1 << dim;
  for (int i = 0; i < array_size; ++i) {
    double temp = 1.0;
    for (int j = 0; j < dim; ++j) {
      int index = 2 * j;
      int id_left = indices(j);
      int id_right = id_left + 1;
      if (i & (1 << j)) {
        temp *= linear_basis_each_dir[index + 1];
        ids[j] = id_right;
      } else {
        temp *= linear_basis_each_dir[index];
        ids[j] = id_left;
      }
    }

    int idx = calculateIndex(dimensions, ids);
    double corner_values = values(idx);

    sum += temp * corner_values;
  }
  return sum;
}

class RegularGridInterpolator
{
private:
  const RealMatView parametric_coords;
  const RealVecView values;
  const IntMatView indices;
  const IntVecView dimensions;

public:
  RegularGridInterpolator(const RealMatView& parametric_coords_,
                          const RealVecView& values_,
                          const IntMatView& indices_,
                          const IntVecView& dimensions_)
    : parametric_coords(parametric_coords_),
      values(values_),
      indices(indices_),
      dimensions(dimensions_){};

  RealVecView linear_interpolation()
  {
    int dim = dimensions.extent(0);
    int N = parametric_coords.extent(0);
    RealVecView interpolated_values("approximated values", N);
    auto parametric_coords_ = parametric_coords;
    auto dimensions_ = dimensions;
    auto values_ = values;
    auto indices_ = indices;
    Kokkos::parallel_for(
      "linear interpolation function", N, KOKKOS_LAMBDA(int j) {
        double linear_basis_each_dir[MAX_DIM] = {0.0};
        auto parametric_coord_each_point =
          Kokkos::subview(parametric_coords_, j, Kokkos::ALL());
        auto index_each_point = Kokkos::subview(indices_, j, Kokkos::ALL());
        basis_function(parametric_coord_each_point, linear_basis_each_dir);
        auto approx_value = linear_interpolant(
          dimensions_, linear_basis_each_dir, index_each_point, values_);
        interpolated_values(j) = approx_value;
      });

    return interpolated_values;
  }
};

struct Result
{
  IntMatView indices_pts;
  RealMatView parametric_coords;
};

Result parametric_indices(const RealMatView& points, const IntVecView& num_bins,
                          const RealVecView& range)
{

  int dim = num_bins.extent(0);
  Result result;
  int N = points.extent(0);
  result.indices_pts = IntMatView("indices", N, dim);
  result.parametric_coords = RealMatView("parametric_coordinates", N, dim);

  Kokkos::parallel_for(
    N, KOKKOS_LAMBDA(int j) {
      int indcs[MAX_DIM] = {0};
      double limits[MAX_DIM] = {0.0};
      double parametric_coord_each[MAX_DIM] = {0.0};
      auto point_coord = Kokkos::subview(points, j, Kokkos::ALL());
      find_indices(num_bins, range, point_coord, indcs);
      find_limits(num_bins, range, indcs, limits);
      evaluateParametricCoord(point_coord, limits, parametric_coord_each);
      for (int i = 0; i < dim; ++i) {
        result.indices_pts(j, i) = indcs[i];
        result.parametric_coords(j, i) = parametric_coord_each[i];
      }
    });

  return result;
}

KOKKOS_INLINE_FUNCTION
double test_function(double* coord)
{
  double fun_value = 0;
  int dim = 5;
  for (int i = 0; i < dim; ++i) {
    fun_value += coord[i];
  }
  return fun_value;
}

#endif
