#ifndef PCMS_COUPLING_UNIFORM_GRID_H
#define PCMS_COUPLING_UNIFORM_GRID_H
#include "pcms/bounding_box.h"
#include "Omega_h_vector.hpp"
#include "Omega_h_bbox.hpp"
#include "Omega_h_mesh.hpp"
#include <numeric>
namespace pcms
{

template <unsigned dim>
struct UniformGrid
{
  // Make private?
  std::array<Real, dim> edge_length;
  std::array<Real, dim> bot_left;
  std::array<LO, dim> divisions;

public:
  [[nodiscard]] LO GetNumCells() const
  {
    return std::accumulate(divisions.begin(), divisions.end(), 1,
                           std::multiplies<LO>{});
  }
  /// return the grid cell ID that the input point is inside or closest to if
  /// the point lies outside
  // take the view as a template because it might be a subview type
  // template <typename T>
  //[[nodiscard]] KOKKOS_INLINE_FUNCTION LO ClosestCellID(const T& point) const
  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO
  ClosestCellID(const Omega_h::Vector<dim>& point) const
  {
    std::array<Real, dim> distance_within_grid;

    for (size_t i = 0; i < dim; ++i) {
      distance_within_grid[i] = point[i] - bot_left[i];
    }

    std::array<LO, dim> indexes;

    for (auto& index : indexes) {
      index = -1;
    }

    for (int i = 0; i < dim; ++i) {
      auto index = static_cast<LO>(
        std::floor(distance_within_grid[i] * divisions[i] / edge_length[i]));
      indexes[i] = Kokkos::clamp(index, 0, divisions[i] - 1);
    }
    // note that the indexes refer to row/columns which have the opposite order
    // of the coordinates i.e. x,y
    reverse(indexes);
    return GetCellIndex(indexes);
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION AABBox<dim> GetCellBBOX(LO idx) const
  {
    auto index = GetDimensionedIndex(idx);
    reverse(index);

    std::array<Real, dim> half_width, center;

    for (size_t i = 0; i < dim; ++i) {
      half_width[i] = edge_length[i] / divisions[i] / 2;
    }

    for (size_t i = 0; i < dim; ++i) {
      center[i] = (2.0 * index[i] + 1.0) * half_width[i] + bot_left[i];
    }

    return {.center = center, .half_width = half_width};
  }

  /// Check if a point is within the grid bounds
  [[nodiscard]] KOKKOS_INLINE_FUNCTION bool IsPointInBounds(
    const Omega_h::Vector<dim>& point) const
  {
    for (size_t i = 0; i < dim; ++i) {
      Real coord = point[i];
      Real grid_min = bot_left[i];
      Real grid_max = bot_left[i] + edge_length[i];
      if (coord < grid_min || coord > grid_max) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION std::array<LO, dim> GetDimensionedIndex(
    LO idx) const
  {
    LO stride = 1;
    for (std::size_t i = 0; i < divisions.size() - 1; ++i) {
      stride *= divisions[i];
    }
    std::array<LO, dim> result;

    for (size_t i = 0; i < dim; ++i) {
      result[i] = idx / stride;
      idx -= result[i] * stride;
      stride /= divisions[i];
    }

    return result;
  }

  [[nodiscard]] KOKKOS_INLINE_FUNCTION LO
  GetCellIndex(std::array<LO, dim> dimensionedIndex) const
  {
    // note that the indexes refer to row/columns which have the opposite order
    // of the coordinates i.e. x,y
    reverse(dimensionedIndex);

    LO idx = 0;
    LO stride = 1;

    for (size_t i = 0; i < dim; ++i) {
      idx += dimensionedIndex[i] * stride;
      stride *= divisions[i];
    }

    return idx;
  }

private:
  template <typename T, std::size_t N>
  KOKKOS_INLINE_FUNCTION static void reverse(std::array<T, N>& arr)
  {
    for (size_t i = 0, j = arr.size() - 1; i < j; ++i, --j) {
      auto temp = arr[i];
      arr[i] = arr[j];
      arr[j] = temp;
    }
  }
};

using Uniform2DGrid = UniformGrid<2>;
using Uniform3DGrid = UniformGrid<3>;

/**
 * \brief Create a uniform grid layout from an Omega_h mesh as a bounding box
 *
 * This function computes the bounding box of the given mesh and creates a
 * uniform Cartesian grid that covers the entire mesh domain. The grid cells
 * are uniformly distributed across each dimension.
 *
 * \tparam dim Spatial dimension of the mesh (2 or 3)
 * \param mesh The Omega_h mesh to create the grid from
 * \param divisions Array specifying the number of cells in each dimension
 * \return UniformGrid<dim> A uniform grid covering the mesh bounding box
 *
 */
template <unsigned dim = 2>
UniformGrid<dim> CreateUniformGridFromMesh(Omega_h::Mesh& mesh,
                                           const std::array<LO, dim>& divisions)
{
  // Get the bounding box of the mesh
  auto bbox = Omega_h::get_bounding_box<dim>(&mesh);

  // Calculate edge lengths and bottom-left corner
  std::array<Real, dim> edge_length;
  std::array<Real, dim> bot_left;

  for (unsigned i = 0; i < dim; ++i) {
    bot_left[i] = bbox.min[i];
    edge_length[i] = bbox.max[i] - bbox.min[i];
  }

  return UniformGrid<dim>{
    .edge_length = edge_length, .bot_left = bot_left, .divisions = divisions};
}

/**
 * \brief Create a uniform grid with equal divisions in all dimensions
 *
 * This is a convenience function that creates a uniform grid with the same
 * number of cells in each dimension.
 *
 * \tparam dim Spatial dimension of the mesh (2 or 3)
 * \param mesh The Omega_h mesh to create the grid from
 * \param cells_per_dim Number of cells per dimension (same for all dimensions)
 * \return UniformGrid<dim> A uniform grid covering the mesh bounding box
 *
 */
template <unsigned dim = 2>
UniformGrid<dim> CreateUniformGridFromMesh(Omega_h::Mesh& mesh,
                                           LO cells_per_dim)
{
  std::array<LO, dim> divisions;
  divisions.fill(cells_per_dim);
  return CreateUniformGridFromMesh<dim>(mesh, divisions);
}

} // namespace pcms

#endif // PCMS_COUPLING_UNIFORM_GRID_H
