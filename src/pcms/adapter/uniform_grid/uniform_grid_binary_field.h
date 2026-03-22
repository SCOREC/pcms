#ifndef PCMS_UNIFORM_GRID_BINARY_FIELD_H
#define PCMS_UNIFORM_GRID_BINARY_FIELD_H

#include "pcms/adapter/uniform_grid/uniform_grid_field.h"
#include "pcms/adapter/uniform_grid/uniform_grid_field_layout.h"
#include "pcms/localization/point_search.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/types.h"
#include "pcms/utility/uniform_grid.h"

#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>
#include <array>
#include <memory>
#include <utility>

namespace pcms
{

/**
 * \brief Create a binary (inside/outside) mask field on a uniform grid.
 *
 * Each DOF in the returned order-1 (vertex-centered) UniformGridField is set
 * to 1.0 if the vertex lies inside the mesh and 0.0 otherwise.  A single
 * point search pass over the grid vertices is performed; no intermediate
 * source field is required.
 *
 * \tparam Dim  Spatial dimension (2 or 3).
 * \param mesh       The Omega_h mesh defining the domain.
 * \param grid       Uniform grid whose vertices are tested.
 * \return Pair of (layout, field) with binary mask values.
 */
template <unsigned Dim>
std::pair<std::shared_ptr<UniformGridFieldLayout<Dim>>,
          std::unique_ptr<UniformGridField<Dim>>>
CreateUniformGridBinaryField(Omega_h::Mesh& mesh, const UniformGrid<Dim>& grid)
{
  auto layout = std::make_shared<UniformGridFieldLayout<Dim>>(
    grid, 1, CoordinateSystem::Cartesian);
  auto field = std::make_unique<UniformGridField<Dim>>(layout);

  auto coord_view = layout->GetDOFHolderCoordinates();
  auto coords     = coord_view.GetCoordinates();
  LO n            = layout->GetNumOwnedDofHolder();

  // Copy DOF coordinates from host mdspan to a Kokkos device view for search
  Kokkos::View<Real* [Dim]> coords_d("coords_d", n);
  auto coords_h = Kokkos::create_mirror_view(coords_d);
  for (LO i = 0; i < n; ++i)
    for (unsigned d = 0; d < Dim; ++d)
      coords_h(i, d) = coords(i, d);
  Kokkos::deep_copy(coords_d, coords_h);

  // Run point-in-mesh search
  Kokkos::View<typename PointLocalizationSearch<Dim>::Result*> results_d;
  if constexpr (Dim == 2) {
    GridPointSearch2D search(mesh, grid.divisions[0], grid.divisions[1]);
    results_d = search(coords_d);
  } else {
    GridPointSearch3D search(mesh, grid.divisions[0], grid.divisions[1],
                              grid.divisions[2]);
    results_d = search(coords_d);
  }
  auto results_h =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, results_d);

  Kokkos::View<Real*, HostMemorySpace> data("binary_mask", n);
  for (LO i = 0; i < n; ++i)
    data(i) = (results_h(i).element_id >= 0) ? 1.0 : 0.0;

  field->SetDOFHolderData(
    Rank1View<const Real, HostMemorySpace>(data.data(), n));

  return {std::move(layout), std::move(field)};
}

/**
 * \brief Convenience overload: derives the grid from the mesh bounding box.
 */
template <unsigned Dim>
std::pair<std::shared_ptr<UniformGridFieldLayout<Dim>>,
          std::unique_ptr<UniformGridField<Dim>>>
CreateUniformGridBinaryField(Omega_h::Mesh& mesh,
                              const std::array<LO, Dim>& divisions)
{
  return CreateUniformGridBinaryField<Dim>(
    mesh, CreateUniformGridFromMesh<Dim>(mesh, divisions));
}

} // namespace pcms

#endif // PCMS_UNIFORM_GRID_BINARY_FIELD_H
