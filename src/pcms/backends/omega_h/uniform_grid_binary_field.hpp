#ifndef PCMS_UNIFORM_GRID_BINARY_FIELD_H
#define PCMS_UNIFORM_GRID_BINARY_FIELD_H

#include "pcms/backends/uniform_grid/layout.hpp"
#include "pcms/field/simple_field_data.hpp"
#include "pcms/field/function_space/lagrange.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/localization/point_search.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"

#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>
#include <array>
#include <memory>
#include <utility>
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/backends/omega_h/omega_h.hpp"

namespace pcms::backend::omega_h
{

/**
 * \brief Create a binary (inside/outside) mask field on a uniform grid.
 *
 * Each DOF in the returned order-1 (vertex-centered) field is set to 1.0
 * if the vertex lies inside the mesh and 0.0 otherwise. A single point
 * search pass over the grid vertices is performed; no intermediate source
 * field is required.
 *
 * \tparam Dim  Spatial dimension (2 or 3).
 * \param mesh       The Omega_h mesh defining the domain.
 * \param grid       Uniform grid whose vertices are tested.
 * \return Pair of (layout, field) with binary mask values.
 */
template <unsigned Dim>
std::pair<std::shared_ptr<const pcms::backend::uniform_grid::UniformGridFieldLayout<Dim>>, Field<Real>>
CreateUniformGridBinaryField(Omega_h::Mesh& mesh, const UniformGrid<Dim>& grid)
{
  auto function_space = pcms::field::LagrangeFunctionSpace::FromUniformGrid(
    grid, 1, csys::Cartesian::Deferred());
  auto layout = std::dynamic_pointer_cast<const pcms::backend::uniform_grid::UniformGridFieldLayout<Dim>>(
    function_space->GetLayout());
  PCMS_ALWAYS_ASSERT(layout != nullptr);
  Field<Real> field = function_space->template CreateFunction<Real>("");

  auto coord_view = layout->GetDOFHolderCoordinates();
  auto coords = coord_view.GetValues();
  LO n = layout->GetNumOwnedDofHolder();

  Kokkos::View<Real* [Dim]> coords_d("coords_d", n);
  Kokkos::parallel_for(
    "CopyCoords", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(int i) {
      for (int d = 0; d < Dim; ++d) {
        coords_d(i, d) = coords(i, d);
      }
    });

  // Run point-in-mesh search
  Kokkos::View<typename pcms::localization::PointLocalizationSearch<Dim>::Result*> results_d;
  if constexpr (Dim == 2) {
    pcms::localization::GridPointSearch2D search(mesh, grid.divisions[0], grid.divisions[1]);
    results_d = search(coords_d);
  } else {
    pcms::localization::GridPointSearch3D search(mesh, grid.divisions[0], grid.divisions[1],
                             grid.divisions[2]);
    results_d = search(coords_d);
  }
  auto results_h =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, results_d);

  Kokkos::View<Real*, HostMemorySpace> data("binary_mask", n);
  for (LO i = 0; i < n; ++i)
    data(i) = (results_h(i).element_id >= 0) ? 1.0 : 0.0;

  field.SetDOFHolderDataUncheckedHost(
    Rank2View<const Real, HostMemorySpace>(data.data(), n, 1));

  return {std::move(layout), std::move(field)};
}

/**
 * \brief Convenience overload: derives the grid from the mesh bounding box.
 */
template <unsigned Dim>
std::pair<std::shared_ptr<const pcms::backend::uniform_grid::UniformGridFieldLayout<Dim>>, Field<Real>>
CreateUniformGridBinaryField(Omega_h::Mesh& mesh,
                             const std::array<LO, Dim>& divisions)
{
  return CreateUniformGridBinaryField<Dim>(
    mesh, CreateUniformGridFromMesh<Dim>(mesh, divisions));
}

} // namespace pcms::backend::omega_h

#endif // PCMS_UNIFORM_GRID_BINARY_FIELD_H
