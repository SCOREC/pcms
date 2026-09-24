#ifndef PCMS_UTILITY_OMEGA_H_ARRAY_UTILS_H
#define PCMS_UTILITY_OMEGA_H_ARRAY_UTILS_H

#include <Omega_h_array.hpp>
#include <Omega_h_mesh.hpp>
#include "pcms/utility/memory_spaces.h"

#include <type_traits>

namespace pcms
{

template <typename View2D>
Omega_h::Reals flatten_to_omega_h_reals_host(const View2D& coords,
                                             const char* label = "flat_coords")
{
  const int n = static_cast<int>(coords.extent(0));
  const int dim = static_cast<int>(coords.extent(1));

  Omega_h::HostWrite<Omega_h::Real> flat(n * dim, label);
  for (int i = 0; i < n; ++i)
    for (int d = 0; d < dim; ++d)
      flat[i * dim + d] = coords(i, d);

  return Omega_h::Reals(flat);
}

template <typename View2D>
Omega_h::Reals flatten_to_omega_h_reals(const View2D& coords,
                                        const char* label = "flat_coords")
{
  const int n = static_cast<int>(coords.extent(0));
  const int dim = static_cast<int>(coords.extent(1));

  Omega_h::Write<Omega_h::Real> flat(n * dim, label);
  Kokkos::parallel_for(
    "flatten_to_omega_h_reals", Kokkos::RangePolicy<>(0, n),
    KOKKOS_LAMBDA(int i) {
      for (int d = 0; d < dim; ++d) {
        flat[i * dim + d] = coords(i, d);
      }
    });

  return Omega_h::Reals(flat);
}

// Convert 1D Omega_h coordinate array to 2D Kokkos view
// This copy is needed because Omega_h provides coordinates as a 1D array of
// length nents*dim, which is always layout_right. Directly creating a 2D mdspan
// with a different layout would access the 1D array with the wrong stride and
// lead to incorrect coordinates. By copying to a new 2D view, we ensure correct
// memory access regardless of the layout of the destination view.
template <typename T, typename IntType1, typename IntType2>
inline Kokkos::View<T**, pcms::DeviceMemorySpace> ConvertCoordsTo2D(
  const Omega_h::Read<T>& coords_1d, IntType1 nents, IntType2 dim)
{
  const int n = static_cast<int>(nents);
  const int d = static_cast<int>(dim);
  Kokkos::View<T**, pcms::DeviceMemorySpace> coords_2d("coords_2d", n, d);

  Kokkos::parallel_for(
    "convert_coords_to_2d", n * d, KOKKOS_LAMBDA(int i) {
      int e = i / d;
      int dd = i % d;
      coords_2d(e, dd) = coords_1d[i];
    });

  return coords_2d;
}

// Synchronize a block of mesh entity data with Omega_h's internal arrays.
template <typename T, typename BlockView>
void SynchronizeOmegaHBlock(Omega_h::Mesh& mesh, int dim, int nc,
                            BlockView block)
{
  if constexpr (std::is_same_v<T, float>) {
    Omega_h::Write<Omega_h::Real> arr(static_cast<Omega_h::LO>(block.size()),
                                      "field_data_block");
    Kokkos::deep_copy(arr.view(), block);
    auto synced = mesh.sync_array(dim, Omega_h::Read<Omega_h::Real>(arr), nc);
    Kokkos::deep_copy(block, synced.view());
  } else {
    auto synced =
      mesh.sync_array(dim, Omega_h::Read<T>(Omega_h::Write<T>(block)), nc);
    Kokkos::deep_copy(block, synced.view());
  }
}

} // namespace pcms

#endif // PCMS_UTILITY_OMEGA_H_ARRAY_UTILS_H
