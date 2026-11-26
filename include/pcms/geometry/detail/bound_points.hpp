#ifndef PCMS_GEOMETRY_DETAIL_BOUND_POINTS_HPP
#define PCMS_GEOMETRY_DETAIL_BOUND_POINTS_HPP

#include <Kokkos_Core.hpp>
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <cstddef>
#include <string>
#include <type_traits>

namespace pcms::detail
{

/// @brief Makes `storage` hold at least `rows` entries, reallocating it,
///        uninitialized, when it does not; the first `rows` are then the ones
///        in use.
/// @param label label of a reallocated view
template <typename T>
void ReserveRows(Kokkos::View<T*, DeviceMemorySpace>& storage,
                 const std::string& label, std::size_t rows)
{
  if (storage.extent(0) < rows) {
    storage = Kokkos::View<T*, DeviceMemorySpace>(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, label), rows);
  }
}

/// @brief Makes `storage` hold at least `rows` rows of `cols` columns, as the
///        rank-1 overload.
///
/// The first rows of a taller view are a padded view on a LayoutLeft memory
/// space (the default on GPUs), which the conversions to mdspan do not
/// honour, so there a view is reused only for as many rows as it has.
template <typename T>
void ReserveRows(Kokkos::View<T**, DeviceMemorySpace>& storage,
                 const std::string& label, std::size_t rows, std::size_t cols)
{
  constexpr bool rows_contiguous =
    std::is_same_v<typename Kokkos::View<T**, DeviceMemorySpace>::array_layout,
                   Kokkos::LayoutRight>;
  if (storage.extent(0) < rows || storage.extent(1) != cols ||
      (!rows_contiguous && storage.extent(0) != rows)) {
    storage = Kokkos::View<T**, DeviceMemorySpace>(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, label), rows, cols);
  }
}

} // namespace pcms::detail

#endif // PCMS_GEOMETRY_DETAIL_BOUND_POINTS_HPP
