#ifndef PCMS_PYTHONAPI_NUMPY_VIEW_H
#define PCMS_PYTHONAPI_NUMPY_VIEW_H

// Conversions between numpy arrays and pcms mdspan views / Kokkos views.
//
// numpy array storage is host memory. Inputs are views on the host 
// and copied when they are needed in another memory space.
// outputs are always returned as owning numpy arrays so that
// numpy controls the lifetime

#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <Kokkos_Core.hpp>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"

namespace pcms
{

namespace py = pybind11;

// Parameter type for numpy inputs. pybind11 converts the caller's array to a
// C-contiguous array of T, copying only when the array is strided (e.g. a
// column slice or a transpose) or has a different dtype.
template <typename T>
using NumpyInput = py::array_t<T, py::array::c_style | py::array::forcecast>;

namespace detail
{

template <typename T, std::size_t Rank>
struct KokkosDataType
{
  using type = typename KokkosDataType<T, Rank - 1>::type*;
};
template <typename T>
struct KokkosDataType<T, 0>
{
  using type = T;
};
// KokkosDataType<T, 2>::type is T**
template <typename T, std::size_t Rank>
using kokkos_data_type_t = typename KokkosDataType<T, Rank>::type;

template <typename Layout>
struct KokkosLayout;
template <>
struct KokkosLayout<Kokkos::layout_right>
{
  using type = Kokkos::LayoutRight;
};
template <>
struct KokkosLayout<Kokkos::layout_left>
{
  using type = Kokkos::LayoutLeft;
};

template <std::size_t Rank>
std::array<LO, Rank> NumpyExtents(const py::array& arr)
{
  if (arr.ndim() != static_cast<py::ssize_t>(Rank)) {
    throw std::invalid_argument("expected a " + std::to_string(Rank) +
                                "D array, got a " + std::to_string(arr.ndim()) +
                                "D array");
  }
  std::array<LO, Rank> extents{};
  for (std::size_t i = 0; i < Rank; ++i) {
    extents[i] = static_cast<LO>(arr.shape(i));
  }
  return extents;
}

template <typename ViewType, typename Pointer, typename Extents,
          std::size_t... Is>
ViewType MakeKokkosViewImpl(Pointer ptr, const Extents& extents,
                            std::index_sequence<Is...>)
{
  return ViewType(ptr, extents[Is]...);
}

template <typename ViewType, typename Extents, std::size_t... Is>
ViewType AllocateKokkosViewImpl(const std::string& label,
                                const Extents& extents,
                                std::index_sequence<Is...>)
{
  return ViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, label),
                  extents[Is]...);
}

} // namespace detail

// Unmanaged Kokkos::View over the same memory as an mdspan, so Kokkos
// deep_copy and mirror utilities can be used on it.
template <typename ElementType, typename Extents, typename LayoutPolicy,
          typename AccessorPolicy>
auto MakeUnmanagedKokkosView(
  const Kokkos::mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>&
    view)
{
  constexpr std::size_t rank = Extents::rank();
  using MemorySpace = detail::memory_space_selector_t<AccessorPolicy>;
  using KokkosViewType =
    Kokkos::View<detail::kokkos_data_type_t<ElementType, rank>,
                 typename detail::KokkosLayout<LayoutPolicy>::type, MemorySpace,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  std::array<std::size_t, rank> extents{};
  for (std::size_t i = 0; i < rank; ++i) {
    extents[i] = view.extent(i);
  }
  return detail::MakeKokkosViewImpl<KokkosViewType>(
    view.data_handle(), extents, std::make_index_sequence<rank>{});
}

// Read-only host view of a numpy input, without copying. The view aliases
// `arr`, so `arr` must outlive it.
template <std::size_t Rank, typename T>
View<Rank, const T, HostMemorySpace, Kokkos::layout_right> NumpyToHostView(
  const NumpyInput<T>& arr)
{
  return View<Rank, const T, HostMemorySpace, Kokkos::layout_right>(
    arr.data(), detail::NumpyExtents<Rank>(arr));
}

// Writable host view of a caller-provided numpy output buffer. Results written
// to a converted copy would be lost, so instead of converting, require the
// buffer to already be a writeable, C-contiguous array of T.
template <std::size_t Rank, typename T>
View<Rank, T, HostMemorySpace, Kokkos::layout_right> NumpyOutputToHostView(
  py::array& arr)
{
  if (!py::isinstance<py::array_t<T, py::array::c_style>>(arr)) {
    throw std::invalid_argument(
      "output array must be C-contiguous with dtype " +
      std::string(py::str(py::dtype::of<T>())));
  }
  if (!arr.writeable()) {
    throw std::invalid_argument("output array must be writeable");
  }
  return View<Rank, T, HostMemorySpace, Kokkos::layout_right>(
    static_cast<T*>(arr.mutable_data()), detail::NumpyExtents<Rank>(arr));
}

// Copy a numpy input into a new Kokkos::View in MemorySpace (with that
// space's default layout).
template <std::size_t Rank, typename T,
          typename MemorySpace = DeviceMemorySpace>
Kokkos::View<detail::kokkos_data_type_t<T, Rank>, MemorySpace>
NumpyToKokkosView(const NumpyInput<T>& arr,
                  const std::string& label = "numpy_copy")
{
  using KokkosViewType =
    Kokkos::View<detail::kokkos_data_type_t<T, Rank>, MemorySpace>;
  const auto src = MakeUnmanagedKokkosView(NumpyToHostView<Rank, T>(arr));
  auto dest = detail::AllocateKokkosViewImpl<KokkosViewType>(
    label, detail::NumpyExtents<Rank>(arr), std::make_index_sequence<Rank>{});
  DeepCopyMismatchLayouts(dest, src);
  return dest;
}

// Copy an mdspan (in any memory space, layout_left or layout_right) into a
// new C-ordered numpy array of the same shape.
template <typename ElementType, typename Extents, typename LayoutPolicy,
          typename AccessorPolicy>
py::array_t<std::remove_const_t<ElementType>> ToNumpy(
  const Kokkos::mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>&
    view)
{
  using T = std::remove_const_t<ElementType>;
  constexpr std::size_t rank = Extents::rank();
  static_assert(rank > 0, "ToNumpy requires a view of rank >= 1");
  std::vector<py::ssize_t> shape(rank);
  for (std::size_t i = 0; i < rank; ++i) {
    shape[i] = static_cast<py::ssize_t>(view.extent(i));
  }
  py::array_t<T> result(shape);
  if (result.size() == 0) {
    return result;
  }
  auto dest = MakeUnmanagedKokkosView(
    View<rank, T, HostMemorySpace, Kokkos::layout_right>(
      result.mutable_data(), detail::NumpyExtents<rank>(result)));
  DeepCopyMismatchLayouts(dest, MakeUnmanagedKokkosView(view));
  return result;
}

} // namespace pcms

#endif // PCMS_PYTHONAPI_NUMPY_VIEW_H
