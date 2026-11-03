#ifndef PCMS_PYTHONAPI_NUMPY_VIEW_H
#define PCMS_PYTHONAPI_NUMPY_VIEW_H

// Conversions between numpy arrays and pcms mdspan views / Kokkos views.
//
// Array semantics for python conversions
//
// Arrays returned for pcms-owned data (ViewToNumpy, ToNumpy of Omega_h data):
//   - are host numpy arrays shaped like the C++ view, and read-only
//   - are valid until the object they came from is modified
//   - never dangle: they keep their owner alive or are copies. On host builds
//     they share memory with pcms (zero-copy); otherwise they are copies
// Array inputs (NumpyInput):
//   - accept any array-like convertible to the element type
//   - may be read in place when only used during the call (NumpyToHostView)
//   - are copied when pcms keeps the data (CopyNumpyToKokkosView), so later
//     changes to the caller's array have no effect
// Output buffers (NumpyOutputToHostView):
//   - must already be writeable, C-contiguous arrays of the element type
// Per-entity data (field DOF values, evaluation results):
//   - is shaped (n, num_components)
//   - returned arrays always have the component axis, so their rank does not
//     depend on the number of components
//   - inputs and output buffers may omit the component axis when
//     num_components is 1
//
// TODO: zero-copy for device memory via DLPack, with the same semantics.

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
// column slice or a transpose) or has a different dtype. Inside the binding the
// parameter refers to the converted array, so it keeps the memory seen by
// NumpyToHostView alive for the duration of the call.
template <typename T>
using NumpyInput = py::array_t<T, py::array::c_style | py::array::forcecast>;

// True when memory in MemorySpace can be read directly by host code.
template <typename MemorySpace>
inline constexpr bool is_host_accessible_v =
  Kokkos::SpaceAccessibility<Kokkos::HostSpace, MemorySpace>::accessible;

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

template <typename ElementType, typename Extents, typename LayoutPolicy,
          typename AccessorPolicy>
std::vector<py::ssize_t> NumpyShape(
  const Kokkos::mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>&
    view)
{
  std::vector<py::ssize_t> shape(Extents::rank());
  for (std::size_t i = 0; i < Extents::rank(); ++i) {
    shape[i] = static_cast<py::ssize_t>(view.extent(i));
  }
  return shape;
}

inline void SetReadOnly(py::array& arr)
{
  py::setattr(arr.attr("flags"), "writeable", py::bool_(false));
}

inline std::string ShapeString(const py::array& arr)
{
  std::string result = "(";
  for (py::ssize_t i = 0; i < arr.ndim(); ++i) {
    result += (i > 0 ? ", " : "") + std::to_string(arr.shape(i));
  }
  return result + (arr.ndim() == 1 ? ",)" : ")");
}

// Extents of per-entity data shaped (n, num_components). A 1D array of length
// n is read as (n, 1) when num_components is 1.
inline std::array<LO, 2> ComponentExtents(const py::array& arr,
                                          LO num_components)
{
  if (arr.ndim() == 1 && num_components == 1) {
    return {static_cast<LO>(arr.shape(0)), 1};
  }
  if (arr.ndim() == 2 && arr.shape(1) == num_components) {
    return {static_cast<LO>(arr.shape(0)), num_components};
  }
  throw std::invalid_argument(
    "expected an array of shape (n, " + std::to_string(num_components) + ")" +
    (num_components == 1 ? " or (n,)" : "") + ", got " + ShapeString(arr));
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
namespace detail
{
template <typename T>
T* OutputBufferData(py::array& arr)
{
  if (!py::isinstance<py::array_t<T, py::array::c_style>>(arr)) {
    throw std::invalid_argument(
      "output array must be C-contiguous with dtype " +
      std::string(py::str(py::dtype::of<T>())));
  }
  if (!arr.writeable()) {
    throw std::invalid_argument("output array must be writeable");
  }
  return static_cast<T*>(arr.mutable_data());
}
} // namespace detail

template <std::size_t Rank, typename T>
View<Rank, T, HostMemorySpace, Kokkos::layout_right> NumpyOutputToHostView(
  py::array& arr)
{
  return View<Rank, T, HostMemorySpace, Kokkos::layout_right>(
    detail::OutputBufferData<T>(arr), detail::NumpyExtents<Rank>(arr));
}

// NumpyToHostView for per-entity data of shape (n, num_components), where the
// component axis may be omitted when num_components is 1.
template <typename T>
Rank2View<const T, HostMemorySpace, Kokkos::layout_right>
NumpyToHostComponentView(const NumpyInput<T>& arr, LO num_components)
{
  return Rank2View<const T, HostMemorySpace, Kokkos::layout_right>(
    arr.data(), detail::ComponentExtents(arr, num_components));
}

// NumpyOutputToHostView for per-entity data of shape (n, num_components),
// where the component axis may be omitted when num_components is 1.
template <typename T>
Rank2View<T, HostMemorySpace, Kokkos::layout_right>
NumpyOutputToHostComponentView(py::array& arr, LO num_components)
{
  return Rank2View<T, HostMemorySpace, Kokkos::layout_right>(
    detail::OutputBufferData<T>(arr),
    detail::ComponentExtents(arr, num_components));
}

// Copy a numpy input into a new, owning Kokkos::View in MemorySpace (with that
// space's default layout).
template <std::size_t Rank, typename T,
          typename MemorySpace = DeviceMemorySpace>
Kokkos::View<detail::kokkos_data_type_t<T, Rank>, MemorySpace>
CopyNumpyToKokkosView(const NumpyInput<T>& arr,
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
  py::array_t<T> result(detail::NumpyShape(view));
  if (result.size() == 0) {
    return result;
  }
  auto dest = MakeUnmanagedKokkosView(
    View<rank, T, HostMemorySpace, Kokkos::layout_right>(
      result.mutable_data(), detail::NumpyExtents<rank>(result)));
  DeepCopyMismatchLayouts(dest, MakeUnmanagedKokkosView(view));
  return result;
}

// numpy array for an mdspan. When the view's memory is host-accessible the
// array aliases it (zero-copy, with the view's strides) and holds a reference
// to `owner`, which must keep the memory alive. Otherwise the data is copied.
// Either way the array is read-only when ElementType is const, so callers see
// the same behavior on host and device builds.
template <typename ElementType, typename Extents, typename LayoutPolicy,
          typename AccessorPolicy>
py::array_t<std::remove_const_t<ElementType>> ViewToNumpy(
  const Kokkos::mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>&
    view,
  py::handle owner)
{
  using T = std::remove_const_t<ElementType>;
  using MemorySpace = detail::memory_space_selector_t<AccessorPolicy>;
  if constexpr (is_host_accessible_v<MemorySpace>) {
    constexpr std::size_t rank = Extents::rank();
    std::vector<py::ssize_t> strides(rank);
    for (std::size_t i = 0; i < rank; ++i) {
      strides[i] = static_cast<py::ssize_t>(view.stride(i) * sizeof(T));
    }
    py::array_t<T> result(detail::NumpyShape(view), strides, view.data_handle(),
                          owner);
    if constexpr (std::is_const_v<ElementType>) {
      detail::SetReadOnly(result);
    }
    return result;
  } else {
    auto result = ToNumpy(view);
    if constexpr (std::is_const_v<ElementType>) {
      detail::SetReadOnly(result);
    }
    return result;
  }
}

} // namespace pcms

#endif // PCMS_PYTHONAPI_NUMPY_VIEW_H
