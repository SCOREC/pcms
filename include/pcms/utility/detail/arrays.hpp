#ifndef PCMS_UTILITY_DETAIL_ARRAYS_HPP
#define PCMS_UTILITY_DETAIL_ARRAYS_HPP

#include "pcms/utility/arrays.hpp"

#include <vector>

namespace pcms::detail
{

template <typename ContainerType, typename ElementType, typename Extents,
          typename LayoutPolicy, typename AccessorPolicy>
auto make_mdspan(const ContainerType& /* unused */)
  -> Kokkos::mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy>
{
  static_assert(detail::dependent_always_false<ContainerType>::type,
                "creating mdspan is not implemented for type");
}

// TODO make_mdspan

template <typename ElementType, typename MemorySpace,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<MemorySpace>>
using Rank3View = View<3, ElementType, MemorySpace, LayoutPolicy>;

template <typename ElementType, typename MemorySpace,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<MemorySpace>>
using Rank4View = View<4, ElementType, MemorySpace, LayoutPolicy>;

// default implementation of make_array_view
template <typename T, typename MemorySpace = detail::memory_space_selector_t<T>,
          typename ElementType = detail::element_type_t<T>>
auto make_array_view(const T& array)
  -> Rank1View<const ElementType, MemorySpace>
{
  using std::data;
  using std::size;
  return Rank1View<const ElementType, MemorySpace>{data(array), size(array)};
}
template <typename T, typename MemorySpace = detail::memory_space_selector_t<T>,
          typename ElementType = detail::element_type_t<T>>
auto make_array_view(T& array) -> Rank1View<ElementType, MemorySpace>
{
  using std::data;
  using std::size;
  return Rank1View<ElementType, MemorySpace>{data(array), size(array)};
}

template <typename T, typename MemorySpace = detail::memory_space_selector_t<T>,
          typename ElementType = detail::element_type_t<T>>
auto make_const_array_view(T& array)
  -> Rank1View<const ElementType, MemorySpace>
{
  using std::data;
  using std::size;
  return Rank1View<const ElementType, MemorySpace>{data(array), size(array)};
}

// Flatten a contiguous (exhaustive) 2D [dof][comp] view into a 1D view over the
// same storage. Intended for boundary code (serialization, element-wise
// comparison) that wants a flat [num_dof_holder * num_components] view rather
// than indexing pointers directly.
template <typename ElementType, typename MemorySpace, typename LayoutPolicy>
Rank1View<ElementType, MemorySpace> FlattenToRank1View(
  Rank2View<ElementType, MemorySpace, LayoutPolicy> view)
{
  return Rank1View<ElementType, MemorySpace>(view.data_handle(),
                                             static_cast<LO>(view.size()));
}

// utility function to fill a view with sequentially increasing values
template <typename T>
void iota_view(Kokkos::View<T*> view, T start = 0)
{
  Kokkos::parallel_for(
    "iota_view", view.extent(0), KOKKOS_LAMBDA(LO i) { view[i] = start + i; });
}

// utility function to copy device Rank1View to device Kokkos View
template <typename T>
void CopyDeviceRank1ViewToDeviceView(Kokkos::View<T*, DeviceMemorySpace> dest,
                                     Rank1View<const T, DeviceMemorySpace> src)
{
  if (dest.extent(0) != src.size()) {
    throw pcms_error("CopyDeviceRank1ViewToDeviceView: size mismatch");
  }
  Kokkos::parallel_for(
    "CopyDeviceRank1ViewToDeviceView", dest.extent(0),
    KOKKOS_LAMBDA(LO i) { dest(i) = src(i); });
}

// utility function to copy device Rank1View to host Kokkos View
template <typename T>
void CopyDeviceRank1ViewToHostView(Kokkos::View<T*, HostMemorySpace> dest,
                                   Rank1View<const T, DeviceMemorySpace> src)
{
  if (dest.extent(0) != src.size()) {
    throw pcms_error("CopyDeviceRank1ViewToHostView: size mismatch");
  }
  Kokkos::View<T*, DeviceMemorySpace> src_tmp(
    "CopyDeviceRank1ViewToHostView_tmp", src.size());
  Kokkos::parallel_for(
    "CopyDeviceRank1ViewToHostView", src.size(),
    KOKKOS_LAMBDA(LO i) { src_tmp(i) = src(i); });
  Kokkos::deep_copy(dest, src_tmp);
}

template <typename T>
void CopyHostRank1ViewToDeviceView(Kokkos::View<T*, DeviceMemorySpace> dest,
                                   Rank1View<const T, HostMemorySpace> src)
{
  if (dest.extent(0) != src.size()) {
    throw pcms_error("CopyHostRank1ViewToDeviceView: size mismatch");
  }
  Kokkos::View<T*, HostMemorySpace> src_tmp("CopyHostRank1ViewToDeviceView_tmp",
                                            src.size());
  for (size_t i = 0; i < src.size(); ++i) {
    src_tmp(i) = src(i);
  }
  Kokkos::deep_copy(dest, src_tmp);
}

// Copy a [dof][comp] device Rank2View into a flat node-major device View.
// Uses the 2D indexing src(i, c) so the result is correct regardless of the
// source view's layout (unlike flattening through data_handle(), which assumes
// node-major storage and silently reorders wide/multi-component data on GPU).
template <typename T>
void CopyDeviceRank2ViewToDeviceView(Kokkos::View<T*, DeviceMemorySpace> dest,
                                     Rank2View<const T, DeviceMemorySpace> src)
{
  const size_t num_dof = src.extent(0);
  const size_t num_comp = src.extent(1);
  if (dest.extent(0) != num_dof * num_comp) {
    throw pcms_error("CopyDeviceRank2ViewToDeviceView: size mismatch");
  }
  Kokkos::parallel_for(
    "CopyDeviceRank2ViewToDeviceView", num_dof, KOKKOS_LAMBDA(LO i) {
      for (size_t c = 0; c < num_comp; ++c) {
        dest(i * num_comp + c) = src(i, c);
      }
    });
}

// Rank-2 destination overload, so callers can migrate their internal storage
// from a flat View<T*> to a View<T**> without changing the call site.
template <typename T>
void CopyDeviceRank2ViewToDeviceView(Kokkos::View<T**, DeviceMemorySpace> dest,
                                     Rank2View<const T, DeviceMemorySpace> src)
{
  const size_t num_dof = src.extent(0);
  const size_t num_comp = src.extent(1);
  if (dest.extent(0) != num_dof || dest.extent(1) != num_comp) {
    throw pcms_error("CopyDeviceRank2ViewToDeviceView: size mismatch");
  }
  Kokkos::parallel_for(
    "CopyDeviceRank2ViewToDeviceView", num_dof, KOKKOS_LAMBDA(LO i) {
      for (size_t c = 0; c < num_comp; ++c) {
        dest(i, c) = src(i, c);
      }
    });
}

// Inverse of CopyDeviceRank2ViewToDeviceView(Kokkos::View<T*,...>, Rank2View):
// scatter a flat node-major device View back into a [dof][comp] device
// Kokkos::View<T**,...>. T is deduced from the (non-alias) Kokkos::View
// destination parameter, since alias-template argument deduction through
// Rank2View's multi-level alias chain is not reliable. Uses 2D indexing on
// the destination so the result is correct regardless of its layout.
template <typename T>
void CopyDeviceViewToRank2View(Kokkos::View<T**, DeviceMemorySpace> dest,
                               Kokkos::View<const T*, DeviceMemorySpace> src)
{
  const size_t num_dof = dest.extent(0);
  const size_t num_comp = dest.extent(1);
  if (src.extent(0) != num_dof * num_comp) {
    throw pcms_error("CopyDeviceViewToRank2View: size mismatch");
  }
  Kokkos::parallel_for(
    "CopyDeviceViewToRank2View", num_dof, KOKKOS_LAMBDA(LO i) {
      for (size_t c = 0; c < num_comp; ++c) {
        dest(i, c) = src(i * num_comp + c);
      }
    });
}

// Host analogue: copy a [dof][comp] host Rank2View into a flat node-major
// device View. The host Rank2View is layout_right (contiguous node-major), so
// we alias its storage with an unmanaged view and deep_copy straight to device.
template <typename T>
void CopyHostRank2ViewToDeviceView(Kokkos::View<T*, DeviceMemorySpace> dest,
                                   Rank2View<const T, HostMemorySpace> src)
{
  const size_t num_dof = src.extent(0);
  const size_t num_comp = src.extent(1);
  if (dest.extent(0) != num_dof * num_comp) {
    throw pcms_error("CopyHostRank2ViewToDeviceView: size mismatch");
  }
  Kokkos::View<const T*, HostMemorySpace,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>>
    src_view(src.data_handle(), num_dof * num_comp);
  Kokkos::deep_copy(dest, src_view);
}

// Rank-2 destination overload for host-to-device copies. The device View<T**>
// may use a different layout than the host source, so stage through a mirror
// that matches the destination layout.
template <typename T>
void CopyHostRank2ViewToDeviceView(Kokkos::View<T**, DeviceMemorySpace> dest,
                                   Rank2View<const T, HostMemorySpace> src)
{
  const size_t num_dof = src.extent(0);
  const size_t num_comp = src.extent(1);
  if (dest.extent(0) != num_dof || dest.extent(1) != num_comp) {
    throw pcms_error("CopyHostRank2ViewToDeviceView: size mismatch");
  }
  auto src_mirror = Kokkos::create_mirror_view(dest);
  for (size_t i = 0; i < num_dof; ++i) {
    for (size_t c = 0; c < num_comp; ++c) {
      src_mirror(i, c) = src(i, c);
    }
  }
  Kokkos::deep_copy(dest, src_mirror);
}

// Element-wise copy between device Rank2Views. Indexing through (i, c) makes
// the copy correct whatever layouts the two views carry, unlike a flat copy
// through data_handle(), which assumes both element sequences agree.
template <typename T>
void CopyDeviceRank2ViewToRank2View(Rank2View<T, DeviceMemorySpace> dest,
                                    Rank2View<const T, DeviceMemorySpace> src)
{
  const size_t num_dof = src.extent(0);
  const size_t num_comp = src.extent(1);
  if (dest.extent(0) != num_dof || dest.extent(1) != num_comp) {
    throw pcms_error("CopyDeviceRank2ViewToRank2View: extent mismatch");
  }
  Kokkos::parallel_for(
    "CopyDeviceRank2ViewToRank2View",
    Kokkos::RangePolicy<typename DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(num_dof)),
    KOKKOS_LAMBDA(LO i) {
      for (size_t c = 0; c < num_comp; ++c) {
        dest(i, c) = src(i, c);
      }
    });
}

// utility function to copy from Rank1View to Rank1View
template <typename T>
void CopyRank1ViewToHost(Rank1View<T, HostMemorySpace> dest,
                         Rank1View<const T, DeviceMemorySpace> src)
{
  if (dest.size() != src.size()) {
    throw pcms_error("CopyRank1ViewToHost: size mismatch");
  }
  Kokkos::View<T*, DeviceMemorySpace> src_tmp("CopyRank1ViewToHost_tmp",
                                              src.size());
  Kokkos::parallel_for(
    "CopyRank1ViewToHost", src.size(),
    KOKKOS_LAMBDA(LO i) { src_tmp(i) = src(i); });
  Kokkos::View<T*, HostMemorySpace> src_tmp_h("CopyRank1ViewToHost_tmp_h",
                                              src.size());
  Kokkos::deep_copy(src_tmp_h, src_tmp);
  for (size_t i = 0; i < src.size(); ++i) {
    dest(i) = src_tmp_h(i);
  }
}

template <typename T>
void CopyRank1ViewToDevice(Rank1View<T, DeviceMemorySpace> dest,
                           Rank1View<const T, HostMemorySpace> src)
{
  if (dest.size() != src.size()) {
    throw pcms_error("CopyRank1ViewToDevice: size mismatch");
  }
  Kokkos::View<T*, HostMemorySpace> src_tmp("CopyRank1ViewToDevice_tmp",
                                            src.size());
  for (size_t i = 0; i < src.size(); ++i) {
    src_tmp(i) = src(i);
  }
  Kokkos::View<T*, DeviceMemorySpace> src_tmp_d("CopyRank1ViewToDevice_tmp_d",
                                                src.size());
  Kokkos::deep_copy(src_tmp_d, src_tmp);
  Kokkos::parallel_for(
    "CopyRank1ViewToDevice", src.size(),
    KOKKOS_LAMBDA(LO i) { dest(i) = src_tmp_d(i); });
}

// dest(rows(i), c) = src(i, c) for every row i of src, on the device.
template <typename T>
void ScatterRows(Rank2View<const T, DeviceMemorySpace> src,
                 Kokkos::View<const LO*, DeviceMemorySpace> rows,
                 Rank2View<T, DeviceMemorySpace> dest)
{
  const size_t num_rows = src.extent(0);
  const size_t num_comp = src.extent(1);
  Kokkos::parallel_for(
    "pcms::ScatterRows",
    Kokkos::MDRangePolicy<DeviceMemorySpace::execution_space, Kokkos::Rank<2>>(
      {0, 0}, {num_rows, num_comp}),
    KOKKOS_LAMBDA(size_t i, size_t c) { dest(rows(i), c) = src(i, c); });
}

// ---- std::vector <-> device copies ----
//
// Each performs at most one host/device transfer. A row-major `std::vector`
// is the layout MPI buffers use; the device side uses the default layout for
// the device memory space, which may be column-major.

// `host` copied to a new device view.
template <typename T>
Kokkos::View<T*, DeviceMemorySpace> CopyToDevice(const std::vector<T>& host)
{
  Kokkos::View<T*, DeviceMemorySpace> device("CopyToDevice", host.size());
  const Kokkos::View<const T*, HostMemorySpace, Kokkos::MemoryUnmanaged>
    host_view(host.data(), host.size());
  Kokkos::deep_copy(device, host_view);
  return device;
}

// A device matrix as a row-major host vector, `extent(0) * extent(1)` long.
template <typename T, typename... Properties>
std::vector<std::remove_const_t<T>> CopyToHostRowMajor(
  const Kokkos::View<T**, Properties...>& device)
{
  using U = std::remove_const_t<T>;
  std::vector<U> host(device.extent(0) * device.extent(1));
  Kokkos::View<U**, Kokkos::LayoutRight, HostMemorySpace,
               Kokkos::MemoryUnmanaged>
    host_view(host.data(), device.extent(0), device.extent(1));
  DeepCopyMismatchLayouts(host_view, device);
  return host;
}

// A row-major host vector of `n_cols` columns as a new device matrix.
template <typename T>
Kokkos::View<T**, DeviceMemorySpace> CopyRowMajorToDevice(
  const std::vector<T>& host, int n_cols)
{
  const size_t n_rows =
    n_cols > 0 ? host.size() / static_cast<size_t>(n_cols) : 0;
  const Kokkos::View<const T**, Kokkos::LayoutRight, HostMemorySpace,
                     Kokkos::MemoryUnmanaged>
    host_view(host.data(), n_rows, n_cols);
  Kokkos::View<T**, DeviceMemorySpace> device("CopyRowMajorToDevice", n_rows,
                                              n_cols);
  DeepCopyMismatchLayouts(device, host_view);
  return device;
}

// The rows `rows(k)` of a rank-2 device view, gathered into a new device
// matrix on the device (no transfer).
template <typename SrcView>
Kokkos::View<typename SrcView::non_const_value_type**, DeviceMemorySpace>
GatherRows(const SrcView& src, Kokkos::View<const LO*, DeviceMemorySpace> rows)
{
  static_assert(Kokkos::is_view<SrcView>::value && SrcView::rank() == 2,
                "GatherRows takes a rank-2 Kokkos::View");
  using T = typename SrcView::non_const_value_type;
  const int n_cols = static_cast<int>(src.extent(1));
  Kokkos::View<T**, DeviceMemorySpace> out("GatherRows", rows.extent(0),
                                           n_cols);
  Kokkos::parallel_for(
    "GatherRows",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, rows.extent(0)),
    KOKKOS_LAMBDA(LO k) {
      for (int c = 0; c < n_cols; ++c) {
        out(k, c) = src(rows(k), c);
      }
    });
  return out;
}

} // namespace pcms::detail

#endif // PCMS_UTILITY_DETAIL_ARRAYS_HPP
