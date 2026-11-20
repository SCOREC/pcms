#ifndef PCMS_FIELD_DETAIL_FIELD_DATA_OPS_HPP
#define PCMS_FIELD_DETAIL_FIELD_DATA_OPS_HPP

#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include <Kokkos_Core.hpp>

namespace pcms::detail
{

// The kernels of the generic FieldData operations, on the owned values.

template <typename T>
using OwnedValues = Rank2View<const T, DeviceMemorySpace>;

template <typename T>
void CheckSameOwnedShape(const OwnedValues<T>& a, const OwnedValues<T>& b)
{
  if (a.extent(0) != b.extent(0) || a.extent(1) != b.extent(1)) {
    throw pcms_error("FieldData operation: the operands have different "
                     "numbers of owned values");
  }
}

// a * x + y, elementwise (y = x with a = 0 is not what this is for: see
// ScaledOwned)
template <typename T>
Kokkos::View<T**, DeviceMemorySpace> AxpyValues(T a, OwnedValues<T> x,
                                                OwnedValues<T> y)
{
  CheckSameOwnedShape<T>(x, y);
  Kokkos::View<T**, DeviceMemorySpace> out(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "pcms::Axpy"), x.extent(0),
    x.extent(1));
  const auto n = static_cast<int>(out.extent(0));
  const auto nc = static_cast<int>(out.extent(1));
  Kokkos::parallel_for(
    "pcms::FieldData::Axpy",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n * nc),
    KOKKOS_LAMBDA(const int k) {
      out(k / nc, k % nc) = a * x(k / nc, k % nc) + y(k / nc, k % nc);
    });
  return out;
}

template <typename T>
Kokkos::View<T**, DeviceMemorySpace> ScaleValues(T a, OwnedValues<T> x)
{
  Kokkos::View<T**, DeviceMemorySpace> out(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "pcms::Scale"), x.extent(0),
    x.extent(1));
  const auto n = static_cast<int>(out.extent(0));
  const auto nc = static_cast<int>(out.extent(1));
  Kokkos::parallel_for(
    "pcms::FieldData::Scale",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n * nc),
    KOKKOS_LAMBDA(const int k) {
      out(k / nc, k % nc) = a * x(k / nc, k % nc);
    });
  return out;
}

template <typename T>
T DotValues(OwnedValues<T> x, OwnedValues<T> y)
{
  CheckSameOwnedShape<T>(x, y);
  const auto n = static_cast<int>(x.extent(0));
  const auto nc = static_cast<int>(x.extent(1));
  T sum = T{};
  Kokkos::parallel_reduce(
    "pcms::FieldData::Dot",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n * nc),
    KOKKOS_LAMBDA(const int k, T& partial) {
      partial += x(k / nc, k % nc) * y(k / nc, k % nc);
    },
    sum);
  return sum;
}

template <typename T>
T MaxAbsValues(OwnedValues<T> x)
{
  const auto n = static_cast<int>(x.extent(0));
  const auto nc = static_cast<int>(x.extent(1));
  T result = T{};
  Kokkos::parallel_reduce(
    "pcms::FieldData::MaxAbs",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n * nc),
    KOKKOS_LAMBDA(const int k, T& partial) {
      const T v = Kokkos::abs(x(k / nc, k % nc));
      if (v > partial) {
        partial = v;
      }
    },
    Kokkos::Max<T>(result));
  return result;
}

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_FIELD_DATA_OPS_HPP
