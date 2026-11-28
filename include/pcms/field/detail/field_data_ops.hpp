#ifndef PCMS_FIELD_DETAIL_FIELD_DATA_OPS_HPP
#define PCMS_FIELD_DETAIL_FIELD_DATA_OPS_HPP

#include "pcms/field/field_data.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include <Kokkos_Core.hpp>
#include <cstdint>

namespace pcms::detail
{

// The kernels of the FieldData coefficient operations, in place on a kind's
// [dof_holder][component] storage.

template <typename T>
using CoefficientView = Kokkos::View<T**, DeviceMemorySpace>;

/// @return `x` as the kind `Kind`
/// @throws pcms_error if `x` is field data of another kind
template <typename Kind, typename T>
const Kind& OperandOfKind(const FieldData<T>& x)
{
  const auto* kind = dynamic_cast<const Kind*>(&x);
  if (kind == nullptr) {
    throw pcms_error("FieldData: the operand is field data of another kind");
  }
  return *kind;
}

template <typename T>
void CheckSameShape(const CoefficientView<T>& a, const CoefficientView<T>& b)
{
  if (a.extent(0) != b.extent(0) || a.extent(1) != b.extent(1)) {
    throw pcms_error("FieldData: the operands have different numbers of "
                     "coefficients");
  }
}

template <typename T>
auto CoefficientPolicy(const CoefficientView<T>& out)
{
  return Kokkos::MDRangePolicy<DeviceMemorySpace::execution_space,
                               Kokkos::Rank<2>>(
    {0, 0},
    {static_cast<int64_t>(out.extent(0)), static_cast<int64_t>(out.extent(1))});
}

/// @brief out = a x + b y
template <typename T>
void LinearSumValues(const CoefficientView<T>& out, T a,
                     const CoefficientView<T>& x, T b,
                     const CoefficientView<T>& y)
{
  CheckSameShape(out, x);
  CheckSameShape(out, y);
  Kokkos::parallel_for(
    "pcms::FieldData::LinearSum", CoefficientPolicy(out),
    KOKKOS_LAMBDA(const int64_t i, const int64_t c) {
      out(i, c) = a * x(i, c) + b * y(i, c);
    });
}

/// @brief out = c x
template <typename T>
void ScaleValues(const CoefficientView<T>& out, T c,
                 const CoefficientView<T>& x)
{
  CheckSameShape(out, x);
  Kokkos::parallel_for(
    "pcms::FieldData::Scale", CoefficientPolicy(out),
    KOKKOS_LAMBDA(const int64_t i, const int64_t k) {
      out(i, k) = c * x(i, k);
    });
}

/// @brief out = c
template <typename T>
void ConstValues(const CoefficientView<T>& out, T c)
{
  Kokkos::deep_copy(out, c);
}

/// @brief out = x / y, entry by entry
template <typename T>
void DivValues(const CoefficientView<T>& out, const CoefficientView<T>& x,
               const CoefficientView<T>& y)
{
  CheckSameShape(out, x);
  CheckSameShape(out, y);
  Kokkos::parallel_for(
    "pcms::FieldData::Div", CoefficientPolicy(out),
    KOKKOS_LAMBDA(const int64_t i, const int64_t c) {
      out(i, c) = x(i, c) / y(i, c);
    });
}

/// @return the sum of x_i y_i
template <typename T>
T DotProdValues(const CoefficientView<T>& x, const CoefficientView<T>& y)
{
  CheckSameShape(x, y);
  T sum = T{};
  Kokkos::parallel_reduce(
    "pcms::FieldData::DotProdLocal", CoefficientPolicy(x),
    KOKKOS_LAMBDA(const int64_t i, const int64_t c, T& partial) {
      partial += x(i, c) * y(i, c);
    },
    sum);
  return sum;
}

/// @return the smallest x_i; the largest T for no coefficients
template <typename T>
T MinValues(const CoefficientView<T>& x)
{
  T min{};
  Kokkos::parallel_reduce(
    "pcms::FieldData::MinLocal", CoefficientPolicy(x),
    KOKKOS_LAMBDA(const int64_t i, const int64_t c, T& partial) {
      partial = x(i, c) < partial ? x(i, c) : partial;
    },
    Kokkos::Min<T>(min));
  return min;
}

/// @return the largest |x_i|; 0 for no coefficients
template <typename T>
T MaxNormValues(const CoefficientView<T>& x)
{
  T max{};
  Kokkos::parallel_reduce(
    "pcms::FieldData::MaxNormLocal", CoefficientPolicy(x),
    KOKKOS_LAMBDA(const int64_t i, const int64_t c, T& partial) {
      const T a = Kokkos::abs(x(i, c));
      partial = a > partial ? a : partial;
    },
    Kokkos::Max<T>(max));
  return max;
}

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_FIELD_DATA_OPS_HPP
