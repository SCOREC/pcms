#ifndef PCMS_FIELD_FIELD_OPS_HPP
#define PCMS_FIELD_FIELD_OPS_HPP

#include "pcms/field/detail/field_ops.hpp"
#include "pcms/field/field.hpp"
#include "pcms/utility/assert.hpp"
#include <Kokkos_Core.hpp>
#include <cmath>

namespace pcms
{

// Operations on the DOF vectors of fields: the owned entries of a field's
// data, in the layout's order. They are basis-dependent by design (coupling
// relaxation, error control by weighted norms); the basis-independent
// pairing of a functional with a function is Apply.
//
// The operands of an operation share a layout. The operations that write leave
// the field consistent (ghost copies are synchronized). They are collective
// over the data's ranks in the kinds that override the FieldData operations;
// the generic defaults are single-rank (see FieldData::Scale).

/// y *= a
template <typename T>
void Scale(T a, Field<T>& y)
{
  y.GetData().Scale(a);
  y.Synchronize();
}

/// y += a x
template <typename T>
void Axpy(T a, const Field<T>& x, Field<T>& y)
{
  detail::CheckSameLayout(x, y, "Axpy");
  y.GetData().Axpy(a, x.GetData());
  y.Synchronize();
}

/// @return the sum over the owned entries of x_i y_i
template <typename T>
[[nodiscard]] T Dot(const Field<T>& x, const Field<T>& y)
{
  detail::CheckSameLayout(x, y, "Dot");
  return x.GetData().Dot(y.GetData());
}

/// @return the Euclidean norm of the DOF vector
template <typename T>
[[nodiscard]] T Norm2(const Field<T>& x)
{
  return std::sqrt(x.GetData().Dot(x.GetData()));
}

/// @return the largest absolute entry of the DOF vector
template <typename T>
[[nodiscard]] T NormInf(const Field<T>& x)
{
  return x.GetData().MaxAbs();
}

/// @return sqrt(sum_i (w_i x_i)^2 / N) over the N owned entries, the weighted
///         RMS norm of error control (the weights are a field of the same
///         layout)
template <typename T>
[[nodiscard]] T WeightedRmsNorm(const Field<T>& x, const Field<T>& weights)
{
  detail::CheckSameLayout(x, weights, "WeightedRmsNorm");
  const auto* xa = x.GetData().GetDOFAccess();
  const auto* wa = weights.GetData().GetDOFAccess();
  if (xa == nullptr || wa == nullptr) {
    throw pcms_error("WeightedRmsNorm: the data provides no DOFAccess");
  }
  const auto xv = xa->GetValues();
  const auto wv = wa->GetValues();
  const auto n = static_cast<int>(xv.extent(0));
  const auto nc = static_cast<int>(xv.extent(1));
  T sum = T{};
  Kokkos::parallel_reduce(
    "pcms::WeightedRmsNorm",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n * nc),
    KOKKOS_LAMBDA(const int k, T& partial) {
      const T v = wv(k / nc, k % nc) * xv(k / nc, k % nc);
      partial += v * v;
    },
    sum);
  return n * nc == 0 ? T{} : std::sqrt(sum / static_cast<T>(n * nc));
}

/// @brief The pairing of a functional (e.g. the moments E^T W f of a function
///        f) with a function: l(u) = sum_i b_i u_i. Basis-independent: the
///        functional is dual to the coefficients.
/// @throws pcms_error if the layouts differ
template <typename T>
[[nodiscard]] T Apply(const Field<T>& functional, const Function<T>& u)
{
  detail::CheckSameLayout(functional, u, "Apply");
  return functional.GetData().Dot(u.GetData());
}

} // namespace pcms

#endif // PCMS_FIELD_FIELD_OPS_HPP
