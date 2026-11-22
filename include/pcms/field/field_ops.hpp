#ifndef PCMS_FIELD_FIELD_OPS_HPP
#define PCMS_FIELD_FIELD_OPS_HPP

#include "pcms/field/detail/field_ops.hpp"
#include "pcms/field/field.hpp"
#include "pcms/utility/assert.hpp"
#include <Kokkos_Core.hpp>
#include <cmath>
#include <mpi.h>

namespace pcms
{

// Operations on the DOF vectors of fields: the owned entries of a field's
// data, in the layout's order. They are basis-dependent by design (coupling
// relaxation, error control by weighted norms); the basis-independent
// pairing of a functional with a function is Apply.
//
// The operands of an operation share a layout. The operations that write leave
// the field consistent (ghost copies are synchronized); the reductions are
// collective over the layout's communicator.

/// @brief y *= a
template <typename T>
void Scale(T a, Field<T>& y)
{
  const auto values = y.GetDOFHolderData().GetValues();
  y.SetDOFHolderDataUnchecked(MakeConstRank2View(detail::MapOwned<T>(
    values, values, KOKKOS_LAMBDA(const T v, const T) { return a * v; })));
  y.Synchronize();
}

/// @brief y += a x
template <typename T>
void Axpy(T a, const Field<T>& x, Field<T>& y)
{
  detail::CheckSameLayout(x, y, "Axpy");
  y.SetDOFHolderDataUnchecked(MakeConstRank2View(detail::MapOwned<T>(
    x.GetDOFHolderData().GetValues(), y.GetDOFHolderData().GetValues(),
    KOKKOS_LAMBDA(const T u, const T v) { return a * u + v; })));
  y.Synchronize();
}

/// @brief y = x
/// @throws pcms_error if the layouts or the value bases differ
template <typename T>
void Copy(const Field<T>& x, Field<T>& y)
{
  detail::CheckSameLayout(x, y, "Copy");
  y.SetDOFHolderData(x.GetDOFHolderData());
  y.Synchronize();
}

/// @return the sum over the owned entries of x_i y_i
template <typename T>
[[nodiscard]] T Dot(const Field<T>& x, const Field<T>& y)
{
  detail::CheckSameLayout(x, y, "Dot");
  return detail::SumOwned<T>(
    x.GetLayout(), x.GetDOFHolderData().GetValues(),
    y.GetDOFHolderData().GetValues(),
    KOKKOS_LAMBDA(const T u, const T v) { return u * v; });
}

/// @return the Euclidean norm of the DOF vector
template <typename T>
[[nodiscard]] T Norm2(const Field<T>& x)
{
  return std::sqrt(Dot(x, x));
}

/// @return the largest absolute entry of the DOF vector
template <typename T>
[[nodiscard]] T NormInf(const Field<T>& x)
{
  const auto values = x.GetDOFHolderData().GetValues();
  const auto nc = static_cast<int>(values.extent(1));
  T local{};
  Kokkos::parallel_reduce(
    "pcms::NormInf",
    Kokkos::RangePolicy<DefaultExecutionSpace>(
      0, static_cast<int>(values.extent(0)) * nc),
    KOKKOS_LAMBDA(const int k, T& partial) {
      partial = Kokkos::max(partial, T(Kokkos::abs(values(k / nc, k % nc))));
    },
    Kokkos::Max<T>(local));
  return detail::Allreduce(x.GetLayout(), Kokkos::max(local, T{}), MPI_MAX);
}

/// @return sqrt(sum_i (w_i x_i)^2 / N) over the N entries, the weighted RMS
///         norm of error control (the weights are a field of the same layout)
template <typename T>
[[nodiscard]] T WeightedRmsNorm(const Field<T>& x, const Field<T>& weights)
{
  detail::CheckSameLayout(x, weights, "WeightedRmsNorm");
  const T sum = detail::SumOwned<T>(
    x.GetLayout(), x.GetDOFHolderData().GetValues(),
    weights.GetDOFHolderData().GetValues(),
    KOKKOS_LAMBDA(const T u, const T w) { return (w * u) * (w * u); });
  const auto n = x.GetLayout().GlobalSize();
  return n == 0 ? T{} : std::sqrt(sum / static_cast<T>(n));
}

/// @brief The pairing of a functional (e.g. the moments E^T W f of a function
///        f) with a function: l(u) = sum_i b_i u_i. Basis-independent: the
///        functional is dual to the coefficients.
/// @throws pcms_error if the layouts differ
template <typename T>
[[nodiscard]] T Apply(const Field<T>& functional, const Function<T>& u)
{
  return Dot<T>(functional, u);
}

} // namespace pcms

#endif // PCMS_FIELD_FIELD_OPS_HPP
