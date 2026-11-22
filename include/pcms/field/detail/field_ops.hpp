#ifndef PCMS_FIELD_DETAIL_FIELD_OPS_HPP
#define PCMS_FIELD_DETAIL_FIELD_OPS_HPP

#include "pcms/field/field.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/mpi_type.hpp"
#include <Kokkos_Core.hpp>
#include <mpi.h>
#include <string>

namespace pcms::detail
{

template <typename T>
void CheckSameLayout(const Field<T>& a, const Field<T>& b, const char* what)
{
  if (&a.GetLayout() != &b.GetLayout()) {
    throw pcms_error(std::string(what) +
                     ": the fields must have the same layout");
  }
}

template <typename T>
using OwnedValues = Rank2View<const T, DeviceMemorySpace>;

/// @brief op(x, y) elementwise over values of the same shape
template <typename T, typename Op>
Kokkos::View<T**, DeviceMemorySpace> MapOwned(OwnedValues<T> x,
                                              OwnedValues<T> y, Op op)
{
  Kokkos::View<T**, DeviceMemorySpace> out(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "pcms::MapOwned"),
    x.extent(0), x.extent(1));
  const auto nc = static_cast<int>(x.extent(1));
  Kokkos::parallel_for(
    "pcms::MapOwned",
    Kokkos::RangePolicy<DefaultExecutionSpace>(
      0, static_cast<int>(x.extent(0)) * nc),
    KOKKOS_LAMBDA(const int k) {
      out(k / nc, k % nc) = op(x(k / nc, k % nc), y(k / nc, k % nc));
    });
  return out;
}

/// @brief `local` reduced with `op` over the layout's communicator
template <typename T>
[[nodiscard]] T Allreduce(const FieldLayout& layout, T local, MPI_Op op)
{
  T global{};
  MPI_Allreduce(&local, &global, 1, GetMPIType(local), op,
                layout.GetCommunicator());
  return global;
}

/// @brief The sum over the layout's ranks of op(x, y) over values of the same
///        shape
template <typename T, typename Op>
[[nodiscard]] T SumOwned(const FieldLayout& layout, OwnedValues<T> x,
                         OwnedValues<T> y, Op op)
{
  const auto nc = static_cast<int>(x.extent(1));
  T sum{};
  Kokkos::parallel_reduce(
    "pcms::SumOwned",
    Kokkos::RangePolicy<DefaultExecutionSpace>(
      0, static_cast<int>(x.extent(0)) * nc),
    KOKKOS_LAMBDA(const int k, T& partial) {
      partial += op(x(k / nc, k % nc), y(k / nc, k % nc));
    },
    sum);
  return Allreduce(layout, sum, MPI_SUM);
}

/// @brief field += addend over the owned values, then synchronizes the field.
/// @param addend shaped [owned DOF holders][components]
template <typename T>
void AddToField(Field<T>& field,
                const Kokkos::View<T**, DeviceMemorySpace>& addend)
{
  field.SetDOFHolderDataUnchecked(MakeConstRank2View(MapOwned<T>(
    field.GetDOFHolderData().GetValues(), MakeConstRank2View(addend),
    KOKKOS_LAMBDA(const T a, const T b) { return a + b; })));
  field.Synchronize();
}

/// @brief y = x over the owned values, whatever their value bases
template <typename T>
void CopyValues(const Field<T>& x, Field<T>& y)
{
  CheckSameLayout(x, y, "Copy");
  y.SetDOFHolderDataUnchecked(x.GetDOFHolderData().GetValues());
  y.Synchronize();
}

/// @brief Sets every owned value of the field to zero.
template <typename T>
void SetZero(Field<T>& field)
{
  Kokkos::View<T**, DeviceMemorySpace> zeros(
    "pcms::SetZero", field.GetLayout().GetNumDofHolder(),
    field.GetLayout().GetNumComponents());
  field.SetDOFHolderDataUnchecked(MakeConstRank2View(zeros));
  field.Synchronize();
}

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_FIELD_OPS_HPP
