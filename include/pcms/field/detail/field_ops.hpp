#ifndef PCMS_FIELD_DETAIL_FIELD_OPS_HPP
#define PCMS_FIELD_DETAIL_FIELD_OPS_HPP

#include "pcms/field/field.hpp"
#include "pcms/utility/assert.hpp"
#include <Kokkos_Core.hpp>
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

/// @brief field += addend over the owned values, then synchronizes the field.
/// @param field the field to add to
/// @param addend shaped [owned DOF holders][components]; overwritten with
///        the sum
template <typename T>
void AddToField(Field<T>& field,
                const Kokkos::View<T**, DeviceMemorySpace>& addend)
{
  const auto current = field.GetDOFHolderData().GetValues();
  Kokkos::parallel_for(
    "pcms::AddToField",
    Kokkos::MDRangePolicy<DefaultExecutionSpace, Kokkos::Rank<2>>(
      {0, 0}, {static_cast<int64_t>(addend.extent(0)),
               static_cast<int64_t>(addend.extent(1))}),
    KOKKOS_LAMBDA(const int64_t i, const int64_t c) {
      addend(i, c) += current(i, c);
    });
  field.SetDOFHolderDataUnchecked(MakeConstRank2View(addend));
  field.Synchronize();
}

/// Sets every owned value of the field to zero.
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
