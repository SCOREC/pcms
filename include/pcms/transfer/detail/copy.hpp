#ifndef PCMS_TRANSFER_DETAIL_COPY_HPP
#define PCMS_TRANSFER_DETAIL_COPY_HPP

#include "pcms/field/field.hpp"
#include "pcms/field/detail/field_ops.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/value_transform.hpp"
#include "pcms/transfer/transfer_operator.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/profile.hpp"
#include "pcms/utility/types.hpp"
#include <memory>
#include <optional>
#include <type_traits>

namespace pcms::detail
{

template <typename T>
void CheckCopyCompatible(const Field<T>& source, const Field<T>& target)
{
  if (&source.GetLayout() != &target.GetLayout()) {
    throw pcms_error("Copy: source and target layouts differ");
  }
  const auto& sd = source.GetData();
  const auto& td = target.GetData();
  if (sd.GetValueType() != td.GetValueType() ||
      !detail::SameValueBasis(sd.GetValueBasis(), td.GetValueBasis())) {
    throw pcms_error("Copy: source and target value type/basis differ");
  }
}

// Copies DOF values between fields of spaces with the same layout. Create it
// with CreateCopy (pcms/transfer/copy.hpp).
template <typename T>
class Copy : public TransferOperator<T>
{
public:
  /// @throws pcms_error if a space is null or the spaces have different
  ///         layouts
  Copy(std::shared_ptr<const FunctionSpace> source_space,
       std::shared_ptr<const FunctionSpace> target_space)
    : source_space_(std::move(source_space)),
      target_space_(std::move(target_space))
  {
    if (source_space_ == nullptr || target_space_ == nullptr) {
      throw pcms_error("Copy: the function spaces must not be null");
    }
    if (source_space_->GetLayout().get() != target_space_->GetLayout().get()) {
      throw pcms_error("Copy: source and target function spaces have "
                       "different layouts");
    }
  }

private:
  void ApplyImpl(const Field<T>& source, Field<T>& target) const override
  {
    PCMS_FUNCTION_TIMER;
    if constexpr (std::is_same_v<T, Real>) {
      const auto& stored = source.GetData().GetValueBasis();
      const auto& declared = target.GetData().GetValueBasis();
      const auto* scaling =
        SameValueBasisExceptScaling(stored, declared)
          ? CachedScalingTransform(
              scaling_, stored, declared,
              target_space_->GetLayout()->GetDOFHolderCoordinates())
          : nullptr;
      if (scaling != nullptr) {
        detail::CheckSameLayout(source, target, "Copy");
        const auto values = source.GetDOFHolderData();
        Kokkos::View<Real**, DeviceMemorySpace> converted(
          Kokkos::view_alloc(Kokkos::WithoutInitializing, "copy_converted"),
          values.extent(0), values.extent(1));
        scaling->Apply(values, ValueView<Real, DeviceMemorySpace>(
                                 declared, MakeRank2View(converted)));
        target.SetDOFHolderData(ValueView<const Real, DeviceMemorySpace>(
          declared, MakeConstRank2View(converted)));
        target.Synchronize();
        return;
      }
    }
    detail::CheckCopyCompatible(source, target);
    target.SetDOFHolderData(source.GetDOFHolderData());
    target.Synchronize();
  }

  std::shared_ptr<const FunctionSpace> source_space_;
  std::shared_ptr<const FunctionSpace> target_space_;
  mutable std::optional<ValueTransform> scaling_;
};

} // namespace pcms::detail

#endif // PCMS_TRANSFER_DETAIL_COPY_HPP
