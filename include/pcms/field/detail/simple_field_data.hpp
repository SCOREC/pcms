#ifndef PCMS_FIELD_DETAIL_SIMPLE_FIELD_DATA_HPP
#define PCMS_FIELD_DETAIL_SIMPLE_FIELD_DATA_HPP

#include "pcms/field/detail/field_data_ops.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <string>

namespace pcms::detail
{

// SimpleFieldData<T> is a generic concrete FieldData<T> backed by shaped
// [dof_holder][component] Kokkos Views.
//
// Ownership of the layout is shared — the layout is typically held by the
// factory that created this field data object.
template <typename T>
class SimpleFieldData
  : public FieldData<T>
  , public DOFAccess<T>
{
public:
  SimpleFieldData(std::shared_ptr<const FieldLayout> layout, ValueBasis basis)
    : layout_(std::move(layout)),
      basis_(std::move(basis)),
      device_data_("simple_field_data_device",
                   static_cast<size_t>(layout_->GetNumDofHolder()),
                   static_cast<size_t>(layout_->GetNumComponents()))
  {
  }

  FieldValueType GetValueType() const override
  {
    return detail::ValueTypeOfRank(basis_.Rank());
  }
  const ValueBasis& GetValueBasis() const override { return basis_; }

  const DOFAccess<T>* GetDOFAccess() const noexcept override { return this; }
  DOFAccess<T>* GetDOFAccess() noexcept override { return this; }

  Rank2View<const T, DeviceMemorySpace> GetValues() const override
  {
    return MakeConstRank2View(device_data_);
  }

  void SetValues(Rank2View<const T, DeviceMemorySpace> values) override
  {
    detail::CopyDeviceRank2ViewToDeviceView(device_data_, values);
  }

protected:
  void LinearSum(T a, const FieldData<T>& x, T b,
                 const FieldData<T>& y) override
  {
    LinearSumValues(device_data_, a, Data(x), b, Data(y));
  }

  void Scale(T c, const FieldData<T>& x) override
  {
    ScaleValues(device_data_, c, Data(x));
  }

  void Const(T c) override { ConstValues(device_data_, c); }

  void Div(const FieldData<T>& x, const FieldData<T>& y) override
  {
    DivValues(device_data_, Data(x), Data(y));
  }

  [[nodiscard]] T DotProdLocal(const FieldData<T>& x) const override
  {
    return DotProdValues(device_data_, Data(x));
  }

  [[nodiscard]] T MinLocal() const override { return MinValues(device_data_); }

  [[nodiscard]] T MaxNormLocal() const override
  {
    return MaxNormValues(device_data_);
  }

private:
  static const Kokkos::View<T**, DeviceMemorySpace>& Data(const FieldData<T>& x)
  {
    return OperandOfKind<SimpleFieldData>(x).device_data_;
  }

  std::shared_ptr<const FieldLayout> layout_;
  ValueBasis basis_;
  Kokkos::View<T**, DeviceMemorySpace> device_data_;
};

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_SIMPLE_FIELD_DATA_HPP
