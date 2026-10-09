#ifndef PCMS_COUPLING_FIELD_H
#define PCMS_COUPLING_FIELD_H

#include "field_data.h"
#include "field_layout.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"
#include <memory>
#include <string>
#include <utility>
#include <variant>

namespace pcms
{

class FunctionSpace;

// Field<T> is the (layout + data) bundle: it pairs coefficient data with the
// FieldLayout that gives the DOF-holder ordering its meaning. It is the
// currency of the communication path and of TransferOperator::Apply; it is NOT
// evaluatable on its own (evaluation semantics live on FunctionSpace).
//
// Fields are created via FunctionSpace::CreateField(). They are move-only
// (unique_ptr member).
template <typename T>
class Field
{
public:
  Field(Field&&) = default;
  Field& operator=(Field&&) = default;
  Field(const Field&) = delete;
  Field& operator=(const Field&) = delete;
  ~Field() = default;

  FieldData<T>& GetData() noexcept { return *data_; }
  const FieldData<T>& GetData() const noexcept { return *data_; }

  const FieldLayout& GetLayout() const { return *layout_; }

  [[nodiscard]] const std::string& GetName() const noexcept { return name_; }

  [[nodiscard]] ValueView<const T, HostMemorySpace> GetDOFHolderDataHost() const
  {
    return ValueView<const T, HostMemorySpace>(data_->GetValueBasis(),
                                               data_->GetDOFHolderDataHost());
  }

  [[nodiscard]] ValueView<const T, DeviceMemorySpace> GetDOFHolderData() const
  {
    return ValueView<const T, DeviceMemorySpace>(data_->GetValueBasis(),
                                                 data_->GetDOFHolderData());
  }

  void SetDOFHolderDataHost(ValueView<const T, HostMemorySpace> v)
  {
    CheckWriteTag(v.GetBasis());
    CheckWriteExtents(v.GetValues().extent(0), v.GetValues().extent(1));
    data_->SetDOFHolderDataHost(v.GetValues());
  }

  void SetDOFHolderData(ValueView<const T, DeviceMemorySpace> v)
  {
    CheckWriteTag(v.GetBasis());
    CheckWriteExtents(v.GetValues().extent(0), v.GetValues().extent(1));
    data_->SetDOFHolderData(v.GetValues());
  }

  // unchecked is needed for optimization in some places, but should generally
  // not be used unless you know what you are doing
  void SetDOFHolderDataUncheckedHost(Rank2View<const T, HostMemorySpace> v)
  {
    data_->SetDOFHolderDataHost(v);
  }

  // unchecked is needed for optimization in some places, but should generally
  // not be used unless you know what you are doing
  void SetDOFHolderDataUnchecked(Rank2View<const T, DeviceMemorySpace> v)
  {
    data_->SetDOFHolderData(v);
  }

protected:
  Field(std::string name, std::shared_ptr<const FieldLayout> layout,
        std::unique_ptr<FieldData<T>> data)
    : name_(std::move(name)), layout_(std::move(layout)), data_(std::move(data))
  {
  }

  void CheckWriteTag(const ValueBasis& basis) const
  {
    if (!SameValueBasis(basis, data_->GetValueBasis())) {
      throw pcms_error(
        "Field::SetDOFHolderData: the written view's declared value "
        "type/basis does not match this field's declaration (use the "
        "Unchecked overload only for untagged wire/raw data)");
    }
  }

  // The basis tag cannot distinguish shapes (all scalar bases compare equal),
  // so the checked setters also verify the view's extents against the layout
  void CheckWriteExtents(size_t num_holders, size_t num_components) const
  {
    if (num_holders != static_cast<size_t>(layout_->GetNumOwnedDofHolder()) ||
        num_components != static_cast<size_t>(layout_->GetNumComponents())) {
      throw pcms_error("Field::SetDOFHolderData: the written view is shaped [" +
                       std::to_string(num_holders) + "][" +
                       std::to_string(num_components) +
                       "] but this field's layout is [" +
                       std::to_string(layout_->GetNumOwnedDofHolder()) + "][" +
                       std::to_string(layout_->GetNumComponents()) + "]");
    }
  }

  // FieldFactory constructs Fields (via WrapField); FunctionSpace constructs
  // Fields the same way and moves a Field's name, layout and data out to build
  // a Function (via CreateFunction). Both stamp the field name.
  friend class FieldFactory;
  friend class FunctionSpace;

private:
  std::string name_;
  std::shared_ptr<const FieldLayout> layout_;
  std::unique_ptr<FieldData<T>> data_;
};

// Function<T> is a Field that additionally knows its FunctionSpace, so it can
// be evaluated and used to build transfer operators. It is-a Field so the
// communication path (FieldCommunicator holds Field&) and
// TransferOperator::Apply(const Field&) bind a Function slice-free. It adds no
// virtuals: moving a Function into a Field slot (e.g. AddField(Field&&))
// intentionally slices off the space reference, leaving a valid comm-only
// Field.
template <typename T>
class Function : public Field<T>
{
public:
  Function(Function&&) = default;
  Function& operator=(Function&&) = default;
  Function(const Function&) = delete;
  Function& operator=(const Function&) = delete;

  const FunctionSpace& GetSpace() const noexcept { return *space_; }
  std::shared_ptr<const FunctionSpace> GetSpacePtr() const noexcept
  {
    return space_;
  }

protected:
  Function(std::string name, std::shared_ptr<const FieldLayout> layout,
           std::unique_ptr<FieldData<T>> data,
           std::shared_ptr<const FunctionSpace> space)
    : Field<T>(std::move(name), std::move(layout), std::move(data)),
      space_(std::move(space))
  {
  }

  friend class FunctionSpace;

private:
  std::shared_ptr<const FunctionSpace> space_;
};

using FieldVariant = std::variant<Field<int8_t>, Field<int32_t>, Field<int64_t>,
                                  Field<float>, Field<double>>;

using FunctionVariant =
  std::variant<Function<int8_t>, Function<int32_t>, Function<int64_t>,
               Function<float>, Function<double>>;

} // namespace pcms

#endif // PCMS_COUPLING_FIELD_H
