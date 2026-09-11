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

  // Optional identifier for this field. Empty by default; set at creation.
  [[nodiscard]] const std::string& GetName() const noexcept { return name_; }

  Rank2View<const T, HostMemorySpace> GetDOFHolderDataHost() const
  {
    return data_->GetDOFHolderDataHost();
  }

  Rank2View<const T, HostMemorySpace> GetOwnedDOFHolderDataHost() const
  {
    return data_->GetOwnedDOFHolderDataHost();
  }

  void SetDOFHolderDataHost(Rank2View<const T, HostMemorySpace> v)
  {
    data_->SetDOFHolderDataHost(v);
  }

  Rank2View<const T, DeviceMemorySpace> GetDOFHolderData() const
  {
    return data_->GetDOFHolderData();
  }

  Rank2View<const T, DeviceMemorySpace> GetOwnedDOFHolderData() const
  {
    return data_->GetOwnedDOFHolderData();
  }

  void SetDOFHolderData(Rank2View<const T, DeviceMemorySpace> v)
  {
    data_->SetDOFHolderData(v);
  }

  // Synchronize ghost DOF-holder values from their owning ranks.
  void SynchronizeGhosts() { data_->SynchronizeGhosts(); }

protected:
  Field(std::string name, std::shared_ptr<const FieldLayout> layout,
        std::unique_ptr<FieldData<T>> data)
    : name_(std::move(name)), layout_(std::move(layout)), data_(std::move(data))
  {
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
