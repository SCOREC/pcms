#ifndef PCMS_COUPLING_FIELD_H
#define PCMS_COUPLING_FIELD_H

#include "field_data.h"
#include "field_layout.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"
#include <memory>
#include <string>
#include <type_traits>
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

  // The DOF-holder accessors below read and write the values of the DOF
  // holders this rank owns through the data's DOFAccess, and throw if the data
  // has none (e.g. a field of an evaluation-only source). Writes synchronize
  // the backend's other copies and are collective; see DOFAccess.

  // Copies to the host unless device memory is host memory. The returned view
  // is valid until the next call or until the field is written or destroyed.
  [[nodiscard]] ValueView<const T, HostMemorySpace> GetDOFHolderDataHost() const
  {
    return ValueView<const T, HostMemorySpace>(data_->GetValueBasis(),
                                               OwnedValuesHost());
  }

  [[nodiscard]] ValueView<const T, DeviceMemorySpace> GetDOFHolderData() const
  {
    return ValueView<const T, DeviceMemorySpace>(
      data_->GetValueBasis(), RequireDOFAccess().GetOwnedValues());
  }

  void SetDOFHolderDataHost(ValueView<const T, HostMemorySpace> v)
  {
    CheckWriteTag(v.GetBasis());
    CheckWriteExtents(v.GetValues().extent(0), v.GetValues().extent(1));
    SetOwnedValuesFromHost(v.GetValues());
  }

  void SetDOFHolderData(ValueView<const T, DeviceMemorySpace> v)
  {
    CheckWriteTag(v.GetBasis());
    CheckWriteExtents(v.GetValues().extent(0), v.GetValues().extent(1));
    RequireDOFAccess().SetOwnedValues(v.GetValues());
  }

  // unchecked is needed for optimization in some places, but should generally
  // not be used unless you know what you are doing
  void SetDOFHolderDataUncheckedHost(Rank2View<const T, HostMemorySpace> v)
  {
    SetOwnedValuesFromHost(v);
  }

  // unchecked is needed for optimization in some places, but should generally
  // not be used unless you know what you are doing
  void SetDOFHolderDataUnchecked(Rank2View<const T, DeviceMemorySpace> v)
  {
    RequireDOFAccess().SetOwnedValues(v);
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

  const DOFAccess<T>& RequireDOFAccess() const
  {
    const auto* access = data_->GetDOFAccess();
    if (access == nullptr) {
      throw pcms_error("Field '" + name_ +
                       "': its values are not accessible through PCMS (the "
                       "field data provides no DOFAccess)");
    }
    return *access;
  }

  DOFAccess<T>& RequireDOFAccess()
  {
    return const_cast<DOFAccess<T>&>(std::as_const(*this).RequireDOFAccess());
  }

  Rank2View<const T, HostMemorySpace> OwnedValuesHost() const
  {
    const auto values = RequireDOFAccess().GetOwnedValues();
    if constexpr (std::is_same_v<HostMemorySpace, DeviceMemorySpace>) {
      return values;
    } else {
      // Stage in a device view of the default layout so the copy to the
      // host handles whatever layout the backend's view carries.
      Kokkos::View<T**, DeviceMemorySpace> staged(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "owned_values_staged"),
        values.extent(0), values.extent(1));
      CopyDeviceRank2ViewToRank2View(MakeRank2View(staged), values);
      host_values_ = Kokkos::View<T**, HostMemorySpace>(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "owned_values_host"),
        values.extent(0), values.extent(1));
      DeepCopyMismatchLayouts(host_values_, staged);
      return MakeConstRank2View(host_values_);
    }
  }

  void SetOwnedValuesFromHost(Rank2View<const T, HostMemorySpace> values)
  {
    if constexpr (std::is_same_v<HostMemorySpace, DeviceMemorySpace>) {
      RequireDOFAccess().SetOwnedValues(values);
    } else {
      Kokkos::View<T**, DeviceMemorySpace> staged(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "owned_values_staged"),
        values.extent(0), values.extent(1));
      CopyHostRank2ViewToDeviceView(staged, values);
      RequireDOFAccess().SetOwnedValues(MakeConstRank2View(staged));
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
  // Host copy backing GetDOFHolderDataHost when device memory is not host
  // memory.
  mutable Kokkos::View<T**, HostMemorySpace> host_values_;
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
