#ifndef PCMS_FIELD_DATA_H
#define PCMS_FIELD_DATA_H

#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include <memory>
#include <variant>

namespace pcms
{

// DOFAccess<T> is the optional capability of field data whose values PCMS can
// read and write directly: the values of the DOF holders this rank owns, shaped
// [owned DOF holder][component] in the layout's order. Coupler exchanges and
// transfer targets use it; evaluation never needs it.
template <typename T>
class DOFAccess
{
public:
  // Valid until the field data is next written or destroyed.
  virtual Rank2View<const T, DeviceMemorySpace> GetValues() const = 0;

  // Writes the owned values, then brings any internal copies the backend
  // keeps (e.g. replicas on other ranks) up to date. Ghost copies keep their
  // values until FieldData::SynchronizeGhosts. Collective over the ranks that
  // share those copies.
  virtual void SetValues(Rank2View<const T, DeviceMemorySpace> values) = 0;

  virtual ~DOFAccess() noexcept = default;
};

// FieldData<T> holds whatever a backend needs to represent a field. It is
// opaque. It should only be directly used by the type that created it.
template <typename T>
class FieldData
{
public:
  using value_type = T;

  virtual FieldValueType GetValueType() const = 0;
  virtual const ValueBasis& GetValueBasis() const = 0;

  virtual const DOFAccess<T>* GetDOFAccess() const noexcept { return nullptr; }
  virtual DOFAccess<T>* GetDOFAccess() noexcept { return nullptr; }

  // Synchronize ghost DOF-holder values from their owning ranks. A no-op for
  // backends without ghosts.
  virtual void SynchronizeGhosts() {}

  virtual ~FieldData() noexcept = default;
};

using FieldDataVariant = std::variant<
  std::unique_ptr<FieldData<int8_t>>, std::unique_ptr<FieldData<int32_t>>,
  std::unique_ptr<FieldData<int64_t>>, std::unique_ptr<FieldData<float>>,
  std::unique_ptr<FieldData<double>>>;

} // namespace pcms

#endif // PCMS_FIELD_DATA_H
