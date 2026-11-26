#ifndef PCMS_FIELD_FIELD_DATA_HPP
#define PCMS_FIELD_FIELD_DATA_HPP

#include "pcms/field/detail/field_data_ops.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include <memory>
#include <utility>
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
  // values until FieldData::Synchronize. Collective over the ranks that
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

  // Raw vector operations on the owned entries (the DOFAccess values), which
  // the public free functions of pcms/field/field_ops.hpp call. The defaults
  // work through DOFAccess with Kokkos kernels on the device and throw for
  // data without one; a kind overrides them to work in place and to reduce
  // over its own communicator. The defaults are SINGLE-RANK: layouts carry no
  // communicator yet, so Dot and MaxAbs sum/maximize over this rank's owned
  // entries only. The operands share a layout (the free functions check it).
  // They write the owned values and leave ghost copies stale (Synchronize).
  virtual void Scale(T a) // this *= a
  {
    Access().SetValues(
      MakeConstRank2View(detail::ScaleValues<T>(a, Access().GetValues())));
  }
  virtual void Axpy(T a, const FieldData<T>& x) // this += a x
  {
    Access().SetValues(MakeConstRank2View(
      detail::AxpyValues<T>(a, x.Access().GetValues(), Access().GetValues())));
  }
  [[nodiscard]] virtual T Dot(const FieldData<T>& x) const // sum this_i x_i
  {
    return detail::DotValues<T>(Access().GetValues(), x.Access().GetValues());
  }
  [[nodiscard]] virtual T MaxAbs() const // max |this_i|
  {
    return detail::MaxAbsValues<T>(Access().GetValues());
  }
  // Bring any copies this rank keeps of DOF holders other ranks own (e.g. the
  // ghost copies of an Omega_h mesh) up to date with their owners' values.
  // Collective; a no-op for backends that keep no such copies.
  virtual void Synchronize() {}

  virtual ~FieldData() noexcept = default;

private:
  // the DOFAccess the generic operations work through
  const DOFAccess<T>& Access() const
  {
    const auto* access = GetDOFAccess();
    if (access == nullptr) {
      throw pcms_error("FieldData operation: the data provides no DOFAccess; "
                       "its kind must override Scale/Axpy/Dot/MaxAbs");
    }
    return *access;
  }
  DOFAccess<T>& Access()
  {
    return const_cast<DOFAccess<T>&>(std::as_const(*this).Access());
  }
};

using FieldDataVariant = std::variant<
  std::unique_ptr<FieldData<int8_t>>, std::unique_ptr<FieldData<int32_t>>,
  std::unique_ptr<FieldData<int64_t>>, std::unique_ptr<FieldData<float>>,
  std::unique_ptr<FieldData<double>>>;

} // namespace pcms

#endif // PCMS_FIELD_FIELD_DATA_HPP
