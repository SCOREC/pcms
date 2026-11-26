#ifndef PCMS_FIELD_FIELD_DATA_HPP
#define PCMS_FIELD_FIELD_DATA_HPP

#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include <memory>
#include <span>
#include <variant>

namespace pcms
{
template <typename T>
class FieldData;
} // namespace pcms

#include "pcms/field/detail/field_data.hpp"

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

  // Bring any copies this rank keeps of DOF holders other ranks own (e.g. the
  // ghost copies of an Omega_h mesh) up to date with their owners' values.
  // Collective; a no-op for backends that keep no such copies.
  virtual void Synchronize() {}

  virtual ~FieldData() noexcept = default;

protected:
  // Coefficient operations: the extension API through which PCMS computes
  // with the coefficients of fields (the linear solvers and projections),
  // reached through the pcms::detail functions of the same names. The set
  // follows the SUNDIALS N_Vector naming and grows by adding virtuals whose
  // default throws. They act in place on the coefficients this rank owns, on
  // the kind's own storage, and leave the ghost copies stale: PCMS
  // synchronizes the field (Synchronize) before anything reads them. The
  // operands are field data of this kind, layout and value basis (PCMS checks
  // them); all but the rank reductions throw unless the kind overrides them.

  /// @brief this = a x + b y
  virtual void LinearSum(T /*a*/, const FieldData& /*x*/, T /*b*/,
                         const FieldData& /*y*/)
  {
    NoCoefficientOperations();
  }

  /// @brief this = c x
  virtual void Scale(T /*c*/, const FieldData& /*x*/)
  {
    NoCoefficientOperations();
  }

  /// @brief this = c
  virtual void Const(T /*c*/) { NoCoefficientOperations(); }

  /// @brief this = x / y, entry by entry
  virtual void Div(const FieldData& /*x*/, const FieldData& /*y*/)
  {
    NoCoefficientOperations();
  }

  /// @return the sum of this_i x_i over the coefficients this rank owns
  [[nodiscard]] virtual T DotProdLocal(const FieldData& /*x*/) const
  {
    NoCoefficientOperations();
  }

  /// @return the smallest coefficient this rank owns
  [[nodiscard]] virtual T MinLocal() const { NoCoefficientOperations(); }

  /// @return the largest coefficient magnitude this rank owns
  [[nodiscard]] virtual T MaxNormLocal() const { NoCoefficientOperations(); }

  /// @brief Replaces each value by its sum over the ranks the data is
  ///        partitioned over. Collective; a single rank by default.
  virtual void SumOverRanks(std::span<T> /*values*/) const {}

  /// @brief Replaces each value by its maximum over the ranks the data is
  ///        partitioned over. Collective; a single rank by default.
  virtual void MaxOverRanks(std::span<T> /*values*/) const {}

private:
  template <typename U>
  friend void detail::LinearSum(U a, const FieldData<U>& x, U b,
                                const FieldData<U>& y, FieldData<U>& z);
  template <typename U>
  friend void detail::Scale(U c, const FieldData<U>& x, FieldData<U>& z);
  template <typename U>
  friend void detail::Const(U c, FieldData<U>& z);
  template <typename U>
  friend void detail::Div(const FieldData<U>& x, const FieldData<U>& y,
                          FieldData<U>& z);
  template <typename U>
  friend U detail::DotProdLocal(const FieldData<U>& x, const FieldData<U>& y);
  template <typename U>
  friend U detail::MinLocal(const FieldData<U>& x);
  template <typename U>
  friend U detail::MaxNormLocal(const FieldData<U>& x);
  template <typename U>
  friend void detail::SumOverRanks(const FieldData<U>& data,
                                   std::span<U> values);
  template <typename U>
  friend void detail::MaxOverRanks(const FieldData<U>& data,
                                   std::span<U> values);

  [[noreturn]] static void NoCoefficientOperations()
  {
    throw pcms_error("FieldData: this field data provides no coefficient "
                     "operations");
  }
};

using FieldDataVariant = std::variant<
  std::unique_ptr<FieldData<int8_t>>, std::unique_ptr<FieldData<int32_t>>,
  std::unique_ptr<FieldData<int64_t>>, std::unique_ptr<FieldData<float>>,
  std::unique_ptr<FieldData<double>>>;

} // namespace pcms

#endif // PCMS_FIELD_FIELD_DATA_HPP
