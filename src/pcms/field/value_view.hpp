#ifndef PCMS_FIELD_VALUE_VIEW_HPP
#define PCMS_FIELD_VALUE_VIEW_HPP

#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include <array>
#include <cstdint>
#include <initializer_list>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>

namespace pcms
{

enum class FieldValueType
{
  Scalar,
  Vector,
  Tensor
};

// Index variance: Contravariant == components on the basis dx/dq^i,
// Covariant == components on the dual grad q^i.
enum class Variance : std::uint8_t
{
  Contravariant,
  Covariant
};

// How the components are scaled relative to the coordinate-induced basis.
//
// each component physical coordinates all have the same units
// the natural coordinates are the standard, unscaled coordinate bases
enum class ComponentScaling : std::uint8_t
{
  Physical,
  Natural
};

class VarianceSignature
{
public:
  static constexpr int max_rank = 4;

  constexpr VarianceSignature() = default;
  constexpr VarianceSignature(std::initializer_list<Variance> entries)
  {
    if (entries.size() > static_cast<std::size_t>(max_rank)) {
      throw pcms_error("VarianceSignature: rank must be 0.." +
                       std::to_string(max_rank) + ", got " +
                       std::to_string(entries.size()));
    }
    for (Variance v : entries) {
      entries_[static_cast<std::size_t>(rank_++)] = v;
    }
  }

  [[nodiscard]] constexpr int Rank() const noexcept { return rank_; }
  [[nodiscard]] constexpr Variance operator[](int i) const noexcept
  {
    return entries_[static_cast<std::size_t>(i)];
  }

  /// @return a copy with a covariant index appended
  /// @throws pcms_error if the rank is already max_rank
  [[nodiscard]] constexpr VarianceSignature WithCovariantIndexAppended() const
  {
    if (rank_ >= max_rank) {
      throw pcms_error("VarianceSignature::WithCovariantIndexAppended: rank " +
                       std::to_string(rank_) + " is already the maximum");
    }
    VarianceSignature appended = *this;
    appended.entries_[static_cast<std::size_t>(appended.rank_++)] =
      Variance::Covariant;
    return appended;
  }

  friend constexpr bool operator==(const VarianceSignature& a,
                                   const VarianceSignature& b) noexcept
  {
    if (a.rank_ != b.rank_) {
      return false;
    }
    for (int i = 0; i < a.rank_; ++i) {
      if (a.entries_[static_cast<std::size_t>(i)] !=
          b.entries_[static_cast<std::size_t>(i)]) {
        return false;
      }
    }
    return true;
  }

private:
  std::array<Variance, max_rank> entries_{};
  int rank_ = 0;
};

/// Sign under an orientation reversal: Odd values (pseudo-tensors, e.g. a
/// curl) change sign; Even values do not.
enum class Parity
{
  Even,
  Odd
};

/// Basis the components of a value are expressed in: coordinate system,
/// component scaling, index variance and parity. Curl produces Odd output;
/// reflections negate Odd values.
struct ValueBasis
{
  std::shared_ptr<const CoordinateSystem> system = nullptr;
  ComponentScaling scaling = ComponentScaling::Physical;
  VarianceSignature variance{};
  Parity parity = Parity::Even;

  [[nodiscard]] int Rank() const noexcept { return variance.Rank(); }

  friend bool operator==(const ValueBasis&, const ValueBasis&) = default;
};

[[nodiscard]] inline FieldValueType ValueTypeOfRank(int rank)
{
  switch (rank) {
    case 0: return FieldValueType::Scalar;
    case 1: return FieldValueType::Vector;
    case 2: return FieldValueType::Tensor;
    default: break;
  }
  throw pcms_error("ValueTypeOfRank: no FieldValueType label for rank " +
                   std::to_string(rank));
}

[[nodiscard]] inline std::optional<ComponentScaling> DerivedComponentScaling(
  const CoordinateSystem& system) noexcept
{
  if (!system.HasOrthogonalBasis()) {
    return ComponentScaling::Natural;
  }
  if (system.HasUnitScaleFactors()) {
    return ComponentScaling::Physical;
  }
  return std::nullopt;
}

[[nodiscard]] inline bool SameValueBasis(const ValueBasis& a,
                                         const ValueBasis& b) noexcept
{
  if (a.Rank() == 0 && b.Rank() == 0) {
    return true;
  }
  if (a.Rank() != b.Rank() || !SameCoordinateSystem(a.system, b.system) ||
      a.scaling != b.scaling) {
    return false;
  }
  return a.scaling == ComponentScaling::Physical || a.variance == b.variance;
}

namespace detail
{

/// Gate for a tagged write: throws pcms_error unless the basis the caller
/// claimed on the output view is the basis the callee is about to write.
/// @param context callee name used to prefix the error message
/// @param claimed basis tagged onto the caller-supplied output view
/// @param written basis the callee actually writes
inline void CheckWrittenValueBasis(const char* context,
                                   const ValueBasis& claimed,
                                   const ValueBasis& written)
{
  if (!SameValueBasis(claimed, written)) {
    throw pcms_error(std::string(context) +
                     ": the output view's declared value type/basis does not "
                     "match the basis this call writes");
  }
}

inline void ValidateValueBasis(const ValueBasis& basis)
{
  if (basis.Rank() == 0) {
    return;
  }
  if (basis.Rank() > 2) {
    throw pcms_error("ValidateValueBasis: rank-" +
                     std::to_string(basis.Rank()) +
                     " values are not yet supported (maximum supported rank "
                     "is 2)");
  }
  if (basis.system == nullptr) {
    throw pcms_error(
      "ValidateValueBasis: Vector/Tensor values require a basis coordinate "
      "system");
  }
  if (basis.scaling == ComponentScaling::Physical &&
      !basis.system->HasOrthogonalBasis()) {
    throw pcms_error(
      "ValidateValueBasis: coordinate system '" +
      std::string(basis.system->Kind()) +
      "' has a non-orthogonal basis, so physical components are undefined; "
      "use ComponentScaling::Natural, or borrow an orthogonal system's basis");
  }
}

inline void ValidateValueSemantics(const ValueBasis& basis, int num_components)
{
  ValidateValueBasis(basis);
  const int rank = basis.Rank();
  if (rank == 0) {
    return;
  }
  int expected = 1;
  for (int i = 0; i < rank; ++i) {
    expected *= basis.system->Dimension();
  }
  if (num_components != expected) {
    throw pcms_error("ValidateValueSemantics: rank-" + std::to_string(rank) +
                     " values on coordinate system '" +
                     std::string(basis.system->Kind()) + "' (dimension " +
                     std::to_string(basis.system->Dimension()) + ") require " +
                     std::to_string(expected) + " components, got " +
                     std::to_string(num_components));
  }
}

inline ValueBasis MakeValueBasis(
  const std::shared_ptr<const CoordinateSystem>& space_system,
  VarianceSignature variance,
  std::shared_ptr<const CoordinateSystem> system, // null == the space's own
  std::optional<ComponentScaling> scaling)
{
  if (variance.Rank() == 0) {
    return ValueBasis{};
  }
  if (system == nullptr) {
    system = space_system;
  }
  ComponentScaling resolved;
  if (scaling.has_value()) {
    resolved = *scaling;
  } else if (auto derived = DerivedComponentScaling(*system)) {
    resolved = *derived;
  } else {
    throw pcms_error(
      "CreateFunction: coordinate system '" + std::string(system->Kind()) +
      "' admits both physical and natural components and they differ, so "
      "the scaling must be stated (pass ComponentScaling::Physical or "
      "ComponentScaling::Natural)");
  }
  return ValueBasis{std::move(system), resolved, variance};
}

} // namespace detail

namespace values
{

inline constexpr VarianceSignature Scalar{};
inline constexpr VarianceSignature Vector{Variance::Contravariant};
inline constexpr VarianceSignature Covector{Variance::Covariant};
inline constexpr VarianceSignature Tensor{Variance::Contravariant,
                                          Variance::Contravariant};

[[nodiscard]] inline VarianceSignature Of(
  std::initializer_list<Variance> variance)
{
  return VarianceSignature{variance};
}

} // namespace values

/// @brief Basis of the gradient of values stored in `stored`.
///
/// Component (j..., i) is the derivative along coordinate i of stored
/// component (j...). Components are physical on orthogonal systems and
/// natural otherwise.
/// @param stored basis of the differentiated values
/// @param system coordinate system of the derivatives
/// @return the stored basis with one covariant index appended; the parity is
///         the stored basis's
/// @throws pcms_error if `system` is null, or if `stored` has rank above 0
///         and `stored` and `system` are not one Cartesian system (the
///         gradient of a rank-N field needs the coordinate system's
///         connection)
[[nodiscard]] inline ValueBasis GradientBasis(
  const ValueBasis& stored, std::shared_ptr<const CoordinateSystem> system)
{
  if (system == nullptr) {
    throw pcms_error("GradientBasis: derivatives require a coordinate system");
  }
  const ComponentScaling scaling = system->HasOrthogonalBasis()
                                     ? ComponentScaling::Physical
                                     : ComponentScaling::Natural;
  if (stored.Rank() == 0) {
    return ValueBasis{std::move(system), scaling, values::Covector,
                      stored.parity};
  }
  const bool cartesian =
    dynamic_cast<const csys::Cartesian*>(system.get()) != nullptr &&
    SameCoordinateSystem(stored.system, system);
  if (!cartesian) {
    throw pcms_error("GradientBasis: the gradient of a rank-" +
                     std::to_string(stored.Rank()) +
                     " field needs the coordinate system's connection, which "
                     "is not implemented");
  }
  return ValueBasis{std::move(system), scaling,
                    stored.variance.WithCovariantIndexAppended(),
                    stored.parity};
}

// view for field values that is tagged with the coordinate system
template <typename ElementType, typename MemorySpace,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<MemorySpace>>
class ValueView
{
public:
  using element_type = ElementType;

  ValueView(ValueBasis basis,
            Rank2View<ElementType, MemorySpace, LayoutPolicy> values)
    : basis_(std::move(basis)), values_(values)
  {
    detail::ValidateValueSemantics(basis_, static_cast<int>(values_.extent(1)));
  }

  template <typename OtherElement,
            typename =
              std::enable_if_t<!std::is_same_v<OtherElement, ElementType> &&
                               std::is_same_v<ElementType, const OtherElement>>>
  ValueView(const ValueView<OtherElement, MemorySpace, LayoutPolicy>& other)
    : basis_(other.GetBasis()), values_(other.GetValues())
  {
  }

  [[nodiscard]] int Rank() const noexcept { return basis_.Rank(); }
  [[nodiscard]] FieldValueType GetValueType() const
  {
    return ValueTypeOfRank(basis_.Rank());
  }

  [[nodiscard]] const ValueBasis& GetBasis() const noexcept { return basis_; }

  [[nodiscard]] Rank2View<ElementType, MemorySpace, LayoutPolicy> GetValues()
    const noexcept
  {
    return values_;
  }

  [[nodiscard]] size_t extent(size_t r) const noexcept
  {
    return values_.extent(r);
  }
  [[nodiscard]] size_t size() const noexcept { return values_.size(); }
  // Index through GetValues(): a ValueView carries a shared_ptr and must not
  // be captured by a device kernel.
  decltype(auto) operator()(size_t i, size_t c) const = delete;

private:
  ValueBasis basis_;
  Rank2View<ElementType, MemorySpace, LayoutPolicy> values_;
};

} // namespace pcms

#endif // PCMS_FIELD_VALUE_VIEW_HPP
