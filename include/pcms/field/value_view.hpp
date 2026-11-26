#ifndef PCMS_FIELD_VALUE_VIEW_HPP
#define PCMS_FIELD_VALUE_VIEW_HPP

#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
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

/// @brief How the components are scaled: Natural on the coordinate basis,
///        Physical on its unit vectors (natural times sqrt(g_ii) per
///        contravariant index, times sqrt(g^ii) per covariant index)
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

/// @brief What the components of a value are.
enum class ValueKind : std::uint8_t
{
  /// components of a tensor of the variance's rank on the coordinates of
  /// `system`, transformed by J or J^-T per index; rank 0 is one component
  Tensor,
  /// a point of `system`, one component per coordinate, moved by the
  /// coordinate map
  Coordinates,
  /// any number of components, never transformed
  Array
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
  ValueKind kind = ValueKind::Tensor;

  [[nodiscard]] int Rank() const noexcept { return variance.Rank(); }

  /// @brief The basis of values that are points of `system`: the position of
  ///        an embedded geometry.
  ///
  /// The Gradient of a Function with this basis is the Jacobian of the map
  /// it describes, [component][derivative]: column j * (parameter dimension)
  /// + i is d x_j / d xi_i.
  [[nodiscard]] static ValueBasis Coordinates(
    std::shared_ptr<const CoordinateSystem> system)
  {
    return ValueBasis{std::move(system), ComponentScaling::Natural,
                      VarianceSignature{}, Parity::Even,
                      ValueKind::Coordinates};
  }

  friend bool operator==(const ValueBasis&, const ValueBasis&) = default;
};

namespace values
{

inline constexpr VarianceSignature Scalar{};
inline constexpr VarianceSignature Vector{Variance::Contravariant};
inline constexpr VarianceSignature Covector{Variance::Covariant};
inline constexpr VarianceSignature Tensor{Variance::Contravariant,
                                          Variance::Contravariant};

} // namespace values

} // namespace pcms

// the checks and conversions of value bases, which ValueView uses
#include "pcms/field/detail/value_view.hpp"

namespace pcms
{

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
    detail::ValidateValueSemantics(basis_, static_cast<int>(values_.extent(1)),
                                   /*blocks=*/true);
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
    return detail::ValueTypeOfRank(basis_.Rank());
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
