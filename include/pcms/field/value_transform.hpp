#ifndef PCMS_FIELD_VALUE_TRANSFORM_HPP
#define PCMS_FIELD_VALUE_TRANSFORM_HPP

#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <memory>

namespace pcms
{

/// @brief A linear operator on values at bound points: one matrix per point.
///
/// Users obtain it from a bound map (BoundCoordinateMap::ValueTransform), or
/// never see it: the evaluators apply it. Apply picks its matrix by the
/// variance of the values: scalars (rank 0, including coordinates) pass
/// through unchanged; vectors and covectors (rank 1) are transformed;
/// values already in the output coordinate system pass through. ApplyTranspose
/// is the dual map under the plain pairing of values (sum over components, no
/// metric): <T v, w> = <v, T^T w>; it is what the transpose of an evaluation
/// (EvaluateTranspose) applies.
///
/// Values may hold several blocks per point (e.g. the gradient of a field
/// with several scalar components): each block of the input basis's
/// dimension is transformed on its own.
///
/// At a point where a scale factor of an orthogonal coordinate system is zero
/// (the axis R = 0 of the cylindrical coordinate system), the factor is taken
/// as 1: a map whose Jacobian there uses the same convention (the maps between
/// coordinate systems do) gives the components in the frame of the angle it
/// assigned to the point.
class ValueTransform
{
public:
  /// @param in basis of the input values
  /// @return the basis Apply writes for `in`
  /// @throws pcms_error if the transform cannot re-express `in`: rank 2 or
  ///         more, a coordinate system other than its input coordinate system,
  ///         or a variance whose matrix needs the inverse of a Jacobian that is
  ///         not square
  [[nodiscard]] ValueBasis OutputBasis(const ValueBasis& in) const;

  /// @brief out = T in, at every point.
  /// @param in [NumPoints()][components] values
  /// @param out [NumPoints()][components of OutputBasis(in)], tagged with
  ///        OutputBasis(in.GetBasis()); may be `in`'s storage when the input
  ///        and output have as many components
  /// @throws pcms_error if the bases or shapes do not match
  void Apply(ValueView<const Real, DeviceMemorySpace> in,
             ValueView<Real, DeviceMemorySpace> out) const;

  /// @brief out = T^T in, at every point: the transpose of Apply for values
  ///        in the basis of `out`.
  /// @param in [NumPoints()][components], tagged with
  ///        OutputBasis(out.GetBasis())
  /// @param out [NumPoints()][components], tagged with the input basis of
  ///        the transform it is the transpose of; may be `in`'s storage when
  ///        both have as many components
  /// @throws pcms_error if the bases or shapes do not match
  void ApplyTranspose(ValueView<const Real, DeviceMemorySpace> in,
                      ValueView<Real, DeviceMemorySpace> out) const;

  /// @return number of points
  [[nodiscard]] LO NumPoints() const noexcept { return num_points_; }

private:
  // how the matrix relates to J: along the map (Pushforward), against it
  // (Pullback), or the identity (Scale)
  enum class Direction : std::uint8_t
  {
    With,
    Against,
    None
  };

  ValueTransform(
    Direction direction,
    View<3, const Real, DeviceMemorySpace, Kokkos::layout_stride> jacobian,
    std::shared_ptr<const CoordinateSystem> from,
    Kokkos::View<const Real**, DeviceMemorySpace> from_scale,
    std::shared_ptr<const CoordinateSystem> to,
    Kokkos::View<const Real**, DeviceMemorySpace> to_scale,
    ComponentScaling scale_to);

  // whether the matrix of `variance` needs an inverse of J
  [[nodiscard]] bool NeedsInverse(Variance variance) const noexcept;

  // Apply (transpose false) or ApplyTranspose of values whose input-side
  // basis is `in_side` and whose output side has `out_scaling`
  void Run(const ValueBasis& in_side, ComponentScaling out_scaling,
           bool transpose, Rank2View<const Real, DeviceMemorySpace> in,
           Rank2View<Real, DeviceMemorySpace> out) const;

  // the factories (detail::Pullback, detail::Pushforward, detail::Scale)
  friend struct ValueTransformAccess;

  Direction direction_;
  View<3, const Real, DeviceMemorySpace, Kokkos::layout_stride> jacobian_;
  std::shared_ptr<const CoordinateSystem> from_;
  std::shared_ptr<const CoordinateSystem> to_;
  // the scale factors of the input and output coordinate systems at the
  // points, [points][2 dimension]: sqrt(g_ii), then sqrt(g^ii); empty where
  // physical and natural components coincide or there is no metric
  Kokkos::View<const Real**, DeviceMemorySpace> from_scale_;
  Kokkos::View<const Real**, DeviceMemorySpace> to_scale_;
  // the output scaling of a Scale transform
  ComponentScaling scale_to_;
  LO num_points_;
};

} // namespace pcms

#include "pcms/field/detail/value_transform.hpp"

#endif // PCMS_FIELD_VALUE_TRANSFORM_HPP
