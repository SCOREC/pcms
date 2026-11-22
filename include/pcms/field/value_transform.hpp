#ifndef PCMS_FIELD_VALUE_TRANSFORM_HPP
#define PCMS_FIELD_VALUE_TRANSFORM_HPP

#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <optional>

namespace pcms
{

/// @brief A linear operator on values at bound points: one matrix per point
///        and variance.
///
/// Users obtain it from a bound map (BoundCoordinateMap::ValueTransforms), or
/// never see it: the evaluators apply it. Apply picks its matrix by the
/// variance of the values: vectors and covectors (rank 1) in its input
/// coordinate system are transformed; scalars (rank 0, including coordinates)
/// and values in another coordinate system pass through unchanged.
/// ApplyTranspose is the dual map under the plain pairing of values (sum over
/// components, no metric): <T v, w> = <v, T^T w>; it is what the transpose of
/// an evaluation (EvaluateTranspose) applies.
///
/// Values may hold several blocks per point (e.g. the gradient of a field
/// with several scalar components): each block of the input basis's
/// dimension is transformed on its own.
///
/// At a point where a scale factor sqrt(g_ii) or sqrt(g^ii) is zero (the axis
/// R = 0 of the cylindrical coordinate system), the factor is taken as 1: a
/// map whose Jacobian there uses the same convention (the maps between
/// coordinate systems do) gives the components in the frame of the angle it
/// assigned to the point.
class ValueTransform
{
public:
  /// @param in basis of the input values
  /// @return the basis Apply writes for `in`
  /// @throws pcms_error if the transform cannot re-express `in`: rank 2 or
  ///         more, physical components it does not convert (convert them
  ///         to natural first), or a variance whose matrix needs the inverse
  ///         of a Jacobian that is not square
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
  using Matrices = Kokkos::View<const Real***, DeviceMemorySpace>;

  ValueTransform(Matrices vector, Matrices covector,
                 std::shared_ptr<const CoordinateSystem> from,
                 std::shared_ptr<const CoordinateSystem> to,
                 std::optional<ComponentScaling> scaling);

  // whether values in `in` leave the transform unchanged
  [[nodiscard]] bool PassesThrough(const ValueBasis& in) const;

  // Apply (transpose false) or ApplyTranspose of values whose input-side basis
  // is `in_side`
  void Run(const char* who, const ValueBasis& in_side, bool transpose,
           Rank2View<const Real, DeviceMemorySpace> in,
           Rank2View<Real, DeviceMemorySpace> out) const;

  // the factories (detail::Pullback, Pushforward, ToNatural, ToPhysical)
  friend struct ValueTransformAccess;

  // the matrices of vectors and covectors, [points][dim to][dim from]; empty
  // for a variance that needs the inverse of a Jacobian that is not square
  Matrices vector_;
  Matrices covector_;
  std::shared_ptr<const CoordinateSystem> from_;
  std::shared_ptr<const CoordinateSystem> to_;
  // the scaling a ToNatural or ToPhysical writes; unset for the transform of
  // a map, which takes and writes natural components
  std::optional<ComponentScaling> scaling_;
  LO num_points_;
};

} // namespace pcms

#include "pcms/field/detail/value_transform.hpp"

#endif // PCMS_FIELD_VALUE_TRANSFORM_HPP
