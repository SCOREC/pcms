#ifndef PCMS_FIELD_COORDINATE_MAP_HPP
#define PCMS_FIELD_COORDINATE_MAP_HPP

#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/point_status.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include <Kokkos_Core.hpp>
#include <memory>

namespace pcms
{

// A CoordinateMap bound to one query point set. Bind performs every per-point
// computation the map needs: the mapped coordinates, their status, and the
// state that re-expressing component values at those points requires. That
// state is the map's differential, which does not depend on the value basis --
// the basis only selects which arithmetic is applied to it -- so one bound map
// serves every field evaluated at that point set.
//
// Rank-0 values never reach OutputBasis or TransformValues: scalars are basis
// invariant, so callers forward them without consulting the map.
class BoundCoordinateMap
{
public:
  BoundCoordinateMap(std::shared_ptr<const CoordinateSystem> system,
                     Kokkos::View<Real**, DeviceMemorySpace> mapped,
                     PointStatusView status = {});

  /// Query points expressed in the map's target coordinate system.
  [[nodiscard]] CoordinateView<DeviceMemorySpace> MappedPoints() const;

  /// Per-point outcome of the mapping; empty when every point mapped cleanly.
  [[nodiscard]] const PointStatusView& Status() const noexcept
  {
    return status_;
  }

  [[nodiscard]] LO NumPoints() const noexcept
  {
    return static_cast<LO>(mapped_.extent(0));
  }

  /// Value basis produced for rank >= 1 values stored in `stored`. Throws when
  /// this map cannot re-express that basis.
  [[nodiscard]] virtual ValueBasis OutputBasis(const ValueBasis& stored) const;

  /// Re-expresses `in` into `out`, which must be tagged with
  /// `OutputBasis(in.GetBasis())`.
  virtual void TransformValues(ValueView<const Real, DeviceMemorySpace> in,
                               ValueView<Real, DeviceMemorySpace> out) const;

  virtual ~BoundCoordinateMap() = default;

  BoundCoordinateMap(const BoundCoordinateMap&) = delete;
  BoundCoordinateMap& operator=(const BoundCoordinateMap&) = delete;

protected:
  /// Throws unless both views hold NumPoints() rows of equal width and `out`
  /// is tagged OutputBasis(in.GetBasis()).
  void ValidateTransformViews(
    const char* who, const ValueView<const Real, DeviceMemorySpace>& in,
    const ValueView<Real, DeviceMemorySpace>& out) const;

private:
  std::shared_ptr<const CoordinateSystem> system_;
  Kokkos::View<Real**, DeviceMemorySpace> mapped_;
  PointStatusView status_;
};

class CoordinateMap
{
public:
  [[nodiscard]] virtual std::shared_ptr<const CoordinateSystem>
  GetSourceCoordinateSystem() const noexcept = 0;
  [[nodiscard]] virtual std::shared_ptr<const CoordinateSystem>
  GetTargetCoordinateSystem() const noexcept = 0;

  /// Performs every per-point computation for `query_points`, which must be in
  /// this map's source coordinate system. The returned object owns its
  /// results; `query_points` may be released once Bind returns.
  [[nodiscard]] virtual std::unique_ptr<BoundCoordinateMap> Bind(
    const CoordinateView<DeviceMemorySpace>& query_points) const = 0;

  virtual ~CoordinateMap() = default;

protected:
  void ValidateSourcePoints(const CoordinateView<DeviceMemorySpace>& points,
                            const char* who) const;
};

// Maps points Cartesian -> cylindrical, so values run the other way: a field
// stored in cylindrical components is re-expressed in Cartesian ones.
class CartesianToCylindrical final : public CoordinateMap
{
public:
  [[nodiscard]] std::shared_ptr<const CoordinateSystem>
  GetSourceCoordinateSystem() const noexcept override;
  [[nodiscard]] std::shared_ptr<const CoordinateSystem>
  GetTargetCoordinateSystem() const noexcept override;
  [[nodiscard]] std::unique_ptr<BoundCoordinateMap> Bind(
    const CoordinateView<DeviceMemorySpace>& query_points) const override;
};

class CylindricalToCartesian final : public CoordinateMap
{
public:
  [[nodiscard]] std::shared_ptr<const CoordinateSystem>
  GetSourceCoordinateSystem() const noexcept override;
  [[nodiscard]] std::shared_ptr<const CoordinateSystem>
  GetTargetCoordinateSystem() const noexcept override;
  [[nodiscard]] std::unique_ptr<BoundCoordinateMap> Bind(
    const CoordinateView<DeviceMemorySpace>& query_points) const override;
};

} // namespace pcms

#endif // PCMS_FIELD_COORDINATE_MAP_HPP
