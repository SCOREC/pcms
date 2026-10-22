#ifndef PCMS_FIELD_COORDINATE_MAP_HPP
#define PCMS_FIELD_COORDINATE_MAP_HPP

#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include <Kokkos_Core.hpp>
#include <memory>

namespace pcms
{

// A CoordinateMap bound to one query point set. Bind performs every per-point
// computation the map needs: the mapped coordinates and the state that
// re-expressing component values at those points requires. That state is the
// map's differential, which does not depend on the value basis -- the basis
// only selects which arithmetic is applied to it -- so one bound map serves
// every field evaluated at that point set.
//
// Rank-0 values never reach OutputBasis or TransformValues: scalars are basis
// invariant, so callers forward them without consulting the map.
class BoundCoordinateMap
{
public:
  BoundCoordinateMap(std::shared_ptr<const CoordinateSystem> system,
                     Kokkos::View<Real**, DeviceMemorySpace> mapped);

  /// Query points expressed in the map's target coordinate system.
  [[nodiscard]] CoordinateView<DeviceMemorySpace> MappedPoints() const;

  [[nodiscard]] LO NumPoints() const noexcept
  {
    return static_cast<LO>(mapped_.extent(0));
  }

  /// Value basis produced for rank >= 1 values stored in `stored`. Throws when
  /// this map cannot re-express that basis.
  [[nodiscard]] virtual ValueBasis OutputBasis(const ValueBasis& stored) const;

  /// Re-expresses `values` in place and returns the same buffer tagged
  /// `OutputBasis(values.GetBasis())`.
  virtual ValueView<Real, DeviceMemorySpace> TransformValues(
    ValueView<Real, DeviceMemorySpace> values) const;

  virtual ~BoundCoordinateMap() = default;

  BoundCoordinateMap(const BoundCoordinateMap&) = delete;
  BoundCoordinateMap& operator=(const BoundCoordinateMap&) = delete;

protected:
  /// Throws unless `values` holds NumPoints() rows in a basis this map
  /// re-expresses; returns the basis the transformed values are in.
  ValueBasis ValidateTransformValues(
    const char* who, const ValueView<Real, DeviceMemorySpace>& values) const;

private:
  std::shared_ptr<const CoordinateSystem> system_;
  Kokkos::View<Real**, DeviceMemorySpace> mapped_;
};

class CoordinateMap
{
public:
  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetSourceCoordinateSystem() const noexcept
  {
    return source_;
  }
  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetTargetCoordinateSystem() const noexcept
  {
    return target_;
  }

  /// Performs every per-point computation for `query_points`, which must be in
  /// this map's source coordinate system. The returned object owns its
  /// results and is complete on return, so `query_points` may be released.
  [[nodiscard]] std::unique_ptr<BoundCoordinateMap> Bind(
    const CoordinateView<DeviceMemorySpace>& query_points) const;

  virtual ~CoordinateMap() = default;

protected:
  CoordinateMap(std::shared_ptr<const CoordinateSystem> source,
                std::shared_ptr<const CoordinateSystem> target);

  /// Binds `query_points`, already checked to be in the source coordinate
  /// system. May leave kernels in flight; Bind fences after it returns.
  [[nodiscard]] virtual std::unique_ptr<BoundCoordinateMap> BindImpl(
    const CoordinateView<DeviceMemorySpace>& query_points) const = 0;

private:
  std::shared_ptr<const CoordinateSystem> source_;
  std::shared_ptr<const CoordinateSystem> target_;
};

// Maps points Cartesian -> cylindrical, so values run the other way: a field
// stored in cylindrical components is re-expressed in Cartesian ones.
class CartesianToCylindrical final : public CoordinateMap
{
public:
  CartesianToCylindrical();

protected:
  [[nodiscard]] std::unique_ptr<BoundCoordinateMap> BindImpl(
    const CoordinateView<DeviceMemorySpace>& query_points) const override;
};

class CylindricalToCartesian final : public CoordinateMap
{
public:
  CylindricalToCartesian();

protected:
  [[nodiscard]] std::unique_ptr<BoundCoordinateMap> BindImpl(
    const CoordinateView<DeviceMemorySpace>& query_points) const override;
};

} // namespace pcms

#endif // PCMS_FIELD_COORDINATE_MAP_HPP
