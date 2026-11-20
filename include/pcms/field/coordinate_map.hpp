#ifndef PCMS_FIELD_COORDINATE_MAP_HPP
#define PCMS_FIELD_COORDINATE_MAP_HPP

#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/field/value_transform.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <optional>

namespace pcms
{

// A CoordinateMap bound to one query point set. Bind performs every per-point
// computation the map needs: the mapped coordinates and the map's Jacobian
// there. Values re-expressed at those points flow against the map, from the
// target coordinate system where the points were located back to the source
// system of the query points, through ValueTransform(): the Pullback of the
// Jacobian, which does not depend on the value basis -- the basis only
// selects the matrix -- so one bound map serves every field evaluated at
// that point set. Scalars pass through it unchanged.
class BoundCoordinateMap
{
public:
  BoundCoordinateMap(std::shared_ptr<const CoordinateSystem> system,
                     Kokkos::View<Real**, DeviceMemorySpace> mapped);

  /// Query points expressed in the map's target coordinate system.
  [[nodiscard]] CoordinateView<DeviceMemorySpace> MappedPoints() const;

  /// @return the mapped points, shared, not copied
  [[nodiscard]] Kokkos::View<const Real**, DeviceMemorySpace> Mapped()
    const noexcept
  {
    return mapped_;
  }

  [[nodiscard]] LO NumPoints() const noexcept
  {
    return static_cast<LO>(mapped_.extent(0));
  }

  /// @return J(p, j, i) = d x_j / d xi_i, the Jacobian of the map at the
  ///         bound points (xi: source coordinates, x: mapped coordinates),
  ///         [points][target dimension][source dimension]
  [[nodiscard]] virtual View<3, const Real, DeviceMemorySpace> Jacobian()
    const = 0;

  /// @return the transform of values evaluated at the mapped points into the
  ///         source coordinate system of the query points: the Pullback of
  ///         Jacobian() (covectors by J^T; vectors by J^-1, which needs a
  ///         square J); built on the first call (like the lazy factors of a
  ///         BoundPoints, not thread-safe)
  [[nodiscard]] const pcms::ValueTransform& ValueTransform() const;

  virtual ~BoundCoordinateMap() = default;

  BoundCoordinateMap(const BoundCoordinateMap&) = delete;
  BoundCoordinateMap& operator=(const BoundCoordinateMap&) = delete;

private:
  friend class CoordinateMap;

  std::shared_ptr<const CoordinateSystem> system_;
  Kokkos::View<Real**, DeviceMemorySpace> mapped_;
  // the source system, set by CoordinateMap::Bind
  std::shared_ptr<const CoordinateSystem> source_;
  mutable std::optional<pcms::ValueTransform> value_transform_;
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
  /// this map's source coordinate system: the mapped points and the Jacobian,
  /// from which the value transform is built on first use. The returned
  /// object owns its results and is complete on return, so `query_points`
  /// may be released.
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

/// @brief A coordinate map given by a kernel that writes the image of a point
///        and the Jacobian there: the way to write a closed-form map between
///        coordinate systems.
///
/// `kernel` is a copyable functor with, for one source dimension S and target
/// dimension T (at most 3),
/// `KOKKOS_FUNCTION void operator()(const Kokkos::Array<Real, S>& x,
/// Kokkos::Array<Real, T>& out, Kokkos::Array<Real, T * S>& jacobian) const`,
/// with jacobian[S * j + i] = d out_j / d x_i.
/// @throws pcms_error if a system is null or its dimension is not the
///         kernel's
template <typename Kernel>
[[nodiscard]] std::shared_ptr<const CoordinateMap> CreateCoordinateMap(
  std::shared_ptr<const CoordinateSystem> source,
  std::shared_ptr<const CoordinateSystem> target, Kernel kernel);

/// Map from Cartesian (x, y, z) to cylindrical (R, theta, Z) points, so values
/// run the other way: a field stored in cylindrical components is re-expressed
/// in Cartesian ones. On the axis (R = 0) theta is atan2(0, 0) and the
/// components are those of the frame of that angle.
[[nodiscard]] std::shared_ptr<const CoordinateMap>
CreateCartesianToCylindrical();

/// Map from cylindrical (R, theta, Z) to Cartesian (x, y, z) points.
[[nodiscard]] std::shared_ptr<const CoordinateMap>
CreateCylindricalToCartesian();

/// @brief The map that inserts coordinate `axis` of `target`, fixed at
///        `value`, among the coordinates of `source`, e.g. (R, Z) ->
///        (R, theta0, Z) from CylindricalRZ to CylindricalRThetaZ with axis 1.
///
/// Scalars and covectors pass (covectors by J^T); vectors would need the
/// inverse of its rectangular Jacobian and are refused.
/// @throws pcms_error if a system is null, dim(target) != dim(source) + 1,
///         dim(target) > 3, or `axis` is not a coordinate of `target`
[[nodiscard]] std::shared_ptr<const CoordinateMap> CreateCoordinateEmbedding(
  std::shared_ptr<const CoordinateSystem> source,
  std::shared_ptr<const CoordinateSystem> target, int axis, Real value);

} // namespace pcms

#include "pcms/field/detail/coordinate_map.hpp"

namespace pcms
{

template <typename Kernel>
std::shared_ptr<const CoordinateMap> CreateCoordinateMap(
  std::shared_ptr<const CoordinateSystem> source,
  std::shared_ptr<const CoordinateSystem> target, Kernel kernel)
{
  constexpr auto dims = detail::KernelDimensions<Kernel>();
  static_assert(dims[2] == 1,
                "CreateCoordinateMap: the kernel must take (const "
                "Kokkos::Array<Real, S>&, Kokkos::Array<Real, T>&, "
                "Kokkos::Array<Real, T * S>&) for one S and T in 1..3");
  return std::make_shared<const detail::KernelMap<Kernel, dims[0], dims[1]>>(
    std::move(source), std::move(target), std::move(kernel));
}

} // namespace pcms

#endif // PCMS_FIELD_COORDINATE_MAP_HPP
