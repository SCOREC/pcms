#ifndef PCMS_GEOMETRY_DISTRIBUTED_BINDING_HPP
#define PCMS_GEOMETRY_DISTRIBUTED_BINDING_HPP

#include "pcms/geometry/bound_points.hpp"
#include "pcms/geometry/distributed_geometry.hpp"
#include "pcms/geometry/distributed_routing_plan.hpp"

#include <cstddef>
#include <memory>

namespace pcms
{

/**
 * \brief Bound points of a DistributedGeometry.
 *
 * Holds the routing plan of the query points, and two bound points of the
 * local geometry: one of this rank's points it owns, one of the points other
 * ranks sent it. A query point is found if the rank that owns it found it.
 * Rebinding locates the new points over the decomposition again, so it is
 * collective.
 */
class DistributedBinding final : public BoundPoints
{
public:
  // Collective.
  DistributedBinding(std::shared_ptr<const DistributedGeometry> geometry,
                     const BindRequest& request);

  [[nodiscard]] const DistributedRoutingPlan& GetRoutingPlan() const noexcept
  {
    return plan_;
  }

  /** Local bound points of the plan's own points, in their order. */
  [[nodiscard]] const std::shared_ptr<const BoundPoints>& OwnBound()
    const noexcept
  {
    return own_bound_;
  }

  /** Local bound points of the plan's incoming points, in their order. */
  [[nodiscard]] const std::shared_ptr<const BoundPoints>& IncomingBound()
    const noexcept
  {
    return incoming_bound_;
  }

  /** How many query points get their value from another rank. */
  [[nodiscard]] std::size_t RoutedPointCount() const noexcept
  {
    return plan_.RoutedPointCount();
  }

protected:
  void RebindImpl() override { Locate({}); }

private:
  // Locates LocatedPointValues() over the decomposition, starting from
  // `owners`, and sets the found points. Collective.
  void Locate(const OwnerHints& owners);

  const DistributedGeometry& distributed_;
  DistributedRoutingPlan plan_;
  std::shared_ptr<const BoundPoints> own_bound_;
  std::shared_ptr<const BoundPoints> incoming_bound_;
};

} // namespace pcms

#endif // PCMS_GEOMETRY_DISTRIBUTED_BINDING_HPP
