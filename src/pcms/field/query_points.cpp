#include "pcms/field/query_points.hpp"
#include "pcms/field/coordinate_map.hpp"

namespace pcms
{

QueryPoints::QueryPoints(std::shared_ptr<const Geometry> geometry,
                         OutOfBoundsPolicy policy,
                         std::shared_ptr<const CoordinateMap> map)
  : geometry_(std::move(geometry)), policy_(policy), map_(std::move(map))
{
  if (geometry_ == nullptr) {
    throw pcms_error("QueryPoints: the geometry must not be null (a space "
                     "without a geometry does not evaluate at points)");
  }
}

void QueryPoints::Locate(CoordinateView<DeviceMemorySpace> coords)
{
  // the evaluators hold the bound source; released, Rebind can reuse it
  evaluators_.clear();
  if (bound_ == nullptr) {
    bound_ = geometry_->Bind(std::move(coords), policy_, map_);
  } else {
    bound_ = geometry_->Rebind(std::move(bound_), std::move(coords));
  }
}

} // namespace pcms
