#ifndef PCMS_LOCALIZATION_LOCALIZATION_FACTORY_HPP
#define PCMS_LOCALIZATION_LOCALIZATION_FACTORY_HPP

#include "pcms/localization/support_results.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/utility/memory_spaces.hpp"

namespace pcms::localization
{

// Abstract interface for building MLS support structures.
//
// A LocalizationFactory encapsulates the source geometry and search strategy
// (N^2 point-cloud search or mesh-adjacency BFS). Concrete implementations
// are created at FunctionSpace construction time and reused across repeated
// CreatePointEvaluator calls with different target coordinate sets.
class LocalizationFactory
{
public:
  virtual ~LocalizationFactory() = default;

  virtual SupportResults Build(
    CoordinateView<DeviceMemorySpace> target_coords) const = 0;
};

} // namespace pcms::localization

#endif // PCMS_LOCALIZATION_LOCALIZATION_FACTORY_HPP
