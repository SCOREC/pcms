#ifndef PCMS_LOCALIZATION_LOCALIZATION_FACTORY_H
#define PCMS_LOCALIZATION_LOCALIZATION_FACTORY_H

#include "pcms/localization/adj_search.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "pcms/utility/memory_spaces.h"

namespace pcms
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

} // namespace pcms

#endif // PCMS_LOCALIZATION_LOCALIZATION_FACTORY_H
