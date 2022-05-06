#ifndef WDM_COUPLING_BOUNDING_BOX_H
#define WDM_COUPLING_BOUNDING_BOX_H
#include "wdmcpl/types.h"
#include <array>
namespace wdmcpl
{

template <int DIM = 2>
struct AABBox
{
  static constexpr int dim = DIM;
  std::array<Real, dim> center;
  // half length of bounding box
  std::array<Real, dim> half_width;
};

template <int dim>
bool intersects(const AABBox<dim>& a, const AABBox<dim>& b)
{
  for (int i = 0; i < dim; ++i) {
    if (std::abs(a.center[i] - b.center[i]) >
        (a.half_width[i] + b.half_width[i]))
      return false;
  }
  return true;
}
} // namespace wdmcpl

#endif // WDM_COUPLING_BOUNDING_BOX_H
