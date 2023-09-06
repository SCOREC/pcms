#ifndef WDM_COUPLING_BOUNDING_BOX_H
#define WDM_COUPLING_BOUNDING_BOX_H
#include "pcms/types.h"
#include <Kokkos_Core.hpp>
#include <array>
namespace pcms
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
KOKKOS_INLINE_FUNCTION
bool intersects(const AABBox<dim>& a, const AABBox<dim>& b)
{
  for (int i = 0; i < dim; ++i) {
    if (std::abs(a.center[i] - b.center[i]) >
        (a.half_width[i] + b.half_width[i]))
      return false;
  }
  return true;
}
} // namespace pcms

#endif // WDM_COUPLING_BOUNDING_BOX_H
