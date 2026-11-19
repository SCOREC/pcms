#ifndef PCMS_UTILITY_BOUNDING_BOX_HPP
#define PCMS_UTILITY_BOUNDING_BOX_HPP
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <array>
namespace pcms
{

template <int DIM = 2>
struct AABBox
{
  static constexpr int dim = DIM;
  Kokkos::Array<Real, dim> center;
  // half length of bounding box
  Kokkos::Array<Real, dim> half_width;
};

} // namespace pcms

#include "pcms/utility/detail/bounding_box.hpp"

#endif // PCMS_UTILITY_BOUNDING_BOX_HPP
