#include "pcms/transfer/interpolator.hpp"
#include "pcms/field/coordinate_map.hpp"

namespace pcms
{

namespace detail
{

bool SameDOFHolderSites(const BoundPoints& bound, const FieldLayout& layout)
{
  const auto sites = layout.GetDOFHolderCoordinates();
  if (!pcms::detail::SameCoordinateSystem(bound.GetCoordinateSystem(),
                                          sites.GetCoordinateSystem())) {
    return false;
  }
  const auto vs = sites.GetValues();
  const auto num_points = static_cast<LO>(vs.extent(0));
  if (bound.NumPoints() != num_points) {
    return false;
  }
  // a bind with a map located the mapped points
  const auto va = bound.LocatedPointValues();
  const bool mapped = bound.GetCoordinateMap() != nullptr;
  const auto vb = mapped ? bound.GetCoordinateMap()->Bind(sites)->Mapped()
                         : Kokkos::View<const Real**, DeviceMemorySpace>();
  const auto num_dims = static_cast<int>(va.extent(1));
  LO num_different = 0;
  Kokkos::parallel_reduce(
    "same_dof_holder_sites",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_points),
    KOKKOS_LAMBDA(const LO i, LO& count) {
      for (int d = 0; d < num_dims; ++d) {
        if (va(i, d) != (mapped ? vb(i, d) : vs(i, d))) {
          ++count;
          return;
        }
      }
    },
    num_different);
  return num_different == 0;
}

} // namespace detail

} // namespace pcms
