#ifndef PCMS_BACKENDS_POINT_CLOUD_MLS_EVALUATOR_HPP
#define PCMS_BACKENDS_POINT_CLOUD_MLS_EVALUATOR_HPP

#include "pcms/backends/point_cloud/mls_interpolation.hpp"
#include "pcms/geometry/bound_points.hpp"
#include "pcms/geometry/geometry.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"

#include <Kokkos_Core.hpp>

#include <string>

namespace pcms::detail
{

/// Whether each of `n` query points has supports, from the supports' row
/// pointers. (A free function: CUDA does not allow device lambdas in
/// constructors or non-public member functions.)
inline Kokkos::View<bool*, DeviceMemorySpace> HasSupports(
  const Kokkos::View<const LO*, DeviceMemorySpace>& supports_ptr, LO n)
{
  Kokkos::View<bool*, DeviceMemorySpace> found(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "mls_found"), n);
  pcms::detail::ForEach(
    "MLSBinding::found", n, KOKKOS_LAMBDA(const LO i) {
      found(i) = supports_ptr[i + 1] > supports_ptr[i];
    });
  return found;
}

/// The rows of `coords` (n x dim) as one flat array [x0, y0, ..., x1, ...].
template <typename View2D>
Kokkos::View<Real*, DeviceMemorySpace> FlattenRows(const View2D& coords,
                                                   const char* label)
{
  const auto n = static_cast<LO>(coords.extent(0));
  const auto dim = static_cast<LO>(coords.extent(1));
  Kokkos::View<Real*, DeviceMemorySpace> flat(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(label)),
    n * dim);
  Kokkos::parallel_for(
    "FlattenRows",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const LO i) {
      for (LO d = 0; d < dim; ++d) {
        flat(i * dim + d) = coords(i, d);
      }
    });
  return flat;
}

/// @brief MLS BoundPoints: the supports of the query points and the
///        coordinates the solve uses. A query point without supports (e.g.
///        one the adjacency search has no start for) is missing.
class MLSBinding final : public BoundPoints
{
public:
  /// @param geometry geometry the points are localized on
  /// @param request query points and out-of-bounds policy
  /// @param source_coords flattened coordinates of the source points
  /// @param supports MLS supports of the query points
  /// @param site_dim the site dimension and
  /// @param options the options `supports` were found with; a rebind finds
  ///        the supports of its points with the geometry's Neighbourhood
  MLSBinding(std::shared_ptr<const Geometry> geometry,
             const BindRequest& request,
             Kokkos::View<const Real*, DeviceMemorySpace> source_coords,
             SupportResults supports, int site_dim,
             NeighbourhoodOptions options)
    : BoundPoints(std::move(geometry), request),
      source_coords_(std::move(source_coords)),
      target_coords_(
        detail::FlattenRows(LocatedPointValues(), "mls_target_coords")),
      supports_(std::move(supports)),
      site_dim_(site_dim),
      options_(options)
  {
    SetFoundFromSupports();
  }

  /// @return flattened coordinates of the source points
  [[nodiscard]] const Kokkos::View<const Real*, DeviceMemorySpace>&
  SourceCoordinates() const noexcept
  {
    return source_coords_;
  }

  /// @return flattened coordinates of the query points in the geometry's
  ///         coordinate system (LocatedPoints())
  [[nodiscard]] const Kokkos::View<const Real*, DeviceMemorySpace>&
  TargetCoordinates() const noexcept
  {
    return target_coords_;
  }

  /// @return MLS supports of the query points
  [[nodiscard]] const SupportResults& Supports() const noexcept
  {
    return supports_;
  }

protected:
  // the arrays are reassigned, not written, so the default UnsharedImpl holds
  void RebindImpl() override
  {
    supports_ =
      GetGeometry()->Neighbourhood(LocatedPoints(), site_dim_, options_);
    target_coords_ =
      detail::FlattenRows(LocatedPointValues(), "mls_target_coords");
    SetFoundFromSupports();
  }

private:
  /// a point is found when it has supports
  void SetFoundFromSupports()
  {
    SetFound(detail::HasSupports(supports_.supports_ptr, NumPoints()));
  }

  Kokkos::View<const Real*, DeviceMemorySpace> source_coords_;
  Kokkos::View<const Real*, DeviceMemorySpace> target_coords_;
  SupportResults supports_;
  int site_dim_;
  NeighbourhoodOptions options_;
};

} // namespace pcms::detail

#endif // PCMS_BACKENDS_POINT_CLOUD_MLS_EVALUATOR_HPP
