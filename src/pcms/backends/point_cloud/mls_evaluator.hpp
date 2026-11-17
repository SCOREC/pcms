#ifndef PCMS_FIELD_EVALUATOR_MLS_POINT_CLOUD_H
#define PCMS_FIELD_EVALUATOR_MLS_POINT_CLOUD_H

#include "pcms/field/field.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/localization/mls_options.hpp"
#include "pcms/field/evaluator/mls_interpolation.hpp"
#include "pcms/backends/omega_h/adj_search.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/omega_h_array_utils.hpp"

#include <Omega_h_array.hpp>
#include <Omega_h_array_ops.hpp>

namespace pcms
{

namespace detail
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

} // namespace detail

/// @brief MLS bound source: the supports of the query points and the
///        coordinates the solve uses. A query point without supports (e.g.
///        one the adjacency search has no start for) is missing.
class MLSBinding final : public BoundSource
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
             const BindRequest& request, Omega_h::Reals source_coords,
             SupportResults supports, int site_dim,
             NeighbourhoodOptions options)
    : BoundSource(std::move(geometry), request),
      source_coords_(std::move(source_coords)),
      target_coords_(
        flatten_to_omega_h_reals(LocatedPointValues(), "mls_target_coords")),
      supports_(std::move(supports)),
      site_dim_(site_dim),
      options_(options)
  {
    SetFoundFromSupports();
  }

  /// @return flattened coordinates of the source points
  [[nodiscard]] const Omega_h::Reals& SourceCoordinates() const noexcept
  {
    return source_coords_;
  }

  /// @return flattened coordinates of the query points in the geometry's
  ///         coordinate system (LocatedPointValues())
  [[nodiscard]] const Omega_h::Reals& TargetCoordinates() const noexcept
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
    const auto& geometry = *GetGeometry();
    supports_ =
      geometry.Neighbourhood(CoordinateView<DeviceMemorySpace>(
                               geometry.GetCoordinateSystem(),
                               MakeConstRank2View(LocatedPointValues())),
                             site_dim_, options_);
    target_coords_ =
      flatten_to_omega_h_reals(LocatedPointValues(), "mls_target_coords");
    SetFoundFromSupports();
  }

private:
  /// a point is found when it has supports
  void SetFoundFromSupports()
  {
    SetFound(detail::HasSupports(supports_.supports_ptr, NumPoints()));
  }

  Omega_h::Reals source_coords_;
  Omega_h::Reals target_coords_;
  SupportResults supports_;
  int site_dim_;
  NeighbourhoodOptions options_;
};

/// @brief Evaluates a scalar point-cloud field at the bound points by MLS.
template <typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class MLSPointEvaluator : public PointEvaluator<Real, LayoutPolicy>
{
public:
  /// @param binding MLS supports of the query points
  /// @param dim spatial dimension
  /// @param options MLS options
  /// @param layout layout of the evaluated fields
  MLSPointEvaluator(std::shared_ptr<const MLSBinding> binding, int dim,
                    MLSOptions options,
                    std::shared_ptr<const FieldLayout> layout)
    : PointEvaluator<Real, LayoutPolicy>(binding, DifferentialOperator::Value,
                                         std::move(layout)),
      binding_(std::move(binding)),
      dim_(dim),
      options_(options)
  {
  }

  // Public because CUDA forbids extended lambdas in non-public member
  // functions. User call Evaluate.
  // TODO: add passkey or other mechanism to protect this
  void EvaluateImpl(
    const Field<Real>& field,
    ValueView<Real, DeviceMemorySpace, LayoutPolicy> out) const override
  {
    const auto values = out.GetValues();
    if (values.extent(1) != 1) {
      throw pcms_error(
        "MLSPointEvaluator: only scalar (num_components==1) evaluation is "
        "supported in this phase");
    }

    const auto device_data = field.GetDOFHolderData().GetValues();
    const int n_sources = static_cast<int>(device_data.size());

    Omega_h::Write<Omega_h::Real> src_w(n_sources, "mls_source_values");
    Kokkos::parallel_for(
      "CopyDeviceDataToOmegaHWrite",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_sources),
      KOKKOS_LAMBDA(int i) { src_w[i] = device_data(i, 0); });
    Omega_h::Reals source_values(src_w);

    auto result = mls_interpolation(
      source_values, binding_->SourceCoordinates(),
      binding_->TargetCoordinates(), binding_->Supports(), dim_,
      static_cast<Omega_h::LO>(options_.degree), options_.basis,
      options_.lambda, options_.tol, options_.decay_factor);

    const int n_targets = result.size();
    Kokkos::parallel_for(
      "CopyMLSResultToValues",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_targets),
      KOKKOS_LAMBDA(int i) { values(i, 0) = result[i]; });
  }

private:
  std::shared_ptr<const MLSBinding> binding_;
  int dim_;
  MLSOptions options_;
};

} // namespace pcms

#endif // PCMS_FIELD_EVALUATOR_MLS_POINT_CLOUD_H
