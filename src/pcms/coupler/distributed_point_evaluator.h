#ifndef PCMS_DISTRIBUTED_POINT_EVALUATOR_H
#define PCMS_DISTRIBUTED_POINT_EVALUATOR_H
#include <limits>
#include <type_traits>
#include <Kokkos_Core.hpp>
#include "pcms/coupler/distributed_evaluation_channel.h"
#include "pcms/field/evaluation_request.h"
#include "pcms/field/field.h"
#include "pcms/field/field_evaluator_factory.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/profile.h"

namespace pcms
{

namespace detail
{

// Copies a (possibly device-resident) Rank2View into a freshly allocated host
// Kokkos::View, regardless of the view's own mdspan layout.
//
// The result always uses HostMemorySpace's default layout. A device (CUDA)
// source view is layout_left, and Kokkos cannot construct a layout_right host
// view directly from it ("Incompatible View copy construction"), so the
// layouts are bridged with DeepCopyMismatchLayouts rather than by returning
// the mirror's own layout -- see the equivalent copy_coordinates_to_host in
// pcms/test/test_distributed_field.cpp.
template <typename ElementType, typename MemorySpace, typename LayoutPolicy>
Kokkos::View<std::remove_const_t<ElementType>**, HostMemorySpace> CopyToHostView(
  Rank2View<ElementType, MemorySpace, LayoutPolicy> src)
{
  const LO n = static_cast<LO>(src.extent(0));
  const int cols = static_cast<int>(src.extent(1));
  using ViewLayout =
    std::conditional_t<std::is_same_v<LayoutPolicy, Kokkos::layout_left>,
                       Kokkos::LayoutLeft, Kokkos::LayoutRight>;
  Kokkos::View<ElementType**, ViewLayout, MemorySpace, Kokkos::MemoryUnmanaged>
    unmanaged(src.data_handle(), n, cols);
  Kokkos::View<std::remove_const_t<ElementType>**, HostMemorySpace> host(
    "copy_to_host_view", n, cols);
  DeepCopyMismatchLayouts(host, unmanaged);
  return host;
}

// Uploads a host Rank2View into a freshly allocated device Kokkos::View.
template <typename ElementType>
Kokkos::View<ElementType**, DeviceMemorySpace> CopyToDeviceView(
  Rank2View<const ElementType, HostMemorySpace> src)
{
  const LO n = static_cast<LO>(src.extent(0));
  const int cols = static_cast<int>(src.extent(1));
  Kokkos::View<ElementType**, DeviceMemorySpace> dst("copy_to_device", n, cols);
  auto mirror = Kokkos::create_mirror_view(dst);
  for (LO i = 0; i < n; ++i) {
    for (int c = 0; c < cols; ++c) {
      mirror(i, c) = src(i, c);
    }
  }
  Kokkos::deep_copy(dst, mirror);
  return dst;
}

// A value used to mark "not localized" in an internal, diagnostic-only FILL
// pass, distinct from the caller's own OutOfBoundsPolicy::fill_value. Chosen
// to be exactly representable in T and to round-trip exactly through the
// OutOfBoundsPolicy::fill_value Real, since PointEvaluator implementations
// convert it back to T internally (numeric_limits<int64_t>::lowest() is not
// exactly representable as a double, so it is special-cased).
template <typename T>
constexpr T DiagnosticSentinel()
{
  if constexpr (std::is_same_v<T, int64_t>) {
    return static_cast<T>(-(int64_t(1) << 53));
  } else {
    return std::numeric_limits<T>::lowest();
  }
}

} // namespace detail

// DistributedPointEvaluator<T> is a PointEvaluator<T> decorator that adds
// cross-rank fallback for query points a distributed field's local mesh
// partition cannot localize: points the wrapped factory cannot resolve
// locally are routed, through a DistributedEvaluationChannel<T> (and thus
// through a standalone evaluation server -- see
// pcms/src/pcms/coupler/distributed_evaluation_channel.h), to whichever
// rank(s) the server's routing partition names as candidate owners. Points no
// candidate rank can resolve are genuinely outside the global domain and are
// handled by `policy`, exactly as OutOfBoundsPolicy documents.
//
// factory must outlive this object: it is reused every resolution round to
// build a fresh, correctly-localized PointEvaluator for whatever points the
// server redirects to this rank that round (a PointEvaluator is bound to one
// fixed query point set at construction time, so it cannot be reused for a
// different point set -- see point_evaluator.h).
//
// Every rank must construct one of these and call Evaluate collectively, in
// lockstep with every other rank's DistributedEvaluationChannel and with the
// standalone server's RunServerRound calls.
template <typename T>
class DistributedPointEvaluator : public PointEvaluator<T>
{
public:
  DistributedPointEvaluator(const FieldEvaluatorFactory<T>& factory,
                            CoordinateView<DeviceMemorySpace> coords,
                            DistributedEvaluationChannel<T>& channel,
                            OutOfBoundsPolicy policy)
    : factory_(factory),
      channel_(channel),
      policy_(policy),
      num_components_(factory.GetLayout().GetNumComponents()),
      coordinate_system_(coords.GetCoordinateSystem()),
      host_points_(detail::CopyToHostView(coords.GetValues())),
      local_evaluator_(factory.CreatePointEvaluator(
        EvaluationRequest::FromCoordinates(
          coords,
          OutOfBoundsPolicy{OutOfBoundsMode::FILL,
                            static_cast<Real>(detail::DiagnosticSentinel<T>())})))
  {
    if (policy_.mode == OutOfBoundsMode::NEAREST_BOUNDARY) {
      throw pcms_error(
        "DistributedPointEvaluator: NearestBoundary is not supported");
    }
  }

  void Evaluate(const Field<T>& field,
               Rank2View<T, DeviceMemorySpace> values) const override
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(values.extent(1) == static_cast<size_t>(num_components_));
    const LO n = static_cast<LO>(host_points_.extent(0));
    constexpr T sentinel = detail::DiagnosticSentinel<T>();

    // Local pass: fill with a sentinel so "the local evaluator rejected this
    // point" and "the real value happens to equal the policy's fill_value"
    // can never be confused.
    Kokkos::View<T**, DeviceMemorySpace> device_values("dpe_values", n,
                                                       num_components_);
    Kokkos::deep_copy(device_values, sentinel);
    local_evaluator_->Evaluate(field, MakeRank2View(device_values));

    auto host_values = detail::CopyToHostView(MakeRank2View(device_values));
    Kokkos::View<bool*, HostMemorySpace> resolved("dpe_resolved", n);
    for (LO i = 0; i < n; ++i) {
      resolved(i) = host_values(i, 0) != sentinel;
    }

    // Answers queries the server redirects to this rank during the same
    // round, using this rank's own field/factory -- see
    // DistributedEvaluationChannel's class comment for why this is what lets
    // a remote rank answer without ever transmitting field data itself.
    auto local_evaluate = [this, &field, sentinel](
                            Rank2View<const Real, HostMemorySpace> pts,
                            Rank2View<T, HostMemorySpace> out_values,
                            Rank1View<bool, HostMemorySpace> out_resolved) {
      const LO m = static_cast<LO>(pts.extent(0));
      auto device_coords = detail::CopyToDeviceView<Real>(pts);
      CoordinateView<DeviceMemorySpace> query(coordinate_system_,
                                             MakeRank2View(device_coords));
      auto evaluator = factory_.CreatePointEvaluator(
        EvaluationRequest::FromCoordinates(
          query, OutOfBoundsPolicy{OutOfBoundsMode::FILL,
                                  static_cast<Real>(sentinel)}));
      Kokkos::View<T**, DeviceMemorySpace> vals("dpe_redirected_values", m,
                                                num_components_);
      Kokkos::deep_copy(vals, sentinel);
      evaluator->Evaluate(field, MakeRank2View(vals));
      auto host_vals = detail::CopyToHostView(MakeRank2View(vals));
      for (LO i = 0; i < m; ++i) {
        out_resolved(i) = host_vals(i, 0) != sentinel;
        for (int c = 0; c < num_components_; ++c) {
          out_values(i, c) = host_vals(i, c);
        }
      }
    };

    // Publish this rank's owned bounding box so the server can build its
    // routing partition from the field's real distribution. Idempotent, and a
    // no-op after the first call.
    channel_.PublishOwnedBounds(factory_.GetLayout());

    channel_.ResolveUnresolved(MakeRank2View(host_points_),
                              MakeRank2View(host_values),
                              MakeRank1View(resolved), local_evaluate);

    for (LO i = 0; i < n; ++i) {
      if (resolved(i)) {
        continue;
      }
      if (policy_.mode == OutOfBoundsMode::ERROR) {
        throw pcms_error(
          "DistributedPointEvaluator: query point is outside the global "
          "domain (no rank could localize it)");
      }
      // FILL: overwrite the internal diagnostic sentinel with the caller's
      // real fill value.
      for (int c = 0; c < num_components_; ++c) {
        host_values(i, c) = static_cast<T>(policy_.fill_value);
      }
    }

    using ValuesLayout = std::conditional_t<
      std::is_same_v<typename decltype(values)::layout_type,
                     Kokkos::layout_left>,
      Kokkos::LayoutLeft, Kokkos::LayoutRight>;
    Kokkos::View<T**, ValuesLayout, DeviceMemorySpace, Kokkos::MemoryUnmanaged>
      values_device_view(values.data_handle(), n, num_components_);
    DeepCopyMismatchLayouts(values_device_view, host_values);
  }

private:
  const FieldEvaluatorFactory<T>& factory_;
  DistributedEvaluationChannel<T>& channel_;
  OutOfBoundsPolicy policy_;
  int num_components_;
  CoordinateSystem coordinate_system_;
  Kokkos::View<Real**, HostMemorySpace> host_points_;
  std::unique_ptr<PointEvaluator<T>> local_evaluator_;
};

} // namespace pcms

#endif // PCMS_DISTRIBUTED_POINT_EVALUATOR_H
