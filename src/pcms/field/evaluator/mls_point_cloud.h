#ifndef PCMS_FIELD_EVALUATOR_MLS_POINT_CLOUD_H
#define PCMS_FIELD_EVALUATOR_MLS_POINT_CLOUD_H

#include "pcms/field/point_evaluator.h"
#include "pcms/field/field_data.h"
#include "pcms/field/evaluator/mls_options.h"
#include "pcms/field/evaluator/mls_interpolation.hpp"
#include "pcms/localization/adj_search.hpp"
#include "pcms/utility/assert.h"

#include <Omega_h_array.hpp>
#include <Omega_h_array_ops.hpp>

namespace pcms
{

// MLSPointEvaluator is a concrete PointEvaluator<Real> that evaluates a
// point-cloud-backed scalar field at a fixed set of query points using MLS.
//
// Supports only scalar fields (num_components == 1). Throws for any other
// component count.
//
// The support structure is built once at construction and reused across
// repeated Evaluate calls at zero additional localization cost.
template <typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class MLSPointEvaluator : public PointEvaluator<Real, LayoutPolicy>
{
public:
  MLSPointEvaluator(Omega_h::Reals source_coords, Omega_h::Reals target_coords,
                    SupportResults supports, int dim, MLSOptions options)
    : source_coords_(std::move(source_coords)),
      target_coords_(std::move(target_coords)),
      supports_(std::move(supports)),
      dim_(dim),
      options_(options)
  {
  }

  void Evaluate(
    const Field<Real>& field,
    Rank2View<Real, DeviceMemorySpace, LayoutPolicy> values) const override
  {
    if (values.extent(1) != 1) {
      throw pcms_error(
        "MLSPointEvaluator: only scalar (num_components==1) evaluation is "
        "supported in this phase");
    }

    auto device_data = field.GetDOFHolderData();
    const int n_sources = static_cast<int>(device_data.size());

    Omega_h::Write<Omega_h::Real> src_w(n_sources, "mls_source_values");
    Kokkos::parallel_for(
      "CopyDeviceDataToOmegaHWrite",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_sources),
      KOKKOS_LAMBDA(int i) { src_w[i] = device_data(i, 0); });
    Omega_h::Reals source_values(src_w);

    auto result = mls_interpolation(
      source_values, source_coords_, target_coords_, supports_, dim_,
      static_cast<Omega_h::LO>(options_.degree), options_.basis,
      options_.lambda, options_.tol, options_.decay_factor);

    const int n_targets = result.size();
    PCMS_ALWAYS_ASSERT(static_cast<size_t>(n_targets) == values.extent(0));
    Kokkos::parallel_for(
      "CopyMLSResultToValues",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_targets),
      KOKKOS_LAMBDA(int i) { values(i, 0) = result[i]; });
  }

private:
  Omega_h::Reals source_coords_;
  Omega_h::Reals target_coords_;
  SupportResults supports_;
  int dim_;
  MLSOptions options_;
};

} // namespace pcms

#endif // PCMS_FIELD_EVALUATOR_MLS_POINT_CLOUD_H
