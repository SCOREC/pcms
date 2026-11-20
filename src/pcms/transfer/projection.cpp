#include "pcms/transfer/projection.hpp"
#include "pcms/configuration.hpp"
#include "pcms/transfer/galerkin_projection.hpp"
#include "pcms/utility/assert.hpp"
#if defined(PCMS_ENABLE_OMEGA_H_MESHFIELDS)
#include "pcms/transfer/omega_h_form_integrator_utils.hpp"
#include "pcms/transfer/omega_h_intersection_quadrature.hpp"
#endif

namespace pcms
{

std::unique_ptr<TransferOperator<Real>> CreateProjection(
  std::shared_ptr<const PointEvaluatorFactory> source,
  std::shared_ptr<const FunctionSpace> target_space, ProjectionOptions options)
{
  if (source == nullptr || target_space == nullptr) {
    throw pcms_error("CreateProjection: the source and the target space must "
                     "not be null");
  }
  return std::make_unique<transfer::GalerkinProjection>(
    std::move(source), std::move(target_space), std::move(options));
}

std::unique_ptr<TransferOperator<Real>> CreateConservativeProjection(
  std::shared_ptr<const FunctionSpace> source_space,
  std::shared_ptr<const FunctionSpace> target_space,
  ConservativeProjectionOptions options)
{
  if (source_space == nullptr || target_space == nullptr) {
    throw pcms_error(
      "CreateConservativeProjection: the function spaces must not be null");
  }
#if defined(PCMS_ENABLE_OMEGA_H_MESHFIELDS)
  // The exact right-hand side integrates (source basis) x (target basis) over
  // the pieces of the intersection of the meshes, which lie in single source
  // and target elements.
  const auto [source_layout, source_geometry] =
    transfer::detail::ResolveOmegaHLagrangeSpace(
      *source_space, 0, 1, "CreateConservativeProjection", "source");
  const auto [target_layout, target_geometry] =
    transfer::detail::ResolveOmegaHLagrangeSpace(
      *target_space, 0, 1, "CreateConservativeProjection", "target");
  ProjectionOptions projection;
  projection.mass_matrix_type = options.mass_matrix_type;
  projection.quadrature = omega_h::CreateIntersectionQuadrature(
    source_geometry, target_geometry,
    source_layout->GetOrder() + target_layout->GetOrder());
  return CreateProjection(std::move(source_space), std::move(target_space),
                          std::move(projection));
#else
  (void)source_space;
  (void)target_space;
  (void)options;
  throw pcms_error("CreateConservativeProjection: not available in this "
                   "build (needs PCMS_ENABLE_OMEGA_H_MESHFIELDS)");
#endif
}

std::unique_ptr<TransferOperator<Real>> CreateControlVariateProjection(
  std::shared_ptr<const PointEvaluatorFactory> source,
  std::shared_ptr<const FunctionSpace> target_space,
  ControlVariateProjectionOptions options)
{
  if (source == nullptr || target_space == nullptr) {
    throw pcms_error("CreateControlVariateProjection: the source and the "
                     "target space must not be null");
  }
  ProjectionOptions projection;
  projection.quadrature = CreateMonteCarloQuadrature(
    target_space->GetGeometry(), options.samples_per_element, options.seed);
  projection.control_variate = true;
  projection.mass_matrix_type = options.mass_matrix_type;
  projection.to_source = std::move(options.to_source);
  return CreateProjection(std::move(source), std::move(target_space),
                          std::move(projection));
}

} // namespace pcms
