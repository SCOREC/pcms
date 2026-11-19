#include "pcms/transfer/projection.hpp"
#include "pcms/configuration.hpp"
#include "pcms/utility/assert.hpp"
#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
#include "pcms/transfer/omega_h_conservative_projection.hpp"
#include "pcms/transfer/omega_h_control_variate_projection.hpp"
#endif

namespace pcms
{

std::unique_ptr<TransferOperator<Real>> CreateConservativeProjection(
  std::shared_ptr<const FunctionSpace> source_space,
  std::shared_ptr<const FunctionSpace> target_space,
  ConservativeProjectionOptions options)
{
  if (source_space == nullptr || target_space == nullptr) {
    throw pcms_error(
      "CreateConservativeProjection: the function spaces must not be null");
  }
#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
  return std::make_unique<transfer::OmegaHConservativeProjection>(
    std::move(source_space), std::move(target_space), options.mass_matrix_type);
#else
  (void)source_space;
  (void)target_space;
  (void)options;
  throw pcms_error(
    "CreateConservativeProjection: not available in this "
    "build (needs PCMS_ENABLE_PETSC and PCMS_ENABLE_MESHFIELDS)");
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
#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
  return std::make_unique<transfer::OmegaHControlVariateProjection>(
    std::move(source), std::move(target_space), options.samples_per_element,
    options.seed, options.mass_matrix_type, std::move(options.to_source));
#else
  (void)source;
  (void)target_space;
  (void)options;
  throw pcms_error(
    "CreateControlVariateProjection: not available in this "
    "build (needs PCMS_ENABLE_PETSC and PCMS_ENABLE_MESHFIELDS)");
#endif
}

} // namespace pcms
