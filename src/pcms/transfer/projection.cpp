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
  const FunctionSpace& source_space, const FunctionSpace& target_space,
  ConservativeProjectionOptions options)
{
#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
  return std::make_unique<pcms::transfer::OmegaHConservativeProjection>(
    source_space, target_space, options.mass_matrix_type);
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
  const PointEvaluatorFactory& source, const FunctionSpace& target_space,
  ControlVariateProjectionOptions options)
{
#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
  return std::make_unique<pcms::transfer::OmegaHControlVariateProjection>(
    source, target_space, options.samples_per_element, options.sampling,
    options.seed, options.mass_matrix_type);
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
