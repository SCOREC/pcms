#ifndef PCMS_POLYNOMIAL_RECONSTRUCTION_FUNCTION_SPACE_H
#define PCMS_POLYNOMIAL_RECONSTRUCTION_FUNCTION_SPACE_H

#include "pcms/field/evaluator/mls_options.h"
#include "pcms/field/function_space.h"
#include "pcms/geometry/geometry.hpp"

#include <memory>

namespace pcms
{

/// An MLS reconstruction space on a point-cloud geometry (its sites are the
/// DOF holders); `options` are the reconstruction parameters, the support
/// search is the geometry's. Spaces on one geometry accept each other's bound
/// sources.
/// @throws pcms_error if `geometry` is null or not a point-cloud geometry
[[nodiscard]] std::shared_ptr<const FunctionSpace>
CreatePolynomialReconstructionSpace(std::shared_ptr<const Geometry> geometry,
                                    MLSOptions options = {});

} // namespace pcms

#endif // PCMS_POLYNOMIAL_RECONSTRUCTION_FUNCTION_SPACE_H
