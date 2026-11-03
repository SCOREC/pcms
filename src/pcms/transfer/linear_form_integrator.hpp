#ifndef PCMS_TRANSFER_LINEAR_FORM_INTEGRATOR_H
#define PCMS_TRANSFER_LINEAR_FORM_INTEGRATOR_H

#include "pcms/field/coordinate_view.hpp"
#include "pcms/utility/arrays.h"
#include <petscvec.h>

namespace pcms
{

class LinearFormIntegrator
{
public:
  // Returns the ordered integration points. The index ordering must match the
  // first dimension of sampled_values passed to Assemble.
  virtual CoordinateView<DeviceMemorySpace> GetIntegrationPoints()
    const noexcept = 0;

  // Returns the internally-owned vector. The integrator owns the vector;
  // callers receive a non-owning handle. Valid after construction.
  virtual Vec GetVector() const noexcept = 0;

  // Zeros the owned vector, scatter-adds weighted contributions from the
  // sampled source-field values, and finalizes assembly.
  // sampled_values shape: [num_integration_points][num_components]
  virtual void Assemble(
    Rank2View<const Real, DeviceMemorySpace> sampled_values) = 0;

  virtual ~LinearFormIntegrator() noexcept = default;
};

} // namespace pcms

#endif // PCMS_TRANSFER_LINEAR_FORM_INTEGRATOR_H
