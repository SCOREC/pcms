#ifndef PCMS_TRANSFER_INTERPOLATION_HELPERS_H
#define PCMS_TRANSFER_INTERPOLATION_HELPERS_H

#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include <Omega_h_mesh.hpp>
#include <Omega_h_array.hpp>

namespace pcms
{

// FIXME, evaluate if this is needed. If so, move to utility library
void copyHostScalarArrayView2HostWrite(
  pcms::Rank1View<double, pcms::HostMemorySpace> source,
  Omega_h::HostWrite<Omega_h::Real>& target);

void copyHostWrite2ScalarArrayView(
  const Omega_h::HostWrite<Omega_h::Real>& source,
  pcms::Rank1View<double, pcms::HostMemorySpace> target);

} // namespace pcms
#endif // PCMS_TRANSFER_INTERPOLATION_HELPERS_H
