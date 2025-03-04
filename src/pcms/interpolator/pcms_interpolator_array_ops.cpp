#ifndef PCMS_INTERPOLATOR_ARRAY_OPS_HPP
#define PCMS_INTERPOLATOR_ARRAY_OPS_HPP

#include "pcms_interpolator_aliases.hpp"
#include <Omega_h_array_ops.hpp>
#include <cmath>

KOKKOS_INLINE_FUNCTION
void find_sq_root_each(ScratchVecView& array)
{
  int size = array.size();

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, size), [=](int i) {
    OMEGA_CHECK_PRINTF(array(i) < 0,
                       "[Error:] Square root of a negative number is invalid!\n"
                       "value is %12.6f\d",
                       array(i));
    array(i) = sqrt(array(i));
  });
}
