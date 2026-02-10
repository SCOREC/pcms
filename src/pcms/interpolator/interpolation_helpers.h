//
// Created by hasanm4 on 10/10/25.
//

#ifndef PCMS_INTERPOLATION_HELPERS_H
#define PCMS_INTERPOLATION_HELPERS_H
#include <pcms/arrays.h>
#include <Omega_h_mesh.hpp>
#include <Omega_h_array.hpp>

void copyHostScalarArrayView2HostWrite(
  pcms::Rank1View<double, pcms::HostMemorySpace> source,
  Omega_h::HostWrite<Omega_h::Real>& target);

void copyHostWrite2ScalarArrayView(
  const Omega_h::HostWrite<Omega_h::Real>& source,
  pcms::Rank1View<double, pcms::HostMemorySpace> target);

Omega_h::Reals getCentroids(Omega_h::Mesh& mesh);

inline bool within_number_of_support_range(unsigned min_supports_found,
                                           unsigned max_supports_found,
                                           unsigned min_req_supports,
                                           unsigned max_allowed_supports)
{
  return (min_supports_found >= min_req_supports) &&
         (max_supports_found <= max_allowed_supports);
}

void minmax(Omega_h::Read<Omega_h::LO> num_supports,
            unsigned& min_supports_found, unsigned& max_supports_found);
void adapt_radii(unsigned min_req_supports, unsigned max_allowed_supports,
                 Omega_h::LO n_targets, Omega_h::Write<Omega_h::Real> radii2_l,
                 Omega_h::Write<Omega_h::LO> num_supports);

#endif // PCMS_INTERPOLATION_HELPERS_H
