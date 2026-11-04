#ifndef PCMS_LOCALIZATION_MLS_SUPPORT_HELPERS_H
#define PCMS_LOCALIZATION_MLS_SUPPORT_HELPERS_H

#include <Omega_h_array.hpp>

namespace pcms
{

struct SupportResults;

// Returns true when both the minimum and maximum support counts fall within
// the requested range.
inline bool within_number_of_support_range(unsigned min_supports_found,
                                           unsigned max_supports_found,
                                           unsigned min_req_supports,
                                           unsigned max_allowed_supports)
{
  return (min_supports_found >= min_req_supports) &&
         (max_supports_found <= max_allowed_supports);
}

// Compute the element-wise min and max of a support-count array.
void minmax(Omega_h::Read<Omega_h::LO> num_supports,
            unsigned& min_supports_found, unsigned& max_supports_found);

// Scale per-target squared radii to drive support counts toward [min, max].
void adapt_radii(unsigned min_req_supports, unsigned max_allowed_supports,
                 Omega_h::LO n_targets, Omega_h::Write<Omega_h::Real> radii2_l,
                 Omega_h::Write<Omega_h::LO> num_supports);

// Build SupportResults for a point-cloud source/target pair using a
// distance-based N² search with optional adaptive radius adjustment.
// source_coords and target_coords are flat arrays in [x0,y0,..., x1,y1,...]
// order.
SupportResults BuildPointCloudSupports(const Omega_h::Reals& source_coords,
                                       const Omega_h::Reals& target_coords,
                                       int dim, double radius,
                                       unsigned min_req_supports,
                                       bool adapt_radius,
                                       unsigned max_iterations = 100);

} // namespace pcms
#endif // PCMS_LOCALIZATION_MLS_SUPPORT_HELPERS_H
