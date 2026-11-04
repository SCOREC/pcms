#ifndef PCMS_LOCALIZATION_SUPPORT_RESULTS_HPP
#define PCMS_LOCALIZATION_SUPPORT_RESULTS_HPP

#include <Omega_h_array.hpp>

namespace pcms
{

struct SupportResults
{
  Omega_h::LOs supports_ptr;
  Omega_h::LOs supports_idx;
  Omega_h::Write<Omega_h::Real> radii2;
};

} // namespace pcms

#endif // PCMS_LOCALIZATION_SUPPORT_RESULTS_HPP
