#ifndef PCMS_GEOMETRY_SUPPORTS_HPP
#define PCMS_GEOMETRY_SUPPORTS_HPP

#include <Omega_h_array.hpp>

namespace pcms
{

/// @brief Supports of query points among a geometry's sites, the result of
///        Geometry::Neighbourhood, in CSR form: the sites of query point i
///        are supports_idx[k] for supports_ptr[i] <= k < supports_ptr[i + 1],
///        found within the squared radius radii2[i].
struct SupportResults
{
  Omega_h::LOs supports_ptr;
  Omega_h::LOs supports_idx;
  Omega_h::Write<Omega_h::Real> radii2;
};

} // namespace pcms

#endif // PCMS_GEOMETRY_SUPPORTS_HPP
