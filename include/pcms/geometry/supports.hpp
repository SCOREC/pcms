#ifndef PCMS_GEOMETRY_SUPPORTS_HPP
#define PCMS_GEOMETRY_SUPPORTS_HPP

#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"

#include <Kokkos_Core.hpp>

namespace pcms
{

/// @brief Supports of query points among a geometry's sites, the result of
///        Geometry::Neighbourhood, in CSR form: the sites of query point i
///        are supports_idx[k] for supports_ptr[i] <= k < supports_ptr[i + 1],
///        found within the squared radius radii2[i].
struct SupportResults
{
  Kokkos::View<const LO*, DeviceMemorySpace> supports_ptr;
  Kokkos::View<const LO*, DeviceMemorySpace> supports_idx;
  Kokkos::View<Real*, DeviceMemorySpace> radii2;
};

} // namespace pcms

#endif // PCMS_GEOMETRY_SUPPORTS_HPP
