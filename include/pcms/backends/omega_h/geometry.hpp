#ifndef PCMS_BACKENDS_OMEGA_H_GEOMETRY_HPP
#define PCMS_BACKENDS_OMEGA_H_GEOMETRY_HPP

#include "pcms/geometry/geometry.hpp"
#include <Omega_h_array.hpp>
#include <Omega_h_defines.hpp>
#include <string>

namespace pcms::omega_h
{

/// Options of an Omega_h geometry.
struct OmegaHGeometryOptions
{
  /// divisions and tolerance of the containing-element search
  SearchOptions search;
  /// the tag holding global entity ids
  std::string global_id_name = "global";
  /// per-vertex ownership override: a vertex the mask is 0 for is not owned,
  /// whatever the mesh says; empty means the mesh's own ownership. Layouts
  /// whose DOF holders are the mesh vertices apply it; layouts of other
  /// entities do not.
  Omega_h::Read<Omega_h::I8> owned_mask;
};

} // namespace pcms::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_GEOMETRY_HPP
