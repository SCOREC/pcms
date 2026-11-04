#ifndef PCMS_BACKENDS_OMEGA_H_OMEGA_H_HPP
#define PCMS_BACKENDS_OMEGA_H_OMEGA_H_HPP

// Public entry points of the Omega_h backend: build PCMS objects from an
// Omega_h mesh. This header (and only backend headers) may use Omega_h types.

#include "pcms/discretization/discretization.hpp"

#include <Omega_h_array.hpp>
#include <Omega_h_mesh.hpp>

#include <memory>
#include <string>

namespace pcms::backend::omega_h
{

/// Discretization on an Omega_h mesh. The mesh is borrowed and must outlive
/// the discretization and everything created from it. global_id_name names
/// the tag holding global entity ids.
std::shared_ptr<const Discretization> CreateDiscretization(
  Omega_h::Mesh& mesh, std::string global_id_name = "global");

/// As above, with per-vertex ownership given by owned_mask (1 = owned)
/// instead of the mesh's partition.
std::shared_ptr<const Discretization> CreateDiscretization(
  Omega_h::Mesh& mesh, Omega_h::Read<Omega_h::I8> owned_mask,
  std::string global_id_name = "global");

/// The Omega_h mesh of a discretization created by CreateDiscretization.
/// Throws pcms_error if the discretization is not Omega_h-based.
Omega_h::Mesh& GetMesh(const Discretization& discretization);

} // namespace pcms::backend::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_OMEGA_H_HPP
