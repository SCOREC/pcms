#ifndef PCMS_UTILITY_MESH_GEOMETRY_H
#define PCMS_UTILITY_MESH_GEOMETRY_H

#include <Omega_h_array.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>

namespace pcms
{

KOKKOS_INLINE_FUNCTION
Omega_h::Real distance_squared(const Omega_h::Real* p1, const Omega_h::Real* p2,
                               int dim)
{
  Omega_h::Real dx = p1[0] - p2[0];
  Omega_h::Real dy = p1[1] - p2[1];
  Omega_h::Real dz = (dim == 3) ? (p1[2] - p2[2]) : 0.0;
  return dx * dx + dy * dy + dz * dz;
}

Omega_h::Reals get_entity_centroids(Omega_h::Mesh& mesh,
                                    Omega_h::Int entity_dim);

// Some meshes (e.g. certain XGC meshes) are stored without any geometric
// classification (no class_id/class_dim tags). Layouts read these tags on
// every mesh dimension, so derive a default classification when it is
// missing. classify_elements + finalize_classification populate a
// geometrically meaningful class_dim on all dimensions but do not create
// class_id (that normally comes from a geometric model we don't have here),
// so fill any missing class_id with 0.
// Classification is only consulted for coupling DOF matching
// (field_exchange_planner); standalone field transfer never reads it, so a
// default id is sufficient to keep such meshes usable.
inline void EnsureClassification(Omega_h::Mesh& mesh)
{
  if (mesh.has_tag(mesh.dim(), "class_id") &&
      mesh.has_tag(mesh.dim(), "class_dim")) {
    return;
  }
  if (!mesh.has_tag(mesh.dim(), "class_dim")) {
    Omega_h::classify_elements(&mesh);
    Omega_h::finalize_classification(&mesh);
  }
  for (int d = 0; d <= mesh.dim(); ++d) {
    if (!mesh.has_tag(d, "class_id")) {
      mesh.add_tag<Omega_h::ClassId>(
        d, "class_id", 1,
        Omega_h::Read<Omega_h::ClassId>(mesh.nents(d), 0, "class_id"));
    }
  }
}

} // namespace pcms

#endif // PCMS_UTILITY_MESH_GEOMETRY_H
