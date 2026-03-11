#ifndef PCMS_UTILITY_MESH_GEOMETRY_H
#define PCMS_UTILITY_MESH_GEOMETRY_H

#include <Omega_h_array.hpp>
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

inline Omega_h::Reals get_entity_centroids(Omega_h::Mesh& mesh,
                                           Omega_h::Int entity_dim)
{
  OMEGA_H_CHECK_PRINTF(entity_dim >= Omega_h::VERT && entity_dim <= mesh.dim(),
                       "Unsupported entity_dim=%d for mesh dim=%d\n",
                       entity_dim, mesh.dim());

  if (entity_dim == Omega_h::VERT) {
    return mesh.coords();
  }

  const auto dim = mesh.dim();
  const auto nents = mesh.nents(entity_dim);
  const auto ent2verts = mesh.ask_down(entity_dim, Omega_h::VERT).ab2b;
  const auto coords = mesh.coords();
  Omega_h::Write<Omega_h::Real> centroids(nents * dim, 0.0, "entity centroids");

  if (dim == 2 && entity_dim == Omega_h::EDGE) {
    Omega_h::parallel_for(
      "entity_centroids_edge_2d", nents, OMEGA_H_LAMBDA(const Omega_h::LO ent) {
        const auto verts = Omega_h::gather_verts<2>(ent2verts, ent);
        const auto vert_coords = Omega_h::gather_vectors<2, 2>(coords, verts);
        const auto centroid = Omega_h::average(vert_coords);
        centroids[2 * ent + 0] = centroid[0];
        centroids[2 * ent + 1] = centroid[1];
      });
  } else if (dim == 2 && entity_dim == Omega_h::FACE) {
    Omega_h::parallel_for(
      "entity_centroids_face_2d", nents, OMEGA_H_LAMBDA(const Omega_h::LO ent) {
        const auto verts = Omega_h::gather_verts<3>(ent2verts, ent);
        const auto vert_coords = Omega_h::gather_vectors<3, 2>(coords, verts);
        const auto centroid = Omega_h::average(vert_coords);
        centroids[2 * ent + 0] = centroid[0];
        centroids[2 * ent + 1] = centroid[1];
      });
  } else if (dim == 3 && entity_dim == Omega_h::EDGE) {
    Omega_h::parallel_for(
      "entity_centroids_edge_3d", nents, OMEGA_H_LAMBDA(const Omega_h::LO ent) {
        const auto verts = Omega_h::gather_verts<2>(ent2verts, ent);
        const auto vert_coords = Omega_h::gather_vectors<2, 3>(coords, verts);
        const auto centroid = Omega_h::average(vert_coords);
        centroids[3 * ent + 0] = centroid[0];
        centroids[3 * ent + 1] = centroid[1];
        centroids[3 * ent + 2] = centroid[2];
      });
  } else if (dim == 3 && entity_dim == Omega_h::FACE) {
    Omega_h::parallel_for(
      "entity_centroids_face_3d", nents, OMEGA_H_LAMBDA(const Omega_h::LO ent) {
        const auto verts = Omega_h::gather_verts<3>(ent2verts, ent);
        const auto vert_coords = Omega_h::gather_vectors<3, 3>(coords, verts);
        const auto centroid = Omega_h::average(vert_coords);
        centroids[3 * ent + 0] = centroid[0];
        centroids[3 * ent + 1] = centroid[1];
        centroids[3 * ent + 2] = centroid[2];
      });
  } else if (dim == 3 && entity_dim == Omega_h::REGION) {
    Omega_h::parallel_for(
      "entity_centroids_region_3d", nents,
      OMEGA_H_LAMBDA(const Omega_h::LO ent) {
        const auto verts = Omega_h::gather_verts<4>(ent2verts, ent);
        const auto vert_coords = Omega_h::gather_vectors<4, 3>(coords, verts);
        const auto centroid = Omega_h::average(vert_coords);
        centroids[3 * ent + 0] = centroid[0];
        centroids[3 * ent + 1] = centroid[1];
        centroids[3 * ent + 2] = centroid[2];
      });
  } else {
    OMEGA_H_CHECK_PRINTF(
      false, "Centroid computation not implemented for entity_dim=%d dim=%d\n",
      entity_dim, dim);
  }

  return Omega_h::Reals(centroids);
}

} // namespace pcms

#endif // PCMS_UTILITY_MESH_GEOMETRY_H
