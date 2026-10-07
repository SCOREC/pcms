#include "pcms/localization/localization.h"

namespace pcms
{

LO GetOwningElementId(Omega_h::Mesh& mesh, int mesh_dim, int entity_dim,
                      LO element_id)
{
  if (element_id < 0)
    return -1;

  const int target_dim = mesh_dim; // faces for 2D, regions for 3D

  // If entity is already at target dimension, return it directly
  if (entity_dim == target_dim)
    return element_id;

  // Get the upward adjacency from entity_dim to target_dim
  auto upward_adj = mesh.ask_up(entity_dim, target_dim);
  auto a2ab_h = Omega_h::HostRead<LO>(upward_adj.a2ab);
  auto ab2b_h = Omega_h::HostRead<LO>(upward_adj.ab2b);

  const auto begin = a2ab_h[element_id];
  const auto end = a2ab_h[element_id + 1];
  if (begin >= end)
    return -1;

  // Find the smallest owning element ID
  LO owner = ab2b_h[begin];
  for (auto i = begin + 1; i < end; ++i) {
    const LO candidate = ab2b_h[i];
    if (candidate < owner)
      owner = candidate;
  }
  return owner;
}
} // namespace pcms
