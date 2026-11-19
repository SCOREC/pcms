#ifndef PCMS_FIELD_DETAIL_FIELD_LAYOUT_HPP
#define PCMS_FIELD_DETAIL_FIELD_LAYOUT_HPP

#include "pcms/field/field_layout.hpp"

namespace pcms::detail
{

// Returns the mesh entity dimension for the DOF holder at local_index, based on
// the entity-offset array: ent_offsets[d]..ent_offsets[d+1] is the range of DOF
// indices belonging to mesh entity dimension d.
inline int GetMeshEntityDim(LO local_index, const EntOffsetsArray& ent_offsets)
{
  for (int d = 0; d < ent_offsets_len - 1; ++d) {
    if (local_index >= static_cast<LO>(ent_offsets[d]) &&
        local_index < static_cast<LO>(ent_offsets[d + 1])) {
      return d;
    }
  }
  return ent_offsets_len - 1;
}

// Entity offsets over the owned DOF holders, in owned order, from the entity
// offsets of the DOF holders with ghosts and the owned indices.
inline EntOffsetsArray OwnedEntOffsets(
  const EntOffsetsArray& local_offsets,
  Kokkos::View<const LO*, HostMemorySpace> owned_indices)
{
  EntOffsetsArray offsets{};
  offsets.fill(0);
  for (size_t o = 0; o < owned_indices.size(); ++o) {
    const int d = detail::GetMeshEntityDim(owned_indices(o), local_offsets);
    for (size_t e = static_cast<size_t>(d) + 1; e < ent_offsets_len; ++e) {
      offsets[e] += 1;
    }
  }
  return offsets;
}

// The geometric classification (dimension and id) of the owned DOF holders,
// as host arrays in owned order. Converted and gathered on the device from the
// per-holder (ghosts included) device arrays; only the owned entries are
// copied to the host.
template <typename DimsView, typename IdsView>
void BuildOwnedClassificationHost(
  const DimsView& class_dims, const IdsView& class_ids,
  Kokkos::View<const LO*, DeviceMemorySpace> owned_indices,
  Kokkos::View<LO*, HostMemorySpace>& owned_dims_host,
  Kokkos::View<LO*, HostMemorySpace>& owned_ids_host)
{
  using exec = DeviceMemorySpace::execution_space;
  const auto num_owned = owned_indices.extent(0);
  Kokkos::View<LO*, DeviceMemorySpace> owned_dims("owned_class_dims",
                                                  num_owned);
  Kokkos::View<LO*, DeviceMemorySpace> owned_ids("owned_class_ids", num_owned);
  Kokkos::parallel_for(
    "pcms::GatherOwnedClassification", Kokkos::RangePolicy<exec>(0, num_owned),
    KOKKOS_LAMBDA(LO o) {
      const LO i = owned_indices(o);
      owned_dims(o) =
        static_cast<LO>(static_cast<unsigned char>(class_dims[i]));
      owned_ids(o) = static_cast<LO>(class_ids[i]);
    });
  owned_dims_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), owned_dims);
  owned_ids_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), owned_ids);
}

// Compact "owned" (rank-exclusive) views derived from a layout's arrays with
// ghosts.
struct OwnedLayoutData
{
  LO num_owned = 0;
  Kokkos::View<LO*, HostMemorySpace> owned_indices_host;
  Kokkos::View<GO*, HostMemorySpace> owned_gids_host;
  Kokkos::View<Real**, DeviceMemorySpace> owned_coords_2d;
  Kokkos::View<LO*, DeviceMemorySpace> owned_indices;
};

// Built on the device from the device-resident ownership mask, GIDs and
// coordinates; only the owned indices and owned GIDs are copied to the host,
// for their host-side users.
template <typename CoordsView>
OwnedLayoutData BuildOwnedLayoutData(
  Kokkos::View<const bool*, DeviceMemorySpace> owned,
  Kokkos::View<const GO*, DeviceMemorySpace> gids, const CoordsView& coords,
  int dim)
{
  using exec = DeviceMemorySpace::execution_space;
  const LO n = static_cast<LO>(owned.extent(0));
  OwnedLayoutData out;

  Kokkos::parallel_reduce(
    "pcms::CountOwned", Kokkos::RangePolicy<exec>(0, n),
    KOKKOS_LAMBDA(LO i, LO & count) { count += owned(i) ? 1 : 0; },
    out.num_owned);

  // Owned holders in increasing order: an exclusive scan over the mask gives
  // each owned holder its owned index.
  out.owned_indices =
    Kokkos::View<LO*, DeviceMemorySpace>("owned_indices", out.num_owned);
  const auto owned_indices = out.owned_indices;
  Kokkos::parallel_scan(
    "pcms::OwnedIndices", Kokkos::RangePolicy<exec>(0, n),
    KOKKOS_LAMBDA(LO i, LO & next, bool final) {
      if (owned(i)) {
        if (final) {
          owned_indices(next) = i;
        }
        ++next;
      }
    });

  out.owned_coords_2d =
    Kokkos::View<Real**, DeviceMemorySpace>("owned_coords", out.num_owned, dim);
  const auto owned_coords = out.owned_coords_2d;
  Kokkos::View<GO*, DeviceMemorySpace> owned_gids("owned_gids", out.num_owned);
  Kokkos::parallel_for(
    "pcms::GatherOwnedGidsAndCoords",
    Kokkos::RangePolicy<exec>(0, out.num_owned), KOKKOS_LAMBDA(LO o) {
      const LO i = owned_indices(o);
      owned_gids(o) = gids(i);
      for (int d = 0; d < dim; ++d) {
        owned_coords(o, d) = coords(i, d);
      }
    });

  out.owned_indices_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), out.owned_indices);
  out.owned_gids_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), owned_gids);
  return out;
}

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_FIELD_LAYOUT_HPP
