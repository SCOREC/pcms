#ifndef PCMS_FIELD_LAYOUT_H
#define PCMS_FIELD_LAYOUT_H
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "pcms/discretization/discretization.h"
#include "pcms/utility/types.h"
#include "pcms/utility/arrays.h"
#include "coordinate_system.h"

namespace pcms
{

constexpr int ent_offsets_len = 5;
using EntOffsetsArray = std::array<size_t, ent_offsets_len>;

using ReversePartitionMap = std::map<pcms::LO, std::vector<pcms::LO>>;

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

class FieldLayout
{
public:
  // Optional identifier for this layout. Set once at construction (via a From*
  // factory) and empty by default. Consumers that need to identify a layout by
  // a stable name use it; it carries no meaning within the field layer itself.
  [[nodiscard]] const std::string& GetName() const noexcept { return name_; }
  void SetName(std::string name) { name_ = std::move(name); }

  virtual std::shared_ptr<const Discretization> GetDiscretization()
    const noexcept = 0;

  // number of components
  int virtual GetNumComponents() const = 0;

  // number of local DOF holders resident on this rank (owned + ghost)
  LO virtual GetNumLocalDofHolder() const = 0;

  // number of owned DOF holders (the subset this rank exclusively owns). For
  // non-distributed layouts this equals the local count.
  LO virtual GetNumOwnedDofHolder() const { return GetNumLocalDofHolder(); }

  GO virtual GetNumGlobalDofHolder() const = 0;

  // size of buffer needed to hold all local coefficients
  // # components * NumLocalDofHolder
  LO LocalSize() const { return GetNumComponents() * GetNumLocalDofHolder(); }

  // size of buffer needed to hold owned coefficients
  // # components * NumOwnedDofHolder
  LO OwnedSize() const { return GetNumComponents() * GetNumOwnedDofHolder(); }

  GO GlobalSize() const
  {
    return GetNumComponents() * GetNumGlobalDofHolder();
  };

  virtual Rank1View<const bool, HostMemorySpace> GetOwnedHost() const = 0;
  virtual GlobalIDView<HostMemorySpace> GetGidsHost() const = 0;

  // Device-resident global IDs for all local DOF holders.
  virtual GlobalIDView<DeviceMemorySpace> GetGids() const = 0;

  // Maps each local DOF holder to its contiguous active index, ordered by GID
  // within each entity block. Components are not included in the permutation.
  // For local GIDs [102, 7, 41, 19], the permutation is [3, 0, 2, 1]. Thus
  // holder i, component c is read from values(permutation(i), c), or from
  // flat_values[permutation(i) * num_components + c].
  Kokkos::View<const LO*, HostMemorySpace> GetGlobalToLocalPermutationHost()
    const;
  Kokkos::View<const LO*, DeviceMemorySpace> GetGlobalToLocalPermutation()
    const;

  // returns true if the field layout is distributed (holds ghost DOF holders in
  // addition to the owned ones); owned and local counts then differ
  [[nodiscard]] virtual bool IsDistributed() const = 0;

  virtual EntOffsetsArray GetEntOffsets() const = 0;

  // Entity offsets over the OWNED holders (owned-index order), derived from the
  // local entity offsets and the owned->local permutation.
  EntOffsetsArray GetOwnedEntOffsets() const
  {
    const auto owned_to_local = GetOwnedToLocalHost();
    if (owned_to_local.size() == 0) {
      return GetEntOffsets(); // non-distributed: owned == local
    }
    EntOffsetsArray offsets{};
    offsets.fill(0);
    for (size_t o = 0; o < owned_to_local.size(); ++o) {
      const int d = GetMeshEntityDim(owned_to_local(o), GetEntOffsets());
      for (size_t e = static_cast<size_t>(d) + 1; e < ent_offsets_len; ++e) {
        offsets[e] += 1;
      }
    }
    return offsets;
  }

  virtual CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const = 0;

  // Owned (rank-exclusive) global IDs, compact and owned-indexed. Defaults to
  // the local array (owned == local for non-distributed layouts).
  virtual GlobalIDView<HostMemorySpace> GetOwnedGidsHost() const
  {
    return GetGidsHost();
  }

  // Owned (rank-exclusive) global IDs, compact and owned-indexed, on device.
  // Defaults to the local device array (owned == local for non-distributed
  // layouts).
  virtual GlobalIDView<DeviceMemorySpace> GetOwnedGids() const
  {
    return GetGids();
  }

  // Owned DOF-holder coordinates, compact and owned-indexed. Defaults to the
  // local coordinates.
  virtual CoordinateView<DeviceMemorySpace> GetOwnedDOFHolderCoordinates() const
  {
    return GetDOFHolderCoordinates();
  }

  // Maps owned index (0..GetNumOwnedDofHolder()-1) to its local index. An empty
  // view means the identity map (owned == local).
  virtual Kokkos::View<const LO*, HostMemorySpace> GetOwnedToLocalHost() const
  {
    return {};
  }

  // Device analogue of GetOwnedToLocalHost(): maps owned index to its local
  // index. An empty view means the identity map (owned == local).
  virtual Kokkos::View<const LO*, DeviceMemorySpace> GetOwnedToLocal() const
  {
    return {};
  }

  virtual int GetDimension() const = 0;

  // Entity dimension of the DOF holders (0 = vertex, ..., mesh dim = element),
  // used for ghost synchronization. -1 means the layout cannot synchronize
  // ghost values via a single mesh entity dimension (e.g. non-mesh or
  // multi-dimension layouts).
  virtual int GetDOFHolderEntityDim() const { return -1; }

  virtual Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const = 0;

  virtual Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationIdsHost() const = 0;

  virtual ~FieldLayout() noexcept = default;

protected:
  // Builds the global-to-local permutation from GetGidsHost()/GetEntOffsets().
  void BuildGlobalToLocalPermutation();

private:
  std::string name_;
  Kokkos::View<LO*, HostMemorySpace> global_to_local_host_;
  Kokkos::View<LO*, DeviceMemorySpace> global_to_local_;
};

// Compact "owned" (rank-exclusive) views derived from a layout's local arrays.
struct OwnedLayoutData
{
  LO num_owned = 0;
  Kokkos::View<LO*, HostMemorySpace> owned_to_local_host;
  Kokkos::View<GO*, HostMemorySpace> owned_gids_host;
  Kokkos::View<Real**, DeviceMemorySpace> owned_coords_2d;
  Kokkos::View<LO*, DeviceMemorySpace> owned_to_local;
  Kokkos::View<GO*, DeviceMemorySpace> owned_gids;
};

template <typename CoordsView>
OwnedLayoutData BuildOwnedLayoutData(
  Kokkos::View<const bool*, HostMemorySpace> owned,
  GlobalIDView<HostMemorySpace> gids, const CoordsView& coords, int dim)
{
  const LO n_local = static_cast<LO>(owned.size());
  OwnedLayoutData out;

  out.num_owned = 0;
  for (LO i = 0; i < n_local; ++i) {
    if (owned(i)) {
      ++out.num_owned;
    }
  }

  out.owned_to_local_host =
    Kokkos::View<LO*, HostMemorySpace>("owned_to_local", out.num_owned);
  out.owned_gids_host =
    Kokkos::View<GO*, HostMemorySpace>("owned_gids", out.num_owned);

  LO o = 0;
  for (LO i = 0; i < n_local; ++i) {
    if (owned(i)) {
      out.owned_to_local_host(o) = i;
      out.owned_gids_host(o) = static_cast<GO>(gids(i));
      ++o;
    }
  }

  // Gather owned coordinates through a same-layout host mirror of the device
  // destination, so the host/device layout mismatch doesn't break deep_copy.
  out.owned_coords_2d =
    Kokkos::View<Real**, DeviceMemorySpace>("owned_coords", out.num_owned, dim);
  auto owned_coords_mirror = Kokkos::create_mirror_view(out.owned_coords_2d);
  auto coords_host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), coords);
  for (LO j = 0; j < out.num_owned; ++j) {
    const LO i = out.owned_to_local_host(j);
    for (int d = 0; d < dim; ++d) {
      owned_coords_mirror(j, d) = coords_host(i, d);
    }
  }
  Kokkos::deep_copy(out.owned_coords_2d, owned_coords_mirror);

  // Device-resident copies of the owned->local map and owned GIDs, for the
  // device-side owned accessors.
  out.owned_to_local = Kokkos::create_mirror_view_and_copy(
    DeviceMemorySpace(), out.owned_to_local_host);
  out.owned_gids = Kokkos::create_mirror_view_and_copy(DeviceMemorySpace(),
                                                       out.owned_gids_host);

  return out;
}

template <typename T>
Rank2View<const T, HostMemorySpace> GatherOwnedHostData(
  const FieldLayout& layout,
  const Kokkos::View<T*, DeviceMemorySpace>& device_data,
  Kokkos::View<T*, HostMemorySpace>& host_data,
  Kokkos::View<T*, HostMemorySpace>& owned_host_data)
{
  const auto owned_to_local = layout.GetOwnedToLocalHost();
  const LO num_comp = layout.GetNumComponents();
  if (owned_to_local.size() == 0) {
    // Non-distributed layout: owned == local.
    Kokkos::deep_copy(host_data, device_data);
    return Rank2View<const T, HostMemorySpace>(
      host_data.data(), layout.GetNumLocalDofHolder(), num_comp);
  }
  const LO num_owned = layout.GetNumOwnedDofHolder();
  Kokkos::deep_copy(host_data, device_data);
  owned_host_data = Kokkos::View<T*, HostMemorySpace>(
    "owned_host_data",
    static_cast<size_t>(num_owned) * static_cast<size_t>(num_comp));
  for (LO o = 0; o < num_owned; ++o) {
    const LO local = owned_to_local(o);
    for (LO c = 0; c < num_comp; ++c) {
      owned_host_data(static_cast<size_t>(o) * static_cast<size_t>(num_comp) +
                      static_cast<size_t>(c)) =
        host_data(static_cast<size_t>(local) * static_cast<size_t>(num_comp) +
                  static_cast<size_t>(c));
    }
  }
  return Rank2View<const T, HostMemorySpace>(owned_host_data.data(), num_owned,
                                             num_comp);
}

template <typename T>
Rank2View<const T, DeviceMemorySpace> GatherOwnedDeviceData(
  const FieldLayout& layout,
  const Kokkos::View<T*, DeviceMemorySpace>& device_data,
  Kokkos::View<T*, DeviceMemorySpace>& owned_device_data)
{
  const auto owned_to_local = layout.GetOwnedToLocal();
  const LO num_comp = layout.GetNumComponents();
  if (owned_to_local.size() == 0) {
    // Non-distributed layout: owned == local.
    return Rank2View<const T, DeviceMemorySpace>(
      device_data.data(), layout.GetNumLocalDofHolder(), num_comp);
  }
  const LO num_owned = layout.GetNumOwnedDofHolder();
  owned_device_data = Kokkos::View<T*, DeviceMemorySpace>(
    "owned_device_data",
    static_cast<size_t>(num_owned) * static_cast<size_t>(num_comp));
  Kokkos::parallel_for(
    "gather_owned_device_data",
    static_cast<size_t>(num_owned) * static_cast<size_t>(num_comp),
    KOKKOS_LAMBDA(LO i) {
      const LO o = i / num_comp;
      const LO c = i % num_comp;
      const LO local = owned_to_local(o);
      owned_device_data(i) = device_data(local * num_comp + c);
    });
  return Rank2View<const T, DeviceMemorySpace>(owned_device_data.data(),
                                               num_owned, num_comp);
}

} // namespace pcms
#endif // PCMS_FIELD_LAYOUT_H
