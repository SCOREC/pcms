#ifndef PCMS_FIELD_LAYOUT_H
#define PCMS_FIELD_LAYOUT_H
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "pcms/discretization/discretization.h"
#include "pcms/utility/types.h"
#include "pcms/utility/arrays.h"
#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_view.hpp"

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

// FieldLayout describes how a field's DOF holders are laid out on this rank.
//
// The accessors describe the DOF holders this rank owns (rank-exclusive),
// compact and in owned order. A distributed layout also holds ghost copies of
// DOF holders other ranks own, for internal use (field storage, ghost
// synchronization, element assembly): together they are the DOF holders with
// ghosts, in mesh order. Only their ownership and the owned indices are exposed
// here. For a layout without ghosts the two coincide (!HasGhosts()).
// Entity offsets over the owned DOF holders, in owned order, from the entity
// offsets of the DOF holders with ghosts and the owned indices.
inline EntOffsetsArray OwnedEntOffsets(
  const EntOffsetsArray& local_offsets,
  Kokkos::View<const LO*, HostMemorySpace> owned_indices)
{
  EntOffsetsArray offsets{};
  offsets.fill(0);
  for (size_t o = 0; o < owned_indices.size(); ++o) {
    const int d = GetMeshEntityDim(owned_indices(o), local_offsets);
    for (size_t e = static_cast<size_t>(d) + 1; e < ent_offsets_len; ++e) {
      offsets[e] += 1;
    }
  }
  return offsets;
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

  // number of DOF holders this rank owns
  LO virtual GetNumDofHolder() const = 0;

  GO virtual GetNumGlobalDofHolder() const = 0;

  // size of buffer needed to hold the owned coefficients
  // # components * NumDofHolder
  LO Size() const { return GetNumComponents() * GetNumDofHolder(); }

  GO GlobalSize() const
  {
    return GetNumComponents() * GetNumGlobalDofHolder();
  };

  // True when the layout holds ghost DOF holders, i.e. owned order differs
  // from the order with ghosts. Layouts that hold ghosts override it.
  [[nodiscard]] virtual bool HasGhosts() const { return false; }

  // returns true if the field layout is distributed: on several ranks it may
  // hold ghost DOF holders in addition to the owned ones (see HasGhosts)
  [[nodiscard]] virtual bool IsDistributed() const = 0;

  // Per DOF holder, ghosts included: whether this rank owns it.
  virtual Rank1View<const bool, HostMemorySpace> GetOwnedHost() const = 0;

  // Maps owned index (0..GetNumDofHolder()-1) to its index among the DOF
  // holders with ghosts. An empty view means the identity map (no ghosts). A
  // distributed layout on a rank that owns nothing also returns an empty view,
  // so callers that need to know whether the layout has no ghosts must ask
  // `!HasGhosts()`, not test for empty.
  virtual Kokkos::View<const LO*, HostMemorySpace> GetOwnedIndicesHost() const
  {
    return {};
  }

  // Device analogue of GetOwnedIndicesHost().
  virtual Kokkos::View<const LO*, DeviceMemorySpace> GetOwnedIndices() const
  {
    return {};
  }

  virtual GlobalIDView<HostMemorySpace> GetGidsHost() const = 0;

  // Maps each local DOF holder to its contiguous active index, ordered by GID
  // within each entity block. Components are not included in the permutation.
  // For local GIDs [102, 7, 41, 19], the permutation is [3, 0, 2, 1]. Thus
  // holder i, component c is read from values(permutation(i), c), or from
  // flat_values[permutation(i) * num_components + c].
  Kokkos::View<const LO*, HostMemorySpace> GetGlobalToLocalPermutationHost()
    const;
  Kokkos::View<const LO*, DeviceMemorySpace> GetGlobalToLocalPermutation()
    const;

  virtual EntOffsetsArray GetEntOffsets() const = 0;

  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const
  {
    if (system_ == nullptr) {
      throw pcms_error("FieldLayout::GetCoordinateSystem: layout has no "
                       "coordinate system set");
    }
    return system_;
  }

  virtual CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const = 0;

  virtual int GetDimension() const = 0;

  // Entity dimension of the DOF holders (0 = vertex, ..., mesh dim = element)
  // when the layout has exactly one DOF holder per local mesh entity of that
  // dimension, ghosts included; used for ghost synchronization and to pick
  // mesh-based localization. -1 otherwise (e.g. non-mesh or multi-dimension
  // layouts).
  virtual int GetDOFHolderEntityDim() const { return -1; }

  virtual Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const = 0;

  virtual Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationIdsHost() const = 0;

  virtual ~FieldLayout() noexcept = default;

protected:
  // Builds the global-to-local permutation from the GIDs and entity offsets of
  // the DOF holders with ghosts (the permutation is indexed with ghosts).
  // Derived layouts that support GetGlobalToLocalPermutation() must call this
  // from their constructor, once their GIDs and entity offsets are available;
  // building eagerly keeps first access free of the data race that lazy
  // construction of the mutable cache would introduce.
  void BuildGlobalToLocalPermutation(
    Kokkos::View<const GO*, DeviceMemorySpace> local_gids,
    const EntOffsetsArray& local_ent_offsets);

  void SetCoordinateSystem(
    std::shared_ptr<const CoordinateSystem> coordinate_system)
  {
    if (coordinate_system == nullptr) {
      throw pcms_error(
        "FieldLayout::SetCoordinateSystem: coordinate system must not be null");
    }
    system_ = std::move(coordinate_system);
  }

private:
  std::shared_ptr<const CoordinateSystem> system_;
  std::string name_;
  Kokkos::View<LO*, HostMemorySpace> global_to_local_host_;
  Kokkos::View<LO*, DeviceMemorySpace> global_to_local_;
};

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

} // namespace pcms
#endif // PCMS_FIELD_LAYOUT_H
