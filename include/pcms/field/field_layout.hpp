#ifndef PCMS_FIELD_FIELD_LAYOUT_HPP
#define PCMS_FIELD_FIELD_LAYOUT_HPP
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "pcms/utility/types.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/geometry/bound_points.hpp"
#include <array>
#include <mpi.h>

namespace pcms
{

constexpr int ent_offsets_len = 5;
using EntOffsetsArray = std::array<size_t, ent_offsets_len>;

// FieldLayout describes how a field's DOF holders are laid out on this rank.
//
// The accessors describe the DOF holders this rank owns (rank-exclusive),
// compact and in owned order. A distributed layout also holds ghost copies of
// DOF holders other ranks own, for internal use (field storage, ghost
// synchronization, element assembly): together they are the DOF holders with
// ghosts, in mesh order. Only their ownership and the owned indices are exposed
// here. For a layout without ghosts the two coincide (!HasGhosts()).
class FieldLayout
{
public:
  // Optional identifier for this layout. Given at construction (the layouts a
  // function space or factory creates take it from their options) and empty by
  // default. Consumers that need to identify a layout by a stable name use it;
  // it carries no meaning within the field layer itself.
  [[nodiscard]] const std::string& GetName() const noexcept { return name_; }

  /// The geometry the layout's entities are defined on: spaces and layouts on
  /// one geometry object share their entities (OwnerHints name entities of a
  /// geometry). Null for a layout without one: XGC, MFEM, empty, and a bare
  /// point cloud.
  ///
  /// A layout never holds a point-cloud geometry built on it (that geometry
  /// holds the layout as its site layout, which would be a cycle of
  /// shared_ptr); the geometry of a layout is the one whose entities the DOF
  /// holders are.
  [[nodiscard]] const std::shared_ptr<const Geometry>& GetGeometry()
    const noexcept
  {
    return geometry_;
  }

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

  /// @brief The communicator of the ranks the layout is distributed over
  [[nodiscard]] virtual MPI_Comm GetCommunicator() const = 0;

  /// @return whether the layout spans more than one rank
  [[nodiscard]] bool IsDistributed() const;

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

  /// @return the coordinate system of the layout's geometry; of a layout
  ///         without a geometry (an empty layout, a bare point cloud, XGC),
  ///         the one it was given (SetCoordinateSystem)
  /// @throws pcms_error if the layout has neither
  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const;

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
  FieldLayout() = default;
  explicit FieldLayout(std::string name) : name_(std::move(name)) {}

  // Builds the global-to-local permutation from the GIDs and entity offsets of
  // the DOF holders with ghosts (the permutation is indexed with ghosts).
  // Derived layouts that support GetGlobalToLocalPermutation() must call this
  // from their constructor, once their GIDs and entity offsets are available;
  // building eagerly keeps first access free of the data race that lazy
  // construction of the mutable cache would introduce.
  void BuildGlobalToLocalPermutation(
    Kokkos::View<const GO*, DeviceMemorySpace> local_gids,
    const EntOffsetsArray& local_ent_offsets);

  /// Sets the geometry GetGeometry() returns; the layout must not hold a
  /// point-cloud geometry built on it.
  void SetGeometry(std::shared_ptr<const Geometry> geometry)
  {
    geometry_ = std::move(geometry);
  }

  /// Sets the coordinate system of a layout that has no geometry.
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
  std::shared_ptr<const Geometry> geometry_;
  std::string name_;
  Kokkos::View<LO*, HostMemorySpace> global_to_local_host_;
  Kokkos::View<LO*, DeviceMemorySpace> global_to_local_;
};

/// @brief The entities of the layout's geometry the DOF holders of `layout`
///        are, as owner hints for binding or searching at the DOF-holder
///        coordinates (the hint names the geometry the layout is on, which
///        verifies it, see OwnerHints).
///
/// DOF holder i is the i-th entity of the layout's one populated entity
/// dimension, so the hinted entity ids are i - (the first DOF holder of that
/// dimension), i.e. 0, 1, ...
/// @return the hint when the layout has a geometry, one DOF holder per entity
///         of one dimension of it and no ghosts; an empty one otherwise
[[nodiscard]] OwnerHints DOFHolderOwners(const FieldLayout& layout);

} // namespace pcms
#include "pcms/field/detail/field_layout.hpp"

#endif // PCMS_FIELD_FIELD_LAYOUT_HPP
