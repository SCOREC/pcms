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

  // nodes for standard lagrange FEM
  LO virtual GetNumOwnedDofHolder() const = 0;
  GO virtual GetNumGlobalDofHolder() const = 0;

  // size of buffer that needs to be allocated to represent the field
  // # components * NumDOFHolder
  LO OwnedSize() const { return GetNumComponents() * GetNumOwnedDofHolder(); };
  GO GlobalSize() const
  {
    return GetNumComponents() * GetNumGlobalDofHolder();
  };

  virtual Rank1View<const bool, HostMemorySpace> GetOwnedHost() const = 0;
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

  // returns true if the field layout is distributed
  // if the field layout is distributed, the owned and global dofs are the same
  [[nodiscard]] virtual bool IsDistributed() const = 0;

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

  virtual Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const = 0;

  virtual Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationIdsHost() const = 0;

  virtual ~FieldLayout() noexcept = default;

protected:
  // Builds the global-to-local permutation from GetGidsHost()/GetEntOffsets().
  // Derived layouts that support GetGlobalToLocalPermutation() must call this
  // from their constructor, once their GIDs and entity offsets are available;
  // building eagerly keeps first access free of the data race that lazy
  // construction of the mutable cache would introduce.
  void BuildGlobalToLocalPermutation();

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

} // namespace pcms
#endif // PCMS_FIELD_LAYOUT_H
