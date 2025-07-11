#ifndef PCMS_FIELD_LAYOUT_H
#define PCMS_FIELD_LAYOUT_H
#include <redev.h>
#include "pcms/arrays.h"

namespace pcms
{
struct PartitionMapping
{
  std::vector<LO> indices;
  std::array<int, 4> ent_offsets;

  PartitionMapping() { ent_offsets.fill(0); }
};

using ReversePartitionMap = std::map<pcms::LO, std::vector<pcms::LO>>;
using ReversePartitionMap2 = std::map<pcms::LO, PartitionMapping>;

class FieldLayout
{
public:
  // number of components
  int virtual GetNumComponents() const = 0;

  // nodes for standard lagrange FEM
  LO virtual GetNumOwnedDofHolder() const = 0;
  GO virtual GetNumGlobalDofHolder() const = 0;

  // size of buffer that needs to be allocated to represent the field
  // # components * NumDOFHolder
  LO OwnedSize() const { return GetNumComponents() * GetNumOwnedDofHolder(); };
  GO GlobalSize() const { return GetNumComponents() * GetNumGlobalDofHolder(); };

  virtual GlobalIDView<HostMemorySpace> GetOwnedGids() const = 0;
  virtual GlobalIDView<HostMemorySpace> GetGids() const = 0;

  // returns true if the field layout is distributed
  // if the field layout is distributed, the owned and global dofs are the same
  virtual bool IsDistributed() = 0;

  // This class should construct the permutation arrays that are needed
  // for serialization / deserialization
  //

  virtual std::array<size_t, 4> GetEntOffsets() const = 0;

  virtual ReversePartitionMap2 GetReversePartitionMap(
    const redev::Partition& partition) const = 0;

  // Serialize, Derserialize, ReversePartitionMap?
  // GetOwnedDofHolderCoordinates(CoordinateSystem);

  // Serialize(FieldDataView, SerializationBuffer);
  // Deserialize(SerializationBuffer, FieldDataView);

  // Adjacency information
  // TODO: Need a GraphView class (simply two rank1 arrays CSR matrix w/o
  // values) virtual bool HasAdjacency() = 0; virtual GraphView GetAdjacency(LO
  // dim) = 0;

  virtual ~FieldLayout() noexcept = default;
};

} // namespace pcms
#endif // PCMS_FIELD_LAYOUT_H
