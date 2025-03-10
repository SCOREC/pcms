#ifndef PCMS_FIELD_LAYOUT_H
#define PCMS_FIELD_LAYOUT_H
#include "pcms/arrays.h"

namespace pcms
{
class FieldLayout
{
public:
  // number of components
  int virtual GetNumComponents() = 0;

  // nodes for standard lagrange FEM
  LO virtual GetNumOwnedDofHolder() = 0;
  GO virtual GetNumGlobalDofHolder() = 0;

  // size of buffer that needs to be allocated to represent the field
  // # components * NumDOFHolder
  LO OwnedSize() { return GetNumComponents() * GetNumOwnedDofHolder(); };
  GO GlobalSize() { return GetNumComponents() * GetNumGlobalDofHolder(); };

  virtual GlobalIDView<HostMemorySpace> GetOwnedGids() = 0;

  // returns true if the field layout is distributed
  // if the field layout is distributed, the owned and global dofs are the same
  virtual bool IsDistributed() = 0;

  // This class should construct the permutation arrays that are needed
  // for serialization / deserialization

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
