#ifndef POINT_CLOUD_LAYOUT_H_
#define POINT_CLOUD_LAYOUT_H_

#include "pcms/field.h"

namespace pcms
{

class PointCloudLayout : public FieldLayout
{
public:
  PointCloudLayout(int dim, Kokkos::View<Real**> coords,
                   CoordinateSystem coordinate_system);

  std::unique_ptr<FieldT<Real>> CreateField() const override;

  int GetNumComponents() const override;
  // nodes for standard lagrange FEM
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwned() const override;
  GlobalIDView<HostMemorySpace> GetGids() const override;
  CoordinateView<HostMemorySpace> GetDOFHolderCoordinates() const override;

  bool IsDistributed() override;
  size_t GetNumEnts() const;
  EntOffsetsArray GetEntOffsets() const override;

  ReversePartitionMap2 GetReversePartitionMap(
    const redev::Partition& partition) const override;

  std::array<int, 4> GetNodesPerDim() const;

private:
  int dim_;
  int components_;
  CoordinateSystem coordinate_system_;
  Kokkos::View<Real**> coords_;
  Kokkos::View<bool*> owned_;
  Kokkos::View<GO*> gids_;
};
} // namespace pcms

#endif // POINT_CLOUD_LAYOUT_H_
