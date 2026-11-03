#ifndef POINT_CLOUD_LAYOUT_H_
#define POINT_CLOUD_LAYOUT_H_

#include "pcms/field/field.h"
#include "pcms/discretization/discretization/point_cloud.hpp"

namespace pcms
{

class PointCloudLayout : public FieldLayout
{
public:
  PointCloudLayout(int dim, Kokkos::View<Real**> coords,
                   std::shared_ptr<const CoordinateSystem> coordinate_system);
  PointCloudLayout(int dim, Kokkos::View<Real**> coords,
                   std::shared_ptr<const CoordinateSystem> coordinate_system,
                   std::shared_ptr<const Discretization> discretization,
                   int classification_entity_dim);

  std::shared_ptr<const Discretization> GetDiscretization()
    const noexcept override;

  int GetNumComponents() const override;
  // nodes for standard lagrange FEM
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override;
  GlobalIDView<HostMemorySpace> GetGidsHost() const override;
  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override;

  [[nodiscard]] bool IsDistributed() const override;
  size_t GetNumEnts() const;
  EntOffsetsArray GetEntOffsets() const override;

  int GetDimension() const override;

  Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override;

  Rank1View<const LO, HostMemorySpace> GetDOFHolderClassificationIdsHost()
    const override;

  std::array<int, 4> GetNodesPerDim() const;

  Kokkos::View<const Real**, DeviceMemorySpace> GetCoordinates() const;

private:
  int dim_;
  int components_;
  Kokkos::View<Real**> coords_;
  Kokkos::View<bool*> owned_;
  Kokkos::View<GO*> gids_;
  Kokkos::View<bool*, HostMemorySpace> owned_host_;
  Kokkos::View<GO*, HostMemorySpace> gids_host_;
  Kokkos::View<LO*, DeviceMemorySpace> classification_dims_;
  Kokkos::View<LO*, DeviceMemorySpace> classification_ids_;
  Kokkos::View<LO*, HostMemorySpace> classification_dims_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_ids_host_;
  std::shared_ptr<const Discretization> discretization_;
};
} // namespace pcms

#endif // POINT_CLOUD_LAYOUT_H_
