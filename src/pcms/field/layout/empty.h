#ifndef PCMS_EMPTY_FIELD_LAYOUT_H
#define PCMS_EMPTY_FIELD_LAYOUT_H

#include "pcms/discretization/discretization/empty.hpp"
#include "pcms/field/field_layout.h"
#include <Kokkos_Core.hpp>

namespace pcms
{

class EmptyFieldLayout : public FieldLayout
{
public:
  explicit EmptyFieldLayout(std::shared_ptr<const CoordinateSystem> system);

  std::shared_ptr<const Discretization> GetDiscretization()
    const noexcept override;

  int GetNumComponents() const override;
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;
  Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override;
  GlobalIDView<HostMemorySpace> GetGidsHost() const override;
  bool IsDistributed() const override;
  EntOffsetsArray GetEntOffsets() const override;
  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override;
  int GetDimension() const override;
  Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override;
  Rank1View<const LO, HostMemorySpace> GetDOFHolderClassificationIdsHost()
    const override;

private:
  Kokkos::View<bool*, DeviceMemorySpace> owned_;
  Kokkos::View<GO*, DeviceMemorySpace> gids_;
  Kokkos::View<LO*, DeviceMemorySpace> class_dims_;
  Kokkos::View<LO*, DeviceMemorySpace> class_ids_;
  Kokkos::View<bool*, HostMemorySpace> owned_host_;
  Kokkos::View<GO*, HostMemorySpace> gids_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_dims_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_ids_host_;
  Kokkos::View<Real**, DeviceMemorySpace> coords_;
  std::shared_ptr<const Discretization> discretization_;
};

} // namespace pcms

#endif // PCMS_EMPTY_FIELD_LAYOUT_H
