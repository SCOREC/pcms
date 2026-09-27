#ifndef PCMS_XGC_FIELD_LAYOUT_H
#define PCMS_XGC_FIELD_LAYOUT_H

#include "pcms/discretization/discretization/xgc.hpp"
#include "pcms/discretization/discretization/xgc_reverse_classification.h"
#include "pcms/field/field_layout.h"
#include <functional>

namespace pcms
{

class XGCFieldLayout : public FieldLayout
{
public:
  XGCFieldLayout(const ReverseClassificationVertex& reverse_classification,
                 std::function<int8_t(int, int)> in_overlap,
                 LO num_plane_nodes);

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

  LO GetFullDataSize() const noexcept;

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
  LO num_plane_nodes_;
  std::shared_ptr<const Discretization> discretization_;
};

} // namespace pcms

#endif // PCMS_XGC_FIELD_LAYOUT_H
