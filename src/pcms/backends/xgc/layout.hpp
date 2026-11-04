#ifndef PCMS_BACKENDS_XGC_LAYOUT_HPP
#define PCMS_BACKENDS_XGC_LAYOUT_HPP

#include "pcms/backends/xgc/discretization.hpp"
#include "pcms/backends/xgc/reverse_classification.hpp"
#include "pcms/field/field_layout.hpp"

namespace pcms::backend::xgc
{

class XGCFieldLayout : public FieldLayout
{
public:
  XGCFieldLayout(const ReverseClassificationVertex& reverse_classification,
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

} // namespace pcms::backend::xgc

#endif // PCMS_BACKENDS_XGC_LAYOUT_HPP
