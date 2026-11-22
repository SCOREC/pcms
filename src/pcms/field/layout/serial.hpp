#ifndef PCMS_FIELD_LAYOUT_SERIAL_HPP
#define PCMS_FIELD_LAYOUT_SERIAL_HPP

#include "pcms/field/field_layout.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <utility>

namespace pcms::detail
{

/// @brief A layout without ghosts whose DOF holders are given as arrays
class SerialLayout final : public FieldLayout
{
public:
  /// @param coords the position of each DOF holder
  /// @param num_global the number of DOF holders of all ranks
  SerialLayout(std::shared_ptr<const Geometry> geometry,
               const Kokkos::View<Real**, HostMemorySpace>& coords,
               Kokkos::View<GO*, HostMemorySpace> gids,
               Kokkos::View<bool*, HostMemorySpace> owned, GO num_global,
               int num_components = 1)
    : num_global_(num_global),
      num_components_(num_components),
      gids_(std::move(gids)),
      owned_(std::move(owned)),
      zeros_("zeros", gids_.extent(0)),
      coords_(ToDeviceView(coords))
  {
    SetCoordinateSystem(geometry->GetCoordinateSystem());
    SetGeometry(std::move(geometry));
  }

  int GetNumComponents() const override { return num_components_; }
  LO GetNumDofHolder() const override
  {
    return static_cast<LO>(gids_.extent(0));
  }
  GO GetNumGlobalDofHolder() const override { return num_global_; }
  bool IsDistributed() const override { return false; }
  Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override
  {
    return MakeRank1View(owned_);
  }
  GlobalIDView<HostMemorySpace> GetGidsHost() const override
  {
    return GlobalIDView<HostMemorySpace>(gids_.data(), gids_.extent(0));
  }
  EntOffsetsArray GetEntOffsets() const override
  {
    const auto n = gids_.extent(0);
    return {0, n, n, n, n};
  }
  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override
  {
    return CoordinateView<DeviceMemorySpace>(GetCoordinateSystem(),
                                             MakeConstRank2View(coords_));
  }
  int GetDimension() const override
  {
    return static_cast<int>(coords_.extent(1));
  }
  Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override
  {
    return MakeRank1View(zeros_);
  }
  Rank1View<const LO, HostMemorySpace> GetDOFHolderClassificationIdsHost()
    const override
  {
    return MakeRank1View(zeros_);
  }

private:
  GO num_global_;
  int num_components_;
  Kokkos::View<GO*, HostMemorySpace> gids_;
  Kokkos::View<bool*, HostMemorySpace> owned_;
  Kokkos::View<LO*, HostMemorySpace> zeros_;
  Kokkos::View<Real**, DeviceMemorySpace> coords_;
};

} // namespace pcms::detail

#endif // PCMS_FIELD_LAYOUT_SERIAL_HPP
