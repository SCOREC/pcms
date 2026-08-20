#include "pcms/field/layout/empty.h"
#include "pcms/utility/assert.h"

namespace pcms
{

namespace
{

// helper/indirection to ensure that we don't end up with a Cartesian system
// with dimension 0
int RequireDimension(const std::shared_ptr<const CoordinateSystem>& system)
{
  if (system == nullptr) {
    throw pcms_error("EmptyFieldLayout: coordinate system must not be null");
  }
  if (system->Dimension() == 0) {
    throw pcms_error(
      "EmptyFieldLayout: requires a concrete coordinate system; the "
      "dimension-deferred Cartesian placeholder has no coordinate data "
      "here to resolve against");
  }
  return system->Dimension();
}

} // namespace

EmptyFieldLayout::EmptyFieldLayout(
  std::shared_ptr<const CoordinateSystem> system)
  : owned_("null_owned", 0),
    gids_("null_gids", 0),
    class_dims_("null_class_dims", 0),
    class_ids_("null_class_ids", 0),
    owned_host_("null_owned_host", 0),
    gids_host_("null_gids_host", 0),
    classification_dims_host_("null_classification_dims_host", 0),
    classification_ids_host_("null_classification_ids_host", 0),
    coords_("null_coords", 0, RequireDimension(system))
{
  SetCoordinateSystem(std::move(system));
  discretization_ = std::make_shared<EmptyDiscretization>();
}

std::shared_ptr<const Discretization> EmptyFieldLayout::GetDiscretization()
  const noexcept
{
  return discretization_;
}

int EmptyFieldLayout::GetNumComponents() const
{
  return 1;
}

LO EmptyFieldLayout::GetNumOwnedDofHolder() const
{
  return 0;
}

GO EmptyFieldLayout::GetNumGlobalDofHolder() const
{
  return 0;
}

Rank1View<const bool, HostMemorySpace> EmptyFieldLayout::GetOwnedHost() const
{
  return make_const_array_view(owned_host_);
}

GlobalIDView<HostMemorySpace> EmptyFieldLayout::GetGidsHost() const
{
  return make_const_array_view(gids_host_);
}

bool EmptyFieldLayout::IsDistributed() const
{
  return false;
}

EntOffsetsArray EmptyFieldLayout::GetEntOffsets() const
{
  return {0, 0, 0, 0, 0};
}

CoordinateView<DeviceMemorySpace> EmptyFieldLayout::GetDOFHolderCoordinates()
  const
{
  auto coords_view = MakeConstRank2View(coords_);
  return CoordinateView<DeviceMemorySpace>{GetCoordinateSystem(), coords_view};
}

int EmptyFieldLayout::GetDimension() const
{
  return GetCoordinateSystem()->Dimension();
}

Rank1View<const LO, HostMemorySpace>
EmptyFieldLayout::GetDOFHolderClassificationDimensionsHost() const
{
  return make_const_array_view(classification_dims_host_);
}

Rank1View<const LO, HostMemorySpace>
EmptyFieldLayout::GetDOFHolderClassificationIdsHost() const
{
  return make_const_array_view(classification_ids_host_);
}

} // namespace pcms
