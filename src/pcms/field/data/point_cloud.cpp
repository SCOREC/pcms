#include "pcms/field/data/point_cloud.h"
#include "pcms/field/layout/point_cloud.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/arrays.h"

namespace pcms
{

PointCloud::PointCloud(std::shared_ptr<const PointCloudLayout> layout)
  : layout_(std::move(layout)),
    metadata_{},
    device_data_("", layout_->GetDOFHolderCoordinates().GetValues().extent(0)),
    data_host_("", layout_->GetDOFHolderCoordinates().GetValues().extent(0))
{
}

const FieldMetadata& PointCloud::GetMetadata() const
{
  return metadata_;
}

Rank2View<const Real, HostMemorySpace> PointCloud::GetDOFHolderDataHost() const
{
  Kokkos::deep_copy(data_host_, device_data_);
  const auto nc = layout_->GetNumComponents();
  return Rank2View<const Real, HostMemorySpace>(
    data_host_.data(), static_cast<LO>(data_host_.size()) / nc, nc);
}

void PointCloud::SetDOFHolderDataHost(
  Rank2View<const Real, HostMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;
  PCMS_ALWAYS_ASSERT(data.size() == device_data_.size());
  CopyHostRank2ViewToDeviceView(device_data_, data);
}

Rank2View<const Real, DeviceMemorySpace> PointCloud::GetDOFHolderData() const
{
  const auto nc = layout_->GetNumComponents();
  return Rank2View<const Real, DeviceMemorySpace>(
    device_data_.data(), static_cast<LO>(device_data_.size()) / nc, nc);
}

void PointCloud::SetDOFHolderData(Rank2View<const Real, DeviceMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;
  PCMS_ALWAYS_ASSERT(data.size() == device_data_.size());
  CopyDeviceRank2ViewToDeviceView(device_data_, data);
}

} // namespace pcms
