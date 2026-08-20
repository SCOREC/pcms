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
    device_data_("", layout_->GetDOFHolderCoordinates().GetValues().extent(0),
                 static_cast<size_t>(layout_->GetNumComponents())),
    data_host_("", layout_->GetDOFHolderCoordinates().GetValues().extent(0),
               static_cast<size_t>(layout_->GetNumComponents()))
{
}

const FieldMetadata& PointCloud::GetMetadata() const
{
  return metadata_;
}

Rank2View<const Real, HostMemorySpace> PointCloud::GetDOFHolderDataHost() const
{
  DeepCopyMismatchLayouts(data_host_, device_data_);
  return MakeConstRank2View(data_host_);
}

void PointCloud::SetDOFHolderDataHost(
  Rank2View<const Real, HostMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;
  CopyHostRank2ViewToDeviceView(device_data_, data);
}

Rank2View<const Real, DeviceMemorySpace> PointCloud::GetDOFHolderData() const
{
  return MakeConstRank2View(device_data_);
}

void PointCloud::SetDOFHolderData(Rank2View<const Real, DeviceMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;
  CopyDeviceRank2ViewToDeviceView(device_data_, data);
}

} // namespace pcms
