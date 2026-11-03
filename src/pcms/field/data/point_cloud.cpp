#include "pcms/field/data/point_cloud.h"
#include "pcms/field/layout/point_cloud.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/arrays.h"

namespace pcms
{

PointCloud::PointCloud(std::shared_ptr<const PointCloudLayout> layout)
  : layout_(std::move(layout)),
    basis_{},
    device_data_("", layout_->GetDOFHolderCoordinates().GetValues().extent(0),
                 static_cast<size_t>(layout_->GetNumComponents()))
{
}

FieldValueType PointCloud::GetValueType() const
{
  return ValueTypeOfRank(basis_.Rank());
}

const ValueBasis& PointCloud::GetValueBasis() const
{
  return basis_;
}

Rank2View<const Real, DeviceMemorySpace> PointCloud::GetOwnedValues() const
{
  return MakeConstRank2View(device_data_);
}

void PointCloud::SetOwnedValues(Rank2View<const Real, DeviceMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;
  CopyDeviceRank2ViewToDeviceView(device_data_, data);
}

} // namespace pcms
