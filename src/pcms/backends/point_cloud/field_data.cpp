#include "pcms/backends/point_cloud/field_data.hpp"
#include "pcms/backends/point_cloud/layout.hpp"
#include "pcms/field/detail/field_data_ops.hpp"
#include "pcms/utility/detail/profile.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/arrays.hpp"

namespace pcms::point_cloud
{

PointCloud::PointCloud(
  std::shared_ptr<const pcms::detail::PointCloudLayout> layout)
  : layout_(std::move(layout)),
    basis_{},
    device_data_("", layout_->GetDOFHolderCoordinates().GetValues().extent(0),
                 static_cast<size_t>(layout_->GetNumComponents()))
{
}

FieldValueType PointCloud::GetValueType() const
{
  return pcms::detail::ValueTypeOfRank(basis_.Rank());
}

const ValueBasis& PointCloud::GetValueBasis() const
{
  return basis_;
}

Rank2View<const Real, DeviceMemorySpace> PointCloud::GetValues() const
{
  return MakeConstRank2View(device_data_);
}

void PointCloud::SetValues(Rank2View<const Real, DeviceMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;
  pcms::detail::CopyDeviceRank2ViewToDeviceView(device_data_, data);
}

const Kokkos::View<Real**, DeviceMemorySpace>& PointCloud::Data(
  const FieldData<Real>& x)
{
  return pcms::detail::OperandOfKind<PointCloud>(x).device_data_;
}

void PointCloud::LinearSum(Real a, const FieldData<Real>& x, Real b,
                           const FieldData<Real>& y)
{
  pcms::detail::LinearSumValues(device_data_, a, Data(x), b, Data(y));
}

void PointCloud::Scale(Real c, const FieldData<Real>& x)
{
  pcms::detail::ScaleValues(device_data_, c, Data(x));
}

void PointCloud::Const(Real c)
{
  pcms::detail::ConstValues(device_data_, c);
}

void PointCloud::Div(const FieldData<Real>& x, const FieldData<Real>& y)
{
  pcms::detail::DivValues(device_data_, Data(x), Data(y));
}

Real PointCloud::DotProdLocal(const FieldData<Real>& x) const
{
  return pcms::detail::DotProdValues(device_data_, Data(x));
}

Real PointCloud::MinLocal() const
{
  return pcms::detail::MinValues(device_data_);
}

Real PointCloud::MaxNormLocal() const
{
  return pcms::detail::MaxNormValues(device_data_);
}

} // namespace pcms::point_cloud
