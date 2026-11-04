#ifndef POINT_CLOUD_H_
#define POINT_CLOUD_H_

#include "pcms/field/field_data.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/field/layout/point_cloud.hpp"
#include <memory>

namespace pcms
{
class PointCloud
  : public FieldData<Real>
  , public DOFAccess<Real>
{
public:
  PointCloud(std::shared_ptr<const PointCloudLayout> layout);

  FieldValueType GetValueType() const override;
  const ValueBasis& GetValueBasis() const override;

  const DOFAccess<Real>* GetDOFAccess() const noexcept override { return this; }
  DOFAccess<Real>* GetDOFAccess() noexcept override { return this; }

  Rank2View<const Real, DeviceMemorySpace> GetOwnedValues() const override;
  void SetOwnedValues(Rank2View<const Real, DeviceMemorySpace> data) override;

private:
  std::shared_ptr<const PointCloudLayout> layout_;
  ValueBasis basis_;
  Kokkos::View<Real**, DeviceMemorySpace> device_data_;
};
} // namespace pcms

#endif // POINT_CLOUD_H_
