#ifndef POINT_CLOUD_H_
#define POINT_CLOUD_H_

#include "pcms/field/field_data.h"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/field/layout/point_cloud.h"
#include <memory>

namespace pcms
{
class PointCloud : public FieldData<Real>
{
public:
  PointCloud(std::shared_ptr<const PointCloudLayout> layout);

  FieldValueType GetValueType() const override;
  const ValueBasis& GetValueBasis() const override;

  Rank2View<const Real, HostMemorySpace> GetDOFHolderDataHost() const override;
  void SetDOFHolderDataHost(
    Rank2View<const Real, HostMemorySpace> data) override;

  Rank2View<const Real, DeviceMemorySpace> GetDOFHolderData() const override;
  void SetDOFHolderData(Rank2View<const Real, DeviceMemorySpace> data) override;

private:
  std::shared_ptr<const PointCloudLayout> layout_;
  ValueBasis basis_;
  Kokkos::View<Real**, DeviceMemorySpace> device_data_;
  mutable Kokkos::View<Real**, HostMemorySpace> data_host_;
};
} // namespace pcms

#endif // POINT_CLOUD_H_
