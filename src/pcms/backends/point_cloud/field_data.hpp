#ifndef PCMS_BACKENDS_POINT_CLOUD_FIELD_DATA_HPP
#define PCMS_BACKENDS_POINT_CLOUD_FIELD_DATA_HPP

#include "pcms/field/field_data.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/backends/point_cloud/layout.hpp"
#include <memory>

namespace pcms::point_cloud
{
class PointCloud
  : public FieldData<Real>
  , public DOFAccess<Real>
{
public:
  PointCloud(std::shared_ptr<const pcms::detail::PointCloudLayout> layout);

  FieldValueType GetValueType() const override;
  const ValueBasis& GetValueBasis() const override;

  const DOFAccess<Real>* GetDOFAccess() const noexcept override { return this; }
  DOFAccess<Real>* GetDOFAccess() noexcept override { return this; }

  Rank2View<const Real, DeviceMemorySpace> GetValues() const override;
  void SetValues(Rank2View<const Real, DeviceMemorySpace> data) override;

protected:
  void LinearSum(Real a, const FieldData<Real>& x, Real b,
                 const FieldData<Real>& y) override;
  void Scale(Real c, const FieldData<Real>& x) override;
  void Const(Real c) override;
  void Div(const FieldData<Real>& x, const FieldData<Real>& y) override;
  [[nodiscard]] Real DotProdLocal(const FieldData<Real>& x) const override;
  [[nodiscard]] Real MinLocal() const override;
  [[nodiscard]] Real MaxNormLocal() const override;

private:
  static const Kokkos::View<Real**, DeviceMemorySpace>& Data(
    const FieldData<Real>& x);

  std::shared_ptr<const pcms::detail::PointCloudLayout> layout_;
  ValueBasis basis_;
  Kokkos::View<Real**, DeviceMemorySpace> device_data_;
};
} // namespace pcms::point_cloud

#endif // PCMS_BACKENDS_POINT_CLOUD_FIELD_DATA_HPP
