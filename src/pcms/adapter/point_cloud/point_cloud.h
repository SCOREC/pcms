#ifndef POINT_CLOUD_H_
#define POINT_CLOUD_H_

#include "pcms/field.h"
#include "pcms/arrays.h"
#include "point_cloud_layout.h"

namespace pcms
{
class PointCloud : public FieldT<Real>
{
public:
  PointCloud(const PointCloudLayout& layout);

  LocalizationHint GetLocalizationHint(
    CoordinateView<HostMemorySpace> coordinate_view) const override;

  void Evaluate(LocalizationHint location,
                FieldDataView<Real, HostMemorySpace> results) const override;

  void EvaluateGradient(FieldDataView<Real, HostMemorySpace> results) override;

  const FieldLayout& GetLayout() const override;

  bool CanEvaluateGradient() override;

  int Serialize(Rank1View<Real, pcms::HostMemorySpace> buffer,
                Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
    const override;

  void Deserialize(
    Rank1View<const Real, pcms::HostMemorySpace> buffer,
    Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) override;

  Rank1View<const Real, HostMemorySpace> GetDOFHolderData() const override;
  void SetDOFHolderData(Rank1View<const Real, HostMemorySpace> data) override;

private:
  const PointCloudLayout& layout_;
  Kokkos::View<Real*> data_;
};
} // namespace pcms

#endif // POINT_CLOUD_H_
