#ifndef POINT_CLOUD_H_
#define POINT_CLOUD_H_

#include "pcms/field.h"
#include "pcms/arrays.h"
#include "pcms/coordinate_system.h"
#include "point_cloud_layout.h"

namespace pcms {
class PointCloud : public FieldT<Real> {
public:
  PointCloud(std::string name, CoordinateSystem coordinate_system,
             const PointCloudLayout& layout);

  const std::string& GetName() const override;

  CoordinateSystem GetCoordinateSystem() const override;

  LocalizationHint GetLocalizationHint(
    CoordinateView<HostMemorySpace> coordinate_view) const override;

  void Evaluate(LocalizationHint location, FieldDataView<double, HostMemorySpace> results) const override;

  void EvaluateGradient(
    FieldDataView<double, HostMemorySpace> results) override;

  const FieldLayout& GetLayout() const override;

  bool CanEvaluateGradient() override;

  int Serialize(Rank1View<double, pcms::HostMemorySpace> buffer,
                Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
    const override;

  void Deserialize(Rank1View<const double, pcms::HostMemorySpace> buffer,
                   Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
    override;

  FieldDataView<const Real, HostMemorySpace> GetDOFHolderData() const override;
  void SetDOFHolderData(FieldDataView<const Real, HostMemorySpace> data) override;
  CoordinateView<HostMemorySpace> GetDOFHolderCoordinates() const;
private:
  std::string name_;
  CoordinateSystem coordinate_system_;
  const PointCloudLayout& layout_;
  Kokkos::View<Real*> data_;
};
}

#endif // POINT_CLOUD_H_
