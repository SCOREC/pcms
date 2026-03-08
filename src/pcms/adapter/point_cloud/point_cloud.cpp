#include "point_cloud.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/assert.h"

namespace pcms
{

struct CopyCoordinatesFunctor
{
  Kokkos::View<Real**> coordinates_;
  Rank2View<const Real, HostMemorySpace> coords_;

  CopyCoordinatesFunctor(Kokkos::View<Real**> coordinates,
                         Rank2View<const Real, HostMemorySpace> coords)
    : coordinates_(coordinates), coords_(coords)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const
  {
    for (int j = 0; j < coords_.extent(1); ++j) {
      coordinates_(i, j) = coords_(i, j);
    }
  }
};

struct PointCloudLocalizationHint
{
  PointCloudLocalizationHint(CoordinateView<HostMemorySpace> coordinate_view)
    : coordinates_("", coordinate_view.GetCoordinates().extent(0),
                   coordinate_view.GetCoordinates().extent(1))
  {
    auto coords = coordinate_view.GetCoordinates();
    int n = coords.extent(0);
    Kokkos::parallel_for(n, CopyCoordinatesFunctor(coordinates_, coords));
  }

  Kokkos::View<Real**> coordinates_;
};

PointCloud::PointCloud(const PointCloudLayout& layout)
  : layout_(layout),
    data_("", layout_.GetDOFHolderCoordinates().GetCoordinates().extent(0)),
    data_host_("", layout_.GetDOFHolderCoordinates().GetCoordinates().extent(0))
{
}

Rank1View<const Real, HostMemorySpace> PointCloud::GetDOFHolderData() const
{
  Kokkos::deep_copy(data_host_, data_);
  return make_const_array_view(data_host_);
}

void PointCloud::SetDOFHolderData(Rank1View<const Real, HostMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;
  PCMS_ALWAYS_ASSERT(data.size() == data_.size());
  Kokkos::parallel_for(
    Kokkos::RangePolicy<pcms::HostMemorySpace::execution_space>(0, data.size()),
    KOKKOS_CLASS_LAMBDA(int i) { data_host_(i) = data[i]; });
  Kokkos::deep_copy(data_, data_host_);
}

LocalizationHint PointCloud::GetLocalizationHint(
  CoordinateView<HostMemorySpace> coordinate_view) const
{
  auto hint = std::make_shared<PointCloudLocalizationHint>(coordinate_view);
  return LocalizationHint{hint};
}

void PointCloud::Evaluate(
  LocalizationHint /* unused */,
  FieldDataView<Real, HostMemorySpace> /* unused */) const
{
  throw std::runtime_error("Not implemented");
}

void PointCloud::EvaluateGradient(
  FieldDataView<Real, HostMemorySpace> /* unused */)
{
  throw std::runtime_error("Not implemented");
}

const FieldLayout& PointCloud::GetLayout() const
{
  return layout_;
}

bool PointCloud::CanEvaluateGradient()
{
  return false;
}

int PointCloud::Serialize(
  Rank1View<Real, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const
{
  PCMS_FUNCTION_TIMER;
  if (buffer.size() > 0) {
    Kokkos::parallel_for(
      data_.size(),
      KOKKOS_CLASS_LAMBDA(int i) { buffer[permutation[i]] = data_(i); });
  }
  return data_.size();
}

void PointCloud::Deserialize(
  Rank1View<const Real, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
{
  PCMS_FUNCTION_TIMER;
  Kokkos::parallel_for(
    data_.size(),
    KOKKOS_CLASS_LAMBDA(int i) { data_(i) = buffer[permutation[i]]; });
}

} // namespace pcms
