#include "pcms/adapter/omega_h/omega_h_field_layout.h"
#include "pcms/adapter/omega_h/omega_h_field2.h"

namespace pcms
{
/*
 * Field
 */
OmegaHField2::OmegaHField2(std::string name, CoordinateSystem coordinate_system,
                         const OmegaHFieldLayout& layout, Omega_h::Mesh& mesh)
  : coordinate_system_(coordinate_system), layout_(layout), mesh_(mesh)
{
}

void OmegaHField2::SetEvaluationCoordinates(CoordinateView<HostMemorySpace> coordinates,
                                           LocalizationHint hint )
{
}

LocalizationHint OmegaHField2::GetLocalizationHint(CoordinateView<HostMemorySpace> coordinates) {}

void OmegaHField2::Evaluate(FieldDataView<double, HostMemorySpace> results)
{
  // TODO decide if we want to implicitly perform the coordinate transformations
  // when possible
  if (results.GetCoordinateSystem() != coordinate_system_) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Coordinate system mismatch");
  }
  // if field is piecewise
  // 1. use localization to get the appropriate element
  // 2. get the face/region value in the element
  // OPEN Question: what to do if point is on the interface between two regions?
  // a) return mean?, b) return 2 values?, c) pick one at random

  // if the field is linear:
  // 1. use localization to get the appropriate element
  // 2. get vertexes of the element (e2v map)
  // 3. get the values of the vertexes
  // 4. get barycentric coordinates of the point in the element
  // 5. interpolate the values of the vertexes to the point
}

void OmegaHField2::EvaluateGradient(FieldDataView<double, HostMemorySpace> results)
{
  // TODO when moved to PCMS throw PCMS exception
  throw std::runtime_error("Not implemented");
}

const FieldLayout& OmegaHField2::GetLayout()
{
  return layout_;
}

bool OmegaHField2::CanEvaluateGradient()
{
  // TODO compute the gradient field using element shape functions
  return false;
}

int OmegaHField2::Serialize(
  Rank1View<double, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const
{
  return 1;
}

void OmegaHField2::Deserialize(
  Rank1View<const double, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const
{
}

} // namespace pcms
