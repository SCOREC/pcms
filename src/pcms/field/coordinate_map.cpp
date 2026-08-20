#include "pcms/field/coordinate_map.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/profile.h"
#include <Kokkos_Core.hpp>
#include <string>
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"

namespace pcms
{

namespace
{

using ExecutionSpace = DeviceMemorySpace::execution_space;

Kokkos::View<Real**, DeviceMemorySpace> AllocateOutput(size_t num_points,
                                                       int dim)
{
  return Kokkos::View<Real**, DeviceMemorySpace>(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "mapped_points"),
    num_points, static_cast<size_t>(dim));
}

// Mapping kernels as free functions: CUDA extended lambdas may not be defined
// inside member functions of classes with private members.
void MapCartesianToCylindrical(Rank2View<const Real, DeviceMemorySpace> in,
                               Kokkos::View<Real**, DeviceMemorySpace> out)
{
  Kokkos::parallel_for(
    "cartesian_to_cylindrical_points",
    Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<LO>(in.extent(0))),
    KOKKOS_LAMBDA(const LO i) {
      const Real x = in(i, 0);
      const Real y = in(i, 1);
      out(i, 0) = Kokkos::sqrt(x * x + y * y);
      out(i, 1) = Kokkos::atan2(y, x);
      out(i, 2) = in(i, 2);
    });
}

void MapCylindricalToCartesian(Rank2View<const Real, DeviceMemorySpace> in,
                               Kokkos::View<Real**, DeviceMemorySpace> out)
{
  Kokkos::parallel_for(
    "cylindrical_to_cartesian_points",
    Kokkos::RangePolicy<ExecutionSpace>(0, static_cast<LO>(in.extent(0))),
    KOKKOS_LAMBDA(const LO i) {
      const Real r = in(i, 0);
      const Real theta = in(i, 1);
      out(i, 0) = r * Kokkos::cos(theta);
      out(i, 1) = r * Kokkos::sin(theta);
      out(i, 2) = in(i, 2);
    });
}

} // namespace

void CoordinateMap::ValidateSourcePoints(
  const CoordinateView<DeviceMemorySpace>& points, const char* who) const
{
  if (!SameCoordinateSystem(points.GetCoordinateSystem(),
                            GetSourceCoordinateSystem())) {
    throw pcms_error(std::string(who) +
                     ": points are tagged with coordinate system '" +
                     std::string(points.GetCoordinateSystem()->Kind()) +
                     "' but this map's source coordinate system is '" +
                     std::string(GetSourceCoordinateSystem()->Kind()) +
                     "' (coordinate systems compare by object identity)");
  }
}

std::unique_ptr<BoundBasisTransformation>
CoordinateMap::MakeBasisTransformation(
  const CoordinateView<DeviceMemorySpace>& /*source_points*/,
  const MappedPoints& /*mapped*/) const
{
  return nullptr;
}

std::shared_ptr<const CoordinateSystem>
CartesianToCylindrical::GetSourceCoordinateSystem() const noexcept
{
  return csys::Cartesian::Create(3);
}

std::shared_ptr<const CoordinateSystem>
CartesianToCylindrical::GetTargetCoordinateSystem() const noexcept
{
  return csys::CylindricalRThetaZ::Create();
}

MappedPoints CartesianToCylindrical::Map(
  const CoordinateView<DeviceMemorySpace>& points) const
{
  PCMS_FUNCTION_TIMER;
  ValidateSourcePoints(points, "CartesianToCylindrical");
  const auto in = points.GetValues();
  auto out = AllocateOutput(in.extent(0), 3);
  MapCartesianToCylindrical(in, out);
  // The kernel reads the caller's point buffer; Map's contract is that the
  // caller may free it as soon as Map returns.
  ExecutionSpace().fence("CartesianToCylindrical::Map: complete mapping");
  return MappedPoints{GetTargetCoordinateSystem(), out, PointStatusView{}};
}

std::unique_ptr<BoundBasisTransformation>
CartesianToCylindrical::MakeBasisTransformation(
  const CoordinateView<DeviceMemorySpace>& /*source_points*/,
  const MappedPoints& mapped) const
{
  // Bound at the cylindrical outputs: theta is already computed there.
  return CylindricalToCartesianBasis{}.Bind(mapped.View());
}

std::shared_ptr<const CoordinateSystem>
CylindricalToCartesian::GetSourceCoordinateSystem() const noexcept
{
  return csys::CylindricalRThetaZ::Create();
}

std::shared_ptr<const CoordinateSystem>
CylindricalToCartesian::GetTargetCoordinateSystem() const noexcept
{
  return csys::Cartesian::Create(3);
}

MappedPoints CylindricalToCartesian::Map(
  const CoordinateView<DeviceMemorySpace>& points) const
{
  PCMS_FUNCTION_TIMER;
  ValidateSourcePoints(points, "CylindricalToCartesian");
  const auto in = points.GetValues();
  auto out = AllocateOutput(in.extent(0), 3);
  MapCylindricalToCartesian(in, out);
  ExecutionSpace().fence("CylindricalToCartesian::Map: complete mapping");
  return MappedPoints{GetTargetCoordinateSystem(), out, PointStatusView{}};
}

std::unique_ptr<BoundBasisTransformation>
CylindricalToCartesian::MakeBasisTransformation(
  const CoordinateView<DeviceMemorySpace>& source_points,
  const MappedPoints& /*mapped*/) const
{
  return CartesianToCylindricalBasis{}.Bind(source_points);
}

} // namespace pcms
