#include "pcms/field/basis_transformation.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/profile.h"
#include <string>
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"

namespace pcms
{

namespace
{

using ExecutionSpace = DeviceMemorySpace::execution_space;

const ValueBasis& CartesianPhysicalVector()
{
  static const ValueBasis basis{csys::Cartesian::Create(3),
                                ComponentScaling::Physical,
                                VarianceSignature{Variance::Contravariant}};
  return basis;
}

const ValueBasis& CylindricalPhysicalVector()
{
  static const ValueBasis basis{csys::CylindricalRThetaZ::Create(),
                                ComponentScaling::Physical,
                                VarianceSignature{Variance::Contravariant}};
  return basis;
}

basis::CartesianToCylindrical::State ComputeCosSinTheta(
  const char* who, const CoordinateView<DeviceMemorySpace>& points)
{
  Kokkos::View<Real*, DeviceMemorySpace> cos_theta;
  Kokkos::View<Real*, DeviceMemorySpace> sin_theta;
  const bool cartesian = SameCoordinateSystem(points.GetCoordinateSystem(),
                                              csys::Cartesian::Create(3));
  const bool cylindrical = SameCoordinateSystem(
    points.GetCoordinateSystem(), csys::CylindricalRThetaZ::Create());
  if (!cartesian && !cylindrical) {
    throw pcms_error(std::string(who) +
                     ": bound points must be tagged Cartesian-3 or "
                     "CylindricalRThetaZ, got '" +
                     std::string(points.GetCoordinateSystem()->Kind()) + "'");
  }
  const auto pts = points.GetValues();
  const auto n = static_cast<LO>(pts.extent(0));
  cos_theta = Kokkos::View<Real*, DeviceMemorySpace>(
    Kokkos::view_alloc(Kokkos::WithoutInitializing,
                       "basis_transformation_cos_theta"),
    n);
  sin_theta = Kokkos::View<Real*, DeviceMemorySpace>(
    Kokkos::view_alloc(Kokkos::WithoutInitializing,
                       "basis_transformation_sin_theta"),
    n);
  auto c = cos_theta;
  auto s = sin_theta;
  Kokkos::parallel_for(
    "basis_transformation_cos_sin_theta",
    Kokkos::RangePolicy<ExecutionSpace>(0, n), KOKKOS_LAMBDA(const LO i) {
      const Real theta =
        cartesian ? Kokkos::atan2(pts(i, 1), pts(i, 0)) : pts(i, 1);
      c(i) = Kokkos::cos(theta);
      s(i) = Kokkos::sin(theta);
    });
  // The kernel reads the caller's point buffer and only the derived cos/sin
  // arrays are cached; the caller may free the points once Bind returns, so
  // the read must complete first.
  ExecutionSpace().fence(
    "basis transformation bind: complete cos/sin derivation");
  return basis::CartesianToCylindrical::State{cos_theta, sin_theta};
}

const char* ValueTypeName(FieldValueType t)
{
  switch (t) {
    case FieldValueType::Scalar: return "Scalar";
    case FieldValueType::Vector: return "Vector";
    case FieldValueType::Tensor: return "Tensor";
  }
  return "<unknown>";
}

// Rotation kernels as free functions: CUDA extended lambdas may not be
// defined inside private member functions (the ApplyVectorImpl overrides).
void RotateCartesianToCylindrical(
  const Kokkos::View<Real*, DeviceMemorySpace>& c,
  const Kokkos::View<Real*, DeviceMemorySpace>& s,
  Rank2View<const Real, DeviceMemorySpace> in,
  Rank2View<Real, DeviceMemorySpace> out, LO n)
{
  Kokkos::parallel_for(
    "cartesian_to_cylindrical_basis", Kokkos::RangePolicy<ExecutionSpace>(0, n),
    KOKKOS_LAMBDA(const LO i) {
      const Real vx = in(i, 0);
      const Real vy = in(i, 1);
      out(i, 0) = vx * c(i) + vy * s(i);  // v_r
      out(i, 1) = -vx * s(i) + vy * c(i); // v_theta
      out(i, 2) = in(i, 2);               // v_z
    });
}

void RotateCylindricalToCartesian(
  const Kokkos::View<Real*, DeviceMemorySpace>& c,
  const Kokkos::View<Real*, DeviceMemorySpace>& s,
  Rank2View<const Real, DeviceMemorySpace> in,
  Rank2View<Real, DeviceMemorySpace> out, LO n)
{
  Kokkos::parallel_for(
    "cylindrical_to_cartesian_basis", Kokkos::RangePolicy<ExecutionSpace>(0, n),
    KOKKOS_LAMBDA(const LO i) {
      const Real vr = in(i, 0);
      const Real vt = in(i, 1);
      out(i, 0) = vr * c(i) - vt * s(i); // v_x
      out(i, 1) = vr * s(i) + vt * c(i); // v_y
      out(i, 2) = in(i, 2);              // v_z
    });
}

} // namespace

BoundBasisTransformation::BoundBasisTransformation(ValueBasis source_basis,
                                                   ValueBasis target_basis,
                                                   LO num_points)
  : source_basis_(std::move(source_basis)),
    target_basis_(std::move(target_basis)),
    num_points_(num_points)
{
  detail::ValidateValueBasis(source_basis_);
  detail::ValidateValueBasis(target_basis_);
  if (source_basis_.Rank() != 1 || target_basis_.Rank() != 1) {
    throw pcms_error("BoundBasisTransformation: both bases must be rank 1");
  }
}

void BoundBasisTransformation::Apply(
  ValueView<const Real, DeviceMemorySpace> in,
  ValueView<Real, DeviceMemorySpace> out) const
{
  PCMS_FUNCTION_TIMER;
  if (in.GetValueType() != out.GetValueType()) {
    throw pcms_error(
      std::string("BoundBasisTransformation::Apply: input value type (") +
      ValueTypeName(in.GetValueType()) + ") differs from output (" +
      ValueTypeName(out.GetValueType()) + ")");
  }
  const auto in_values = in.GetValues();
  const auto out_values = out.GetValues();
  if (static_cast<LO>(in_values.extent(0)) != num_points_ ||
      static_cast<LO>(out_values.extent(0)) != num_points_) {
    throw pcms_error("BoundBasisTransformation::Apply: value counts do not "
                     "match the bound point set");
  }
  if (in_values.extent(1) != out_values.extent(1)) {
    throw pcms_error(
      "BoundBasisTransformation::Apply: input and output widths differ");
  }
  if (in.GetValueType() == FieldValueType::Scalar) {
    // Scalars are basis-invariant: deep copy, any bases accepted.
    CopyDeviceRank2ViewToRank2View(out_values, in_values);
    return;
  }
  if (in.GetValueType() == FieldValueType::Tensor) {
    throw pcms_error("BoundBasisTransformation::Apply: tensor basis "
                     "transformations are not yet supported");
  }
  if (!SameValueBasis(in.GetBasis(), source_basis_)) {
    throw pcms_error("BoundBasisTransformation::Apply: input basis does not "
                     "match this transformation's source basis");
  }
  if (!SameValueBasis(out.GetBasis(), target_basis_)) {
    throw pcms_error("BoundBasisTransformation::Apply: output basis does not "
                     "match this transformation's target basis");
  }
  ApplyVectorImpl(in_values, out_values);
}

const ValueBasis& basis::CartesianToCylindrical::Source() const noexcept
{
  return CartesianPhysicalVector();
}

const ValueBasis& basis::CartesianToCylindrical::Target() const noexcept
{
  return CylindricalPhysicalVector();
}

basis::CartesianToCylindrical::State basis::CartesianToCylindrical::Prepare(
  const CoordinateView<DeviceMemorySpace>& points, PointStatusView&) const
{
  PCMS_FUNCTION_TIMER;
  return ComputeCosSinTheta("CartesianToCylindricalBasis", points);
}

void basis::CartesianToCylindrical::ApplyVector(
  const State& state, Rank2View<const Real, DeviceMemorySpace> in,
  Rank2View<Real, DeviceMemorySpace> out, LO n) const
{
  RotateCartesianToCylindrical(state.cos_theta, state.sin_theta, in, out, n);
}

const ValueBasis& basis::CylindricalToCartesian::Source() const noexcept
{
  return CylindricalPhysicalVector();
}

const ValueBasis& basis::CylindricalToCartesian::Target() const noexcept
{
  return CartesianPhysicalVector();
}

basis::CylindricalToCartesian::State basis::CylindricalToCartesian::Prepare(
  const CoordinateView<DeviceMemorySpace>& points, PointStatusView&) const
{
  PCMS_FUNCTION_TIMER;
  return ComputeCosSinTheta("CylindricalToCartesianBasis", points);
}

void basis::CylindricalToCartesian::ApplyVector(
  const State& state, Rank2View<const Real, DeviceMemorySpace> in,
  Rank2View<Real, DeviceMemorySpace> out, LO n) const
{
  RotateCylindricalToCartesian(state.cos_theta, state.sin_theta, in, out, n);
}

template class BoundBasis<basis::CartesianToCylindrical>;
template class UnboundBasis<basis::CartesianToCylindrical>;
template class BoundBasis<basis::CylindricalToCartesian>;
template class UnboundBasis<basis::CylindricalToCartesian>;

} // namespace pcms
