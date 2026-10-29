#include "pcms/field/field_aligned/magnetic_field.hpp"
#include "pcms/field/evaluator/spline_interpolator.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/profile.h"
#include <Kokkos_Core.hpp>
#include <cmath>
#include <initializer_list>
#include <string>
#include <vector>

namespace pcms
{

namespace
{

using ExecutionSpace = DeviceMemorySpace::execution_space;
using BiCubicSpline = CompactBiCubicSplineInterpolator<Real, DeviceMemorySpace>;
using CubicSpline = CompactCubicSplineInterpolator<Real, DeviceMemorySpace>;

Rank1View<Real, DeviceMemorySpace> ToRank1View(
  const Kokkos::View<Real*, DeviceMemorySpace>& v)
{
  return Rank1View<Real, DeviceMemorySpace>(v.data(), v.extent(0));
}

Kokkos::View<Real*, DeviceMemorySpace> CopyToDevice(
  const std::string& label, const std::vector<Real>& host)
{
  Kokkos::View<Real*, DeviceMemorySpace> device(label, host.size());
  auto mirror = Kokkos::create_mirror_view(device);
  for (size_t i = 0; i < host.size(); ++i) {
    mirror(i) = host[i];
  }
  Kokkos::deep_copy(device, mirror);
  return device;
}

/// @brief Spline derivative selector: 1 marks an entry to compute
Kokkos::View<LO*, DeviceMemorySpace> MakeSelector(
  const std::string& label, std::initializer_list<LO> flags)
{
  Kokkos::View<LO*, DeviceMemorySpace> device(label, flags.size());
  auto mirror = Kokkos::create_mirror_view(device);
  LO i = 0;
  for (const LO flag : flags) {
    mirror(i++) = flag;
  }
  Kokkos::deep_copy(device, mirror);
  return device;
}

// Kernels as free functions: CUDA extended lambdas may not be defined inside
// member functions of classes with private members.

/// @brief Grid node coordinates start + delta * i, i = 0, ..., n - 1
Kokkos::View<Real*, DeviceMemorySpace> UniformCoordinates(
  const std::string& label, LO n, Real start, Real delta)
{
  Kokkos::View<Real*, DeviceMemorySpace> coords(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, label), n);
  Kokkos::parallel_for(
    label, Kokkos::RangePolicy<ExecutionSpace>(0, n),
    KOKKOS_LAMBDA(const LO i) { coords(i) = start + delta * i; });
  return coords;
}

/// @brief Reorders PSIZR (R fastest: i_R + i_Z * nR) into the spline's
///        layout (i_R * nZ + i_Z)
Kokkos::View<Real*, DeviceMemorySpace> TransposePsi(
  const Kokkos::View<Real*, DeviceMemorySpace>& psizr, LO nR, LO nZ)
{
  Kokkos::View<Real*, DeviceMemorySpace> values(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "magnetic_field_psi"),
    static_cast<size_t>(nR) * nZ);
  Kokkos::parallel_for(
    "magnetic_field_transpose_psi",
    Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<2>>({0, 0}, {nZ, nR}),
    KOKKOS_LAMBDA(const LO iz, const LO ir) {
      values(ir * nZ + iz) = psizr(ir + iz * nR);
    });
  return values;
}

/// @brief Flags points inside the R-Z box and moves the others onto a grid node
///
/// @details The spline prints from the kernel and leaves its output unset for
/// points off its grid, so it is only ever given points inside the box.
void LocatePoints(const Kokkos::View<const Real*, DeviceMemorySpace>& R,
                  const Kokkos::View<const Real*, DeviceMemorySpace>& Z,
                  Real r_min, Real r_max, Real z_min, Real z_max,
                  const Kokkos::View<bool*, DeviceMemorySpace>& inside,
                  const Kokkos::View<Real*, DeviceMemorySpace>& safe_R,
                  const Kokkos::View<Real*, DeviceMemorySpace>& safe_Z)
{
  Kokkos::parallel_for(
    "magnetic_field_locate",
    Kokkos::RangePolicy<ExecutionSpace>(0, R.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      const bool in =
        R(i) >= r_min && R(i) <= r_max && Z(i) >= z_min && Z(i) <= z_max;
      inside(i) = in;
      safe_R(i) = in ? R(i) : r_min;
      safe_Z(i) = in ? Z(i) : z_min;
    });
}

/// @brief psi limited to [lo, hi], the psi_grid range of the I(psi) spline
void ClampPsi(
  const Kokkos::View<Real**, Kokkos::LayoutRight, DeviceMemorySpace>& psi_out,
  Real lo, Real hi, const Kokkos::View<Real*, DeviceMemorySpace>& clamped)
{
  Kokkos::parallel_for(
    "magnetic_field_clamp_psi",
    Kokkos::RangePolicy<ExecutionSpace>(0, clamped.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      clamped(i) = Kokkos::min(Kokkos::max(psi_out(i, 0), lo), hi);
    });
}

void WriteOutputs(
  const Kokkos::View<bool*, DeviceMemorySpace>& inside,
  const Kokkos::View<Real**, Kokkos::LayoutRight, DeviceMemorySpace>& psi_out,
  const Kokkos::View<Real**, Kokkos::LayoutRight, DeviceMemorySpace>& I_out,
  const Kokkos::View<Real*, DeviceMemorySpace>& psi,
  const Kokkos::View<Real*, DeviceMemorySpace>& psi_R,
  const Kokkos::View<Real*, DeviceMemorySpace>& psi_Z,
  const Kokkos::View<Real*, DeviceMemorySpace>& I)
{
  Kokkos::parallel_for(
    "magnetic_field_write",
    Kokkos::RangePolicy<ExecutionSpace>(0, psi.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      const bool in = inside(i);
      psi(i) = in ? psi_out(i, 0) : 0.0;
      psi_R(i) = in ? psi_out(i, 1) : 0.0;
      psi_Z(i) = in ? psi_out(i, 2) : 0.0;
      I(i) = in ? I_out(i, 0) : 0.0;
    });
}

} // namespace

struct AxisymmetricMagneticField::Impl
{
  // The interpolators hold non-owning views of their grids.
  Kokkos::View<Real*, DeviceMemorySpace> r_coords;
  Kokkos::View<Real*, DeviceMemorySpace> z_coords;
  Kokkos::View<Real*, DeviceMemorySpace> psi_grid;
  Kokkos::View<LO*, DeviceMemorySpace> psi_select; // psi, psi_R, psi_Z
  Kokkos::View<LO*, DeviceMemorySpace> I_select;   // I
  std::unique_ptr<BiCubicSpline> psi_spline;
  std::unique_ptr<CubicSpline> I_spline;
};

AxisymmetricMagneticField::AxisymmetricMagneticField(const EQDSKData& eqdsk,
                                                     int poloidal_sign)
  : impl_(std::make_shared<Impl>()), poloidal_sign_(poloidal_sign)
{
  PCMS_FUNCTION_TIMER;
  if (poloidal_sign != 1 && poloidal_sign != -1) {
    throw pcms_error(
      "AxisymmetricMagneticField: poloidal_sign must be +1 or -1, got " +
      std::to_string(poloidal_sign));
  }

  const LO nR = eqdsk.GetNW();
  const LO nZ = eqdsk.GetNH();
  const Real dR = eqdsk.GetDeltaR();
  const Real dZ = eqdsk.GetDeltaZ();
  if (nR < 2 || nZ < 2 || !(dR > 0.0) || !(dZ > 0.0)) {
    throw pcms_error("AxisymmetricMagneticField: the R-Z grid needs at least "
                     "2 points and a positive spacing in each direction");
  }
  if (eqdsk.PSIZR.extent(0) != static_cast<size_t>(nR) * nZ) {
    throw pcms_error(
      "AxisymmetricMagneticField: PSIZR size does not match the R-Z grid");
  }
  impl_->r_coords =
    UniformCoordinates("magnetic_field_R", nR, eqdsk.GetRMin(), dR);
  impl_->z_coords =
    UniformCoordinates("magnetic_field_Z", nZ, eqdsk.GetZMin(), dZ);
  r_min_ = eqdsk.GetRMin();
  z_min_ = eqdsk.GetZMin();
  // Read back the last node so the box matches the spline grid bit for bit;
  // the device may fuse start + delta * i differently from the host.
  Kokkos::deep_copy(r_max_, Kokkos::subview(impl_->r_coords, nR - 1));
  Kokkos::deep_copy(z_max_, Kokkos::subview(impl_->z_coords, nZ - 1));

  auto psi_values = TransposePsi(eqdsk.PSIZR, nR, nZ);
  Kokkos::View<Real*, DeviceMemorySpace> no_bc("magnetic_field_no_bc", 0);
  impl_->psi_spline = std::make_unique<BiCubicSpline>(
    ToRank1View(impl_->r_coords), ToRank1View(impl_->z_coords),
    ToRank1View(psi_values), BoundaryCondition::NOT_A_KNOT, ToRank1View(no_bc),
    BoundaryCondition::NOT_A_KNOT, ToRank1View(no_bc),
    BoundaryCondition::NOT_A_KNOT, ToRank1View(no_bc),
    BoundaryCondition::NOT_A_KNOT, ToRank1View(no_bc));

  const LO npsi = eqdsk.GetNPsi();
  if (npsi < 2 || eqdsk.I_psi.extent(0) != static_cast<size_t>(npsi)) {
    throw pcms_error("AxisymmetricMagneticField: psi_grid needs at least 2 "
                     "points and the same size as I_psi");
  }
  auto psi_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, eqdsk.psi_grid);
  auto I_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, eqdsk.I_psi);
  // The spline needs an ascending grid; G-EQDSK runs from the axis to the
  // boundary, which is descending when psi decreases outward.
  const bool descending = psi_host(npsi - 1) < psi_host(0);
  std::vector<Real> psi_sorted(npsi);
  std::vector<Real> I_sorted(npsi);
  for (LO i = 0; i < npsi; ++i) {
    const LO src = descending ? npsi - 1 - i : i;
    psi_sorted[i] = psi_host(src);
    I_sorted[i] = I_host(src);
  }
  // Same uniformity test as the spline (relative tolerance 1e-3), done here so
  // a bad grid gets a clear error instead of the spline's internal assert.
  const Real spacing = (psi_sorted[npsi - 1] - psi_sorted[0]) / (npsi - 1);
  bool uniform = spacing > 0.0;
  for (LO i = 1; uniform && i < npsi; ++i) {
    const Real diff = psi_sorted[i] - psi_sorted[i - 1];
    uniform = diff > 0.0 && std::abs(diff - spacing) <= 1.0e-3 * spacing;
  }
  if (!uniform) {
    throw pcms_error("AxisymmetricMagneticField: psi_grid must be uniform and "
                     "strictly ascending or strictly descending");
  }
  psi_min_ = psi_sorted[0];
  psi_max_ = psi_sorted[npsi - 1];
  impl_->psi_grid = CopyToDevice("magnetic_field_psi_grid", psi_sorted);
  auto I_values = CopyToDevice("magnetic_field_I", I_sorted);
  impl_->I_spline = std::make_unique<CubicSpline>(
    ToRank1View(impl_->psi_grid), ToRank1View(I_values),
    BoundaryCondition::NOT_A_KNOT, 0.0, BoundaryCondition::NOT_A_KNOT, 0.0);

  impl_->psi_select =
    MakeSelector("magnetic_field_psi_select", {1, 1, 1, 0, 0, 0});
  impl_->I_select = MakeSelector("magnetic_field_I_select", {1, 0, 0});
}

void AxisymmetricMagneticField::Evaluate(
  Kokkos::View<const Real*, DeviceMemorySpace> R,
  Kokkos::View<const Real*, DeviceMemorySpace> Z,
  Kokkos::View<Real*, DeviceMemorySpace> psi,
  Kokkos::View<Real*, DeviceMemorySpace> psi_R,
  Kokkos::View<Real*, DeviceMemorySpace> psi_Z,
  Kokkos::View<Real*, DeviceMemorySpace> I,
  Kokkos::View<bool*, DeviceMemorySpace> inside) const
{
  PCMS_FUNCTION_TIMER;
  const size_t n = R.extent(0);
  if (Z.extent(0) != n || psi.extent(0) != n || psi_R.extent(0) != n ||
      psi_Z.extent(0) != n || I.extent(0) != n || inside.extent(0) != n) {
    throw pcms_error("AxisymmetricMagneticField::Evaluate: R, Z, psi, psi_R, "
                     "psi_Z, I and inside must have the same size");
  }
  if (n == 0) {
    return;
  }

  Kokkos::View<Real*, DeviceMemorySpace> safe_R(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "magnetic_field_safe_R"),
    n);
  Kokkos::View<Real*, DeviceMemorySpace> safe_Z(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "magnetic_field_safe_Z"),
    n);
  LocatePoints(R, Z, r_min_, r_max_, z_min_, z_max_, inside, safe_R, safe_Z);

  Kokkos::View<Real**, Kokkos::LayoutRight, DeviceMemorySpace> psi_out(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "magnetic_field_psi_out"),
    n, 3);
  impl_->psi_spline->evaluate(
    impl_->psi_select, ToRank1View(safe_R), ToRank1View(safe_Z),
    Rank2View<Real, DeviceMemorySpace>(psi_out.data(), n, 3));

  Kokkos::View<Real*, DeviceMemorySpace> psi_clamped(
    Kokkos::view_alloc(Kokkos::WithoutInitializing,
                       "magnetic_field_psi_clamped"),
    n);
  ClampPsi(psi_out, psi_min_, psi_max_, psi_clamped);
  Kokkos::View<Real**, Kokkos::LayoutRight, DeviceMemorySpace> I_out(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "magnetic_field_I_out"), n,
    1);
  impl_->I_spline->evaluate(
    impl_->I_select, ToRank1View(psi_clamped),
    Rank2View<Real, DeviceMemorySpace>(I_out.data(), n, 1));

  WriteOutputs(inside, psi_out, I_out, psi, psi_R, psi_Z, I);
}

} // namespace pcms
