#include "pcms/field/field_aligned/field_line_tracer.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/profile.h"
#include <Kokkos_Core.hpp>
#include <cmath>
#include <string>

namespace pcms
{

namespace
{

using ExecutionSpace = DeviceMemorySpace::execution_space;
using Reals = Kokkos::View<Real*, DeviceMemorySpace>;
using Flags = Kokkos::View<bool*, DeviceMemorySpace>;
using Positions = Kokkos::View<Real***, Kokkos::LayoutRight, DeviceMemorySpace>;

Reals Allocate(const std::string& label, size_t n)
{
  return Reals(Kokkos::view_alloc(Kokkos::WithoutInitializing, label), n);
}

/// Output buffers of one batched AxisymmetricMagneticField::Evaluate call.
struct FieldSamples
{
  explicit FieldSamples(size_t n)
    : psi(Allocate("tracer_psi", n)),
      psi_R(Allocate("tracer_psi_R", n)),
      psi_Z(Allocate("tracer_psi_Z", n)),
      I(Allocate("tracer_I", n)),
      inside(Kokkos::view_alloc(Kokkos::WithoutInitializing, "tracer_inside"),
             n)
  {
  }

  void Evaluate(const AxisymmetricMagneticField& field, const Reals& R,
                const Reals& Z) const
  {
    field.Evaluate(R, Z, psi, psi_R, psi_Z, I, inside);
  }

  Reals psi;
  Reals psi_R;
  Reals psi_Z;
  Reals I;
  Flags inside;
};

/// Axis-aligned R-Z box of the equilibrium grid, in m.
struct Box
{
  Real r_min, r_max, z_min, z_max;

  explicit Box(const AxisymmetricMagneticField& field)
    : r_min(field.GetRMin()),
      r_max(field.GetRMax()),
      z_min(field.GetZMin()),
      z_max(field.GetZMax())
  {
  }

  KOKKOS_INLINE_FUNCTION bool Contains(Real R, Real Z) const
  {
    return R >= r_min && R <= r_max && Z >= z_min && Z <= z_max;
  }
};

// Kernels as free functions: CUDA extended lambdas may not be defined inside
// member functions of classes with private members.

/// @brief Field-line slope (dR/dphi, dZ/dphi) at the stage points (R, Z)
///
/// @details R in the slope is the stage's own R, which keeps the scheme
/// fourth-order accurate. ok(i) is cleared where the stage left the box or
/// I(psi) is zero, and that point's slope is set to zero.
void FieldLineSlope(const FieldSamples& f, Real s, const Reals& R,
                    const Reals& kR, const Reals& kZ, const Flags& ok)
{
  const auto psi_R = f.psi_R;
  const auto psi_Z = f.psi_Z;
  const auto I = f.I;
  const auto inside = f.inside;
  Kokkos::parallel_for(
    "tracer_slope", Kokkos::RangePolicy<ExecutionSpace>(0, R.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      if (ok(i) && inside(i) && I(i) != 0.0) {
        kR(i) = s * R(i) * psi_Z(i) / I(i);
        kZ(i) = -s * R(i) * psi_R(i) / I(i);
      } else {
        ok(i) = false;
        kR(i) = 0.0;
        kZ(i) = 0.0;
      }
    });
}

/// @brief Stage point (Rs, Zs) = (R, Z) + a * (kR, kZ)
void StagePoint(const Reals& R, const Reals& Z, const Reals& kR,
                const Reals& kZ, Real a, const Reals& Rs, const Reals& Zs)
{
  Kokkos::parallel_for(
    "tracer_stage", Kokkos::RangePolicy<ExecutionSpace>(0, R.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      Rs(i) = R(i) + a * kR(i);
      Zs(i) = Z(i) + a * kZ(i);
    });
}

/// @brief Combines the four RK4 slopes and advances the valid points
///
/// @details A point whose stages all succeeded and whose new position is in the
/// box moves; any other valid point is marked invalid and stays where it was.
void RK4Update(Real h, const Box& box, const Reals& k1R, const Reals& k1Z,
               const Reals& k2R, const Reals& k2Z, const Reals& k3R,
               const Reals& k3Z, const Reals& k4R, const Reals& k4Z,
               const Flags& ok, const Reals& R, const Reals& Z,
               const Flags& valid)
{
  Kokkos::parallel_for(
    "tracer_rk4_update", Kokkos::RangePolicy<ExecutionSpace>(0, R.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      if (!valid(i)) {
        return;
      }
      const Real R_new =
        R(i) + h / 6.0 * (k1R(i) + 2.0 * k2R(i) + 2.0 * k3R(i) + k4R(i));
      const Real Z_new =
        Z(i) + h / 6.0 * (k1Z(i) + 2.0 * k2Z(i) + 2.0 * k3Z(i) + k4Z(i));
      if (ok(i) && box.Contains(R_new, Z_new)) {
        R(i) = R_new;
        Z(i) = Z_new;
      } else {
        valid(i) = false;
      }
    });
}

/// @brief One Newton step x -= (psi - psi0) * grad(psi) / |grad(psi)|^2
///
/// @details A point the previous step moved out of the box is restored to
/// (R_keep, Z_keep) and marked invalid.
void NewtonPsiStep(const FieldSamples& f, const Reals& psi0,
                   const Reals& R_keep, const Reals& Z_keep, const Reals& R,
                   const Reals& Z, const Flags& valid)
{
  const auto psi = f.psi;
  const auto psi_R = f.psi_R;
  const auto psi_Z = f.psi_Z;
  const auto inside = f.inside;
  Kokkos::parallel_for(
    "tracer_newton_psi", Kokkos::RangePolicy<ExecutionSpace>(0, R.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      if (!valid(i)) {
        return;
      }
      if (!inside(i)) {
        R(i) = R_keep(i);
        Z(i) = Z_keep(i);
        valid(i) = false;
        return;
      }
      const Real grad2 = psi_R(i) * psi_R(i) + psi_Z(i) * psi_Z(i);
      if (grad2 > 0.0) {
        const Real step = (psi(i) - psi0(i)) / grad2;
        R(i) -= step * psi_R(i);
        Z(i) -= step * psi_Z(i);
      }
    });
}

/// @brief Restores and invalidates valid points the projection left outside
///        the box
void RejectOutsideBox(const Box& box, const Reals& R_keep, const Reals& Z_keep,
                      const Reals& R, const Reals& Z, const Flags& valid)
{
  Kokkos::parallel_for(
    "tracer_reject_outside",
    Kokkos::RangePolicy<ExecutionSpace>(0, R.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      if (valid(i) && !box.Contains(R(i), Z(i))) {
        R(i) = R_keep(i);
        Z(i) = Z_keep(i);
        valid(i) = false;
      }
    });
}

void RecordStation(const Positions& positions, size_t k, const Reals& R,
                   const Reals& Z)
{
  Kokkos::parallel_for(
    "tracer_record_station",
    Kokkos::RangePolicy<ExecutionSpace>(0, R.extent(0)),
    KOKKOS_LAMBDA(const size_t i) {
      positions(k, i, 0) = R(i);
      positions(k, i, 1) = Z(i);
    });
}

void RK4Step(const AxisymmetricMagneticField& field, Real h, FieldSamples& f,
             const Reals& R, const Reals& Z, const Flags& valid)
{
  const size_t n = R.extent(0);
  const Real s = field.GetPoloidalSign();
  const Reals k1R = Allocate("tracer_k1R", n), k1Z = Allocate("tracer_k1Z", n);
  const Reals k2R = Allocate("tracer_k2R", n), k2Z = Allocate("tracer_k2Z", n);
  const Reals k3R = Allocate("tracer_k3R", n), k3Z = Allocate("tracer_k3Z", n);
  const Reals k4R = Allocate("tracer_k4R", n), k4Z = Allocate("tracer_k4Z", n);
  const Reals Rs = Allocate("tracer_Rs", n), Zs = Allocate("tracer_Zs", n);
  Flags ok(Kokkos::view_alloc(Kokkos::WithoutInitializing, "tracer_ok"), n);
  Kokkos::deep_copy(ok, valid);

  f.Evaluate(field, R, Z);
  FieldLineSlope(f, s, R, k1R, k1Z, ok);
  StagePoint(R, Z, k1R, k1Z, 0.5 * h, Rs, Zs);
  f.Evaluate(field, Rs, Zs);
  FieldLineSlope(f, s, Rs, k2R, k2Z, ok);
  StagePoint(R, Z, k2R, k2Z, 0.5 * h, Rs, Zs);
  f.Evaluate(field, Rs, Zs);
  FieldLineSlope(f, s, Rs, k3R, k3Z, ok);
  StagePoint(R, Z, k3R, k3Z, h, Rs, Zs);
  f.Evaluate(field, Rs, Zs);
  FieldLineSlope(f, s, Rs, k4R, k4Z, ok);
  RK4Update(h, Box(field), k1R, k1Z, k2R, k2Z, k3R, k3Z, k4R, k4Z, ok, R, Z,
            valid);
}

void ProjectToPsi(const AxisymmetricMagneticField& field, FieldSamples& f,
                  const Reals& psi0, const Reals& R, const Reals& Z,
                  const Flags& valid)
{
  const size_t n = R.extent(0);
  const Reals R_keep = Allocate("tracer_R_keep", n);
  const Reals Z_keep = Allocate("tracer_Z_keep", n);
  Kokkos::deep_copy(R_keep, R);
  Kokkos::deep_copy(Z_keep, Z);
  constexpr int num_newton_steps = 3;
  for (int step = 0; step < num_newton_steps; ++step) {
    f.Evaluate(field, R, Z);
    NewtonPsiStep(f, psi0, R_keep, Z_keep, R, Z, valid);
  }
  RejectOutsideBox(Box(field), R_keep, Z_keep, R, Z, valid);
}

} // namespace

TracedPoints TraceFieldLines(const AxisymmetricMagneticField& field,
                             Kokkos::View<const Real*, DeviceMemorySpace> R0,
                             Kokkos::View<const Real*, DeviceMemorySpace> Z0,
                             const std::vector<Real>& stations, Real max_dphi,
                             bool project_to_psi)
{
  PCMS_FUNCTION_TIMER;
  const size_t n = R0.extent(0);
  if (Z0.extent(0) != n) {
    throw pcms_error("TraceFieldLines: R0 and Z0 must have the same size");
  }
  if (stations.empty()) {
    throw pcms_error("TraceFieldLines: stations must not be empty");
  }
  if (!(max_dphi > 0.0)) {
    throw pcms_error("TraceFieldLines: max_dphi must be positive");
  }
  const size_t num_stations = stations.size();
  const bool increasing = num_stations < 2 || stations[1] > stations[0];
  for (size_t k = 1; k < num_stations; ++k) {
    const Real delta = stations[k] - stations[k - 1];
    if (!(increasing ? delta > 0.0 : delta < 0.0)) {
      throw pcms_error("TraceFieldLines: stations must be strictly increasing "
                       "or strictly decreasing");
    }
  }

  TracedPoints traced;
  traced.stations = Reals("traced_stations", num_stations);
  auto stations_host = Kokkos::create_mirror_view(traced.stations);
  for (size_t k = 0; k < num_stations; ++k) {
    stations_host(k) = stations[k];
  }
  Kokkos::deep_copy(traced.stations, stations_host);
  traced.positions = Positions("traced_positions", num_stations, n, 2);
  traced.valid = Flags("traced_valid", n);
  if (n == 0) {
    return traced;
  }

  const Reals R = Allocate("tracer_R", n);
  const Reals Z = Allocate("tracer_Z", n);
  Kokkos::deep_copy(R, R0);
  Kokkos::deep_copy(Z, Z0);
  FieldSamples f(n);
  f.Evaluate(field, R, Z);
  const Reals psi0 = Allocate("tracer_psi0", n);
  Kokkos::deep_copy(psi0, f.psi);
  Kokkos::deep_copy(traced.valid, f.inside);
  RecordStation(traced.positions, 0, R, Z);

  for (size_t k = 1; k < num_stations; ++k) {
    const Real delta = stations[k] - stations[k - 1];
    const auto num_steps =
      static_cast<long>(std::ceil(std::abs(delta) / max_dphi));
    // The field is axisymmetric, so only the step size matters and equal
    // substeps land exactly on the next station.
    const Real h = delta / static_cast<Real>(num_steps);
    for (long step = 0; step < num_steps; ++step) {
      RK4Step(field, h, f, R, Z, traced.valid);
    }
    if (project_to_psi) {
      ProjectToPsi(field, f, psi0, R, Z, traced.valid);
    }
    RecordStation(traced.positions, k, R, Z);
  }
  return traced;
}

TracedPoints TraceFieldLines(
  const AxisymmetricMagneticField& field,
  Kokkos::View<const Real*, DeviceMemorySpace> R0,
  Kokkos::View<const Real*, DeviceMemorySpace> Z0,
  Kokkos::View<const Real*, DeviceMemorySpace> stations, Real max_dphi,
  bool project_to_psi)
{
  auto stations_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, stations);
  const std::vector<Real> host(stations_host.data(),
                               stations_host.data() + stations_host.extent(0));
  return TraceFieldLines(field, R0, Z0, host, max_dphi, project_to_psi);
}

std::vector<Real> WedgeStations(Real phi0, Real dphi, int num_layers)
{
  if (num_layers < 0) {
    throw pcms_error("WedgeStations: num_layers must be >= 0, got " +
                     std::to_string(num_layers));
  }
  std::vector<Real> stations;
  stations.reserve(3 * static_cast<size_t>(num_layers) + 1);
  for (int k = 0; k < num_layers; ++k) {
    const Real start = phi0 + k * dphi;
    stations.push_back(start);
    stations.push_back(start + dphi / 3.0);
    stations.push_back(start + 2.0 * dphi / 3.0);
  }
  stations.push_back(phi0 + num_layers * dphi);
  return stations;
}

} // namespace pcms
