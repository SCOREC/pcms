#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <pcms/field/field_aligned/field_line_tracer.hpp>
#include <pcms/field/field_aligned/magnetic_field.hpp>
#include <pcms/utility/assert.h>
#include <pcms/utility/eqdsk.h>
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <unistd.h>
#include <vector>

using Catch::Approx;
using pcms::Real;

namespace
{

// psi = sign * c * ((R - R0)^2 + Z^2) has circular flux surfaces about (R0, 0).
// The not-a-knot splines reproduce a quadratic psi and a linear I(psi)
// exactly, so the field below is exact up to rounding.
constexpr Real kAxisR = 2.0;
constexpr Real kC = 0.0625;
constexpr Real kI = 0.5;
constexpr int kCells = 64; // 65 x 65 points on [1, 3] x [-1, 1]
constexpr Real kPi = Kokkos::numbers::pi_v<Real>;

Real CircularPsi(Real R, Real Z)
{
  return kC * ((R - kAxisR) * (R - kAxisR) + Z * Z);
}

// Field line through (R0 + r, 0) at phi = 0 for psi = c*((R - R0)^2 + Z^2),
// with I the value of I(psi) on that line (psi, hence I, is constant along it):
// d(theta)/d(phi) = -s * (2c / I) * (R0 + r*cos(theta)) integrates
// to tan(theta / 2) = sqrt((R0 + r) / (R0 - r)) * tan(x).
struct RZPoint
{
  Real R;
  Real Z;
};

RZPoint AnalyticFieldLine(Real r, Real phi, int s, Real I = kI)
{
  const Real k = std::sqrt(kAxisR * kAxisR - r * r);
  const Real x = 0.5 * k * (-s * (2.0 * kC / I) * phi);
  // atan2 instead of tan: theta passes through pi, where tan(theta / 2) is
  // singular. theta is only known modulo 2*pi, which is enough for R and Z.
  const Real theta =
    2.0 * std::atan2(std::sqrt((kAxisR + r) / (kAxisR - r)) * std::sin(x),
                     std::cos(x));
  return {kAxisR + r * std::cos(theta), r * std::sin(theta)};
}

std::vector<Real> UniformValues(Real first, Real last, int n)
{
  std::vector<Real> values(n);
  for (int i = 0; i < n; ++i) {
    values[i] = first + (last - first) * i / (n - 1);
  }
  return values;
}

pcms::EQDSKData MakeCircularEquilibrium(Real psi_sign,
                                        const std::vector<Real>& psi_grid,
                                        const std::vector<Real>& I_psi)
{
  pcms::Uniform2DGrid grid;
  grid.divisions[0] = kCells;
  grid.divisions[1] = kCells;
  grid.bot_left[0] = 1.0;
  grid.bot_left[1] = -1.0;
  grid.edge_length[0] = 2.0;
  grid.edge_length[1] = 2.0;
  const int n = kCells + 1;
  std::vector<Real> psirz(static_cast<size_t>(n) * n);
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      const Real R = 1.0 + 2.0 * i / kCells;
      const Real Z = -1.0 + 2.0 * j / kCells;
      psirz[i + j * n] = psi_sign * CircularPsi(R, Z);
    }
  }
  return pcms::EQDSKData(grid, psirz, psi_grid, I_psi);
}

// psi_grid covers [0, 0.13], beyond the largest psi on the grid (0.125).
pcms::EQDSKData MakeCircularEquilibrium()
{
  return MakeCircularEquilibrium(1.0, UniformValues(0.0, 0.13, 27),
                                 std::vector<Real>(27, kI));
}

// Same psi on a box that is not square, not centred on the axis and has
// different R and Z spacings, so a PSIZR read with R and Z swapped shows up.
constexpr int kOffsetCellsR = 35; // R in [1.25, 3.0], dR = 0.05
constexpr int kOffsetCellsZ = 50; // Z in [-0.9, 1.1], dZ = 0.04

pcms::EQDSKData MakeOffsetEquilibrium()
{
  pcms::Uniform2DGrid grid;
  grid.divisions[0] = kOffsetCellsR;
  grid.divisions[1] = kOffsetCellsZ;
  grid.bot_left[0] = 1.25;
  grid.bot_left[1] = -0.9;
  grid.edge_length[0] = 1.75;
  grid.edge_length[1] = 2.0;
  const int nR = kOffsetCellsR + 1;
  const int nZ = kOffsetCellsZ + 1;
  std::vector<Real> psirz(static_cast<size_t>(nR) * nZ);
  for (int j = 0; j < nZ; ++j) {
    for (int i = 0; i < nR; ++i) {
      const Real R = 1.25 + 1.75 * i / kOffsetCellsR;
      const Real Z = -0.9 + 2.0 * j / kOffsetCellsZ;
      psirz[i + j * nR] = CircularPsi(R, Z);
    }
  }
  // psi reaches 0.138 at (3.0, 1.1).
  return pcms::EQDSKData(grid, psirz, UniformValues(0.0, 0.14, 29),
                         std::vector<Real>(29, kI));
}

Kokkos::View<Real*, pcms::DeviceMemorySpace> ToDevice(
  const std::vector<Real>& values)
{
  Kokkos::View<Real*, pcms::DeviceMemorySpace> device("test_values",
                                                      values.size());
  auto host = Kokkos::create_mirror_view(device);
  for (size_t i = 0; i < values.size(); ++i) {
    host(i) = values[i];
  }
  Kokkos::deep_copy(device, host);
  return device;
}

struct HostFieldValues
{
  Kokkos::View<Real*, Kokkos::HostSpace> psi, psi_R, psi_Z, I;
  Kokkos::View<bool*, Kokkos::HostSpace> inside;
};

HostFieldValues EvaluateField(const pcms::AxisymmetricMagneticField& field,
                              const std::vector<Real>& R,
                              const std::vector<Real>& Z)
{
  const size_t n = R.size();
  Kokkos::View<Real*, pcms::DeviceMemorySpace> psi("psi", n), psi_R("psi_R", n),
    psi_Z("psi_Z", n), I("I", n);
  Kokkos::View<bool*, pcms::DeviceMemorySpace> inside("inside", n);
  field.Evaluate(ToDevice(R), ToDevice(Z), psi, psi_R, psi_Z, I, inside);
  return {Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, psi),
          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, psi_R),
          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, psi_Z),
          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, I),
          Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, inside)};
}

struct HostTrace
{
  Kokkos::View<Real***, Kokkos::LayoutRight, Kokkos::HostSpace> positions;
  Kokkos::View<bool*, Kokkos::HostSpace> valid;
};

HostTrace Trace(const pcms::AxisymmetricMagneticField& field,
                const std::vector<Real>& R0, const std::vector<Real>& Z0,
                const std::vector<Real>& stations, Real max_dphi,
                bool project_to_psi = false)
{
  const auto traced = pcms::TraceFieldLines(field, ToDevice(R0), ToDevice(Z0),
                                            stations, max_dphi, project_to_psi);
  return {
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, traced.positions),
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, traced.valid)};
}

// Largest distance, over all stations and start points, between the traced
// and the analytic field lines started at (R0 + r, 0).
Real MaxAnalyticError(const pcms::AxisymmetricMagneticField& field,
                      const std::vector<Real>& r,
                      const std::vector<Real>& stations, Real max_dphi)
{
  std::vector<Real> R0, Z0;
  for (const Real ri : r) {
    R0.push_back(kAxisR + ri);
    Z0.push_back(0.0);
  }
  const auto trace = Trace(field, R0, Z0, stations, max_dphi);
  Real max_error = 0.0;
  for (size_t i = 0; i < r.size(); ++i) {
    REQUIRE(trace.valid(i));
    for (size_t k = 0; k < stations.size(); ++k) {
      const auto exact =
        AnalyticFieldLine(r[i], stations[k], field.GetPoloidalSign());
      max_error =
        std::max(max_error, std::hypot(trace.positions(k, i, 0) - exact.R,
                                       trace.positions(k, i, 1) - exact.Z));
    }
  }
  return max_error;
}

// The spline reports out-of-range lookups with printf, so silence is checked by
// capturing stdout at the file-descriptor level.
std::string CaptureStdout(const std::function<void()>& run)
{
  std::fflush(stdout);
  std::FILE* capture = std::tmpfile();
  REQUIRE(capture != nullptr);
  const int saved = dup(fileno(stdout));
  dup2(fileno(capture), fileno(stdout));
  run();
  Kokkos::fence();
  std::fflush(stdout);
  dup2(saved, fileno(stdout));
  close(saved);
  REQUIRE(std::fseek(capture, 0, SEEK_SET) == 0);
  std::string text;
  char buffer[256];
  size_t count;
  while ((count = std::fread(buffer, 1, sizeof(buffer), capture)) > 0) {
    text.append(buffer, count);
  }
  std::fclose(capture);
  return text;
}

} // namespace

TEST_CASE("AxisymmetricMagneticField: psi and its gradient match the analytic "
          "field",
          "[field_line_tracer]")
{
  const pcms::AxisymmetricMagneticField field(MakeCircularEquilibrium(), 1);
  const std::vector<Real> R = {2.3, 1.5, 2.8, 2.0, 3.5};
  const std::vector<Real> Z = {0.1, -0.4, 0.7, 0.0, 0.0};
  const auto values = EvaluateField(field, R, Z);
  for (size_t i = 0; i + 1 < R.size(); ++i) {
    CAPTURE(i, R[i], Z[i]);
    REQUIRE(values.inside(i));
    CHECK(values.psi(i) == Approx(CircularPsi(R[i], Z[i])).margin(1e-12));
    CHECK(values.psi_R(i) == Approx(2.0 * kC * (R[i] - kAxisR)).margin(1e-12));
    CHECK(values.psi_Z(i) == Approx(2.0 * kC * Z[i]).margin(1e-12));
    CHECK(values.I(i) == Approx(kI).margin(1e-12));
  }
  // (3.5, 0) lies outside the R-Z box: flagged, not extrapolated.
  const size_t outside = R.size() - 1;
  CHECK_FALSE(values.inside(outside));
  CHECK(values.psi(outside) == 0.0);
  CHECK(values.I(outside) == 0.0);
}

TEST_CASE("AxisymmetricMagneticField: a non-square box off the axis reads "
          "PSIZR with R fastest",
          "[field_line_tracer]")
{
  const pcms::AxisymmetricMagneticField field(MakeOffsetEquilibrium(), 1);
  CHECK(field.GetRMin() == Approx(1.25));
  CHECK(field.GetRMax() == Approx(3.0));
  CHECK(field.GetZMin() == Approx(-0.9));
  CHECK(field.GetZMax() == Approx(1.1));

  // The last two points lie inside the square test box but outside this one.
  const std::vector<Real> R = {2.3, 1.5, 2.8, 1.4, 2.9, 1.1, 2.0};
  const std::vector<Real> Z = {0.1, -0.4, 0.7, 0.95, -0.85, 0.0, -0.95};
  const auto values = EvaluateField(field, R, Z);
  for (size_t i = 0; i + 2 < R.size(); ++i) {
    CAPTURE(i, R[i], Z[i]);
    REQUIRE(values.inside(i));
    CHECK(values.psi(i) == Approx(CircularPsi(R[i], Z[i])).margin(1e-12));
    CHECK(values.psi_R(i) == Approx(2.0 * kC * (R[i] - kAxisR)).margin(1e-12));
    CHECK(values.psi_Z(i) == Approx(2.0 * kC * Z[i]).margin(1e-12));
  }
  CHECK_FALSE(values.inside(R.size() - 2));
  CHECK_FALSE(values.inside(R.size() - 1));

  const auto stations = pcms::WedgeStations(0.0, 2.0 * kPi / 8, 8);
  for (const Real r : {0.3, 0.7}) {
    const auto trace = Trace(field, {kAxisR + r}, {0.0}, stations, 0.01);
    REQUIRE(trace.valid(0));
    for (size_t k = 0; k < stations.size(); ++k) {
      const auto exact = AnalyticFieldLine(r, stations[k], 1);
      CAPTURE(r, k, stations[k]);
      CHECK(trace.positions(k, 0, 0) == Approx(exact.R).margin(1e-9));
      CHECK(trace.positions(k, 0, 1) == Approx(exact.Z).margin(1e-9));
    }
  }
}

TEST_CASE("AxisymmetricMagneticField: rejects invalid input",
          "[field_line_tracer]")
{
  const auto eqdsk = MakeCircularEquilibrium();
  CHECK_THROWS_AS(pcms::AxisymmetricMagneticField(eqdsk, 0), pcms::pcms_error);
  CHECK_THROWS_AS(pcms::AxisymmetricMagneticField(eqdsk, 2), pcms::pcms_error);

  const std::vector<Real> I(4, kI);
  const auto non_uniform =
    MakeCircularEquilibrium(1.0, {0.0, 0.01, 0.05, 0.13}, I);
  CHECK_THROWS_AS(pcms::AxisymmetricMagneticField(non_uniform, 1),
                  pcms::pcms_error);
  const auto non_monotonic =
    MakeCircularEquilibrium(1.0, {0.0, 0.05, 0.1, 0.0}, I);
  CHECK_THROWS_AS(pcms::AxisymmetricMagneticField(non_monotonic, 1),
                  pcms::pcms_error);
}

TEST_CASE("AxisymmetricMagneticField: a descending psi_grid gives the same "
          "I(psi)",
          "[field_line_tracer]")
{
  // psi decreasing outward with a G-EQDSK style grid from the axis (0) to the
  // edge (-0.13). I varies linearly so a wrong reversal would show.
  const auto descending_psi = UniformValues(0.0, -0.13, 27);
  std::vector<Real> descending_I(descending_psi.size());
  for (size_t i = 0; i < descending_psi.size(); ++i) {
    descending_I[i] = kI + 2.0 * descending_psi[i];
  }
  const std::vector<Real> ascending_psi(descending_psi.rbegin(),
                                        descending_psi.rend());
  const std::vector<Real> ascending_I(descending_I.rbegin(),
                                      descending_I.rend());
  const pcms::AxisymmetricMagneticField descending(
    MakeCircularEquilibrium(-1.0, descending_psi, descending_I), 1);
  const pcms::AxisymmetricMagneticField ascending(
    MakeCircularEquilibrium(-1.0, ascending_psi, ascending_I), 1);

  const std::vector<Real> R = {2.3, 1.5, 2.8, 2.0};
  const std::vector<Real> Z = {0.1, -0.4, 0.7, 0.0};
  const auto from_descending = EvaluateField(descending, R, Z);
  const auto from_ascending = EvaluateField(ascending, R, Z);
  for (size_t i = 0; i < R.size(); ++i) {
    const Real psi = -CircularPsi(R[i], Z[i]);
    CAPTURE(i, psi);
    CHECK(from_descending.psi(i) == Approx(psi).margin(1e-12));
    CHECK(from_descending.I(i) == Approx(from_ascending.I(i)).margin(1e-14));
    CHECK(from_descending.I(i) == Approx(kI + 2.0 * psi).margin(1e-12));
  }
}

TEST_CASE(
  "AxisymmetricMagneticField: I(psi) outside psi_grid is the end value, "
  "with no output",
  "[field_line_tracer]")
{
  // psi_grid covers [0.01, 0.05] only: the axis (psi = 0) lies below it and
  // (2.9, 0) (psi = 0.050625) above it. (3.5, 0) is outside the R-Z box.
  const auto psi_grid = UniformValues(0.01, 0.05, 9);
  std::vector<Real> I_psi(psi_grid.size());
  for (size_t i = 0; i < psi_grid.size(); ++i) {
    I_psi[i] = kI + 2.0 * psi_grid[i];
  }
  const pcms::AxisymmetricMagneticField field(
    MakeCircularEquilibrium(1.0, psi_grid, I_psi), 1);

  const std::vector<Real> R = {2.0, 2.9, 2.5, 3.5};
  const std::vector<Real> Z = {0.0, 0.0, 0.0, 0.0};
  HostFieldValues values;
  const auto output =
    CaptureStdout([&] { values = EvaluateField(field, R, Z); });
  CHECK(output.empty());
  CHECK(values.I(0) == Approx(kI + 2.0 * 0.01).margin(1e-12));
  CHECK(values.I(1) == Approx(kI + 2.0 * 0.05).margin(1e-12));
  CHECK(values.I(2) == Approx(kI + 2.0 * CircularPsi(2.5, 0.0)).margin(1e-12));
  CHECK_FALSE(values.inside(3));
}

TEST_CASE("TraceFieldLines: psi is conserved over a full toroidal turn",
          "[field_line_tracer]")
{
  const pcms::AxisymmetricMagneticField field(MakeCircularEquilibrium(), 1);
  const std::vector<Real> r = {0.3, 0.6, 0.9};
  std::vector<Real> R0, Z0;
  for (const Real ri : r) {
    R0.push_back(kAxisR + ri);
    Z0.push_back(0.0);
  }
  const auto trace = Trace(field, R0, Z0, {0.0, 2.0 * kPi}, 2.0 * kPi / 1000);
  for (size_t i = 0; i < r.size(); ++i) {
    CAPTURE(r[i]);
    REQUIRE(trace.valid(i));
    const Real psi =
      CircularPsi(trace.positions(1, i, 0), trace.positions(1, i, 1));
    CHECK(std::abs(psi - kC * r[i] * r[i]) < 1e-10);
  }
}

TEST_CASE("TraceFieldLines: matches the analytic circular-equilibrium solution",
          "[field_line_tracer]")
{
  const auto eqdsk = MakeCircularEquilibrium();
  const auto stations = pcms::WedgeStations(0.0, 2.0 * kPi / 8, 8);
  REQUIRE(stations.size() == 25);
  for (const int s : {1, -1}) {
    const pcms::AxisymmetricMagneticField field(eqdsk, s);
    for (const Real r : {0.3, 0.6, 0.9}) {
      const auto trace = Trace(field, {kAxisR + r}, {0.0}, stations, 0.01);
      REQUIRE(trace.valid(0));
      for (size_t k = 0; k < stations.size(); ++k) {
        const auto exact = AnalyticFieldLine(r, stations[k], s);
        CAPTURE(s, r, k, stations[k]);
        CHECK(trace.positions(k, 0, 0) == Approx(exact.R).margin(1e-9));
        CHECK(trace.positions(k, 0, 1) == Approx(exact.Z).margin(1e-9));
      }
    }
  }
}

TEST_CASE("TraceFieldLines: fourth-order convergence in phi",
          "[field_line_tracer]")
{
  // Halving the step must cut the error by about 2^4 = 16. Each station
  // interval of 2*pi/24 takes 3 substeps at 0.1 and 6 at 0.05, an exact
  // halving.
  const pcms::AxisymmetricMagneticField field(MakeCircularEquilibrium(), 1);
  const auto stations = pcms::WedgeStations(0.0, 2.0 * kPi / 8, 8);
  const std::vector<Real> r = {0.3, 0.6, 0.9};
  const Real coarse = MaxAnalyticError(field, r, stations, 0.1);
  const Real fine = MaxAnalyticError(field, r, stations, 0.05);
  CAPTURE(coarse, fine);
  REQUIRE(fine > 1e-13); // above rounding, so the ratio is meaningful
  CHECK(coarse / fine > 12.0);
}

TEST_CASE(
  "TraceFieldLines: matches the analytic solution with I varying in psi",
  "[field_line_tracer]")
{
  // psi is constant along a field line, so with I(psi) = I0 + a * psi each line
  // sees the constant I0 + a * c * r^2 and the analytic solution still holds.
  const Real slope = 2.0;
  const auto psi_grid = UniformValues(0.0, 0.13, 27);
  std::vector<Real> I_psi(psi_grid.size());
  for (size_t i = 0; i < psi_grid.size(); ++i) {
    I_psi[i] = kI + slope * psi_grid[i];
  }
  const pcms::AxisymmetricMagneticField field(
    MakeCircularEquilibrium(1.0, psi_grid, I_psi), 1);
  const auto stations = pcms::WedgeStations(0.0, 2.0 * kPi / 8, 8);
  for (const Real r : {0.3, 0.6, 0.9}) {
    const Real I = kI + slope * kC * r * r;
    const auto trace = Trace(field, {kAxisR + r}, {0.0}, stations, 0.01);
    REQUIRE(trace.valid(0));
    for (size_t k = 0; k < stations.size(); ++k) {
      const auto exact = AnalyticFieldLine(r, stations[k], 1, I);
      CAPTURE(r, I, k, stations[k]);
      CHECK(trace.positions(k, 0, 0) == Approx(exact.R).margin(1e-9));
      CHECK(trace.positions(k, 0, 1) == Approx(exact.Z).margin(1e-9));
    }
  }
}

TEST_CASE("TraceFieldLines: traces toward negative phi", "[field_line_tracer]")
{
  // Strictly decreasing stations integrate backwards along the field line.
  const auto eqdsk = MakeCircularEquilibrium();
  const auto stations = pcms::WedgeStations(0.0, -2.0 * kPi / 8, 8);
  REQUIRE(stations.back() == Approx(-2.0 * kPi));
  for (const int s : {1, -1}) {
    const pcms::AxisymmetricMagneticField field(eqdsk, s);
    for (const Real r : {0.3, 0.9}) {
      const auto trace = Trace(field, {kAxisR + r}, {0.0}, stations, 0.01);
      REQUIRE(trace.valid(0));
      for (size_t k = 0; k < stations.size(); ++k) {
        const auto exact = AnalyticFieldLine(r, stations[k], s);
        CAPTURE(s, r, k, stations[k]);
        CHECK(trace.positions(k, 0, 0) == Approx(exact.R).margin(1e-9));
        CHECK(trace.positions(k, 0, 1) == Approx(exact.Z).margin(1e-9));
      }
    }
  }
}

TEST_CASE("TraceFieldLines: stations in a device view give the same result",
          "[field_line_tracer]")
{
  // Both overloads run the same steps, so the results are identical.
  const pcms::AxisymmetricMagneticField field(MakeCircularEquilibrium(), 1);
  const auto stations = pcms::WedgeStations(0.0, 2.0 * kPi / 16, 4);
  const std::vector<Real> R0 = {2.3, 2.6};
  const std::vector<Real> Z0 = {0.0, 0.1};
  const auto from_vector = Trace(field, R0, Z0, stations, 0.05);
  const auto traced = pcms::TraceFieldLines(field, ToDevice(R0), ToDevice(Z0),
                                            ToDevice(stations), 0.05);
  const auto positions =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, traced.positions);
  const auto valid =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, traced.valid);
  const auto stations_out =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, traced.stations);
  REQUIRE(positions.extent(0) == stations.size());
  for (size_t k = 0; k < stations.size(); ++k) {
    CHECK(stations_out(k) == stations[k]);
    for (size_t i = 0; i < R0.size(); ++i) {
      CAPTURE(k, i);
      CHECK(valid(i) == from_vector.valid(i));
      CHECK(positions(k, i, 0) == from_vector.positions(k, i, 0));
      CHECK(positions(k, i, 1) == from_vector.positions(k, i, 1));
    }
  }
}

TEST_CASE("TraceFieldLines: the poloidal sign sets the direction of rotation",
          "[field_line_tracer]")
{
  // From (2.7, 0), dZ/dphi = -s * R * psi_R / I has the sign of -s.
  const auto eqdsk = MakeCircularEquilibrium();
  const pcms::AxisymmetricMagneticField plus(eqdsk, 1);
  const pcms::AxisymmetricMagneticField minus(eqdsk, -1);
  const auto up = Trace(plus, {2.7}, {0.0}, {0.0, 0.1}, 1e-2);
  const auto down = Trace(minus, {2.7}, {0.0}, {0.0, 0.1}, 1e-2);
  CHECK(up.positions(1, 0, 1) < 0.0);
  CHECK(down.positions(1, 0, 1) > 0.0);
  // Flipping s mirrors the trajectory in Z.
  CHECK(up.positions(1, 0, 0) == Approx(down.positions(1, 0, 0)).margin(1e-12));
  CHECK(up.positions(1, 0, 1) ==
        Approx(-down.positions(1, 0, 1)).margin(1e-12));
}

TEST_CASE("WedgeStations: three stations per layer and the closing plane",
          "[field_line_tracer]")
{
  const auto stations = pcms::WedgeStations(0.1, 0.3, 4);
  REQUIRE(stations.size() == 13);
  for (int k = 0; k < 4; ++k) {
    for (int j = 0; j < 3; ++j) {
      CAPTURE(k, j);
      CHECK(stations[3 * k + j] ==
            Approx(0.1 + 0.3 * k + 0.1 * j).margin(1e-15));
    }
  }
  CHECK(stations.back() == Approx(1.3).margin(1e-15));
  CHECK(pcms::WedgeStations(0.5, 0.2, 0) == std::vector<Real>{0.5});
  CHECK_THROWS_AS(pcms::WedgeStations(0.0, 0.1, -1), pcms::pcms_error);
}

TEST_CASE("TraceFieldLines: tracing through all stations equals tracing each "
          "interval",
          "[field_line_tracer]")
{
  // Each interval takes the same substeps either way, so the results agree to
  // rounding.
  const pcms::AxisymmetricMagneticField field(MakeCircularEquilibrium(), 1);
  const auto stations = pcms::WedgeStations(0.0, 2.0 * kPi / 16, 4);
  const std::vector<Real> R0 = {2.3, 2.6};
  const std::vector<Real> Z0 = {0.0, 0.1};
  const Real max_dphi = 0.05;
  const auto full = Trace(field, R0, Z0, stations, max_dphi);

  std::vector<Real> R = R0, Z = Z0;
  for (size_t k = 1; k < stations.size(); ++k) {
    const auto step =
      Trace(field, R, Z, {stations[k - 1], stations[k]}, max_dphi);
    for (size_t i = 0; i < R0.size(); ++i) {
      CAPTURE(k, i);
      REQUIRE(full.valid(i));
      R[i] = step.positions(1, i, 0);
      Z[i] = step.positions(1, i, 1);
      CHECK(full.positions(k, i, 0) == Approx(R[i]).margin(1e-13));
      CHECK(full.positions(k, i, 1) == Approx(Z[i]).margin(1e-13));
    }
  }
}

TEST_CASE("TraceFieldLines: project_to_psi keeps the starting flux surface",
          "[field_line_tracer]")
{
  // A coarse step drifts off the surface; three Newton steps on an exactly
  // quadratic psi bring it back to rounding.
  const pcms::AxisymmetricMagneticField field(MakeCircularEquilibrium(), 1);
  const Real psi0 = CircularPsi(2.6, 0.0);
  const auto plain = Trace(field, {2.6}, {0.0}, {0.0, 2.0 * kPi}, 0.5);
  const auto projected =
    Trace(field, {2.6}, {0.0}, {0.0, 2.0 * kPi}, 0.5, true);
  REQUIRE(projected.valid(0));
  const Real drift = std::abs(
    CircularPsi(plain.positions(1, 0, 0), plain.positions(1, 0, 1)) - psi0);
  const Real projected_drift = std::abs(
    CircularPsi(projected.positions(1, 0, 0), projected.positions(1, 0, 1)) -
    psi0);
  CAPTURE(drift, projected_drift);
  CHECK(projected_drift < 1e-12);
  CHECK(projected_drift < drift);
}

TEST_CASE("TraceFieldLines: a field line leaving the grid is flagged, not "
          "thrown",
          "[field_line_tracer]")
{
  // (2.9, 0.9) lies on a surface of radius 1.27 that leaves the R-Z box,
  // (3.5, 0) starts outside it, and (2.3, 0) stays inside.
  const pcms::AxisymmetricMagneticField field(MakeCircularEquilibrium(), 1);
  const std::vector<Real> R0 = {2.9, 3.5, 2.3};
  const std::vector<Real> Z0 = {0.9, 0.0, 0.0};
  HostTrace trace;
  REQUIRE_NOTHROW(trace = Trace(field, R0, Z0, {0.0, kPi, 2.0 * kPi}, 1e-2));
  CHECK_FALSE(trace.valid(0));
  CHECK_FALSE(trace.valid(1));
  CHECK(trace.valid(2));
  for (size_t k = 0; k < 3; ++k) {
    CAPTURE(k);
    const Real R = trace.positions(k, 0, 0);
    const Real Z = trace.positions(k, 0, 1);
    CHECK(R >= field.GetRMin());
    CHECK(R <= field.GetRMax());
    CHECK(Z >= field.GetZMin());
    CHECK(Z <= field.GetZMax());
    CHECK(trace.positions(k, 1, 0) == 3.5);
    CHECK(trace.positions(k, 1, 1) == 0.0);
  }
}

TEST_CASE("TraceFieldLines: returns to its start after one poloidal turn",
          "[field_line_tracer]")
{
  // d(theta)/d(phi) = -s * (2c / I) * (R0 + r * cos(theta)), so one poloidal
  // turn takes 2*pi / ((2c / I) * sqrt(R0^2 - r^2)) = 8*pi / sqrt(4 - r^2)
  // here, whatever the starting poloidal angle.
  const std::vector<Real> R0 = {2.3, 2.0, 1.6};
  const std::vector<Real> Z0 = {0.0, 0.5, -0.3}; // r = 0.3, 0.5, 0.5
  for (const int s : {1, -1}) {
    const pcms::AxisymmetricMagneticField field(MakeCircularEquilibrium(), s);
    for (size_t i = 0; i < R0.size(); ++i) {
      const Real r = std::hypot(R0[i] - kAxisR, Z0[i]);
      const Real period = 8.0 * kPi / std::sqrt(4.0 - r * r);
      const auto trace = Trace(field, {R0[i]}, {Z0[i]}, {0.0, period}, 0.01);
      CAPTURE(s, i, r, period);
      REQUIRE(trace.valid(0));
      CHECK(trace.positions(1, 0, 0) == Approx(R0[i]).margin(1e-8));
      CHECK(trace.positions(1, 0, 1) == Approx(Z0[i]).margin(1e-8));
    }
  }
}

TEST_CASE("TraceFieldLines: matches an independent reference trace of a mesh "
          "triangle",
          "[field_line_tracer]")
{
  // Triangle 310 of an XGC-style mesh of this equilibrium, traced from phi = 0
  // to pi/4 by an independent Python tracer (SciPy DOP853, rtol 1e-10; it
  // agrees with the analytic field line to 3e-12). expected[k][a] is the
  // (R, Z) of vertex a at phi = k * pi/12.
  const Real expected[4][3][2] = {
    {{2.736144, -0.050356}, {2.736144, 0.050356}, {2.630506, 0.04962}},
    {{2.7155011378747225, -0.1802826590994909},
     {2.7333381574447446, -0.08160185235359393},
     {2.6296808756613186, -0.05917732050701803}},
    {{2.6727471824721314, -0.30307580561280045},
     {2.707267573646819, -0.2102767099317841},
     {2.610331211495885, -0.1658184932688772}},
    {{2.610537244517441, -0.4143524930897984},
     {2.6596092412438246, -0.3306952741392141},
     {2.573550883124193, -0.2665320710600273}}};
  const pcms::AxisymmetricMagneticField field(MakeCircularEquilibrium(), 1);
  std::vector<Real> R0, Z0;
  for (int a = 0; a < 3; ++a) {
    R0.push_back(expected[0][a][0]);
    Z0.push_back(expected[0][a][1]);
  }
  const auto stations = pcms::WedgeStations(0.0, kPi / 4, 1);
  const auto trace = Trace(field, R0, Z0, stations, 0.01);
  for (int a = 0; a < 3; ++a) {
    REQUIRE(trace.valid(a));
    for (int k = 0; k < 4; ++k) {
      CAPTURE(a, k);
      CHECK(trace.positions(k, a, 0) == Approx(expected[k][a][0]).margin(1e-9));
      CHECK(trace.positions(k, a, 1) == Approx(expected[k][a][1]).margin(1e-9));
    }
  }
}
