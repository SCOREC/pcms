// The 1D spaces on a periodic geometry: Lagrange order 1, periodic cubic
// spline of grid values, Hermite cubic (values and derivatives) and real
// Fourier. Each is a kernel that gives the weight of DOF k in the value (or
// derivative) at a point; the evaluation rows are built when points are
// bound (point_rows.hpp).

#include "pcms/field/function_spaces.hpp"
#include "pcms/geometry/periodic_geometry.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include "../../geometry/periodic_binding.hpp"
#include "pcms/field/layout/serial.hpp"
#include "pcms/field/point_rows.hpp"
#include "quadrature_space.hpp"
#include <cmath>
#include <numbers>
#include <set>
#include <vector>

namespace pcms
{

namespace
{

constexpr Real TwoPi = 2 * std::numbers::pi_v<Real>;

/// The serial layout of DOF holders at the given positions (zeros for spaces
/// whose DOFs have none).
std::shared_ptr<const FieldLayout> LayoutAt(
  std::shared_ptr<const Geometry> geometry, const std::vector<Real>& positions,
  int num_components)
{
  const auto n = positions.size();
  Kokkos::View<Real**, HostMemorySpace> coords("coords", n, 1);
  Kokkos::View<GO*, HostMemorySpace> gids("gids", n);
  Kokkos::View<bool*, HostMemorySpace> owned("owned", n);
  for (std::size_t i = 0; i < n; ++i) {
    coords(i, 0) = positions[i];
    gids(i) = static_cast<GO>(i);
    owned(i) = true;
  }
  return std::make_shared<pcms::detail::SerialLayout>(
    std::move(geometry), coords, gids, owned, static_cast<GO>(n),
    num_components);
}

using pcms::detail::PeriodicBinding;
using RowsFactory = pcms::detail::RowsEvaluatorFactory<PeriodicBinding>;

/// The rows of a kernel: entry k of row p is kernel(p, k, derivative, holder,
/// weight); every row has `entries` entries.
template <typename Kernel>
pcms::detail::PointRows BuildRows(const Kernel& kernel, LO num_points,
                                  int entries, bool derivative)
{
  const auto n = static_cast<std::size_t>(num_points) * entries;
  Kokkos::View<LO*, DeviceMemorySpace> holders("rows_holders", n);
  Kokkos::View<Real*, DeviceMemorySpace> values("rows_values", n);
  pcms::detail::ForEach(
    "BuildRows", static_cast<LO>(n), KOKKOS_LAMBDA(const LO i) {
      kernel(i / entries, i % entries, derivative, holders(i), values(i));
    });
  return {pcms::detail::StrideOffsets(num_points, entries), holders, values};
}

// Lagrange order 1: DOFs at the nodes, hat functions.
struct LagrangeKernel
{
  Kokkos::View<const LO*, DeviceMemorySpace> cells;
  Kokkos::View<const Real*, DeviceMemorySpace> xi;
  LO num_cells;
  Real h;
  KOKKOS_INLINE_FUNCTION void operator()(LO p, int k, bool derivative,
                                         LO& holder, Real& weight) const
  {
    const LO c = cells(p);
    holder = k == 0 ? c : (c + 1) % num_cells;
    const Real t = xi(p);
    weight = derivative ? (k == 0 ? -1 / h : 1 / h) : (k == 0 ? 1 - t : t);
  }
};

// Hermite cubic: DOFs (value, derivative) of node i are 2 i and 2 i + 1.
struct HermiteKernel
{
  Kokkos::View<const LO*, DeviceMemorySpace> cells;
  Kokkos::View<const Real*, DeviceMemorySpace> xi;
  LO num_cells;
  Real h;
  KOKKOS_INLINE_FUNCTION void operator()(LO p, int k, bool derivative,
                                         LO& holder, Real& weight) const
  {
    const LO c = cells(p);
    const Real t = xi(p);
    const LO node = k < 2 ? c : (c + 1) % num_cells;
    holder = 2 * node + k % 2;
    if (!derivative) {
      const Real w[4] = {2 * t * t * t - 3 * t * t + 1,
                         h * (t * t * t - 2 * t * t + t),
                         -2 * t * t * t + 3 * t * t, h * (t * t * t - t * t)};
      weight = w[k];
    } else {
      const Real w[4] = {(6 * t * t - 6 * t) / h, 3 * t * t - 4 * t + 1,
                         (-6 * t * t + 6 * t) / h, 3 * t * t - 2 * t};
      weight = w[k];
    }
  }
};

// Periodic cubic spline of the node values: with the moments M = S f, the
// spline on cell c is linear in f and M; S is circulant, S(i, j) =
// circulant((i - j) mod n).
struct SplineKernel
{
  Kokkos::View<const LO*, DeviceMemorySpace> cells;
  Kokkos::View<const Real*, DeviceMemorySpace> xi;
  Kokkos::View<const Real*, DeviceMemorySpace> circulant;
  LO num_cells;
  Real h;
  KOKKOS_INLINE_FUNCTION Real Moment(LO i, LO j) const
  {
    return circulant((i - j + num_cells) % num_cells);
  }
  KOKKOS_INLINE_FUNCTION void operator()(LO p, int k, bool derivative,
                                         LO& holder, Real& weight) const
  {
    const LO c = cells(p);
    const LO c1 = (c + 1) % num_cells;
    const Real t = xi(p);
    const Real s = 1 - t;
    holder = k;
    const Real left = k == c ? 1 : 0;
    const Real right = k == c1 ? 1 : 0;
    if (!derivative) {
      weight =
        s * left + t * right +
        h * h / 6 *
          ((s * s * s - s) * Moment(c, k) + (t * t * t - t) * Moment(c1, k));
    } else {
      weight = (right - left) / h + h / 6 *
                                      ((1 - 3 * s * s) * Moment(c, k) +
                                       (3 * t * t - 1) * Moment(c1, k));
    }
  }
};

// Real Fourier: DOF g is cos or sin of frequency freq(g) of the phase.
struct FourierKernel
{
  Kokkos::View<const Real*, DeviceMemorySpace> phase;
  Kokkos::View<const Real*, DeviceMemorySpace> freq;
  Kokkos::View<const int*, DeviceMemorySpace> is_sin;
  Real omega; // 2 pi / period
  KOKKOS_INLINE_FUNCTION void operator()(LO p, int k, bool derivative,
                                         LO& holder, Real& weight) const
  {
    holder = k;
    const Real theta = freq(k) * TwoPi * phase(p);
    const Real c = Kokkos::cos(theta);
    const Real s = Kokkos::sin(theta);
    if (!derivative) {
      weight = is_sin(k) ? s : c;
    } else {
      weight = freq(k) * omega * (is_sin(k) ? c : -s);
    }
  }
};

const PeriodicGeometry& RequirePeriodic(
  const std::shared_ptr<const Geometry>& g, const char* who)
{
  const auto periodic = std::dynamic_pointer_cast<const PeriodicGeometry>(g);
  if (periodic == nullptr) {
    throw pcms_error(std::string(who) + ": needs a periodic geometry");
  }
  return *periodic;
}

std::vector<Real> NodePositions(const PeriodicGeometry& g, int per_node)
{
  std::vector<Real> x;
  for (LO i = 0; i < g.NumCells(); ++i) {
    for (int j = 0; j < per_node; ++j) {
      x.push_back(g.Origin() + g.Period() * i / g.NumCells());
    }
  }
  return x;
}

/// The space of `layout` whose rows `kernel_rows` makes for a point set.
std::shared_ptr<const FunctionSpace> MakeSpace(
  std::shared_ptr<const Geometry> geometry,
  std::shared_ptr<const FieldLayout> layout, RowsFactory::MakeRows make_rows,
  std::function<Quadrature()> make_quadrature)
{
  auto evaluator = std::make_shared<RowsFactory>(
    geometry, layout, std::set{DifferentialOperator::Gradient},
    std::move(make_rows));
  return std::make_shared<detail::QuadratureSpace>(
    std::move(layout), std::move(evaluator), std::move(geometry),
    std::move(make_quadrature));
}

} // namespace

std::shared_ptr<const FunctionSpace> CreatePeriodicLagrangeSpace(
  std::shared_ptr<const Geometry> geometry)
{
  const auto& g = RequirePeriodic(geometry, "CreatePeriodicLagrangeSpace");
  auto layout = LayoutAt(geometry, NodePositions(g, 1), 1);
  const LO num_cells = g.NumCells();
  const Real h = g.Period() / num_cells;
  return MakeSpace(
    geometry, layout,
    [num_cells, h](const PeriodicBinding& binding, DifferentialOperator op) {
      return BuildRows(LagrangeKernel{binding.CellIndices(),
                                      binding.Parametric(), num_cells, h},
                       binding.NumPoints(), 2,
                       op == DifferentialOperator::Gradient);
    },
    [geometry] { return CreateGaussQuadrature(geometry, 2); });
}

std::shared_ptr<const FunctionSpace> CreatePeriodicHermiteSpace(
  std::shared_ptr<const Geometry> geometry, int num_components)
{
  const auto& g = RequirePeriodic(geometry, "CreatePeriodicHermiteSpace");
  auto layout = LayoutAt(geometry, NodePositions(g, 2), num_components);
  const LO num_cells = g.NumCells();
  const Real h = g.Period() / num_cells;
  return MakeSpace(
    geometry, layout,
    [num_cells, h](const PeriodicBinding& binding, DifferentialOperator op) {
      return BuildRows(HermiteKernel{binding.CellIndices(),
                                     binding.Parametric(), num_cells, h},
                       binding.NumPoints(), 4,
                       op == DifferentialOperator::Gradient);
    },
    [geometry] { return CreateGaussQuadrature(geometry, 6); });
}

std::shared_ptr<const FunctionSpace> CreatePeriodicSplineSpace(
  std::shared_ptr<const Geometry> geometry, int num_components)
{
  const auto& g = RequirePeriodic(geometry, "CreatePeriodicSplineSpace");
  const LO n = g.NumCells();
  if (n < 3) {
    throw pcms_error("CreatePeriodicSplineSpace: needs at least 3 cells");
  }
  auto layout = LayoutAt(geometry, NodePositions(g, 1), num_components);
  const Real h = g.Period() / n;
  // M = S f solves (h / 6) (M[i-1] + 4 M[i] + M[i+1]) = (f[i-1] - 2 f[i] +
  // f[i+1]) / h; both sides are circulant, so S is diagonal in the discrete
  // Fourier basis
  Kokkos::View<Real*, HostMemorySpace> host("spline_circulant", n);
  for (LO d = 0; d < n; ++d) {
    Real sum = 0;
    for (LO k = 0; k < n; ++k) {
      const Real phi = TwoPi * k / n;
      const Real lambda =
        (2 * std::cos(phi) - 2) / h / (h / 6 * (4 + 2 * std::cos(phi)));
      sum += lambda * std::cos(phi * d);
    }
    host(d) = sum / n;
  }
  const Kokkos::View<const Real*, DeviceMemorySpace> circulant =
    Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, host);
  return MakeSpace(
    geometry, layout,
    [n, h, circulant](const PeriodicBinding& binding, DifferentialOperator op) {
      return BuildRows(SplineKernel{binding.CellIndices(), binding.Parametric(),
                                    circulant, n, h},
                       binding.NumPoints(), n,
                       op == DifferentialOperator::Gradient);
    },
    [geometry] { return CreateGaussQuadrature(geometry, 6); });
}

std::shared_ptr<const FunctionSpace> CreateFourierSpace(
  std::shared_ptr<const Geometry> geometry, std::vector<int> modes, int stride,
  int num_components)
{
  const auto& g = RequirePeriodic(geometry, "CreateFourierSpace");
  if (stride < 1) {
    throw pcms_error("CreateFourierSpace: the stride must be positive");
  }
  // mode m is the DOF cos(m stride theta), and for m > 0 also sin(m stride
  // theta), in the order of the list
  std::vector<Real> freq;
  std::vector<int> is_sin;
  int max_mode = 0;
  for (const int m : modes) {
    if (m < 0) {
      throw pcms_error("CreateFourierSpace: the modes must not be negative");
    }
    max_mode = std::max(max_mode, m);
    freq.push_back(static_cast<Real>(m) * stride);
    is_sin.push_back(0);
    if (m > 0) {
      freq.push_back(static_cast<Real>(m) * stride);
      is_sin.push_back(1);
    }
  }
  Kokkos::View<Real*, HostMemorySpace> host_freq("fourier_freq", freq.size());
  Kokkos::View<int*, HostMemorySpace> host_sin("fourier_sin", freq.size());
  for (std::size_t i = 0; i < freq.size(); ++i) {
    host_freq(i) = freq[i];
    host_sin(i) = is_sin[i];
  }
  const Kokkos::View<const Real*, DeviceMemorySpace> device_freq =
    Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, host_freq);
  const Kokkos::View<const int*, DeviceMemorySpace> device_sin =
    Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, host_sin);
  auto layout =
    LayoutAt(geometry, std::vector<Real>(freq.size(), 0.0), num_components);
  const auto entries = static_cast<int>(freq.size());
  const Real omega = TwoPi / g.Period();
  return MakeSpace(
    geometry, layout,
    [entries, omega, device_freq, device_sin](const PeriodicBinding& binding,
                                              DifferentialOperator op) {
      return BuildRows(
        FourierKernel{binding.Phase(), device_freq, device_sin, omega},
        binding.NumPoints(), entries, op == DifferentialOperator::Gradient);
    },
    // the product of two functions has frequencies up to 2 max_mode stride
    [geometry, n = 2 * max_mode * stride + 1] {
      return CreateTrapezoidQuadrature(geometry, n);
    });
}

} // namespace pcms
