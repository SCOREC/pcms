// The 1D spaces on a periodic geometry: Lagrange order 1, periodic cubic
// spline of grid values, Hermite cubic (values and derivatives) and real
// Fourier. Each is a kernel that gives the weight of DOF k in the value (or
// derivative) at a point; the evaluation matrix is built when points are
// bound (sparse_rows.hpp).

#include "pcms/field/function_spaces.hpp"
#include "pcms/geometry/periodic_geometry.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include "../../geometry/periodic_binding.hpp"
#include "sparse_rows.hpp"
#include <cmath>
#include <numbers>
#include <set>
#include <vector>

namespace pcms
{

namespace
{

constexpr Real TwoPi = 2 * std::numbers::pi_v<Real>;

/// A serial layout of `size` DOF holders with one component, at the given
/// positions (zeros for spaces whose DOFs have none).
class PeriodicLayout final : public FieldLayout
{
public:
  PeriodicLayout(std::shared_ptr<const Geometry> geometry,
                 const std::vector<Real>& positions)
    : size_(static_cast<LO>(positions.size())),
      owned_("owned", size_),
      gids_("gids", size_),
      zeros_("zeros", size_)
  {
    Kokkos::View<Real**, HostMemorySpace> host_coords("coords", size_, 1);
    for (LO i = 0; i < size_; ++i) {
      owned_(i) = true;
      gids_(i) = i;
      host_coords(i, 0) = positions[i];
    }
    coords_ = pcms::detail::ToDeviceView(host_coords);
    SetCoordinateSystem(geometry->GetCoordinateSystem());
    SetGeometry(std::move(geometry));
  }
  int GetNumComponents() const override { return 1; }
  LO GetNumDofHolder() const override { return size_; }
  GO GetNumGlobalDofHolder() const override { return size_; }
  bool IsDistributed() const override { return false; }
  Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override
  {
    return MakeRank1View(owned_);
  }
  GlobalIDView<HostMemorySpace> GetGidsHost() const override
  {
    return GlobalIDView<HostMemorySpace>(gids_.data(), size_);
  }
  EntOffsetsArray GetEntOffsets() const override
  {
    const auto n = static_cast<std::size_t>(size_);
    return {0, n, n, n, n};
  }
  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override
  {
    return CoordinateView<DeviceMemorySpace>(GetCoordinateSystem(),
                                             MakeConstRank2View(coords_));
  }
  int GetDimension() const override { return 1; }
  Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override
  {
    return MakeRank1View(zeros_);
  }
  Rank1View<const LO, HostMemorySpace> GetDOFHolderClassificationIdsHost()
    const override
  {
    return MakeRank1View(zeros_);
  }

private:
  LO size_;
  Kokkos::View<bool*, HostMemorySpace> owned_;
  Kokkos::View<GO*, HostMemorySpace> gids_;
  Kokkos::View<LO*, HostMemorySpace> zeros_;
  Kokkos::View<Real**, DeviceMemorySpace> coords_;
};

/// The rows of a kernel: entry k of row p is kernel(p, k, derivative, gid,
/// weight); every row has `entries` entries.
template <typename Kernel>
detail::SparseRows BuildRows(const Kernel& kernel, LO num_points, int entries,
                             bool derivative)
{
  detail::SparseRows rows;
  rows.offsets =
    Kokkos::View<LO*, DeviceMemorySpace>("rows_offsets", num_points + 1);
  const auto n = static_cast<std::size_t>(num_points) * entries;
  rows.columns = Kokkos::View<GO*, DeviceMemorySpace>("rows_columns", n);
  rows.holders = Kokkos::View<LO*, DeviceMemorySpace>("rows_holders", n);
  rows.values = Kokkos::View<Real*, DeviceMemorySpace>("rows_values", n);
  const auto offsets = rows.offsets;
  const auto columns = rows.columns;
  const auto holders = rows.holders;
  const auto values = rows.values;
  pcms::detail::ForEach(
    "BuildRows", num_points + 1,
    KOKKOS_LAMBDA(const LO p) { offsets(p) = p * entries; });
  pcms::detail::ForEach(
    "BuildRows", static_cast<LO>(n), KOKKOS_LAMBDA(const LO i) {
      GO gid = 0;
      Real weight = 0;
      kernel(i / entries, i % entries, derivative, gid, weight);
      columns(i) = gid;
      holders(i) = static_cast<LO>(gid);
      values(i) = weight;
    });
  return rows;
}

/// The periodic binding of a bound set of points of a periodic geometry.
const pcms::detail::PeriodicBinding& AsPeriodic(const BoundPoints& bound)
{
  return dynamic_cast<const pcms::detail::PeriodicBinding&>(bound);
}

// Lagrange order 1: DOFs at the nodes, hat functions.
struct LagrangeKernel
{
  Kokkos::View<const LO*, DeviceMemorySpace> cells;
  Kokkos::View<const Real*, DeviceMemorySpace> xi;
  LO num_cells;
  Real h;
  KOKKOS_INLINE_FUNCTION void operator()(LO p, int k, bool derivative, GO& gid,
                                         Real& weight) const
  {
    const LO c = cells(p);
    gid = k == 0 ? c : (c + 1) % num_cells;
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
  KOKKOS_INLINE_FUNCTION void operator()(LO p, int k, bool derivative, GO& gid,
                                         Real& weight) const
  {
    const LO c = cells(p);
    const Real t = xi(p);
    const LO node = k < 2 ? c : (c + 1) % num_cells;
    gid = 2 * node + k % 2;
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
  KOKKOS_INLINE_FUNCTION void operator()(LO p, int k, bool derivative, GO& gid,
                                         Real& weight) const
  {
    const LO c = cells(p);
    const LO c1 = (c + 1) % num_cells;
    const Real t = xi(p);
    const Real s = 1 - t;
    gid = k;
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
  KOKKOS_INLINE_FUNCTION void operator()(LO p, int k, bool derivative, GO& gid,
                                         Real& weight) const
  {
    gid = k;
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
  std::shared_ptr<const FieldLayout> layout,
  detail::SparseRowEvaluatorFactory::MakeRows make_rows,
  std::function<Quadrature()> make_quadrature)
{
  auto evaluator = std::make_shared<detail::SparseRowEvaluatorFactory>(
    geometry, layout,
    std::set<DifferentialOperator>{DifferentialOperator::Value,
                                   DifferentialOperator::Gradient},
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
  auto layout = std::make_shared<PeriodicLayout>(geometry, NodePositions(g, 1));
  const LO num_cells = g.NumCells();
  const Real h = g.Period() / num_cells;
  return MakeSpace(
    geometry, layout,
    [num_cells, h](const std::shared_ptr<const BoundPoints>& bound,
                   DifferentialOperator op) {
      const auto& binding = AsPeriodic(*bound);
      return BuildRows(LagrangeKernel{binding.CellIndices(),
                                      binding.Parametric(), num_cells, h},
                       bound->NumPoints(), 2,
                       op == DifferentialOperator::Gradient);
    },
    [geometry] { return CreateGaussQuadrature(geometry, 2); });
}

std::shared_ptr<const FunctionSpace> CreatePeriodicHermiteSpace(
  std::shared_ptr<const Geometry> geometry)
{
  const auto& g = RequirePeriodic(geometry, "CreatePeriodicHermiteSpace");
  auto layout = std::make_shared<PeriodicLayout>(geometry, NodePositions(g, 2));
  const LO num_cells = g.NumCells();
  const Real h = g.Period() / num_cells;
  return MakeSpace(
    geometry, layout,
    [num_cells, h](const std::shared_ptr<const BoundPoints>& bound,
                   DifferentialOperator op) {
      const auto& binding = AsPeriodic(*bound);
      return BuildRows(HermiteKernel{binding.CellIndices(),
                                     binding.Parametric(), num_cells, h},
                       bound->NumPoints(), 4,
                       op == DifferentialOperator::Gradient);
    },
    [geometry] { return CreateGaussQuadrature(geometry, 6); });
}

std::shared_ptr<const FunctionSpace> CreatePeriodicSplineSpace(
  std::shared_ptr<const Geometry> geometry)
{
  const auto& g = RequirePeriodic(geometry, "CreatePeriodicSplineSpace");
  const LO n = g.NumCells();
  if (n < 3) {
    throw pcms_error("CreatePeriodicSplineSpace: needs at least 3 cells");
  }
  auto layout = std::make_shared<PeriodicLayout>(geometry, NodePositions(g, 1));
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
    [n, h, circulant](const std::shared_ptr<const BoundPoints>& bound,
                      DifferentialOperator op) {
      const auto& binding = AsPeriodic(*bound);
      return BuildRows(SplineKernel{binding.CellIndices(), binding.Parametric(),
                                    circulant, n, h},
                       bound->NumPoints(), n,
                       op == DifferentialOperator::Gradient);
    },
    [geometry] { return CreateGaussQuadrature(geometry, 6); });
}

std::shared_ptr<const FunctionSpace> CreateFourierSpace(
  std::shared_ptr<const Geometry> geometry, std::vector<int> modes, int stride)
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
  auto layout = std::make_shared<PeriodicLayout>(
    geometry, std::vector<Real>(freq.size(), 0.0));
  const auto entries = static_cast<int>(freq.size());
  const Real omega = TwoPi / g.Period();
  return MakeSpace(
    geometry, layout,
    [entries, omega, device_freq,
     device_sin](const std::shared_ptr<const BoundPoints>& bound,
                 DifferentialOperator op) {
      return BuildRows(FourierKernel{AsPeriodic(*bound).Phase(), device_freq,
                                     device_sin, omega},
                       bound->NumPoints(), entries,
                       op == DifferentialOperator::Gradient);
    },
    // the product of two functions has frequencies up to 2 max_mode stride
    [geometry, n = 2 * max_mode * stride + 1] {
      return CreateTrapezoidQuadrature(geometry, n);
    });
}

} // namespace pcms
