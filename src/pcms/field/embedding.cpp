#include "pcms/field/embedding.hpp"
#include "pcms/backends/uniform_grid/uniform_grid_geometry.hpp"
#include "pcms/field/detail/transformed_point_evaluator.hpp"
#include "pcms/geometry/detail/geometry.hpp"
#include "pcms/localization/candidate_map.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/dense_linear_algebra.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include "pcms/utility/detail/periodic.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace pcms
{

namespace
{

template <typename T>
using DeviceView = Kokkos::View<T, DeviceMemorySpace>;

// jacobian(p, j, i) = d x_j / d xi_i from the gradient of the position, whose
// component j has the partials d x_j / d xi_i at columns j * dp + i.
void ReshapeJacobian(DeviceView<const Real**> gradient,
                     DeviceView<Real***> jacobian)
{
  const int dt = jacobian.extent(1);
  const int dp = jacobian.extent(2);
  detail::ForEach(
    "position_jacobian", static_cast<LO>(jacobian.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      for (int j = 0; j < dt; ++j) {
        for (int i = 0; i < dp; ++i) {
          jacobian(p, j, i) = gradient(p, j * dp + i);
        }
      }
    });
}

// The metric g = J^T G J induced in the parameter system and its volume factor:
// the physical volume element per unit of parameter volume, sqrt(det g) times
// the factor of the coordinate system's volume element relative to its metric.
void InduceMetric(View<3, const Real, DeviceMemorySpace> jacobian,
                  const BoundCoordinateSystem& target,
                  DeviceView<Real***> metric, DeviceView<Real*> volume)
{
  const int dp = metric.extent(1);
  const int dt = jacobian.extent(1);
  const auto g_target = target.Metric();
  const auto volume_target = target.VolumeFactor();
  detail::ForEach(
    "embedding_induced_metric", static_cast<LO>(metric.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      Real g[3][3] = {};
      Real g_t[3][3] = {};
      for (int j = 0; j < dt; ++j) {
        for (int k = 0; k < dt; ++k) {
          g_t[j][k] = g_target(p, j, k);
        }
      }
      for (int a = 0; a < dp; ++a) {
        for (int b = 0; b < dp; ++b) {
          for (int j = 0; j < dt; ++j) {
            for (int k = 0; k < dt; ++k) {
              g[a][b] += jacobian(p, j, a) * g_t[j][k] * jacobian(p, k, b);
            }
          }
          metric(p, a, b) = g[a][b];
        }
      }
      volume(p) = Kokkos::sqrt(detail::Determinant(g, dp)) * volume_target(p) /
                  Kokkos::sqrt(detail::Determinant(g_t, dt));
    });
}

struct Evaluation
{
  DeviceView<Real**> mapped;    // [point][coordinate system dimension]
  DeviceView<Real***> jacobian; // [point][coordinate system][parameter]
};

// A position Function: its value and Jacobian at parameter points.
class Position
{
public:
  /// @param who the caller, for the messages
  Position(std::shared_ptr<const Function<Real>> position, const char* who)
    : position_(std::move(position))
  {
    const std::string prefix = std::string(who) + ": ";
    if (position_ == nullptr) {
      throw pcms_error(prefix + "the position must not be null");
    }
    if (position_->GetData().GetValueBasis().kind != ValueKind::Coordinates) {
      throw pcms_error(prefix + "the position's values must be coordinates "
                                "(ValueBasis::Coordinates)");
    }
    if (!position_->GetSpace().Supports(DifferentialOperator::Gradient)) {
      throw pcms_error(prefix +
                       "the position's space does not evaluate gradients");
    }
    if (ParameterSystem()->Dimension() > 3 || System()->Dimension() > 3) {
      throw pcms_error(prefix + "dimensions above 3 are not supported");
    }
  }

  [[nodiscard]] const std::shared_ptr<const Function<Real>>& Get() const
  {
    return position_;
  }
  /// the coordinate system the position's values are points of
  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>& System() const
  {
    return position_->GetData().GetValueBasis().system;
  }
  /// the coordinate system of the parameter geometry
  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>& ParameterSystem()
    const
  {
    return position_->GetSpace().GetCoordinateSystem();
  }

  /// The position and its Jacobian at parameter points; points outside the
  /// parameter geometry use its nearest cell.
  [[nodiscard]] Evaluation Evaluate(
    Rank2View<const Real, DeviceMemorySpace> xi) const
  {
    const FunctionSpace& space = position_->GetSpace();
    const auto geometry = space.GetGeometry();
    const auto n = xi.extent(0);
    const int dp = ParameterSystem()->Dimension();
    const int dt = System()->Dimension();
    const auto bound = geometry->Bind(
      CoordinateView<DeviceMemorySpace>(geometry->GetCoordinateSystem(), xi),
      OutOfBoundsPolicy{OutOfBoundsMode::NEAREST_BOUNDARY});
    const ValueBasis basis = position_->GetData().GetValueBasis();
    Evaluation e{DeviceView<Real**>("position_mapped", n, dt),
                 DeviceView<Real***>("position_jacobian", n, dt, dp)};
    space.CreatePointEvaluator<Real>(bound)->Evaluate(
      *position_,
      ValueView<Real, DeviceMemorySpace>(basis, MakeRank2View(e.mapped)));
    DeviceView<Real**> gradient("position_gradient", n,
                                static_cast<std::size_t>(dt) * dp);
    const auto evaluator =
      space.CreatePointEvaluator<Real>(bound, DifferentialOperator::Gradient);
    evaluator->Evaluate(
      *position_, ValueView<Real, DeviceMemorySpace>(
                    evaluator->OutputBasis(basis), MakeRank2View(gradient)));
    ReshapeJacobian(gradient, e.jacobian);
    return e;
  }

private:
  std::shared_ptr<const Function<Real>> position_;
};

// The coordinate map of a position: parameter system -> the position's
// coordinate system.
class PositionMap final : public CoordinateMap
{
public:
  explicit PositionMap(Position position)
    : CoordinateMap(position.ParameterSystem(), position.System()),
      position_(std::move(position))
  {
  }

protected:
  [[nodiscard]] std::unique_ptr<BoundCoordinateMap> BindImpl(
    const CoordinateView<DeviceMemorySpace>& query_points) const override
  {
    auto e = position_.Evaluate(query_points.GetValues());
    return std::make_unique<detail::BoundJacobianMap>(
      GetTargetCoordinateSystem(), std::move(e.mapped), std::move(e.jacobian));
  }

private:
  Position position_;
};

// ---------------------------------------------------------------------------
// The embedded geometry
// ---------------------------------------------------------------------------

// Newton iterations of one candidate cell, the residual that counts as
// converged (relative to the extent of the image), the fraction of a cell
// extent the image bounding boxes are padded by (the image of a cell bulges
// beyond its sampled points), and the sample lattice per cell axis.
constexpr int MaxNewtonIterations = 30;
constexpr Real ResidualTol = 1e-10;
constexpr Real BoxPadding = 0.1;
constexpr int SamplesPerAxis = 4;

// Starts candidate round `round` at the points that are not found and have
// that many candidates: they become active at the center of the candidate
// cell.
// @return the number of active points
template <unsigned Dim>
LO StartRound(int round, UniformGrid<Dim> grid,
              detail::CandidateBins candidates, DeviceView<const bool*> found,
              DeviceView<bool*> active, DeviceView<bool*> converged,
              DeviceView<LO*> cell, DeviceView<Real**> xi)
{
  LO started = 0;
  Kokkos::parallel_reduce(
    "embedded_start_round",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(found.extent(0))),
    KOKKOS_LAMBDA(const LO p, LO& count) {
      const LO begin = candidates.row_map(p);
      active(p) = !found(p) && candidates.row_map(p + 1) - begin > round;
      converged(p) = false;
      if (active(p)) {
        cell(p) = candidates.entries(begin + round);
        const auto box = grid.GetCellBBOX(cell(p));
        for (unsigned d = 0; d < Dim; ++d) {
          xi(p, d) = box.center[d];
        }
        ++count;
      }
    },
    started);
  return started;
}

// One Gauss-Newton step of the active points that have not converged: the
// residual r = target - x(xi), to the nearest image along the axes with a
// `period`, is reduced by the step solving
// (J J^T) delta = J r, which is Newton's for a square J. The step is
// regularized where J is singular (a pole of the map): J J^T + eps I.
// @return the number of points that took a step
template <unsigned Dim>
LO NewtonStep(UniformGrid<Dim> grid, Kokkos::Array<Real, 3> period,
              Real residual_tolerance, bool update,
              Rank2View<const Real, DeviceMemorySpace> target,
              DeviceView<const Real**> mapped,
              View<3, const Real, DeviceMemorySpace> jacobian,
              DeviceView<const bool*> active, DeviceView<bool*> converged,
              DeviceView<const LO*> cell, DeviceView<Real**> xi)
{
  const int dt = target.extent(1);
  LO stepped = 0;
  Kokkos::parallel_reduce(
    "embedded_newton_step",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(target.extent(0))),
    KOKKOS_LAMBDA(const LO p, LO& count) {
      if (!active(p) || converged(p)) {
        return;
      }
      Real r[3] = {};
      Real norm2 = 0;
      for (int j = 0; j < dt; ++j) {
        r[j] =
          pcms::detail::NearestImage(target(p, j) - mapped(p, j), period[j]);
        norm2 += r[j] * r[j];
      }
      // a point within the tolerance still takes this step, which puts it
      // within roundoff of the solution (Newton converges quadratically):
      // otherwise a point on the face of a cell can be rejected by its cell
      const bool close = Kokkos::sqrt(norm2) <= residual_tolerance;
      if (!update) {
        converged(p) = close;
        return;
      }
      Real m[3][3] = {};
      Real b[3] = {};
      Real trace = 0;
      for (unsigned a = 0; a < Dim; ++a) {
        for (unsigned c = 0; c < Dim; ++c) {
          for (int j = 0; j < dt; ++j) {
            m[a][c] += jacobian(p, j, a) * jacobian(p, j, c);
          }
        }
        for (int j = 0; j < dt; ++j) {
          b[a] += jacobian(p, j, a) * r[j];
        }
        trace += m[a][a];
      }
      const Real scale = trace / Dim;
      Real reference = 1;
      for (unsigned a = 0; a < Dim; ++a) {
        reference *= scale;
      }
      Real regular[3][3];
      for (unsigned a = 0; a < Dim; ++a) {
        for (unsigned c = 0; c < Dim; ++c) {
          regular[a][c] = m[a][c];
        }
      }
      if (detail::Determinant(m, Dim) <= 1e-12 * reference) {
        for (unsigned a = 0; a < Dim; ++a) {
          regular[a][a] += 1e-6 * scale + 1e-300;
        }
      }
      detail::Solve(regular, Dim, b);
      const auto box = grid.GetCellBBOX(cell(p));
      for (unsigned d = 0; d < Dim; ++d) {
        // stay within a cell of the candidate: the position function is
        // evaluated by extrapolation outside the grid
        xi(p, d) =
          Kokkos::clamp(xi(p, d) + b[d], box.center[d] - 2 * box.half_width[d],
                        box.center[d] + 2 * box.half_width[d]);
      }
      converged(p) = close;
      ++count;
    },
    stepped);
  return stepped;
}

// A converged point whose parameter coordinates are in its candidate cell is
// found there.
template <unsigned Dim>
void AcceptRound(UniformGrid<Dim> grid, Real slack,
                 DeviceView<const bool*> active,
                 DeviceView<const bool*> converged, DeviceView<const LO*> cell,
                 DeviceView<const Real**> xi, DeviceView<bool*> found)
{
  detail::ForEach(
    "embedded_accept_round", static_cast<LO>(active.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      if (!active(p) || !converged(p)) {
        return;
      }
      const auto box = grid.GetCellBBOX(cell(p));
      bool inside = true;
      for (unsigned d = 0; d < Dim; ++d) {
        inside = inside && Kokkos::fabs(xi(p, d) - box.center[d]) <=
                             box.half_width[d] * (1 + slack);
      }
      found(p) = inside;
    });
}

// Marks the listed points found, the others missing.
void MarkFound(DeviceView<const LO*> indices, DeviceView<bool*> found)
{
  Kokkos::deep_copy(found, false);
  detail::ForEach(
    "embedded_mark_found", static_cast<LO>(indices.extent(0)),
    KOKKOS_LAMBDA(const LO k) { found(indices(k)) = true; });
}

// The tensor Gauss-Legendre rule over each cell of the grid, in the parameter
// space: point q of cell c is at the nodes of the digits of q, and its weight
// is the product of the nodes' weights times the cell's half-widths (the
// rule's weights sum to 2).
template <unsigned Dim>
void FillCellRule(UniformGrid<Dim> grid, DeviceView<const Real*> nodes,
                  DeviceView<const Real*> rule_weights, DeviceView<Real**> xi,
                  DeviceView<Real*> weights)
{
  const int n = nodes.extent(0);
  int per_cell = 1;
  for (unsigned d = 0; d < Dim; ++d) {
    per_cell *= n;
  }
  detail::ForEach(
    "embedded_cell_rule", static_cast<LO>(weights.extent(0)),
    KOKKOS_LAMBDA(const LO i) {
      const LO cell = i / per_cell;
      int digits = i % per_cell;
      const auto box = grid.GetCellBBOX(cell);
      Real weight = 1;
      for (unsigned d = 0; d < Dim; ++d) {
        const int k = digits % n;
        digits /= n;
        xi(i, d) = box.center[d] + box.half_width[d] * nodes(k);
        weight *= box.half_width[d] * rule_weights(k);
      }
      weights(i) = weight;
    });
}

// The nodes and weights of the n-point Gauss-Legendre rule on [-1, 1].
std::pair<std::vector<Real>, std::vector<Real>> GaussLegendre(int n)
{
  std::vector<Real> nodes(n);
  std::vector<Real> weights(n);
  for (int i = 0; i < n; ++i) {
    Real x = std::cos(std::numbers::pi * (i + 0.75) / (n + 0.5));
    Real derivative = 1;
    for (int iteration = 0; iteration < 100; ++iteration) {
      Real previous = 1;
      Real current = x;
      for (int k = 2; k <= n; ++k) {
        const Real next = ((2 * k - 1) * x * current - (k - 1) * previous) / k;
        previous = current;
        current = next;
      }
      derivative = n * (x * current - previous) / (x * x - 1);
      const Real step = current / derivative;
      x -= step;
      if (std::abs(step) < 1e-15) {
        break;
      }
    }
    nodes[i] = x;
    weights[i] = 2 / ((1 - x * x) * derivative * derivative);
  }
  return {nodes, weights};
}

// The bounding boxes of the images of the cells under the position: of a
// lattice of sample points of each cell, evaluated in batches of at most
// detail::EmbeddingSampleBatch samples.
template <unsigned Dim>
void SampleCellImages(const Position& position, const UniformGrid<Dim>& grid,
                      std::vector<std::array<Real, 3>>& lo,
                      std::vector<std::array<Real, 3>>& hi)
{
  const int dt = position.System()->Dimension();
  const LO num_cells = grid.GetNumCells();
  std::size_t per_cell = 1;
  for (unsigned d = 0; d < Dim; ++d) {
    per_cell *= SamplesPerAxis + 1;
  }
  const auto batch = static_cast<LO>(
    std::max<std::size_t>(1, detail::EmbeddingSampleBatch / per_cell));
  std::array<Real, 3> empty_lo, empty_hi;
  empty_lo.fill(std::numeric_limits<Real>::max());
  empty_hi.fill(std::numeric_limits<Real>::lowest());
  lo.assign(num_cells, empty_lo);
  hi.assign(num_cells, empty_hi);
  for (LO first = 0; first < num_cells; first += batch) {
    const LO count = std::min(batch, num_cells - first);
    DeviceView<Real**> device("embedded_samples", count * per_cell, Dim);
    const auto samples = Kokkos::create_mirror_view(device);
    for (LO c = 0; c < count; ++c) {
      const auto box = grid.GetCellBBOX(first + c);
      for (std::size_t q = 0; q < per_cell; ++q) {
        std::size_t digits = q;
        for (unsigned d = 0; d < Dim; ++d) {
          const auto digit = static_cast<Real>(digits % (SamplesPerAxis + 1));
          digits /= SamplesPerAxis + 1;
          samples(c * per_cell + q, d) =
            box.center[d] +
            (2 * digit / SamplesPerAxis - 1) * box.half_width[d];
        }
      }
    }
    Kokkos::deep_copy(device, samples);
    const auto images = Kokkos::create_mirror_view_and_copy(
      HostMemorySpace(), position.Evaluate(MakeConstRank2View(device)).mapped);
    for (LO c = 0; c < count; ++c) {
      for (std::size_t q = 0; q < per_cell; ++q) {
        for (int j = 0; j < dt; ++j) {
          const Real x = images(c * per_cell + q, j);
          lo[first + c][j] = std::min(lo[first + c][j], x);
          hi[first + c][j] = std::max(hi[first + c][j], x);
        }
      }
    }
  }
}

// The candidate map of the boxes [lo, hi] of the target system, on `count`
// cells per axis over them, modulo the periods `period` (0: not periodic).
template <int N>
CandidateMap<N> CreateBoxMap(const std::vector<std::array<Real, 3>>& lo,
                             const std::vector<std::array<Real, 3>>& hi,
                             LO count, const Kokkos::Array<Real, 3>& period)
{
  Kokkos::View<AABBox<N>*> boxes("embedded_boxes", lo.size());
  const auto host = Kokkos::create_mirror_view(boxes);
  std::array<Real, 3> grid_lo, grid_hi;
  grid_lo.fill(std::numeric_limits<Real>::max());
  grid_hi.fill(std::numeric_limits<Real>::lowest());
  for (std::size_t c = 0; c < lo.size(); ++c) {
    for (int j = 0; j < N; ++j) {
      host(c).center[j] = (lo[c][j] + hi[c][j]) / 2;
      host(c).half_width[j] = (hi[c][j] - lo[c][j]) / 2;
      grid_lo[j] = std::min(grid_lo[j], lo[c][j]);
      grid_hi[j] = std::max(grid_hi[j], hi[c][j]);
    }
  }
  Kokkos::deep_copy(boxes, host);
  UniformGrid<N> grid{};
  Kokkos::Array<bool, N> periodic{};
  for (int j = 0; j < N; ++j) {
    periodic[j] = period[j] > 0;
    grid.bot_left[j] = grid_lo[j];
    grid.edge_length[j] =
      periodic[j] ? period[j] : std::max(grid_hi[j] - grid_lo[j], 1e-300);
    grid.divisions[j] = count;
  }
  return CreateCandidateMap(Kokkos::View<const AABBox<N>*>(boxes), grid, 0,
                            periodic);
}

// ---------------------------------------------------------------------------

// What the BoundPoints of an embedded geometry hold: the parameter
// coordinates of the points, bound on the parameter geometry, and the
// geometric factors there (computed on first use).
class EmbeddedBound : public BoundPoints
{
public:
  /// @return the points bound on the parameter geometry, at the parameter
  ///         coordinates: what the evaluators of its spaces take
  [[nodiscard]] const std::shared_ptr<const BoundPoints>& Inner() const noexcept
  {
    return inner_;
  }

  /// @return the Pushforward of the Jacobian: parameter-system components to
  ///         components in the geometry's coordinate system, at the points
  [[nodiscard]] const ValueTransform& Pushforward() const
  {
    if (!pushforward_) {
      const auto jacobian = Jacobian();
      pushforward_.emplace(
        detail::Pushforward(jacobian, position_.ParameterSystem(),
                            position_.System(), system_->Metric()));
    }
    return *pushforward_;
  }

protected:
  EmbeddedBound(std::shared_ptr<const Geometry> geometry,
                const BindRequest& request, Position position)
    : BoundPoints(std::move(geometry), request), position_(std::move(position))
  {
  }

  // Takes the parameter coordinates of the points: `xi`, bound on `parameter`.
  void SetParameters(const Geometry& parameter, DeviceView<Real**> xi,
                     OutOfBoundsPolicy policy)
  {
    xi_ = std::move(xi);
    inner_ = parameter.Bind(
      CoordinateView<DeviceMemorySpace>(parameter.GetCoordinateSystem(),
                                        MakeConstRank2View(xi_)),
      policy);
    system_.reset();
    pushforward_.reset();
  }

  // J from the Gradient of the position at the parameter coordinates of the
  // points, and the metric J^T G J it induces
  [[nodiscard]] detail::GeometricFactors ComputeFactors() const override
  {
    auto e = position_.Evaluate(MakeConstRank2View(xi_));
    const auto n = xi_.extent(0);
    const int dp = position_.ParameterSystem()->Dimension();
    const int dt = position_.System()->Dimension();
    system_ = position_.System()->Bind(MakeConstRank2View(e.mapped));
    DeviceView<Real***> metric("induced_metric", n, dp, dp);
    DeviceView<Real*> volume("induced_volume_factor", n);
    InduceMetric(
      View<3, const Real, DeviceMemorySpace>(e.jacobian.data(), n, dt, dp),
      *system_, metric, volume);
    mapped_ = e.mapped;
    return detail::GeometricFactors(
      static_cast<LO>(n), e.jacobian,
      std::make_shared<const BoundCoordinateSystem>(static_cast<LO>(n), metric,
                                                    volume));
  }

private:
  Position position_;
  std::shared_ptr<const BoundPoints> inner_;
  DeviceView<Real**> xi_;
  mutable DeviceView<Real**> mapped_;
  // the coordinate system's metric at the mapped points
  mutable std::shared_ptr<const BoundCoordinateSystem> system_;
  mutable std::optional<ValueTransform> pushforward_;
};

} // namespace

// The search of an embedded geometry: the images of the parameter cells,
// binned modulo the periods of the coordinate system and the symmetry, and the
// position.
struct detail::EmbeddedSearch
{
  Position position;
  std::optional<TranslationSymmetry> symmetry;
  std::variant<CandidateMap<1>, CandidateMap<2>, CandidateMap<3>> map;
  Kokkos::Array<Real, 3> period = {};
  Real residual_tolerance = 0;
  Real tolerance = 0;
};

namespace
{

// The search over the cells of `grid`: the bounding boxes of their images,
// from a lattice of sample points of each cell, padded.
template <unsigned Dim>
std::shared_ptr<const detail::EmbeddedSearch> CreateSearch(
  const UniformGrid<Dim>& grid, Position position,
  const EmbeddingOptions& options)
{
  const auto& system = position.System();
  const int dt = system->Dimension();
  const LO num_cells = grid.GetNumCells();
  std::vector<std::array<Real, 3>> lo, hi;
  SampleCellImages<Dim>(position, grid, lo, hi);
  std::array<Real, 3> image_lo, image_hi;
  image_lo.fill(std::numeric_limits<Real>::max());
  image_hi.fill(std::numeric_limits<Real>::lowest());
  for (LO c = 0; c < num_cells; ++c) {
    for (int j = 0; j < dt; ++j) {
      image_lo[j] = std::min(image_lo[j], lo[c][j]);
      image_hi[j] = std::max(image_hi[j], hi[c][j]);
    }
  }
  Real extent2 = 0;
  for (int j = 0; j < dt; ++j) {
    extent2 += (image_hi[j] - image_lo[j]) * (image_hi[j] - image_lo[j]);
  }
  const Real extent = std::sqrt(extent2);
  // pad: the image of a cell bulges beyond its samples
  for (LO c = 0; c < num_cells; ++c) {
    for (int j = 0; j < dt; ++j) {
      const Real pad = BoxPadding * (hi[c][j] - lo[c][j]) + 1e-9 * extent;
      lo[c][j] -= pad;
      hi[c][j] += pad;
    }
  }
  Kokkos::Array<Real, 3> period = {};
  for (int j = 0; j < dt; ++j) {
    period[j] = system->Period(j).value_or(0);
  }
  if (options.symmetry) {
    period[options.symmetry->axis] = options.symmetry->length;
  }
  const LO count =
    options.search.divisions > 0
      ? options.search.divisions
      : std::max(1, static_cast<int>(std::ceil(std::pow(num_cells, 1.0 / dt))));
  auto search = std::make_shared<detail::EmbeddedSearch>(
    detail::EmbeddedSearch{std::move(position),
                           options.symmetry,
                           {},
                           period,
                           ResidualTol * extent,
                           options.search.tolerance});
  if (dt == 1) {
    search->map = CreateBoxMap<1>(lo, hi, count, period);
  } else if (dt == 2) {
    search->map = CreateBoxMap<2>(lo, hi, count, period);
  } else {
    search->map = CreateBoxMap<3>(lo, hi, count, period);
  }
  return search;
}

// The points with their coordinate on the symmetry's axis reduced into
// [0, length).
DeviceView<Real**> ReduceBySymmetry(
  Kokkos::View<const Real**, DeviceMemorySpace> points,
  TranslationSymmetry symmetry)
{
  DeviceView<Real**> reduced("embedded_reduced", points.extent(0),
                             points.extent(1));
  Kokkos::deep_copy(reduced, points);
  detail::ForEach(
    "embedded_reduce_by_symmetry", static_cast<LO>(points.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      const Real x = reduced(p, symmetry.axis);
      reduced(p, symmetry.axis) =
        x - symmetry.length * Kokkos::floor(x / symmetry.length);
    });
  return reduced;
}

struct Located
{
  DeviceView<Real**> xi;   // [point][parameter dimension]
  DeviceView<bool*> found; // [point]
};

// The parameter coordinates of points of the coordinate system: the candidate
// cells of each point, modulo the periods, are tried in turn by Gauss-Newton
// from the cell center.
template <unsigned Dim>
Located Search(const detail::EmbeddedSearch& search,
               const UniformGrid<Dim>& grid,
               Rank2View<const Real, DeviceMemorySpace> points)
{
  const auto n = points.extent(0);
  DeviceView<Real**> copy("embedded_points", n, points.extent(1));
  detail::CopyDeviceRank2ViewToDeviceView(copy, points);
  const auto candidates = std::visit(
    [&](const auto& map) { return map.GetCandidates(copy); }, search.map);
  DeviceView<LO*> cell("embedded_cell", n);
  DeviceView<bool*> active("embedded_active", n);
  DeviceView<bool*> converged("embedded_converged", n);
  Located located{DeviceView<Real**>("embedded_xi", n, Dim),
                  DeviceView<bool*>("embedded_found", n)};
  for (int round = 0; StartRound<Dim>(round, grid, candidates, located.found,
                                      active, converged, cell, located.xi) > 0;
       ++round) {
    for (int iteration = 0; iteration <= MaxNewtonIterations; ++iteration) {
      const auto e = search.position.Evaluate(MakeConstRank2View(located.xi));
      const LO stepped =
        NewtonStep<Dim>(grid, search.period, search.residual_tolerance,
                        iteration < MaxNewtonIterations, points, e.mapped,
                        View<3, const Real, DeviceMemorySpace>(
                          e.jacobian.data(), e.jacobian.extent(0),
                          e.jacobian.extent(1), e.jacobian.extent(2)),
                        active, converged, cell, located.xi);
      if (stepped == 0) {
        break;
      }
    }
    AcceptRound<Dim>(grid, search.tolerance, active, converged, cell,
                     located.xi, located.found);
  }
  return located;
}

// The tensor Gauss-Legendre rule of the parameter cells, exact for
// polynomials of the degree in the parameters, with the weights of the
// physical measure: the rule's weight times the volume factor at the
// points.
template <unsigned Dim>
Quadrature CellQuadrature(const EmbeddedGeometry& geometry,
                          const UniformGrid<Dim>& grid,
                          const detail::QuadratureRule& rule)
{
  const auto* gauss = std::get_if<detail::GaussRule>(&rule);
  if (gauss == nullptr) {
    throw pcms_error("CreateQuadrature: an embedded geometry has only the "
                     "Gauss rule");
  }
  const int per_axis = (gauss->order + 2) / 2;
  const auto [nodes_host, weights_host] = GaussLegendre(per_axis);
  DeviceView<Real*> nodes("gauss_nodes", per_axis);
  DeviceView<Real*> rule_weights("gauss_weights", per_axis);
  Kokkos::deep_copy(nodes,
                    Kokkos::View<const Real*, HostMemorySpace,
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                      nodes_host.data(), per_axis));
  Kokkos::deep_copy(rule_weights,
                    Kokkos::View<const Real*, HostMemorySpace,
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                      weights_host.data(), per_axis));
  size_t num_points = grid.GetNumCells();
  for (unsigned d = 0; d < Dim; ++d) {
    num_points *= per_axis;
  }
  DeviceView<Real**> xi("quadrature_xi", num_points, Dim);
  DeviceView<Real*> weights("quadrature_weights", num_points);
  FillCellRule<Dim>(grid, nodes, rule_weights, xi, weights);
  auto bound = geometry.Bind(CoordinateView<DeviceMemorySpace>(
    geometry.GetParameterSystem(), MakeConstRank2View(xi)));
  pcms::detail::ScaleByVolumeFactor(weights, *bound);
  return Quadrature{std::move(bound), weights, {}};
}

class EmbeddedBinding final : public EmbeddedBound
{
public:
  EmbeddedBinding(std::shared_ptr<const EmbeddedGeometry> geometry,
                  const detail::EmbeddedSearch& search,
                  const BindRequest& request)
    : EmbeddedBound(geometry, request, search.position),
      geometry_(std::move(geometry)),
      search_(search)
  {
    Locate();
  }

protected:
  void RebindImpl() override { Locate(); }

private:
  // Points in the parameter system are bound on the parameter geometry as
  // they are; points of the target system are searched for first. The
  // parameter geometry gets all the points (a missing one at its last
  // iterate, clamped): the evaluators fill the rows of the missing points.
  void Locate()
  {
    const auto points = LocatedPointValues();
    const auto n = NumPoints();
    const bool parametric =
      GetBoundMap() == nullptr &&
      pcms::detail::SameCoordinateSystem(GetCoordinateSystem(),
                                         geometry_->GetParameterSystem());
    DeviceView<bool*> found("embedded_found", n);
    if (parametric) {
      DeviceView<Real**> xi("embedded_xi", n, points.extent(1));
      Kokkos::deep_copy(xi, points);
      SetParameters(*geometry_->GetParameterGeometry(), xi, Policy());
      MarkFound(Inner()->FoundPoints(), found);
    } else {
      const auto located = uniform_grid::VisitUniformGrid(
        geometry_->GetParameterGeometry(),
        [&]<unsigned Dim>(
          const std::shared_ptr<const uniform_grid::UniformGridGeometry<Dim>>&
            grid) {
          return Search<Dim>(search_, grid->GetGrid(),
                             MakeConstRank2View(points));
        });
      if (search_.symmetry && GetBoundMap() == nullptr) {
        ReplaceLocatedPoints(ReduceBySymmetry(points, *search_.symmetry));
      }
      SetParameters(*geometry_->GetParameterGeometry(), located.xi,
                    OutOfBoundsPolicy{OutOfBoundsMode::NEAREST_BOUNDARY});
      found = located.found;
    }
    SetFound(found);
  }

  std::shared_ptr<const EmbeddedGeometry> geometry_;
  const detail::EmbeddedSearch& search_;
};

// Evaluates the fields of a space on the parameter geometry at the points of
// an embedded geometry: at their parameter coordinates. Vectors and covectors
// in parameter-system components come out in the geometry's coordinate system,
// through the Pushforward of the position's Jacobian.
class EmbeddedEvaluatorFactory final : public PointEvaluatorFactory
{
public:
  EmbeddedEvaluatorFactory(std::shared_ptr<const EmbeddedGeometry> embedded,
                           std::shared_ptr<const FunctionSpace> space)
    : embedded_(std::move(embedded)), space_(std::move(space))
  {
  }

  [[nodiscard]] std::shared_ptr<const Geometry> GetGeometry()
    const noexcept override
  {
    return embedded_;
  }
  [[nodiscard]] bool Supports(DifferentialOperator op) const noexcept override
  {
    return op == DifferentialOperator::Value;
  }

protected:
  PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, std::shared_ptr<const BoundPoints> bound,
    DifferentialOperator op) const override
  {
    const auto embedded = std::dynamic_pointer_cast<const EmbeddedBound>(bound);
    if (value_type != Type::Real || embedded == nullptr) {
      throw pcms_error("CreatePointEvaluator: an embedded space evaluates "
                       "Real fields at the points of its embedded geometry");
    }
    return std::unique_ptr<PointEvaluator<Real>>(
      std::make_unique<detail::TransformedPointEvaluator<Real>>(
        embedded,
        CreateUntransformedPointEvaluator<Real>(*space_, embedded->Inner(), op),
        std::vector<detail::LazyValueTransform>{
          [embedded]() -> const ValueTransform& {
            return embedded->Pushforward();
          }}));
  }

private:
  std::shared_ptr<const EmbeddedGeometry> embedded_;
  std::shared_ptr<const FunctionSpace> space_;
};

} // namespace

EmbeddedGeometry::EmbeddedGeometry(
  std::shared_ptr<const Geometry> parameter,
  std::shared_ptr<const Function<Real>> position,
  std::shared_ptr<const detail::EmbeddedSearch> search)
  : parameter_(std::move(parameter)),
    position_(std::move(position)),
    system_(position_->GetData().GetValueBasis().system),
    search_(std::move(search))
{
}

EmbeddedGeometry::~EmbeddedGeometry() noexcept = default;

int EmbeddedGeometry::TopologicalDimension() const noexcept
{
  return parameter_->GetCoordinateSystem()->Dimension();
}

std::shared_ptr<const CoordinateSystem> EmbeddedGeometry::GetParameterSystem()
  const noexcept
{
  return parameter_->GetCoordinateSystem();
}

std::shared_ptr<BoundPoints> EmbeddedGeometry::BindImpl(
  const BindRequest& request) const
{
  return std::make_shared<EmbeddedBinding>(
    std::static_pointer_cast<const EmbeddedGeometry>(shared_from_this()),
    *search_, request);
}

Quadrature EmbeddedGeometry::CreateQuadrature(
  const detail::QuadratureRule& rule) const
{
  return uniform_grid::VisitUniformGrid(
    parameter_,
    [&]<unsigned Dim>(
      const std::shared_ptr<const uniform_grid::UniformGridGeometry<Dim>>&
        grid) { return CellQuadrature<Dim>(*this, grid->GetGrid(), rule); });
}

std::shared_ptr<const CoordinateMap> CreateCoordinateMap(
  std::shared_ptr<const Function<Real>> position)
{
  return std::make_shared<const PositionMap>(
    Position(std::move(position), "CreateCoordinateMap"));
}

std::shared_ptr<const EmbeddedGeometry> CreateEmbeddedGeometry(
  std::shared_ptr<const Geometry> parameter_geometry,
  std::shared_ptr<const Function<Real>> position, EmbeddingOptions options)
{
  Position checked(std::move(position), "CreateEmbeddedGeometry");
  if (parameter_geometry == nullptr ||
      checked.Get()->GetSpace().GetGeometry() != parameter_geometry) {
    throw pcms_error("CreateEmbeddedGeometry: the position must be a Function "
                     "on the parameter geometry");
  }
  const auto& symmetry = options.symmetry;
  if (symmetry && (symmetry->length <= 0 || symmetry->axis < 0 ||
                   symmetry->axis >= checked.System()->Dimension())) {
    throw pcms_error("CreateEmbeddedGeometry: the symmetry must translate an "
                     "axis of the position's coordinate system by a positive "
                     "length");
  }
  if (detail::SameCoordinateSystem(checked.ParameterSystem(),
                                   checked.System())) {
    throw pcms_error(
      "CreateEmbeddedGeometry: the parameter geometry's coordinate system "
      "must differ from the coordinate system of the position (give the "
      "parameter geometry a parameter system)");
  }
  auto search = uniform_grid::VisitUniformGrid(
    parameter_geometry,
    [&]<unsigned Dim>(
      const std::shared_ptr<const uniform_grid::UniformGridGeometry<Dim>>& grid)
      -> std::shared_ptr<const detail::EmbeddedSearch> {
      if constexpr (Dim == 1) {
        return nullptr;
      } else {
        return CreateSearch<Dim>(grid->GetGrid(), checked, options);
      }
    });
  if (search == nullptr) {
    throw pcms_error("CreateEmbeddedGeometry: only 2D and 3D uniform-grid "
                     "parameter geometries are supported");
  }
  return std::shared_ptr<const EmbeddedGeometry>(new EmbeddedGeometry(
    std::move(parameter_geometry), checked.Get(), std::move(search)));
}

std::shared_ptr<const FunctionSpace> CreateEmbeddedSpace(
  std::shared_ptr<const EmbeddedGeometry> embedded,
  std::shared_ptr<const FunctionSpace> space)
{
  if (space == nullptr || embedded == nullptr ||
      embedded->GetParameterGeometry() != space->GetGeometry()) {
    throw pcms_error("CreateEmbeddedSpace: the geometry must be an embedded "
                     "geometry of the space's geometry");
  }
  auto layout = space->GetLayout();
  return std::make_shared<FunctionSpace>(
    std::move(layout),
    std::make_shared<EmbeddedEvaluatorFactory>(embedded, std::move(space)),
    embedded);
}

} // namespace pcms
