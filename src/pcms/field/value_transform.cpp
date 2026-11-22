#include "pcms/field/value_transform.hpp"
#include "pcms/utility/dense_linear_algebra.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include "pcms/utility/detail/profile.hpp"
#include <algorithm>
#include <string>
#include <utility>

namespace pcms
{

namespace
{

// The largest coordinate system dimension a transform handles.
constexpr int MaxDimension = 3;

using Matrices = Kokkos::View<Real***, DeviceMemorySpace>;

// out = M in, or M^T in, at every point, block by block. A free function: CUDA
// does not allow device lambdas in private members.
void Multiply(Kokkos::View<const Real***, DeviceMemorySpace> m, bool transpose,
              Rank2View<const Real, DeviceMemorySpace> in,
              Rank2View<Real, DeviceMemorySpace> out)
{
  const int rows = static_cast<int>(transpose ? m.extent(2) : m.extent(1));
  const int columns = static_cast<int>(transpose ? m.extent(1) : m.extent(2));
  const int blocks = static_cast<int>(in.extent(1)) / columns;
  detail::ForEach(
    "ValueTransform", static_cast<LO>(in.extent(0)), KOKKOS_LAMBDA(const LO p) {
      for (int b = 0; b < blocks; ++b) {
        Real x[MaxDimension] = {};
        for (int c = 0; c < columns; ++c) {
          x[c] = in(p, b * columns + c);
        }
        for (int r = 0; r < rows; ++r) {
          Real y = 0;
          for (int c = 0; c < columns; ++c) {
            y += (transpose ? m(p, c, r) : m(p, r, c)) * x[c];
          }
          out(p, b * rows + r) = y;
        }
      }
    });
}

} // namespace

ValueTransform::ValueTransform(Matrices vector, Matrices covector,
                               std::shared_ptr<const CoordinateSystem> from,
                               std::shared_ptr<const CoordinateSystem> to,
                               std::optional<ComponentScaling> scaling)
  : vector_(std::move(vector)),
    covector_(std::move(covector)),
    from_(std::move(from)),
    to_(std::move(to)),
    scaling_(scaling),
    num_points_(
      static_cast<LO>(std::max(vector_.extent(0), covector_.extent(0))))
{
}

bool ValueTransform::PassesThrough(const ValueBasis& in) const
{
  return in.Rank() == 0 || !detail::SameCoordinateSystem(in.system, from_) ||
         (scaling_ && in.scaling == scaling_);
}

ValueBasis ValueTransform::OutputBasis(const ValueBasis& in) const
{
  if (in.Rank() > 1) {
    throw pcms_error("ValueTransform::OutputBasis: values of rank " +
                     std::to_string(in.Rank()) + " cannot be transformed yet");
  }
  if (PassesThrough(in)) {
    return in;
  }
  const bool vector = in.variance[0] == Variance::Contravariant;
  if ((vector ? vector_ : covector_).extent(1) == 0) {
    throw pcms_error(std::string("ValueTransform::OutputBasis: ") +
                     (vector ? "vectors" : "covectors") +
                     " need the inverse of the Jacobian here, which is not "
                     "square (" +
                     std::to_string(from_->Dimension()) + " x " +
                     std::to_string(to_->Dimension()) + ")");
  }
  if (!scaling_ && in.scaling == ComponentScaling::Physical &&
      !from_->IsOrthonormal()) {
    throw pcms_error("ValueTransform::OutputBasis: this transform takes "
                     "natural components in '" +
                     std::string(from_->Kind()) + "', not physical ones");
  }
  const ComponentScaling scaling =
    scaling_.value_or(to_->IsOrthonormal() ? ComponentScaling::Physical
                                           : ComponentScaling::Natural);
  return ValueBasis{to_,       scaling,           in.variance,
                    in.parity, ValueKind::Tensor, in.units};
}

void ValueTransform::Run(const char* who, const ValueBasis& in_side,
                         bool transpose,
                         Rank2View<const Real, DeviceMemorySpace> in,
                         Rank2View<Real, DeviceMemorySpace> out) const
{
  if (in.extent(0) != static_cast<size_t>(num_points_) ||
      out.extent(0) != in.extent(0)) {
    throw pcms_error(std::string(who) +
                     ": the values must have one row per point of the "
                     "transform");
  }
  if (PassesThrough(in_side)) {
    detail::CopyDeviceRank2ViewToRank2View<Real>(out, in);
    return;
  }
  const int din = transpose ? to_->Dimension() : from_->Dimension();
  const int dout = transpose ? from_->Dimension() : to_->Dimension();
  if (in.extent(1) % din != 0 || out.extent(1) != in.extent(1) / din * dout) {
    throw pcms_error(std::string(who) +
                     ": the input and output have numbers of components "
                     "that do not match");
  }
  Multiply(in_side.variance[0] == Variance::Contravariant ? vector_ : covector_,
           transpose, in, out);
}

void ValueTransform::Apply(ValueView<const Real, DeviceMemorySpace> in,
                           ValueView<Real, DeviceMemorySpace> out) const
{
  PCMS_FUNCTION_TIMER;
  detail::CheckWrittenValueBasis("ValueTransform::Apply", out.GetBasis(),
                                 OutputBasis(in.GetBasis()));
  Run("ValueTransform::Apply", in.GetBasis(), false, in.GetValues(),
      out.GetValues());
}

void ValueTransform::ApplyTranspose(
  ValueView<const Real, DeviceMemorySpace> in,
  ValueView<Real, DeviceMemorySpace> out) const
{
  PCMS_FUNCTION_TIMER;
  detail::CheckWrittenValueBasis("ValueTransform::ApplyTranspose",
                                 in.GetBasis(), OutputBasis(out.GetBasis()));
  Run("ValueTransform::ApplyTranspose", out.GetBasis(), true, in.GetValues(),
      out.GetValues());
}

// Creates the transforms of the factories once their kernels, which read the
// caller's Jacobian or metric, are complete.
struct ValueTransformAccess
{
  static ValueTransform Create(Matrices vector, Matrices covector,
                               std::shared_ptr<const CoordinateSystem> from,
                               std::shared_ptr<const CoordinateSystem> to,
                               std::optional<ComponentScaling> scaling = {})
  {
    Kokkos::fence("ValueTransform: matrices formed");
    return ValueTransform(std::move(vector), std::move(covector),
                          std::move(from), std::move(to), scaling);
  }
};

namespace detail
{

namespace
{

void CheckDimension(const CoordinateSystem& system)
{
  if (system.Dimension() > MaxDimension) {
    throw pcms_error("ValueTransform: coordinate systems of dimension above 3 "
                     "are not supported");
  }
}

// J or J^T (`transposed`), and J^-1 or J^-T (the other way; empty unless J
// is square).
std::pair<Matrices, Matrices> MapMatrices(const JacobianView& jacobian,
                                          const CoordinateSystem& image,
                                          const CoordinateSystem& domain,
                                          bool transposed)
{
  CheckDimension(image);
  CheckDimension(domain);
  const auto n = jacobian.extent(0);
  const int rows = image.Dimension();
  const int columns = domain.Dimension();
  if (jacobian.extent(1) != static_cast<size_t>(rows) ||
      jacobian.extent(2) != static_cast<size_t>(columns)) {
    throw pcms_error("ValueTransform: the Jacobian must be shaped "
                     "[points][dimension of its image][dimension of its "
                     "domain]");
  }
  Matrices direct("direct", n, transposed ? columns : rows,
                  transposed ? rows : columns);
  Matrices inverse;
  if (rows == columns) {
    inverse = Matrices("inverse", n, rows, rows);
  }
  ForEach(
    "ValueTransform::MapMatrices", static_cast<LO>(n),
    KOKKOS_LAMBDA(const LO p) {
      Real m[3][3] = {};
      for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < columns; ++c) {
          m[r][c] = jacobian(p, r, c);
          (transposed ? direct(p, c, r) : direct(p, r, c)) = m[r][c];
        }
      }
      if (rows == columns) {
        Real inv[3][3] = {};
        Invert(m, rows, inv);
        for (int r = 0; r < rows; ++r) {
          for (int c = 0; c < rows; ++c) {
            inverse(p, r, c) = transposed ? inv[r][c] : inv[c][r];
          }
        }
      }
    });
  return {direct, inverse};
}

// A zero scale factor (a singular point of the coordinate system) counts as 1.
KOKKOS_INLINE_FUNCTION Real SafeFactor(Real h)
{
  return h > 0 ? h : Real(1);
}

// The diagonal transform to `to` components: vectors by sqrt(g_ii), covectors
// by sqrt(g^ii) (or their inverses to natural).
ValueTransform Diagonal(std::shared_ptr<const CoordinateSystem> system,
                        MetricView metric, ComponentScaling to)
{
  CheckDimension(*system);
  const int dim = system->Dimension();
  const auto n = metric.extent(0);
  const bool physical = to == ComponentScaling::Physical;
  Matrices vector("vector", n, dim, dim);
  Matrices covector("covector", n, dim, dim);
  ForEach(
    "ValueTransform::Diagonal", static_cast<LO>(n), KOKKOS_LAMBDA(const LO p) {
      Real g[3][3] = {};
      for (int k = 0; k < dim; ++k) {
        for (int l = 0; l < dim; ++l) {
          g[k][l] = metric(p, k, l);
        }
      }
      Real inverse[3][3] = {};
      const bool singular = !(Determinant(g, dim) > 0);
      if (!singular) {
        Invert(g, dim, inverse);
      }
      for (int k = 0; k < dim; ++k) {
        const Real h = SafeFactor(Kokkos::sqrt(g[k][k]));
        const Real h_dual =
          singular ? Real(1) : SafeFactor(Kokkos::sqrt(inverse[k][k]));
        vector(p, k, k) = physical ? h : 1 / h;
        covector(p, k, k) = physical ? h_dual : 1 / h_dual;
      }
    });
  return ValueTransformAccess::Create(vector, covector, system, system, to);
}

} // namespace

ValueTransform Pullback(JacobianView jacobian,
                        std::shared_ptr<const CoordinateSystem> from,
                        std::shared_ptr<const CoordinateSystem> to)
{
  auto [covector, vector] = MapMatrices(jacobian, *from, *to, true);
  return ValueTransformAccess::Create(vector, covector, std::move(from),
                                      std::move(to));
}

ValueTransform Pushforward(JacobianView jacobian,
                           std::shared_ptr<const CoordinateSystem> from,
                           std::shared_ptr<const CoordinateSystem> to)
{
  auto [vector, covector] = MapMatrices(jacobian, *to, *from, false);
  return ValueTransformAccess::Create(vector, covector, std::move(from),
                                      std::move(to));
}

ValueTransform ToNatural(std::shared_ptr<const CoordinateSystem> system,
                         MetricView metric)
{
  return Diagonal(std::move(system), metric, ComponentScaling::Natural);
}

ValueTransform ToPhysical(std::shared_ptr<const CoordinateSystem> system,
                          MetricView metric)
{
  return Diagonal(std::move(system), metric, ComponentScaling::Physical);
}

Kokkos::View<Real***, DeviceMemorySpace> PullbackMetric(JacobianView jacobian,
                                                        MetricView metric)
{
  const auto n = jacobian.extent(0);
  const int dg = static_cast<int>(jacobian.extent(1));
  const int dim = static_cast<int>(jacobian.extent(2));
  Kokkos::View<Real***, DeviceMemorySpace> pulled("pulled_back_metric", n, dim,
                                                  dim);
  ForEach(
    "PullbackMetric", static_cast<LO>(n), KOKKOS_LAMBDA(const LO p) {
      for (int k = 0; k < dim; ++k) {
        for (int l = 0; l < dim; ++l) {
          for (int a = 0; a < dg; ++a) {
            for (int b = 0; b < dg; ++b) {
              pulled(p, k, l) +=
                jacobian(p, a, k) * metric(p, a, b) * jacobian(p, b, l);
            }
          }
        }
      }
    });
  Kokkos::fence("PullbackMetric");
  return pulled;
}

} // namespace detail

} // namespace pcms
