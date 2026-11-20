#include "pcms/field/value_transform.hpp"
#include "pcms/utility/dense_linear_algebra.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include "pcms/utility/detail/profile.hpp"
#include <string>
#include <utility>

namespace pcms
{

namespace
{

// The largest coordinate system dimension a transform handles.
constexpr int MaxDimension = 3;

// What happens to the components on one side of the matrix: nothing, or a
// multiplication or division by the scale factor of their axis.
enum class Scaling : std::uint8_t
{
  None,
  Multiply,
  Divide
};

// The matrix applied at each point, with Jm(a, c) = J(p, c, a) (a indexes the
// map's domain, c its image): the identity,
// y = Jm x or y = Jm^T x (multiply), or the solution of Jm y = x or
// Jm^T y = x (solve).
struct Operation
{
  bool identity = false;
  bool solve = false;
  bool transposed = false;
  Scaling pre = Scaling::None;
  Scaling post = Scaling::None;
};

// A zero scale factor (a singular point of the coordinate system) counts as 1.
KOKKOS_INLINE_FUNCTION Real SafeFactor(Real h)
{
  return h > 0 ? h : Real(1);
}

KOKKOS_INLINE_FUNCTION void ApplyScaling(
  Scaling scaling, Kokkos::View<const Real**, DeviceMemorySpace> factors, LO p,
  int n, Real (&v)[MaxDimension])
{
  if (scaling == Scaling::Multiply) {
    for (int k = 0; k < n; ++k) {
      v[k] *= SafeFactor(factors(p, k));
    }
  } else if (scaling == Scaling::Divide) {
    for (int k = 0; k < n; ++k) {
      v[k] /= SafeFactor(factors(p, k));
    }
  }
}

// out = op(in) at every point, block by block: `din` input and `dout` output
// components per block. A free function: CUDA does not allow device lambdas
// in private members.
void RunOperation(Operation op, detail::JacobianView jacobian,
                  Kokkos::View<const Real**, DeviceMemorySpace> pre_factors,
                  Kokkos::View<const Real**, DeviceMemorySpace> post_factors,
                  Rank2View<const Real, DeviceMemorySpace> in,
                  Rank2View<Real, DeviceMemorySpace> out, int din, int dout)
{
  const int blocks = static_cast<int>(in.extent(1)) / din;
  const int e1 = op.identity ? 0 : static_cast<int>(jacobian.extent(2));
  const int e2 = op.identity ? 0 : static_cast<int>(jacobian.extent(1));
  detail::ForEach(
    "ValueTransform", static_cast<LO>(in.extent(0)), KOKKOS_LAMBDA(const LO p) {
      for (int b = 0; b < blocks; ++b) {
        Real x[MaxDimension] = {};
        Real y[MaxDimension] = {};
        for (int k = 0; k < din; ++k) {
          x[k] = in(p, b * din + k);
        }
        ApplyScaling(op.pre, pre_factors, p, din, x);
        if (op.identity) {
          for (int k = 0; k < din; ++k) {
            y[k] = x[k];
          }
        } else if (!op.solve && !op.transposed) {
          for (int a = 0; a < e1; ++a) {
            for (int c = 0; c < e2; ++c) {
              y[a] += jacobian(p, c, a) * x[c];
            }
          }
        } else if (!op.solve) {
          for (int c = 0; c < e2; ++c) {
            for (int a = 0; a < e1; ++a) {
              y[c] += jacobian(p, c, a) * x[a];
            }
          }
        } else {
          Real m[3][3] = {};
          for (int a = 0; a < e1; ++a) {
            for (int c = 0; c < e1; ++c) {
              m[a][c] = op.transposed ? jacobian(p, a, c) : jacobian(p, c, a);
            }
          }
          detail::Solve(m, e1, x);
          for (int k = 0; k < e1; ++k) {
            y[k] = x[k];
          }
        }
        ApplyScaling(op.post, post_factors, p, dout, y);
        for (int k = 0; k < dout; ++k) {
          out(p, b * dout + k) = y[k];
        }
      }
    });
}

void CopyValues(Rank2View<const Real, DeviceMemorySpace> in,
                Rank2View<Real, DeviceMemorySpace> out)
{
  if (in.data_handle() == out.data_handle()) {
    return;
  }
  const int columns = static_cast<int>(in.extent(1));
  detail::ForEach(
    "ValueTransform::Copy", static_cast<LO>(in.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      for (int c = 0; c < columns; ++c) {
        out(p, c) = in(p, c);
      }
    });
}

// The scale factors sqrt(g_ii) of a coordinate system with them
// (HasScaleFactors) at the points, from its metric there; empty for another
// coordinate system or without a metric. With the Jacobian J(p, j, i) = d y_j /
// d x_i of a map x -> y of the points, `metric` is the metric at their images y
// and the coordinate system's is J^T g J.
Kokkos::View<const Real**, DeviceMemorySpace> ScaleFactors(
  const CoordinateSystem& system, detail::MetricView metric,
  detail::JacobianView jacobian = {})
{
  if (!detail::HasScaleFactors(system) || metric.extent(0) == 0) {
    return {};
  }
  const int dim = system.Dimension();
  const int dm = static_cast<int>(metric.extent(1));
  const bool pulled = jacobian.extent(0) > 0;
  Kokkos::View<Real**, DeviceMemorySpace> factors(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "scale_factors"),
    metric.extent(0), dim);
  detail::ForEach(
    "ValueTransform::ScaleFactors", static_cast<LO>(factors.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      for (int k = 0; k < dim; ++k) {
        Real g = pulled ? Real(0) : metric(p, k, k);
        for (int a = 0; pulled && a < dm; ++a) {
          for (int b = 0; b < dm; ++b) {
            g += jacobian(p, a, k) * metric(p, a, b) * jacobian(p, b, k);
          }
        }
        factors(p, k) = Kokkos::sqrt(g);
      }
    });
  // the metric is the caller's, which may release it on return
  Kokkos::fence("ValueTransform::ScaleFactors");
  return factors;
}

void CheckRows(const char* who, size_t rows, LO num_points)
{
  if (rows != static_cast<size_t>(num_points)) {
    throw pcms_error(std::string(who) +
                     ": the values must have one row per point of the "
                     "transform");
  }
}

} // namespace

ValueTransform::ValueTransform(
  Direction direction,
  View<3, const Real, DeviceMemorySpace, Kokkos::layout_stride> jacobian,
  std::shared_ptr<const CoordinateSystem> from,
  Kokkos::View<const Real**, DeviceMemorySpace> from_scale,
  std::shared_ptr<const CoordinateSystem> to,
  Kokkos::View<const Real**, DeviceMemorySpace> to_scale,
  ComponentScaling scale_to)
  : direction_(direction),
    jacobian_(jacobian),
    from_(std::move(from)),
    to_(std::move(to)),
    from_scale_(std::move(from_scale)),
    to_scale_(std::move(to_scale)),
    scale_to_(scale_to),
    num_points_(direction == Direction::None
                  ? static_cast<LO>(from_scale_.extent(0))
                  : static_cast<LO>(jacobian.extent(0)))
{
  if (from_->Dimension() > MaxDimension || to_->Dimension() > MaxDimension) {
    throw pcms_error(
      "ValueTransform: coordinate systems of dimension above 3 are not "
      "supported");
  }
}

bool ValueTransform::NeedsInverse(Variance variance) const noexcept
{
  // vectors with the map and covectors against it need no inverse
  return direction_ != Direction::None &&
         (direction_ == Direction::With) !=
           (variance == Variance::Contravariant);
}

namespace
{

// Whether values in `in` leave a transform unchanged: scalars and
// coordinates, values already in the output coordinate system of a map (in
// any coordinate system but the parameter system along a map), values already
// in the output scaling of a Scale.
bool PassesThrough(const ValueBasis& in, bool scale, bool along,
                   const std::shared_ptr<const CoordinateSystem>& from,
                   const std::shared_ptr<const CoordinateSystem>& to,
                   ComponentScaling scale_to)
{
  if (in.Rank() == 0) {
    return true;
  }
  if (scale) {
    return detail::SameCoordinateSystem(in.system, from) &&
           in.scaling == scale_to;
  }
  return (along || detail::SameCoordinateSystem(in.system, to)) &&
         !detail::SameCoordinateSystem(in.system, from);
}

} // namespace

ValueBasis ValueTransform::OutputBasis(const ValueBasis& in) const
{
  const bool scale = direction_ == Direction::None;
  if (in.Rank() > 1) {
    throw pcms_error("ValueTransform::OutputBasis: values of rank " +
                     std::to_string(in.Rank()) +
                     " cannot be transformed yet (scalars, vectors and "
                     "covectors can)");
  }
  if (PassesThrough(in, scale, direction_ == Direction::With, from_, to_,
                    scale_to_)) {
    return in;
  }
  if (!detail::SameCoordinateSystem(in.system, from_)) {
    throw pcms_error("ValueTransform::OutputBasis: the values are components "
                     "in '" +
                     std::string(in.system->Kind()) +
                     "'; this transform takes components in '" +
                     std::string(from_->Kind()) + "'");
  }
  if (NeedsInverse(in.variance[0]) &&
      jacobian_.extent(1) != jacobian_.extent(2)) {
    throw pcms_error(
      std::string("ValueTransform::OutputBasis: ") +
      (in.variance[0] == Variance::Contravariant ? "vectors" : "covectors") +
      " need the inverse of the Jacobian here, which is not square (" +
      std::to_string(jacobian_.extent(1)) + " x " +
      std::to_string(jacobian_.extent(2)) + ")");
  }
  ComponentScaling scaling = ComponentScaling::Natural;
  if (scale) {
    scaling = scale_to_;
  } else if (to_->IsOrthogonal() && (in.scaling == ComponentScaling::Physical ||
                                     to_->IsNormalized())) {
    scaling = ComponentScaling::Physical;
  }
  return ValueBasis{to_, scaling, in.variance, in.parity};
}

void ValueTransform::Run(const ValueBasis& in_side,
                         ComponentScaling out_scaling, bool transpose,
                         Rank2View<const Real, DeviceMemorySpace> in,
                         Rank2View<Real, DeviceMemorySpace> out) const
{
  const bool vector = in_side.variance[0] == Variance::Contravariant;
  Operation op;
  op.identity = direction_ == Direction::None;
  op.solve = NeedsInverse(in_side.variance[0]);
  op.transposed = vector != transpose;
  // physical components of a vector are the natural ones times the scale
  // factor, of a covector divided by it
  const bool in_physical =
    in_side.scaling == ComponentScaling::Physical && from_scale_.extent(1) > 0;
  const bool out_physical =
    out_scaling == ComponentScaling::Physical && to_scale_.extent(1) > 0;
  const Scaling to_natural = vector ? Scaling::Divide : Scaling::Multiply;
  const Scaling to_physical = vector ? Scaling::Multiply : Scaling::Divide;
  const int din = from_->Dimension();
  const int dout = to_->Dimension();
  if (!transpose) {
    // T = S_out M S_in^-1
    op.pre = in_physical ? to_natural : Scaling::None;
    op.post = out_physical ? to_physical : Scaling::None;
    if (in.extent(1) % din != 0 || out.extent(1) != in.extent(1) / din * dout) {
      throw pcms_error("ValueTransform::Apply: the input and output have "
                       "numbers of components that do not match");
    }
    RunOperation(op, jacobian_, from_scale_, to_scale_, in, out, din, dout);
  } else {
    // T^T = S_in^-1 M^T S_out, for the diagonal S of the forward transform
    op.pre = out_physical ? to_physical : Scaling::None;
    op.post = in_physical ? to_natural : Scaling::None;
    if (in.extent(1) % dout != 0 ||
        out.extent(1) != in.extent(1) / dout * din) {
      throw pcms_error("ValueTransform::ApplyTranspose: the input and output "
                       "have numbers of components that do not match");
    }
    RunOperation(op, jacobian_, to_scale_, from_scale_, in, out, dout, din);
  }
}

void ValueTransform::Apply(ValueView<const Real, DeviceMemorySpace> in,
                           ValueView<Real, DeviceMemorySpace> out) const
{
  PCMS_FUNCTION_TIMER;
  detail::CheckWrittenValueBasis("ValueTransform::Apply", out.GetBasis(),
                                 OutputBasis(in.GetBasis()));
  CheckRows("ValueTransform::Apply", in.extent(0), num_points_);
  CheckRows("ValueTransform::Apply", out.extent(0), num_points_);
  if (PassesThrough(in.GetBasis(), direction_ == Direction::None,
                    direction_ == Direction::With, from_, to_, scale_to_)) {
    CopyValues(in.GetValues(), out.GetValues());
    return;
  }
  Run(in.GetBasis(),
      out.GetBasis().scaling.value_or(ComponentScaling::Physical), false,
      in.GetValues(), out.GetValues());
}

void ValueTransform::ApplyTranspose(
  ValueView<const Real, DeviceMemorySpace> in,
  ValueView<Real, DeviceMemorySpace> out) const
{
  PCMS_FUNCTION_TIMER;
  detail::CheckWrittenValueBasis("ValueTransform::ApplyTranspose",
                                 in.GetBasis(), OutputBasis(out.GetBasis()));
  CheckRows("ValueTransform::ApplyTranspose", in.extent(0), num_points_);
  CheckRows("ValueTransform::ApplyTranspose", out.extent(0), num_points_);
  if (PassesThrough(out.GetBasis(), direction_ == Direction::None,
                    direction_ == Direction::With, from_, to_, scale_to_)) {
    CopyValues(in.GetValues(), out.GetValues());
    return;
  }
  Run(out.GetBasis(),
      in.GetBasis().scaling.value_or(ComponentScaling::Physical), true,
      in.GetValues(), out.GetValues());
}

// Creates the transforms of the factories: against a map, along it, or a
// scaling.
struct ValueTransformAccess
{
  static ValueTransform Create(
    bool along, detail::JacobianView jacobian,
    std::shared_ptr<const CoordinateSystem> from,
    Kokkos::View<const Real**, DeviceMemorySpace> from_scale,
    std::shared_ptr<const CoordinateSystem> to,
    Kokkos::View<const Real**, DeviceMemorySpace> to_scale)
  {
    return ValueTransform(along ? ValueTransform::Direction::With
                                : ValueTransform::Direction::Against,
                          jacobian, std::move(from), std::move(from_scale),
                          std::move(to), std::move(to_scale),
                          ComponentScaling::Natural);
  }

  static ValueTransform CreateScale(
    std::shared_ptr<const CoordinateSystem> system,
    Kokkos::View<const Real**, DeviceMemorySpace> factors, ComponentScaling to)
  {
    return ValueTransform(ValueTransform::Direction::None,
                          detail::JacobianView{}, system, factors, system,
                          factors, to);
  }
};

namespace detail
{

namespace
{

void CheckJacobian(const char* who, const JacobianView& jacobian,
                   const CoordinateSystem& domain,
                   const CoordinateSystem& image, const MetricView& metric)
{
  const auto n = jacobian.extent(0);
  const auto dim = static_cast<size_t>(image.Dimension());
  if (jacobian.extent(1) != dim ||
      jacobian.extent(2) != static_cast<size_t>(domain.Dimension()) ||
      (metric.extent(0) != 0 &&
       (metric.extent(0) != n || metric.extent(1) != dim ||
        metric.extent(2) != dim))) {
    throw pcms_error(std::string(who) +
                     ": the Jacobian must be shaped [points][dimension of "
                     "its image][dimension of its domain], and the metric "
                     "[points][dimension of the image][same]");
  }
}

} // namespace

ValueTransform Pullback(JacobianView jacobian,
                        std::shared_ptr<const CoordinateSystem> from,
                        MetricView metric,
                        std::shared_ptr<const CoordinateSystem> to)
{
  CheckJacobian("Pullback", jacobian, *to, *from, metric);
  auto from_scale = ScaleFactors(*from, metric);
  auto to_scale = ScaleFactors(*to, metric, jacobian);
  return ValueTransformAccess::Create(false, jacobian, std::move(from),
                                      std::move(from_scale), std::move(to),
                                      std::move(to_scale));
}

ValueTransform Pushforward(JacobianView jacobian,
                           std::shared_ptr<const CoordinateSystem> from,
                           std::shared_ptr<const CoordinateSystem> to,
                           MetricView metric)
{
  CheckJacobian("Pushforward", jacobian, *from, *to, metric);
  auto from_scale = ScaleFactors(*from, metric, jacobian);
  auto to_scale = ScaleFactors(*to, metric);
  return ValueTransformAccess::Create(true, jacobian, std::move(from),
                                      std::move(from_scale), std::move(to),
                                      std::move(to_scale));
}

ValueTransform Scale(std::shared_ptr<const CoordinateSystem> system,
                     MetricView metric, ComponentScaling to)
{
  if (to == ComponentScaling::Physical && !system->IsOrthogonal()) {
    throw pcms_error(
      "Scale: physical components need an orthogonal coordinate system; "
      "'" +
      std::string(system->Kind()) + "' is not");
  }
  auto factors = ScaleFactors(*system, metric);
  if (factors.extent(0) == 0) {
    // physical and natural components coincide: an identity of as many points
    factors = Kokkos::View<const Real**, DeviceMemorySpace>(
      Kokkos::View<Real**, DeviceMemorySpace>("unit_scale_factors",
                                              metric.extent(0), 0));
  }
  return ValueTransformAccess::CreateScale(std::move(system), factors, to);
}

} // namespace detail

} // namespace pcms
