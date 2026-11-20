#include "pcms/transfer/l2_inner_product.hpp"
#include "pcms/transfer/forms.hpp"
#include "pcms/utility/assert.hpp"
#include <cmath>
#include <string>

namespace pcms
{

namespace
{

Quadrature MakeQuadrature(const FunctionSpace& space,
                          std::optional<Quadrature> quadrature)
{
  if (!quadrature) {
    return space.CreateQuadrature();
  }
  if (quadrature->points->GetGeometry() != space.GetGeometry()) {
    throw pcms_error("L2: the quadrature is not on the space's geometry");
  }
  return std::move(*quadrature);
}

// The inner product sums products of components, which is the metric inner
// product of vectors and tensors only where the metric is the identity; the
// pointwise metric of the other systems is not applied yet. Scalars are
// unaffected: the weights carry the volume factor.
void CheckInnerProductMetric(const Field<Real>& u, const char* caller)
{
  const ValueBasis& basis = u.GetData().GetValueBasis();
  if (basis.Rank() > 0 && !basis.system->IsOrthonormal()) {
    throw pcms_error(std::string(caller) +
                     ": the inner product of vectors in '" +
                     std::string(basis.system->Kind()) +
                     "' needs the metric, which is not applied yet");
  }
}

void Evaluate(const PointEvaluator<Real>& evaluator, const Field<Real>& field,
              const Kokkos::View<Real**, DeviceMemorySpace>& values)
{
  const ValueBasis basis =
    evaluator.OutputBasis(field.GetData().GetValueBasis());
  auto view = values;
  evaluator.Evaluate(
    field, ValueView<Real, DeviceMemorySpace>(basis, MakeRank2View(view)));
}

// sum_p w_p sum_c (a(p, c) - sign b(p, c))^2 style reductions
Real WeightedSum(const Kokkos::View<const Real*, DeviceMemorySpace>& weights,
                 const Kokkos::View<Real**, DeviceMemorySpace>& a,
                 const Kokkos::View<Real**, DeviceMemorySpace>& b,
                 bool difference_squared)
{
  const auto num_points = static_cast<LO>(a.extent(0));
  const auto num_components = static_cast<int>(a.extent(1));
  Real sum = 0;
  Kokkos::parallel_reduce(
    "pcms::L2::WeightedSum",
    Kokkos::RangePolicy<DefaultExecutionSpace>(0, num_points),
    KOKKOS_LAMBDA(const LO p, Real& partial) {
      Real local = 0;
      for (int c = 0; c < num_components; ++c) {
        const Real d = difference_squared ? a(p, c) - b(p, c) : 0;
        local += difference_squared ? d * d : a(p, c) * b(p, c);
      }
      partial += weights(p) * local;
    },
    sum);
  return sum;
}

} // namespace

L2InnerProduct::L2InnerProduct(std::shared_ptr<const FunctionSpace> u_space,
                               std::shared_ptr<const PointEvaluatorFactory> v,
                               std::optional<Quadrature> quadrature)
  : u_space_(std::move(u_space)),
    v_source_(std::move(v)),
    quadrature_(MakeQuadrature(*u_space_, std::move(quadrature)))
{
  if (v_source_ == nullptr || v_source_->GetGeometry() == nullptr) {
    throw pcms_error("L2InnerProduct: v must be an evaluator source with a "
                     "geometry");
  }
  u_at_points_ = u_space_->CreatePointEvaluator<Real>(quadrature_.points);
  v_at_points_ = v_source_->CreatePointEvaluator<Real>(
    transfer::detail::BindQuadrature(*v_source_->GetGeometry(), quadrature_));
  const auto num_points = quadrature_.points->NumPoints();
  const int num_components = u_space_->GetLayout()->GetNumComponents();
  u_values_ =
    Kokkos::View<Real**, DeviceMemorySpace>("l2_u", num_points, num_components);
  v_values_ =
    Kokkos::View<Real**, DeviceMemorySpace>("l2_v", num_points, num_components);
}

Real L2InnerProduct::Apply(const Field<Real>& u, const Field<Real>& v) const
{
  CheckInnerProductMetric(u, "L2InnerProduct::Apply");
  Evaluate(*u_at_points_, u, u_values_);
  Evaluate(*v_at_points_, v, v_values_);
  return WeightedSum(quadrature_.weights, u_values_, v_values_, false);
}

Real L2InnerProduct::Distance(const Field<Real>& u, const Field<Real>& v) const
{
  CheckInnerProductMetric(u, "L2InnerProduct::Distance");
  Evaluate(*u_at_points_, u, u_values_);
  Evaluate(*v_at_points_, v, v_values_);
  return std::sqrt(
    WeightedSum(quadrature_.weights, u_values_, v_values_, true));
}

L2Integral::L2Integral(std::shared_ptr<const FunctionSpace> space,
                       std::optional<Quadrature> quadrature)
  : space_(std::move(space)),
    quadrature_(MakeQuadrature(*space_, std::move(quadrature)))
{
  at_points_ = space_->CreatePointEvaluator<Real>(quadrature_.points);
  values_ = Kokkos::View<Real**, DeviceMemorySpace>(
    "l2_integral", quadrature_.points->NumPoints(),
    space_->GetLayout()->GetNumComponents());
}

Real L2Integral::Apply(const Field<Real>& u) const
{
  if (values_.extent(1) != 1) {
    throw pcms_error("L2Integral::Apply: only for scalar fields");
  }
  Evaluate(*at_points_, u, values_);
  const auto weights = quadrature_.weights;
  const auto values = values_;
  Real sum = 0;
  Kokkos::parallel_reduce(
    "pcms::L2Integral::Apply",
    Kokkos::RangePolicy<DefaultExecutionSpace>(
      0, static_cast<LO>(values.extent(0))),
    KOKKOS_LAMBDA(const LO p, Real& partial) {
      partial += weights(p) * values(p, 0);
    },
    sum);
  return sum;
}

Real L2Integral::Norm(const Field<Real>& u) const
{
  CheckInnerProductMetric(u, "L2Integral::Norm");
  Evaluate(*at_points_, u, values_);
  return std::sqrt(WeightedSum(quadrature_.weights, values_, values_, false));
}

} // namespace pcms
