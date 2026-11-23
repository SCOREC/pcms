#include "pcms/transfer/l2_inner_product.hpp"
#include "pcms/field/detail/field_ops.hpp"
#include "pcms/geometry/detail/geometry.hpp"
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
  return quadrature ? std::move(*quadrature) : space.CreateQuadrature();
}

// The inner product sums products of components, which is the metric inner
// product of vectors and tensors only where the metric is the identity.
// Scalars are unaffected: the weights carry the volume factor.
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

ValueView<const Real, DeviceMemorySpace> Evaluate(
  const PointEvaluator<Real>& evaluator, const Field<Real>& field,
  const Kokkos::View<Real**, DeviceMemorySpace>& values)
{
  const ValueBasis basis =
    evaluator.OutputBasis(field.GetData().GetValueBasis());
  auto view = values;
  evaluator.Evaluate(
    field, ValueView<Real, DeviceMemorySpace>(basis, MakeRank2View(view)));
  return {basis, MakeConstRank2View(values)};
}

} // namespace

class L2Integral::Form : public transfer::detail::LinearForm
{
public:
  using LinearForm::LinearForm;
};

L2Integral::L2Integral(std::shared_ptr<const FunctionSpace> space,
                       std::shared_ptr<const PointEvaluatorFactory> other,
                       std::optional<Quadrature> quadrature)
  : form_(std::make_unique<const Form>(
      space, MakeQuadrature(*space, std::move(quadrature))))
{
  const Quadrature& rule = form_->GetQuadrature();
  if (other != nullptr) {
    if (other->GetGeometry() == nullptr) {
      throw pcms_error("L2Integral: `other` must have a geometry");
    }
    other_ = other->CreatePointEvaluator<Real>(
      pcms::detail::BindOn(rule, *other->GetGeometry()));
  }
  const auto num_points = rule.points->NumPoints();
  const int num_components = space->GetLayout()->GetNumComponents();
  u_values_ =
    Kokkos::View<Real**, DeviceMemorySpace>("l2_u", num_points, num_components);
  v_values_ =
    Kokkos::View<Real**, DeviceMemorySpace>("l2_v", num_points, num_components);
}

L2Integral::~L2Integral() = default;

Real L2Integral::Apply(const Field<Real>& u) const
{
  if (u_values_.extent(1) != 1) {
    throw pcms_error("L2Integral::Apply: only for scalar fields");
  }
  const auto wu =
    form_->Weight(Evaluate(form_->Evaluator(), u, u_values_)).GetValues();
  return pcms::detail::SumOwned<Real>(
    *form_->Space()->GetLayout(), wu, wu,
    KOKKOS_LAMBDA(const Real a, const Real) { return a; });
}

Real L2Integral::Inner(const Field<Real>& u, const Field<Real>& v) const
{
  CheckInnerProductMetric(u, "L2Integral::Inner");
  const auto wu =
    form_->Weight(Evaluate(form_->Evaluator(), u, u_values_)).GetValues();
  const auto vv =
    Evaluate(other_ ? *other_ : form_->Evaluator(), v, v_values_).GetValues();
  return pcms::detail::SumOwned<Real>(
    *form_->Space()->GetLayout(), wu, vv,
    KOKKOS_LAMBDA(const Real a, const Real b) { return a * b; });
}

Real L2Integral::Distance(const Field<Real>& u, const Field<Real>& v) const
{
  CheckInnerProductMetric(u, "L2Integral::Distance");
  const auto uu = Evaluate(form_->Evaluator(), u, u_values_);
  const auto vv = Evaluate(other_ ? *other_ : form_->Evaluator(), v, v_values_);
  pcms::detail::CheckWrittenValueBasis("L2Integral::Distance", uu.GetBasis(),
                                       vv.GetBasis());
  const auto difference = pcms::detail::MapOwned<Real>(
    uu.GetValues(), vv.GetValues(),
    KOKKOS_LAMBDA(const Real a, const Real b) { return a - b; });
  const auto d = MakeConstRank2View(difference);
  const auto wd =
    form_->Weight(ValueView<const Real, DeviceMemorySpace>(uu.GetBasis(), d))
      .GetValues();
  return std::sqrt(pcms::detail::SumOwned<Real>(
    *form_->Space()->GetLayout(), wd, d,
    KOKKOS_LAMBDA(const Real a, const Real b) { return a * b; }));
}

Real L2Integral::Norm(const Field<Real>& u) const
{
  CheckInnerProductMetric(u, "L2Integral::Norm");
  const auto uu = Evaluate(form_->Evaluator(), u, u_values_);
  const auto wu = form_->Weight(uu).GetValues();
  return std::sqrt(pcms::detail::SumOwned<Real>(
    *form_->Space()->GetLayout(), wu, uu.GetValues(),
    KOKKOS_LAMBDA(const Real a, const Real b) { return a * b; }));
}

} // namespace pcms
