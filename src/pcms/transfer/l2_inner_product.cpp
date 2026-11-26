#include "pcms/transfer/l2_inner_product.hpp"
#include "pcms/transfer/forms.hpp"
#include "pcms/utility/assert.hpp"
#include <cmath>
#include <string>

namespace pcms
{

namespace
{

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

// sum over the points of w |a - b|^2
Real SquaredDistance(const transfer::detail::BilinearForm& form,
                     const Field<Real>& u, const Field<Real>& v)
{
  const auto uu = form.EvaluateTest(u);
  const auto vv = form.EvaluateTrial(v);
  const auto a = uu.GetValues();
  const auto b = vv.GetValues();
  Kokkos::View<Real**, DeviceMemorySpace> d("l2_difference", a.extent(0),
                                            a.extent(1));
  Kokkos::parallel_for(
    "pcms::L2::Difference",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
      {0, 0}, {static_cast<LO>(a.extent(0)), static_cast<LO>(a.extent(1))}),
    KOKKOS_LAMBDA(const LO p, const LO c) { d(p, c) = a(p, c) - b(p, c); });
  const ValueView<const Real, DeviceMemorySpace> difference(
    uu.GetBasis(), MakeConstRank2View(d));
  return form.Contract(difference, difference);
}

} // namespace

class L2Integral::Form : public transfer::detail::BilinearForm
{
public:
  using BilinearForm::BilinearForm;
};

L2Integral::L2Integral(std::shared_ptr<const FunctionSpace> space,
                       std::shared_ptr<const PointEvaluatorFactory> other,
                       std::optional<Quadrature> quadrature)
{
  if (other != nullptr && other->GetGeometry() == nullptr) {
    throw pcms_error("L2Integral: `other` must have a geometry");
  }
  auto trial = other != nullptr ? std::move(other) : space;
  form_ = std::make_unique<const Form>(
    std::move(space), std::move(trial),
    transfer::detail::FormOptions{.quadrature = std::move(quadrature)});
}

L2Integral::~L2Integral() = default;

Real L2Integral::Apply(const Field<Real>& u) const
{
  const auto values = form_->EvaluateTest(u);
  if (values.extent(1) != 1) {
    throw pcms_error("L2Integral::Apply: only for scalar fields");
  }
  Kokkos::View<Real**, DeviceMemorySpace> ones("l2_ones", values.extent(0), 1);
  Kokkos::deep_copy(ones, 1.0);
  return form_->Contract(
    values, ValueView<const Real, DeviceMemorySpace>(values.GetBasis(),
                                                     MakeConstRank2View(ones)));
}

Real L2Integral::Inner(const Field<Real>& u, const Field<Real>& v) const
{
  CheckInnerProductMetric(u, "L2Integral::Inner");
  return form_->Pair(u, v);
}

Real L2Integral::Distance(const Field<Real>& u, const Field<Real>& v) const
{
  CheckInnerProductMetric(u, "L2Integral::Distance");
  return std::sqrt(SquaredDistance(*form_, u, v));
}

Real L2Integral::Norm(const Field<Real>& u) const
{
  CheckInnerProductMetric(u, "L2Integral::Norm");
  const auto values = form_->EvaluateTest(u);
  return std::sqrt(form_->Contract(values, values));
}

} // namespace pcms
