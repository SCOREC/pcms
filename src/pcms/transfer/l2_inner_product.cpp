#include "pcms/transfer/l2_inner_product.hpp"
#include "pcms/field/detail/field_ops.hpp"
#include "pcms/transfer/forms.hpp"
#include "pcms/utility/assert.hpp"
#include <cmath>

namespace pcms
{

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
    transfer::detail::FormOptions{.quadrature = std::move(quadrature),
                                  .metric = true});
}

L2Integral::~L2Integral() = default;

Real L2Integral::Apply(const Field<Real>& u) const
{
  const auto values = EvaluateTest(*form_, u);
  if (values.extent(1) != 1) {
    throw pcms_error("L2Integral::Apply: only for scalar fields");
  }
  Kokkos::View<Real**, DeviceMemorySpace> ones("l2_ones", values.extent(0), 1);
  Kokkos::deep_copy(ones, 1.0);
  return Contract(*form_, values,
                  ValueView<const Real, DeviceMemorySpace>(
                    values.GetBasis(), MakeConstRank2View(ones)));
}

Real L2Integral::Inner(const Field<Real>& u, const Field<Real>& v) const
{
  return form_->Pair(u, v);
}

Real L2Integral::Distance(const Field<Real>& u, const Field<Real>& v) const
{
  const auto uu = EvaluateTest(*form_, u);
  const auto vv = EvaluateTrial(*form_, v);
  pcms::detail::CheckWrittenValueBasis("L2Integral::Distance", uu.GetBasis(),
                                       vv.GetBasis());
  const auto difference = pcms::detail::Elementwise<Real>(
    uu.GetValues(), vv.GetValues(),
    KOKKOS_LAMBDA(const Real a, const Real b) { return a - b; });
  const ValueView<const Real, DeviceMemorySpace> d(
    uu.GetBasis(), MakeConstRank2View(difference));
  return std::sqrt(Contract(*form_, d, d));
}

Real L2Integral::Norm(const Field<Real>& u) const
{
  const auto values = EvaluateTest(*form_, u);
  return std::sqrt(Contract(*form_, values, values));
}

} // namespace pcms
