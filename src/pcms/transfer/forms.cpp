#include "pcms/transfer/forms.hpp"
#include "pcms/geometry/detail/geometry.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include <Kokkos_Core.hpp>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <utility>

namespace pcms::transfer::detail
{

namespace
{

void ApplyWeights(const Kokkos::View<const Real*, DeviceMemorySpace>& weights,
                  const Rank2View<const Real, DeviceMemorySpace>& values,
                  const Kokkos::View<Real**, DeviceMemorySpace>& out)
{
  Kokkos::parallel_for(
    "pcms::transfer::ApplyWeights",
    Kokkos::MDRangePolicy<DeviceMemorySpace::execution_space, Kokkos::Rank<2>>(
      {0, 0}, {static_cast<int64_t>(out.extent(0)),
               static_cast<int64_t>(out.extent(1))}),
    KOKKOS_LAMBDA(const int64_t p, const int64_t c) {
      out(p, c) = weights(p) * values(p, c);
    });
}

// sum over the points of w(p) a(p) . b(p)
Real WeightedDot(const Kokkos::View<const Real*, DeviceMemorySpace>& weights,
                 const Rank2View<const Real, DeviceMemorySpace>& a,
                 const Rank2View<const Real, DeviceMemorySpace>& b)
{
  Real sum = 0;
  Kokkos::parallel_reduce(
    "pcms::transfer::WeightedDot",
    Kokkos::MDRangePolicy<DeviceMemorySpace::execution_space, Kokkos::Rank<2>>(
      {0, 0},
      {static_cast<int64_t>(a.extent(0)), static_cast<int64_t>(a.extent(1))}),
    KOKKOS_LAMBDA(const int64_t i, const int64_t c, Real& partial) {
      partial += weights(i) * a(i, c) * b(i, c);
    },
    sum);
  return sum;
}

// the points of `quadrature` bound on `geometry`: its own points, else a bind
// of them that starts from the cells the quadrature names
std::shared_ptr<const BoundPoints> PointsOn(
  const Quadrature& quadrature, const Geometry& geometry,
  const std::shared_ptr<const CoordinateMap>& map)
{
  if (map == nullptr && quadrature.points->GetGeometry().get() == &geometry) {
    return quadrature.points;
  }
  return geometry.Bind(quadrature.points->LocatedPoints(), quadrature.other,
                       OutOfBoundsPolicy{}, map);
}

// a view shaped [rows][columns], reallocated only when its shape changes
Kokkos::View<Real**, DeviceMemorySpace>& Shaped(
  Kokkos::View<Real**, DeviceMemorySpace>& view, std::size_t rows,
  std::size_t columns)
{
  if (view.extent(0) != rows || view.extent(1) != columns) {
    view = Kokkos::View<Real**, DeviceMemorySpace>(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "bilinear_form_values"),
      rows, columns);
  }
  return view;
}

} // namespace

BilinearForm::BilinearForm(std::shared_ptr<const Test> test,
                           std::shared_ptr<const PointEvaluatorFactory> trial,
                           FormOptions options)
  : test_(std::move(test)),
    trial_(std::move(trial)),
    options_(std::move(options))
{
  if (trial_ == nullptr || trial_->GetGeometry() == nullptr) {
    throw pcms_error("BilinearForm: the trial must have a geometry");
  }
  if (trial_ != test_->space || options_.to_trial != nullptr) {
    trial_evaluator_ = trial_->CreatePointEvaluator<Real>(
      PointsOn(test_->quadrature, *trial_->GetGeometry(), options_.to_trial));
  }
}

BilinearForm::BilinearForm(std::shared_ptr<const FunctionSpace> test,
                           std::shared_ptr<const PointEvaluatorFactory> trial,
                           FormOptions options)
  : BilinearForm(
      [&] {
        auto side = std::make_shared<Test>();
        side->quadrature =
          options.quadrature ? *options.quadrature : test->GetQuadrature();
        const auto& q = side->quadrature;
        if (q.points == nullptr ||
            q.points->GetGeometry() != test->GetGeometry() ||
            q.weights.extent(0) != static_cast<size_t>(q.points->NumPoints())) {
          throw pcms_error("BilinearForm: the quadrature must be bound on the "
                           "test space's geometry and have one weight per "
                           "point");
        }
        side->evaluator = test->CreatePointEvaluator<Real>(q.points);
        side->space = std::move(test);
        return std::shared_ptr<const Test>(std::move(side));
      }(),
      std::move(trial), options)
{
}

BilinearForm::BilinearForm(std::shared_ptr<const FunctionSpace> space,
                           FormOptions options)
  : BilinearForm(space, space, std::move(options))
{
}

ValueView<const Real, DeviceMemorySpace> BilinearForm::Weight(
  ValueView<const Real, DeviceMemorySpace> samples) const
{
  if (samples.extent(0) != test_->quadrature.weights.extent(0)) {
    throw std::invalid_argument(
      "BilinearForm: samples must have one row per point");
  }
  auto& weighted = Shaped(weighted_, samples.extent(0), samples.extent(1));
  ApplyWeights(test_->quadrature.weights, samples.GetValues(), weighted);
  return ValueView<const Real, DeviceMemorySpace>(samples.GetBasis(),
                                                  MakeConstRank2View(weighted));
}

void BilinearForm::Assemble(ValueView<const Real, DeviceMemorySpace> samples,
                            Field<Real>& functional) const
{
  test_->evaluator->EvaluateTranspose(Weight(samples), functional);
}

BilinearForm BilinearForm::WithTrial(
  std::shared_ptr<const PointEvaluatorFactory> trial) const
{
  FormOptions options = options_;
  options.to_trial = nullptr;
  return BilinearForm(test_, std::move(trial), std::move(options));
}

ValueView<const Real, DeviceMemorySpace> BilinearForm::EvaluateTest(
  const Field<Real>& v) const
{
  const auto& evaluator = *test_->evaluator;
  const ValueBasis basis = evaluator.OutputBasis(v.GetData().GetValueBasis());
  auto& values = Shaped(test_values_, evaluator.NumPoints(),
                        v.GetLayout().GetNumComponents());
  evaluator.Evaluate(
    v, ValueView<Real, DeviceMemorySpace>(basis, MakeRank2View(values)));
  return {basis, MakeConstRank2View(values)};
}

ValueView<const Real, DeviceMemorySpace> BilinearForm::EvaluateTrial(
  const Field<Real>& u) const
{
  const auto& evaluator = TrialEvaluator();
  const ValueBasis basis = evaluator.OutputBasis(u.GetData().GetValueBasis());
  auto& values = Shaped(trial_values_, evaluator.NumPoints(),
                        u.GetLayout().GetNumComponents());
  evaluator.Evaluate(
    u, ValueView<Real, DeviceMemorySpace>(basis, MakeRank2View(values)));
  return {basis, MakeConstRank2View(values)};
}

Real BilinearForm::Contract(ValueView<const Real, DeviceMemorySpace> a,
                            ValueView<const Real, DeviceMemorySpace> b) const
{
  Real sum = WeightedDot(GetQuadrature().weights, a.GetValues(), b.GetValues());
  pcms::detail::SumOverRanks(*GetQuadrature().points->GetGeometry(),
                             std::span<Real>(&sum, 1));
  return sum;
}

Real BilinearForm::Pair(const Field<Real>& v, const Field<Real>& u) const
{
  return Contract(EvaluateTest(v), EvaluateTrial(u));
}

void BilinearForm::Apply(const Field<Real>& u, Field<Real>& functional) const
{
  Assemble(EvaluateTrial(u), functional);
}

void BilinearForm::ApplyTranspose(const Field<Real>& v, Field<Real>& functional,
                                  Real beta) const
{
  const ValueBasis written = TrialBasis(functional.GetData().GetValueBasis());
  const auto samples = EvaluateTest(v);
  pcms::detail::CheckWrittenValueBasis("BilinearForm::ApplyTranspose",
                                       samples.GetBasis(), written);
  pcms::detail::EvaluateTranspose(TrialEvaluator(), Weight(samples), functional,
                                  beta);
}

} // namespace pcms::transfer::detail
