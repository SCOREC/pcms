#include "pcms/transfer/forms.hpp"
#include "pcms/geometry/detail/geometry.hpp"
#include "pcms/field/detail/field_ops.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include <Kokkos_Core.hpp>
#include <string>
#include <utility>

namespace pcms::transfer::detail
{

namespace
{

constexpr Real RowSumTolerance = 1e-12;

// whether every value exceeds a round-off fraction of the largest magnitude
bool AllPositive(const Rank2View<const Real, DeviceMemorySpace>& values)
{
  const auto nc = static_cast<int>(values.extent(1));
  Real lo = 0;
  Real hi = 0;
  Kokkos::parallel_reduce(
    "pcms::transfer::AllPositive",
    Kokkos::RangePolicy<DefaultExecutionSpace>(
      0, static_cast<int>(values.extent(0)) * nc),
    KOKKOS_LAMBDA(const int k, Real& partial_lo, Real& partial_hi) {
      const Real v = values(k / nc, k % nc);
      partial_lo = v < partial_lo ? v : partial_lo;
      partial_hi = Kokkos::abs(v) > partial_hi ? Kokkos::abs(v) : partial_hi;
    },
    Kokkos::Min<Real>(lo), Kokkos::Max<Real>(hi));
  return lo > RowSumTolerance * hi;
}

void ApplyWeights(const Kokkos::View<const Real*, DeviceMemorySpace>& weights,
                  const Rank2View<const Real, DeviceMemorySpace>& values,
                  const Kokkos::View<Real**, DeviceMemorySpace>& out)
{
  const auto num_points = static_cast<LO>(out.extent(0));
  const auto num_components = static_cast<int>(out.extent(1));
  Kokkos::parallel_for(
    "pcms::transfer::ApplyWeights",
    Kokkos::RangePolicy<DefaultExecutionSpace>(0, num_points),
    KOKKOS_LAMBDA(const LO p) {
      for (int c = 0; c < num_components; ++c) {
        out(p, c) = weights(p) * values(p, c);
      }
    });
}

// sum over the points of w(p) a(p) . b(p)
Real WeightedDot(const Kokkos::View<const Real*, DeviceMemorySpace>& weights,
                 const Rank2View<const Real, DeviceMemorySpace>& a,
                 const Rank2View<const Real, DeviceMemorySpace>& b)
{
  const auto nc = static_cast<int>(a.extent(1));
  Real sum = 0;
  Kokkos::parallel_reduce(
    "pcms::transfer::WeightedDot",
    Kokkos::RangePolicy<DefaultExecutionSpace>(
      0, static_cast<int>(a.extent(0)) * nc),
    KOKKOS_LAMBDA(const int k, Real& partial) {
      partial += weights(k / nc) * a(k / nc, k % nc) * b(k / nc, k % nc);
    },
    sum);
  return sum;
}

// the diagonal of E^T W E over the DOF holders with ghosts, read back at the
// owned ones
Kokkos::View<Real**, DeviceMemorySpace> Diagonal(
  const PointMatrix<Real>& matrix,
  const Kokkos::View<const Real*, DeviceMemorySpace>& weights,
  const FieldLayout& layout)
{
  const int num_components = layout.GetNumComponents();
  const bool ghosted = layout.HasGhosts();
  const auto owned = layout.GetOwnedIndices();
  const LO num_holders = ghosted
                           ? static_cast<LO>(layout.GetOwnedHost().extent(0))
                           : layout.GetNumDofHolder();
  Kokkos::View<Real*, DeviceMemorySpace> all("bilinear_form_diagonal",
                                             num_holders);
  Kokkos::parallel_for(
    "pcms::transfer::AssembleDiagonal",
    Kokkos::RangePolicy<DefaultExecutionSpace>(
      0, static_cast<LO>(weights.extent(0))),
    KOKKOS_LAMBDA(const LO p) {
      for (LO j = matrix.row_offsets(p); j < matrix.row_offsets(p + 1); ++j) {
        Kokkos::atomic_add(&all(matrix.holders(j)),
                           weights(p) * matrix.values(j) * matrix.values(j));
      }
    });
  Kokkos::View<Real**, DeviceMemorySpace> diagonal(
    "bilinear_form_owned_diagonal", layout.GetNumDofHolder(), num_components);
  Kokkos::parallel_for(
    "pcms::transfer::OwnedDiagonal",
    Kokkos::RangePolicy<DefaultExecutionSpace>(0, layout.GetNumDofHolder()),
    KOKKOS_LAMBDA(const LO o) {
      for (int c = 0; c < num_components; ++c) {
        diagonal(o, c) = all(ghosted ? owned(o) : o);
      }
    });
  return diagonal;
}

BilinearForm::Coo CooEntries(
  const PointMatrix<Real>& matrix,
  const Kokkos::View<const Real*, DeviceMemorySpace>& weights)
{
  const auto num_points = static_cast<LO>(weights.extent(0));
  // the entries of point p start at start(p): its row has n entries, the
  // point contributes n x n
  Kokkos::View<LO*, DeviceMemorySpace> start("coo_start", num_points + 1);
  Kokkos::parallel_scan(
    "pcms::transfer::CooOffsets",
    Kokkos::RangePolicy<DefaultExecutionSpace>(0, num_points + 1),
    KOKKOS_LAMBDA(const LO p, LO& partial, const bool final) {
      if (final) {
        start(p) = partial;
      }
      if (p < num_points) {
        const LO n = matrix.row_offsets(p + 1) - matrix.row_offsets(p);
        partial += n * n;
      }
    });
  LO nnz = 0;
  Kokkos::deep_copy(nnz, Kokkos::subview(start, num_points));
  BilinearForm::Coo coo{
    Kokkos::View<GO*, DeviceMemorySpace>("coo_rows", nnz),
    Kokkos::View<GO*, DeviceMemorySpace>("coo_columns", nnz),
    Kokkos::View<Real*, DeviceMemorySpace>("coo_values", nnz)};
  const auto rows = coo.rows;
  const auto columns = coo.columns;
  const auto values = coo.values;
  Kokkos::parallel_for(
    "pcms::transfer::AssembleCoo",
    Kokkos::RangePolicy<DefaultExecutionSpace>(0, num_points),
    KOKKOS_LAMBDA(const LO p) {
      const LO first = matrix.row_offsets(p);
      const LO n = matrix.row_offsets(p + 1) - first;
      for (LO i = 0; i < n; ++i) {
        for (LO j = 0; j < n; ++j) {
          const LO k = start(p) + i * n + j;
          rows(k) = matrix.columns(first + i);
          columns(k) = matrix.columns(first + j);
          values(k) =
            weights(p) * matrix.values(first + i) * matrix.values(first + j);
        }
      }
    });
  return coo;
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

LinearForm::LinearForm(std::shared_ptr<const FunctionSpace> space,
                       Quadrature quadrature)
  : space_(std::move(space)), quadrature_(std::move(quadrature))
{
  if (quadrature_.points == nullptr ||
      quadrature_.points->GetGeometry() != space_->GetGeometry() ||
      quadrature_.weights.extent(0) !=
        static_cast<size_t>(quadrature_.points->NumPoints())) {
    throw pcms_error("LinearForm: the quadrature must be bound on the "
                     "space's geometry and have one weight per point");
  }
  evaluator_ = space_->CreatePointEvaluator<Real>(quadrature_.points);
}

ValueView<const Real, DeviceMemorySpace> LinearForm::Weight(
  ValueView<const Real, DeviceMemorySpace> samples) const
{
  if (samples.extent(0) != quadrature_.weights.extent(0)) {
    throw std::invalid_argument(
      "LinearForm: samples must have one row per point");
  }
  auto& weighted = Shaped(weighted_, samples.extent(0), samples.extent(1));
  ApplyWeights(quadrature_.weights, samples.GetValues(), weighted);
  return ValueView<const Real, DeviceMemorySpace>(samples.GetBasis(),
                                                  MakeConstRank2View(weighted));
}

void LinearForm::Assemble(ValueView<const Real, DeviceMemorySpace> samples,
                          Field<Real>& functional) const
{
  if (!evaluator_->SupportsTranspose()) {
    throw pcms_error("LinearForm: the evaluator of the space does not "
                     "support the transpose");
  }
  pcms::detail::SetZero(functional);
  evaluator_->EvaluateTranspose(Weight(samples), functional);
}

BilinearForm::BilinearForm(std::shared_ptr<const LinearForm> test,
                           std::shared_ptr<const PointEvaluatorFactory> trial,
                           FormOptions options)
  : test_(std::move(test)),
    trial_(std::move(trial)),
    options_(std::move(options))
{
  if (trial_ == nullptr || trial_->GetGeometry() == nullptr) {
    throw pcms_error("BilinearForm: the trial must have a geometry");
  }
  if (trial_ != test_->Space() || options_.to_trial != nullptr) {
    trial_evaluator_ = trial_->CreatePointEvaluator<Real>(PointsOn(
      test_->GetQuadrature(), *trial_->GetGeometry(), options_.to_trial));
  }
}

BilinearForm::BilinearForm(std::shared_ptr<const FunctionSpace> test,
                           std::shared_ptr<const PointEvaluatorFactory> trial,
                           FormOptions options)
  : BilinearForm(
      std::make_shared<const LinearForm>(
        test, options.quadrature ? *options.quadrature : test->GetQuadrature()),
      std::move(trial), options)
{
}

BilinearForm::BilinearForm(std::shared_ptr<const FunctionSpace> space,
                           FormOptions options)
  : BilinearForm(space, space, std::move(options))
{
}

BilinearForm BilinearForm::WithTrial(
  std::shared_ptr<const PointEvaluatorFactory> trial) const
{
  FormOptions options = options_;
  options.to_trial = nullptr;
  return BilinearForm(test_, std::move(trial), std::move(options));
}

const pcms::detail::TransformChain& BilinearForm::G(
  const ValueBasis& basis) const
{
  for (const auto& g : g_) {
    if (g->chain.InputBasis() == basis) {
      return g->chain;
    }
  }
  auto g = std::make_unique<Operator>();
  std::vector<const ValueTransform*> transforms;
  if (options_.metric && basis.Rank() == 1 && !basis.system->IsOrthonormal()) {
    const auto& points = *GetQuadrature().points;
    const auto& geometry = *points.GetGeometry();
    BoundCoordinateSystem::MetricView metric;
    if (pcms::detail::SameCoordinateSystem(basis.system,
                                           geometry.GetCoordinateSystem())) {
      if (metric_ == nullptr) {
        metric_ = basis.system->Bind(points.LocatedPoints().GetValues());
      }
      metric = metric_->Metric();
    } else {
      throw pcms_error("BilinearForm: the metric of values in '" +
                       std::string(basis.system->Kind()) +
                       "' is not known at the points");
    }
    g->lower = pcms::detail::LowerIndex(basis.system, metric, basis.scaling);
    transforms.push_back(&*g->lower);
  }
  g->chain.Reset(basis, std::move(transforms));
  g_.push_back(std::move(g));
  return g_.back()->chain;
}

ValueView<const Real, DeviceMemorySpace> BilinearForm::EvaluateTest(
  const Field<Real>& v) const
{
  const auto& evaluator = test_->Evaluator();
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
  const auto& g = G(b.GetBasis());
  auto gb = b.GetValues();
  if (!g.Empty()) {
    auto& out = Shaped(g_values_, gb.extent(0), gb.extent(1));
    g.Apply(gb, MakeRank2View(out));
    gb = MakeConstRank2View(out);
  }
  return WeightedDot(GetQuadrature().weights, a.GetValues(), gb);
}

Real BilinearForm::Pair(const Field<Real>& v, const Field<Real>& u) const
{
  return Contract(EvaluateTest(v), EvaluateTrial(u));
}

void BilinearForm::Apply(const Field<Real>& u, Field<Real>& functional) const
{
  if (!SupportsApply()) {
    throw pcms_error("BilinearForm::Apply: the test space's evaluator does "
                     "not support the transpose");
  }
  const auto samples = EvaluateTrial(u);
  const auto& g = G(samples.GetBasis());
  auto values = samples.GetValues();
  if (!g.Empty()) {
    auto& out = Shaped(g_values_, values.extent(0), values.extent(1));
    g.Apply(values, MakeRank2View(out));
    values = MakeConstRank2View(out);
  }
  test_->Assemble(
    ValueView<const Real, DeviceMemorySpace>(g.OutputBasis(), values),
    functional);
}

void BilinearForm::ApplyTranspose(const Field<Real>& v,
                                  Field<Real>& functional) const
{
  if (!SupportsTranspose()) {
    throw pcms_error("BilinearForm::ApplyTranspose: the trial's evaluator "
                     "does not support the transpose");
  }
  const ValueBasis written = TrialBasis(functional.GetData().GetValueBasis());
  const auto& g = G(written);
  const auto samples = EvaluateTest(v);
  pcms::detail::CheckWrittenValueBasis("BilinearForm::ApplyTranspose",
                                       samples.GetBasis(), g.OutputBasis());
  const auto weighted = test_->Weight(samples).GetValues();
  TrialEvaluator().EvaluateTranspose(ValueView<const Real, DeviceMemorySpace>(
                                       written, g.ApplyTranspose(weighted)),
                                     functional);
}

bool BilinearForm::SupportsMatrix() const noexcept
{
  return trial_evaluator_ == nullptr && !options_.metric &&
         test_->Evaluator().SupportsMatrix();
}

void BilinearForm::AssembleLumped(Field<Real>& lumped) const
{
  if (trial_evaluator_ != nullptr) {
    throw pcms_error("BilinearForm::AssembleLumped: the form is not square");
  }
  auto ones = ScratchFunction(*TestSpace());
  const auto& layout = *TestSpace()->GetLayout();
  Kokkos::View<Real**, DeviceMemorySpace> values(
    "bilinear_form_ones", layout.GetNumDofHolder(), layout.GetNumComponents());
  Kokkos::deep_copy(values, 1.0);
  ones.SetDOFHolderDataUnchecked(MakeConstRank2View(values));
  ones.Synchronize();
  Apply(ones, lumped);
  // the integral of each basis function, E^T W 1
  const auto samples = EvaluateTest(ones);
  Kokkos::View<Real**, DeviceMemorySpace> unit(
    "bilinear_form_unit", samples.extent(0), samples.extent(1));
  Kokkos::deep_copy(unit, 1.0);
  test_->Assemble(ValueView<const Real, DeviceMemorySpace>(
                    samples.GetBasis(), MakeConstRank2View(unit)),
                  ones);
  if (!AllPositive(ones.GetDOFHolderData().GetValues())) {
    throw pcms_error("BilinearForm::AssembleLumped: a basis function has "
                     "no positive integral; a lumped mass needs a "
                     "non-negative basis");
  }
}

void BilinearForm::AssembleDiagonal(Field<Real>& diagonal) const
{
  const auto& layout = *TestSpace()->GetLayout();
  if (!SupportsMatrix() || layout.IsDistributed()) {
    throw pcms_error("BilinearForm::AssembleDiagonal: needs a square form "
                     "with the evaluation matrix, on a layout on one rank");
  }
  const auto values =
    Diagonal(test_->Evaluator().GetMatrix(), GetQuadrature().weights, layout);
  diagonal.SetDOFHolderDataUnchecked(MakeConstRank2View(values));
  diagonal.Synchronize();
}

BilinearForm::Coo BilinearForm::AssembleCoo() const
{
  if (!SupportsMatrix()) {
    throw pcms_error("BilinearForm::AssembleCoo: needs a square form with "
                     "the evaluation matrix");
  }
  return CooEntries(test_->Evaluator().GetMatrix(), GetQuadrature().weights);
}

} // namespace pcms::transfer::detail
