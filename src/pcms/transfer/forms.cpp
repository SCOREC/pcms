#include "pcms/transfer/forms.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

namespace pcms::transfer::detail
{

// Maps the global id of a DOF holder this rank owns to the holder: a binary
// search in the sorted gids; -1 for another gid.
struct HolderOfGid
{
  Kokkos::View<GO*, DeviceMemorySpace> sorted;
  Kokkos::View<LO*, DeviceMemorySpace> holders;

  KOKKOS_INLINE_FUNCTION LO operator()(GO gid) const
  {
    LO lo = 0;
    LO hi = static_cast<LO>(sorted.extent(0));
    while (lo < hi) {
      const LO mid = lo + (hi - lo) / 2;
      if (sorted(mid) < gid) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }
    return lo < static_cast<LO>(sorted.extent(0)) && sorted(lo) == gid
             ? holders(lo)
             : -1;
  }
};

namespace
{

constexpr Real RowSumTolerance = 1e-12;

HolderOfGid MakeHolderOfGid(const FieldLayout& layout)
{
  const auto gids = layout.GetGidsHost();
  std::vector<LO> order(gids.extent(0));
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(),
            [&](LO a, LO b) { return gids(a) < gids(b); });
  Kokkos::View<GO*, HostMemorySpace> sorted("sorted_gids", order.size());
  Kokkos::View<LO*, HostMemorySpace> holders("sorted_holders", order.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    sorted(i) = gids(order[i]);
    holders(i) = order[i];
  }
  return {Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, sorted),
          Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, holders)};
}

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

// the diagonal of E^T W E at the owned DOF holders
Kokkos::View<Real**, DeviceMemorySpace> Diagonal(
  const PointMatrix<Real>& matrix,
  const Kokkos::View<const Real*, DeviceMemorySpace>& weights,
  const FieldLayout& layout)
{
  const auto holder_of = MakeHolderOfGid(layout);
  const int num_components = layout.GetNumComponents();
  Kokkos::View<Real**, DeviceMemorySpace> diagonal(
    "bilinear_form_diagonal", layout.GetNumDofHolder(), num_components);
  Kokkos::parallel_for(
    "pcms::transfer::AssembleDiagonal",
    Kokkos::RangePolicy<DefaultExecutionSpace>(
      0, static_cast<LO>(weights.extent(0))),
    KOKKOS_LAMBDA(const LO p) {
      for (LO j = matrix.row_offsets(p); j < matrix.row_offsets(p + 1); ++j) {
        const LO holder = holder_of(matrix.columns(j));
        if (holder < 0) {
          continue;
        }
        for (int c = 0; c < num_components; ++c) {
          Kokkos::atomic_add(&diagonal(holder, c),
                             weights(p) * matrix.values(j) * matrix.values(j));
        }
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
  const auto& layout = functional.GetLayout();
  functional.SetDOFHolderDataUnchecked(MakeConstRank2View(
    Kokkos::View<Real**, DeviceMemorySpace>("linear_form_zeros",
                                            layout.GetNumDofHolder(),
                                            layout.GetNumComponents())));
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
  return WeightedDot(GetQuadrature().weights, a.GetValues(), b.GetValues());
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
  test_->Assemble(EvaluateTrial(u), functional);
}

void BilinearForm::ApplyTranspose(const Field<Real>& v,
                                  Field<Real>& functional) const
{
  if (!SupportsTranspose()) {
    throw pcms_error("BilinearForm::ApplyTranspose: the trial's evaluator "
                     "does not support the transpose");
  }
  const ValueBasis written = TrialBasis(functional.GetData().GetValueBasis());
  const auto samples = EvaluateTest(v);
  pcms::detail::CheckWrittenValueBasis("BilinearForm::ApplyTranspose",
                                       samples.GetBasis(), written);
  TrialEvaluator().EvaluateTranspose(test_->Weight(samples), functional);
}

bool BilinearForm::SupportsMatrix() const noexcept
{
  return trial_evaluator_ == nullptr && test_->Evaluator().SupportsMatrix();
}

void BilinearForm::AssembleLumped(Field<Real>& lumped) const
{
  if (trial_evaluator_ != nullptr) {
    throw pcms_error("BilinearForm::AssembleLumped: the form is not square");
  }
  auto ones = TestSpace()->CreateFunction<Real>();
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
