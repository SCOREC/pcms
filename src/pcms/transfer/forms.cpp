#include "pcms/transfer/forms.hpp"
#include "pcms/geometry/detail/geometry.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include <Kokkos_Core.hpp>
#include <algorithm>
#include <numeric>
#include <vector>

namespace pcms::transfer::detail
{

namespace
{

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

} // namespace

std::shared_ptr<const BoundPoints> BindQuadrature(
  const Geometry& geometry, const Quadrature& quadrature,
  std::shared_ptr<const CoordinateMap> map)
{
  if (map == nullptr && quadrature.points->GetGeometry().get() == &geometry) {
    return quadrature.points;
  }
  return geometry.Bind(quadrature.points->LocatedPoints(), quadrature.other,
                       OutOfBoundsPolicy{}, std::move(map));
}

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
  if (!evaluator_->SupportsTranspose()) {
    throw pcms_error("LinearForm: the space's evaluator does not support the "
                     "transpose");
  }
  const auto& layout = *space_->GetLayout();
  zeros_ = Kokkos::View<Real**, DeviceMemorySpace>(
    "linear_form_zeros", layout.GetNumDofHolder(), layout.GetNumComponents());
  weighted_ = Kokkos::View<Real**, DeviceMemorySpace>(
    "linear_form_weighted", quadrature_.points->NumPoints(),
    layout.GetNumComponents());
}

ValueView<const Real, DeviceMemorySpace> LinearForm::Weight(
  ValueView<const Real, DeviceMemorySpace> samples) const
{
  if (samples.extent(0) != weighted_.extent(0) ||
      samples.extent(1) != weighted_.extent(1)) {
    throw std::invalid_argument(
      "LinearForm: samples must be shaped [points][components]");
  }
  ApplyWeights(quadrature_.weights, samples.GetValues(), weighted_);
  return ValueView<const Real, DeviceMemorySpace>(
    samples.GetBasis(), MakeConstRank2View(weighted_));
}

void LinearForm::Assemble(ValueView<const Real, DeviceMemorySpace> samples,
                          Field<Real>& functional) const
{
  const auto weighted = Weight(samples);
  functional.SetDOFHolderDataUnchecked(MakeConstRank2View(zeros_));
  evaluator_->EvaluateTranspose(weighted, functional);
}

BilinearForm::BilinearForm(std::shared_ptr<const FunctionSpace> space,
                           Quadrature quadrature)
  : form_(std::move(space), std::move(quadrature)),
    samples_("bilinear_form_samples", form_.GetQuadrature().points->NumPoints(),
             form_.Space()->GetLayout()->GetNumComponents())
{
}

void BilinearForm::Apply(const Field<Real>& u, Field<Real>& functional) const
{
  const ValueBasis& basis = u.GetData().GetValueBasis();
  form_.Evaluator().Evaluate(
    u, ValueView<Real, DeviceMemorySpace>(basis, MakeRank2View(samples_)));
  form_.Assemble(ValueView<const Real, DeviceMemorySpace>(
                   basis, MakeConstRank2View(samples_)),
                 functional);
}

void BilinearForm::AssembleLumped(Field<Real>& lumped) const
{
  auto ones = ScratchFunction(*form_.Space());
  const auto& layout = *form_.Space()->GetLayout();
  Kokkos::View<Real**, DeviceMemorySpace> values(
    "bilinear_form_ones", layout.GetNumDofHolder(), layout.GetNumComponents());
  Kokkos::deep_copy(values, 1.0);
  ones.SetDOFHolderDataUnchecked(MakeConstRank2View(values));
  ones.Synchronize();
  Apply(ones, lumped);
}

void BilinearForm::AssembleDiagonal(Field<Real>& diagonal) const
{
  const auto& layout = *form_.Space()->GetLayout();
  if (!SupportsMatrix() || layout.IsDistributed()) {
    throw pcms_error("BilinearForm::AssembleDiagonal: needs the evaluation "
                     "matrix and a layout on one rank");
  }
  const auto matrix = form_.Evaluator().GetMatrix();
  const auto weights = form_.GetQuadrature().weights;
  const auto num_points = static_cast<LO>(weights.extent(0));
  const int num_components = layout.GetNumComponents();
  Kokkos::View<Real**, DeviceMemorySpace> diag(
    "bilinear_form_diagonal", layout.GetNumDofHolder(), num_components);
  Kokkos::parallel_for(
    "pcms::transfer::AssembleDiagonal",
    Kokkos::RangePolicy<DefaultExecutionSpace>(0, num_points),
    KOKKOS_LAMBDA(const LO p) {
      for (LO j = matrix.row_offsets(p); j < matrix.row_offsets(p + 1); ++j) {
        const LO holder = matrix.holders(j);
        for (int c = 0; c < num_components; ++c) {
          Kokkos::atomic_add(&diag(holder, c),
                             weights(p) * matrix.values(j) * matrix.values(j));
        }
      }
    });
  diagonal.SetDOFHolderDataUnchecked(MakeConstRank2View(diag));
  diagonal.Synchronize();
}

BilinearForm::Coo BilinearForm::AssembleCoo() const
{
  if (!SupportsMatrix()) {
    throw pcms_error("BilinearForm::AssembleCoo: the evaluator provides no "
                     "matrix");
  }
  const auto matrix = form_.Evaluator().GetMatrix();
  const auto weights = form_.GetQuadrature().weights;
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
  Coo coo{Kokkos::View<GO*, DeviceMemorySpace>("coo_rows", nnz),
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

} // namespace pcms::transfer::detail
