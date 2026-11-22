#ifndef PCMS_BACKENDS_POINT_CLOUD_MLS_EVALUATOR_HPP
#define PCMS_BACKENDS_POINT_CLOUD_MLS_EVALUATOR_HPP

#include "pcms/field/detail/field_ops.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/mls_options.hpp"
#include "pcms/backends/point_cloud/mls_interpolation.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/memory_spaces.hpp"

#include <Kokkos_Core.hpp>

#include <string>

namespace pcms::detail
{

/// Whether each of `n` query points has supports, from the supports' row
/// pointers. (A free function: CUDA does not allow device lambdas in
/// constructors or non-public member functions.)
inline Kokkos::View<bool*, DeviceMemorySpace> HasSupports(
  const Kokkos::View<const LO*, DeviceMemorySpace>& supports_ptr, LO n)
{
  Kokkos::View<bool*, DeviceMemorySpace> found(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "mls_found"), n);
  pcms::detail::ForEach(
    "MLSBinding::found", n, KOKKOS_LAMBDA(const LO i) {
      found(i) = supports_ptr[i + 1] > supports_ptr[i];
    });
  return found;
}

/// The rows of `coords` (n x dim) as one flat array [x0, y0, ..., x1, ...].
template <typename View2D>
Kokkos::View<Real*, DeviceMemorySpace> FlattenRows(const View2D& coords,
                                                   const char* label)
{
  const auto n = static_cast<LO>(coords.extent(0));
  const auto dim = static_cast<LO>(coords.extent(1));
  Kokkos::View<Real*, DeviceMemorySpace> flat(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(label)),
    n * dim);
  Kokkos::parallel_for(
    "FlattenRows",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(const LO i) {
      for (LO d = 0; d < dim; ++d) {
        flat(i * dim + d) = coords(i, d);
      }
    });
  return flat;
}

/// @brief MLS BoundPoints: the supports of the query points and the
///        coordinates the solve uses. A query point without supports (e.g.
///        one the adjacency search has no start for) is missing.
class MLSBinding final : public BoundPoints
{
public:
  /// @param geometry geometry the points are localized on
  /// @param request query points and out-of-bounds policy
  /// @param source_coords flattened coordinates of the source points
  /// @param supports MLS supports of the query points
  /// @param site_dim the site dimension and
  /// @param options the options `supports` were found with; a rebind finds
  ///        the supports of its points with the geometry's Neighbourhood
  /// @param axes the periods of the axes
  MLSBinding(std::shared_ptr<const Geometry> geometry,
             const BindRequest& request,
             Kokkos::View<const Real*, DeviceMemorySpace> source_coords,
             SupportResults supports, int site_dim,
             NeighbourhoodOptions options, AxisPeriods axes)
    : BoundPoints(std::move(geometry), request),
      source_coords_(std::move(source_coords)),
      target_coords_(
        detail::FlattenRows(LocatedPointValues(), "mls_target_coords")),
      supports_(std::move(supports)),
      site_dim_(site_dim),
      options_(options),
      axes_(axes)
  {
    SetFoundFromSupports();
  }

  /// @return flattened coordinates of the source points
  [[nodiscard]] const Kokkos::View<const Real*, DeviceMemorySpace>&
  SourceCoordinates() const noexcept
  {
    return source_coords_;
  }

  /// @return flattened coordinates of the query points in the geometry's
  ///         coordinate system (LocatedPoints())
  [[nodiscard]] const Kokkos::View<const Real*, DeviceMemorySpace>&
  TargetCoordinates() const noexcept
  {
    return target_coords_;
  }

  /// @return MLS supports of the query points
  [[nodiscard]] const SupportResults& Supports() const noexcept
  {
    return supports_;
  }

  /// @brief The periods of the axes.
  [[nodiscard]] const AxisPeriods& Axes() const noexcept { return axes_; }

protected:
  // the arrays are reassigned, not written, so the default UnsharedImpl holds
  void RebindImpl() override
  {
    supports_ =
      GetGeometry()->Neighbourhood(LocatedPoints(), site_dim_, options_);
    target_coords_ =
      detail::FlattenRows(LocatedPointValues(), "mls_target_coords");
    SetFoundFromSupports();
  }

private:
  /// a point is found when it has supports
  void SetFoundFromSupports()
  {
    SetFound(detail::HasSupports(supports_.supports_ptr, NumPoints()));
  }

  Kokkos::View<const Real*, DeviceMemorySpace> source_coords_;
  Kokkos::View<const Real*, DeviceMemorySpace> target_coords_;
  SupportResults supports_;
  int site_dim_;
  NeighbourhoodOptions options_;
  AxisPeriods axes_;
};

/// The values at the targets from the coefficient rows, for every component:
/// values(i, c) = sum_j coefficients(j) data(idx(j), c) over the supports of i.
/// (Free functions: CUDA does not allow device lambdas in the virtual
/// functions of the evaluator.)
template <typename Data, typename Values>
void ApplyMLSRows(const SupportResults& supports,
                  const Kokkos::View<const Real*, DeviceMemorySpace>& rows,
                  const Data& data, const Values& values, int n_comp)
{
  const auto num_points = static_cast<LO>(supports.supports_ptr.extent(0)) - 1;
  pcms::detail::ForEach(
    "MLSPointEvaluator::Evaluate", num_points, KOKKOS_LAMBDA(const LO i) {
      for (int c = 0; c < n_comp; ++c) {
        Real sum = 0;
        for (LO j = supports.supports_ptr[i]; j < supports.supports_ptr[i + 1];
             ++j) {
          sum += rows(j) * data(supports.supports_idx[j], c);
        }
        values(i, c) = sum;
      }
    });
}

/// The transpose of ApplyMLSRows: adds coefficients(j) values(i, c) to the
/// support sites of each target i (atomically: supports overlap).
template <typename Values, typename Dofs>
void ApplyMLSRowsTranspose(
  const SupportResults& supports,
  const Kokkos::View<const Real*, DeviceMemorySpace>& rows,
  const Values& values, const Dofs& dofs, int n_comp)
{
  const auto num_points = static_cast<LO>(supports.supports_ptr.extent(0)) - 1;
  pcms::detail::ForEach(
    "MLSPointEvaluator::EvaluateTranspose", num_points,
    KOKKOS_LAMBDA(const LO i) {
      for (int c = 0; c < n_comp; ++c) {
        for (LO j = supports.supports_ptr[i]; j < supports.supports_ptr[i + 1];
             ++j) {
          Kokkos::atomic_add(&dofs(supports.supports_idx[j], c),
                             rows(j) * values(i, c));
        }
      }
    });
}

/// @brief Evaluates a point-cloud field at the bound points by MLS. The
///        coefficient rows of the targets (each value is a weighted sum of
///        the values at its supports, the same weights for every component)
///        are computed once, when the evaluator is created, so Evaluate is a
///        sparse matrix-vector product, EvaluateTranspose its transpose and
///        GetMatrix the rows themselves.
template <typename LayoutPolicy =
            pcms::detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class MLSPointEvaluator : public PointEvaluator<Real, LayoutPolicy>
{
public:
  /// @param binding MLS supports of the query points
  /// @param dim spatial dimension
  /// @param options MLS options
  /// @param layout layout of the evaluated fields
  MLSPointEvaluator(std::shared_ptr<const MLSBinding> binding, int dim,
                    MLSOptions options,
                    std::shared_ptr<const FieldLayout> layout)
    : PointEvaluator<Real, LayoutPolicy>(binding, DifferentialOperator::Value,
                                         layout),
      binding_(std::move(binding)),
      layout_(std::move(layout)),
      rows_(pcms::detail::mls_coefficients(
        binding_->SourceCoordinates(), binding_->TargetCoordinates(),
        binding_->Supports(), dim, static_cast<LO>(options.degree),
        options.basis, options.lambda, options.tol, options.decay_factor,
        binding_->Axes()))
  {
  }

  // Public because CUDA forbids extended lambdas in non-public member
  // functions. User call Evaluate.
  // TODO: add passkey or other mechanism to protect this
  void EvaluateImpl(
    const Field<Real>& field,
    ValueView<Real, DeviceMemorySpace, LayoutPolicy> out) const override
  {
    ApplyMLSRows(binding_->Supports(), rows_,
                 field.GetDOFHolderData().GetValues(), out.GetValues(),
                 layout_->GetNumComponents());
  }

  [[nodiscard]] bool SupportsTranspose() const noexcept override
  {
    return true;
  }

  // Public for the same reason as EvaluateImpl. Single rank: the supports are
  // sites of this rank's layout, whose DOFs are all owned.
  void EvaluateTransposeImpl(
    Rank2View<const Real, DeviceMemorySpace, LayoutPolicy> values,
    Field<Real>& functional) const override
  {
    const int n_comp = layout_->GetNumComponents();
    Kokkos::View<Real**, DeviceMemorySpace> dofs(
      "mls_transpose", layout_->GetNumDofHolder(), n_comp);
    ApplyMLSRowsTranspose(binding_->Supports(), rows_, values, dofs, n_comp);
    pcms::detail::AddToField(functional, dofs);
  }

  [[nodiscard]] bool SupportsMatrix() const noexcept override { return true; }

  /// The rows, with the layout's global ids as columns (single rank).
  PointMatrix<Real> GetMatrixImpl() const override
  {
    const auto& supports = binding_->Supports();
    const auto gids_host = layout_->GetGidsHost();
    Kokkos::View<GO*, HostMemorySpace> gids_h("mls_gids", gids_host.extent(0));
    for (std::size_t i = 0; i < gids_h.extent(0); ++i) {
      gids_h(i) = gids_host(i);
    }
    const auto gids =
      Kokkos::create_mirror_view_and_copy(DeviceMemorySpace{}, gids_h);
    Kokkos::View<GO*, DeviceMemorySpace> columns(
      "mls_columns", supports.supports_idx.extent(0));
    const auto idx = supports.supports_idx;
    pcms::detail::ForEach(
      "MLSPointEvaluator::GetMatrix", static_cast<LO>(columns.extent(0)),
      KOKKOS_LAMBDA(const LO j) { columns(j) = gids(idx(j)); });
    return PointMatrix<Real>{supports.supports_ptr, columns, rows_};
  }

private:
  std::shared_ptr<const MLSBinding> binding_;
  std::shared_ptr<const FieldLayout> layout_;
  // one weight per support entry, aligned with the binding's supports
  Kokkos::View<const Real*, DeviceMemorySpace> rows_;
};

} // namespace pcms::detail

#endif // PCMS_BACKENDS_POINT_CLOUD_MLS_EVALUATOR_HPP
