#ifndef PCMS_BACKENDS_OMEGA_H_LAGRANGE_EVALUATOR_HPP
#define PCMS_BACKENDS_OMEGA_H_LAGRANGE_EVALUATOR_HPP

#include "pcms/backends/omega_h/omega_h_geometry.hpp"
#include "pcms/backends/omega_h/lagrange_layout.hpp"
#include "pcms/backends/omega_h/field_data.hpp"
#include "pcms/field/point_rows.hpp"
#include "pcms/utility/detail/profile.hpp"
#include "pcms/utility/types.hpp"

#include <Omega_h_mesh.hpp>
#include <Kokkos_Core.hpp>
#include <memory>
#include <utility>

namespace pcms::omega_h
{

// ---------------------------------------------------------------------------
// PointEvaluator
// ---------------------------------------------------------------------------

namespace detail
{
// DOF-holder values read from an Omega_h tag: one value per component per
// mesh entity, entity-major.
template <typename T, typename TagValue>
struct OmegaHLagrangeTagValues
{
  Omega_h::Read<TagValue> values;
  int num_components;
  KOKKOS_INLINE_FUNCTION T operator()(LO holder, int c) const
  {
    return static_cast<T>(values[holder * num_components + c]);
  }
};

/// @brief Fills the evaluation rows of an order 0 or 1 Lagrange field at the
///        found points of a binding: order 0 takes the lowest element
///        adjacent to the entity of the point's record (1 entry per point);
///        order 1 the vertices of that entity, with the weights its reduced
///        coordinates give (Dim + 1 entries per point).
template <int Dim>
void FillOmegaHLagrangeRows(
  Omega_h::Mesh& mesh, int order,
  const pcms::detail::OmegaHMeshBinding& binding,
  const Kokkos::View<LO*, DeviceMemorySpace>& holders,
  const Kokkos::View<Real*, DeviceMemorySpace>& weights)
{
  static_assert(Dim == 2 || Dim == 3, "Omega_h meshes are 2D or 3D");
  const auto found = binding.FoundPoints();
  const auto results = binding.Results<Dim>();
  const auto up = UpAdjacencies<Dim>(mesh);
  const auto verts_of = VerticesOfEntities<Dim>(mesh);
  pcms::detail::ForEach(
    "OmegaHLagrangePointEvaluator::Rows", static_cast<LO>(found.extent(0)),
    KOKKOS_LAMBDA(const LO k) {
      const LO p = found(k);
      if (order == 0) {
        holders(p) = ElementOfResult<Dim>(up, results(p));
        weights(p) = 1.0;
      } else {
        EntityVertexWeights<Dim>(verts_of, results(p), &holders(p * (Dim + 1)),
                                 &weights(p * (Dim + 1)));
      }
    });
}

/// @brief The rows of an order 0 or 1 Lagrange field at the points of a
///        binding; zero at the points that are not found.
inline pcms::detail::PointRows OmegaHLagrangeRows(
  const OmegaHLagrangeLayout& layout,
  const pcms::detail::OmegaHMeshBinding& binding)
{
  auto& mesh = const_cast<Omega_h::Mesh&>(layout.GetMesh());
  const int stride = layout.GetOrder() == 0 ? 1 : mesh.dim() + 1;
  Kokkos::View<LO*, DeviceMemorySpace> holders(
    "omega_h_lagrange_holders",
    static_cast<std::size_t>(binding.NumPoints()) * stride);
  Kokkos::View<Real*, DeviceMemorySpace> weights("omega_h_lagrange_weights",
                                                 holders.extent(0));
  VisitMeshDimension(mesh, [&]<int Dim>() {
    FillOmegaHLagrangeRows<Dim>(mesh, layout.GetOrder(), binding, holders,
                                weights);
  });
  return {pcms::detail::StrideOffsets(binding.NumPoints(), stride), holders,
          weights};
}

} // namespace detail

// OmegaHLagrangePointEvaluator implements the evaluation of Real fields on
// simplex meshes backed by Omega_h, using the search records of the
// BoundPoints.
//
// Output shape: [num_query_points][num_components].
//
// The values are read from the field's Omega_h tag
// (OmegaHFieldData::GetTagName). Its owned entries are always current; its
// ghost entries, which elements on the partition boundary read, are correct
// only after the field's last write was followed by Synchronize.
class OmegaHLagrangePointEvaluator final
  : public pcms::detail::RowsPointEvaluator
{
public:
  /// @param layout Omega_h Lagrange layout of the evaluated fields
  /// @param binding the search records of the query points
  OmegaHLagrangePointEvaluator(
    const std::shared_ptr<const OmegaHLagrangeLayout>& layout,
    const std::shared_ptr<const pcms::detail::OmegaHMeshBinding>& binding)
    : RowsPointEvaluator(binding, layout, DifferentialOperator::Value,
                         detail::OmegaHLagrangeRows(*layout, *binding))
  {
  }

  // Public because CUDA forbids extended lambdas in non-public member
  // functions. User call Evaluate.
  void EvaluateImpl(const Field<Real>& field,
                    ValueView<Real, DeviceMemorySpace> out) const override
  {
    PCMS_FUNCTION_TIMER;
    // Elements read the values of all their DOF holders, ghosts included, so
    // read the field's Omega_h tag rather than the owned values.
    const auto& data =
      dynamic_cast<const OmegaHFieldData<Real>&>(field.GetData());
    using TagValue = typename OmegaHFieldData<Real>::tag_value_type;
    pcms::detail::Apply(
      Rows(),
      detail::OmegaHLagrangeTagValues<Real, TagValue>{
        Layout().GetMesh().get_array<TagValue>(Layout().GetDOFHolderEntityDim(),
                                               data.GetTagName()),
        Layout().GetNumComponents()},
      out.GetValues());
  }

  // Public for the same reason as EvaluateImpl.
  void EvaluateTransposeImpl(Rank2View<const Real, DeviceMemorySpace> values,
                             Field<Real>& functional) const override
  {
    PCMS_FUNCTION_TIMER;
    using TagValue = utility::OmegaHValue_t<Real>;
    const auto& layout = Layout();
    const int n_comp = layout.GetNumComponents();
    Omega_h::Write<TagValue> sums(layout.GetNumDofHolderWithGhosts() * n_comp,
                                  0, "transpose_sums");
    pcms::detail::ApplyTranspose(
      Rows(), values,
      Kokkos::View<TagValue**, Kokkos::LayoutRight, DeviceMemorySpace,
                   Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
        sums.data(), layout.GetNumDofHolderWithGhosts(), n_comp));
    const bool ghosted = layout.HasGhosts();
    // the copies of a holder on other ranks add into its owner (collective)
    const auto reduced =
      layout.GetMesh().comm()->size() > 1
        ? layout.GetMesh().reduce_array(layout.GetDOFHolderEntityDim(),
                                        Omega_h::Read<TagValue>(sums), n_comp,
                                        OMEGA_H_SUM)
        : Omega_h::Read<TagValue>(sums);
    const auto owned = layout.GetOwnedIndices();
    Kokkos::View<Real**, DeviceMemorySpace> addend(
      "transpose_addend", layout.GetNumDofHolder(), n_comp);
    Kokkos::parallel_for(
      "OmegaHLagrangePointEvaluator::EvaluateTranspose::owned",
      Kokkos::RangePolicy<DefaultExecutionSpace>(0, layout.GetNumDofHolder()),
      KOKKOS_LAMBDA(const LO o) {
        const LO h = ghosted ? owned(o) : o;
        for (int c = 0; c < n_comp; ++c) {
          addend(o, c) = static_cast<Real>(reduced[h * n_comp + c]);
        }
      });
    pcms::detail::AddToField(functional, addend);
  }

  [[nodiscard]] PointMatrix<Real> GetMatrixImpl() const override
  {
    return pcms::detail::ToPointMatrix(Rows(), Layout().GetGids());
  }

private:
  [[nodiscard]] const OmegaHLagrangeLayout& Layout() const
  {
    return static_cast<const OmegaHLagrangeLayout&>(*AcceptedLayout());
  }
};

// ---------------------------------------------------------------------------
// EvaluatorFactory
// ---------------------------------------------------------------------------

/// @brief Point evaluation of Real Lagrange fields on an Omega_h simplex
///        mesh.
class OmegaHLagrangeEvaluatorFactory
  : public pcms::detail::TypedPointEvaluatorFactory<
      Real, pcms::detail::OmegaHMeshBinding>
{
public:
  /// @param geometry geometry of the layout's mesh
  /// @param layout Omega_h Lagrange layout of the evaluated fields
  OmegaHLagrangeEvaluatorFactory(
    std::shared_ptr<const OmegaHGeometry> geometry,
    std::shared_ptr<const OmegaHLagrangeLayout> layout)
    : pcms::detail::TypedPointEvaluatorFactory<Real,
                                               pcms::detail::OmegaHMeshBinding>(
        std::move(geometry)),
      layout_(std::move(layout))
  {
  }

protected:
  std::unique_ptr<PointEvaluator<Real>> CreateEvaluator(
    std::shared_ptr<const pcms::detail::OmegaHMeshBinding> bound,
    DifferentialOperator /*op*/) const override
  {
    return std::make_unique<OmegaHLagrangePointEvaluator>(layout_, bound);
  }

private:
  std::shared_ptr<const OmegaHLagrangeLayout> layout_;
};

} // namespace pcms::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_LAGRANGE_EVALUATOR_HPP
