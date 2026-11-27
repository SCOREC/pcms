#ifndef PCMS_BACKENDS_OMEGA_H_LAGRANGE_EVALUATOR_HPP
#define PCMS_BACKENDS_OMEGA_H_LAGRANGE_EVALUATOR_HPP

#include "pcms/backends/omega_h/omega_h_geometry.hpp"
#include "pcms/backends/omega_h/lagrange_layout.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/backends/omega_h/field_data.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/arrays.hpp"
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

} // namespace detail

// OmegaHLagrangePointEvaluator<T> implements PointEvaluator<T> for simplex
// meshes backed by Omega_h, using the search records of the BoundPoints.
//
// Output shape: [num_query_points][num_components].
//
// The values are read from the field's Omega_h tag
// (OmegaHFieldData::GetTagName). Its owned entries are always current; its
// ghost entries, which elements on the partition boundary read, are correct
// only after the field's last write was followed by Synchronize.
template <typename T,
          typename LayoutPolicy =
            pcms::detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class OmegaHLagrangePointEvaluator : public PointEvaluator<T, LayoutPolicy>
{
public:
  /// @param layout Omega_h Lagrange layout of the evaluated fields
  /// @param binding the search records of the query points
  OmegaHLagrangePointEvaluator(
    std::shared_ptr<const OmegaHLagrangeLayout> layout,
    std::shared_ptr<const pcms::detail::OmegaHMeshBinding> binding)
    : PointEvaluator<T, LayoutPolicy>(binding, DifferentialOperator::Value,
                                      layout),
      layout_(std::move(layout)),
      binding_(std::move(binding)),
      stride_(layout_->GetOrder() == 0 ? 1 : layout_->GetMesh().dim() + 1),
      holders_("omega_h_lagrange_holders",
               static_cast<std::size_t>(binding_->NumPoints()) * stride_),
      weights_("omega_h_lagrange_weights", holders_.extent(0))
  {
    Omega_h::Mesh& mesh = const_cast<Omega_h::Mesh&>(layout_->GetMesh());
    if (mesh.dim() == 2) {
      detail::FillOmegaHLagrangeRows<2>(mesh, layout_->GetOrder(), *binding_,
                                        holders_, weights_);
    } else {
      detail::FillOmegaHLagrangeRows<3>(mesh, layout_->GetOrder(), *binding_,
                                        holders_, weights_);
    }
  }

  // Public because CUDA forbids extended lambdas in non-public member
  // functions. User call Evaluate.
  // TODO: add passkey or other mechanism to protect this
  void EvaluateImpl(
    const Field<T>& field,
    ValueView<T, DeviceMemorySpace, LayoutPolicy> out) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto values = out.GetValues();
    // Elements read the values of all their DOF holders, ghosts included, so
    // read the field's Omega_h tag rather than the owned values.
    const auto& data = dynamic_cast<const OmegaHFieldData<T>&>(field.GetData());
    using TagValue = typename OmegaHFieldData<T>::tag_value_type;
    const auto dof_data = detail::OmegaHLagrangeTagValues<T, TagValue>{
      layout_->GetMesh().get_array<TagValue>(layout_->GetDOFHolderEntityDim(),
                                             data.GetTagName()),
      layout_->GetNumComponents()};
    const auto found = binding_->FoundPoints();
    const int n_comp = layout_->GetNumComponents();
    const int stride = stride_;
    const auto holders = holders_;
    const auto weights = weights_;
    pcms::detail::ForEach(
      "OmegaHLagrangePointEvaluator::Evaluate",
      static_cast<LO>(found.extent(0)), KOKKOS_LAMBDA(const LO k) {
        const LO p = found(k);
        for (int c = 0; c < n_comp; ++c) {
          T val = T{};
          for (int v = p * stride; v < (p + 1) * stride; ++v) {
            val += static_cast<T>(weights(v)) * dof_data(holders(v), c);
          }
          values(p, c) = val;
        }
      });
  }

  [[nodiscard]] bool SupportsTranspose() const noexcept override
  {
    return true;
  }

  // Public for the same reason as EvaluateImpl.
  void EvaluateTransposeImpl(
    Rank2View<const T, DeviceMemorySpace, LayoutPolicy> values,
    Field<T>& functional, T beta) const override
  {
    PCMS_FUNCTION_TIMER;
    using TagValue = utility::OmegaHValue_t<T>;
    const auto found = binding_->FoundPoints();
    const int n_comp = layout_->GetNumComponents();
    const int stride = stride_;
    const auto holders = holders_;
    const auto weights = weights_;
    Omega_h::Write<TagValue> sums(layout_->GetNumDofHolderWithGhosts() * n_comp,
                                  0, "transpose_sums");
    pcms::detail::ForEach(
      "OmegaHLagrangePointEvaluator::EvaluateTranspose",
      static_cast<LO>(found.extent(0)), KOKKOS_LAMBDA(const LO k) {
        const LO p = found(k);
        for (int v = p * stride; v < (p + 1) * stride; ++v) {
          for (int c = 0; c < n_comp; ++c) {
            Kokkos::atomic_add(
              &sums[holders(v) * n_comp + c],
              static_cast<TagValue>(weights(v) * values(p, c)));
          }
        }
      });
    // the copies of a holder on other ranks add into its owner
    const auto reduced = layout_->GetMesh().reduce_array(
      layout_->GetDOFHolderEntityDim(), Omega_h::Read<TagValue>(sums), n_comp,
      OMEGA_H_SUM);
    dynamic_cast<OmegaHFieldData<T>&>(functional.GetData())
      .SetFromOwners(reduced, beta);
  }

  [[nodiscard]] bool SupportsMatrix() const noexcept override { return true; }

  PointMatrix<T> GetMatrixImpl() const override
  {
    const LO n = binding_->NumPoints();
    const int stride = stride_;
    const auto holders = holders_;
    const auto weights = weights_;
    const auto gids = layout_->GetGidsWithGhosts();
    Kokkos::View<LO*, DeviceMemorySpace> offsets("offsets", n + 1);
    Kokkos::View<GO*, DeviceMemorySpace> columns("columns", holders.extent(0));
    Kokkos::View<T*, DeviceMemorySpace> entries("entries", holders.extent(0));
    Kokkos::parallel_for(
      "OmegaHLagrangePointEvaluator::GetMatrix",
      Kokkos::RangePolicy<DefaultExecutionSpace>(0, n + 1),
      KOKKOS_LAMBDA(const LO p) {
        offsets(p) = p * stride;
        for (int v = p * stride; p < n && v < (p + 1) * stride; ++v) {
          columns(v) = gids[holders(v)];
          entries(v) = static_cast<T>(weights(v));
        }
      });
    return PointMatrix<T>{offsets, columns, entries};
  }

private:
  std::shared_ptr<const OmegaHLagrangeLayout> layout_;
  std::shared_ptr<const pcms::detail::OmegaHMeshBinding> binding_;
  // the evaluation matrix: stride_ entries per point
  int stride_;
  Kokkos::View<LO*, DeviceMemorySpace> holders_;
  Kokkos::View<Real*, DeviceMemorySpace> weights_;
};

// ---------------------------------------------------------------------------
// EvaluatorFactory
// ---------------------------------------------------------------------------

/// @brief Point evaluation of Lagrange fields on an Omega_h simplex mesh.
/// @tparam T field value type
template <typename T>
class OmegaHLagrangeEvaluatorFactory
  : public pcms::detail::TypedPointEvaluatorFactory<
      T, pcms::detail::OmegaHMeshBinding>
{
public:
  /// @param geometry geometry of the layout's mesh
  /// @param layout Omega_h Lagrange layout of the evaluated fields
  OmegaHLagrangeEvaluatorFactory(
    std::shared_ptr<const OmegaHGeometry> geometry,
    std::shared_ptr<const OmegaHLagrangeLayout> layout)
    : pcms::detail::TypedPointEvaluatorFactory<T,
                                               pcms::detail::OmegaHMeshBinding>(
        std::move(geometry)),
      layout_(std::move(layout))
  {
  }

protected:
  std::unique_ptr<PointEvaluator<T>> CreateEvaluator(
    std::shared_ptr<const pcms::detail::OmegaHMeshBinding> bound,
    DifferentialOperator /*op*/) const override
  {
    return std::make_unique<OmegaHLagrangePointEvaluator<T>>(layout_,
                                                             std::move(bound));
  }

private:
  std::shared_ptr<const OmegaHLagrangeLayout> layout_;
};

} // namespace pcms::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_LAGRANGE_EVALUATOR_HPP
