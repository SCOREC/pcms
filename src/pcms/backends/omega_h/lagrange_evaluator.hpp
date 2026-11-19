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

/// Evaluates an order 0 or 1 Lagrange field at the found points of a binding
/// of a mesh of dimension Dim. Order 0 takes the value of the lowest element
/// adjacent to the entity of the point's record; order 1 interpolates on the
/// vertices of that entity itself, with the weights its reduced coordinates
/// give, and picks no element. (A free function: CUDA does not allow device
/// lambdas in the public virtual function that dispatches on the dimension.)
template <int Dim, typename T, typename Dof, typename Values>
void EvaluateOmegaHLagrange(Omega_h::Mesh& mesh, int order,
                            const pcms::detail::OmegaHMeshBinding& binding,
                            const Dof& dof_data, int n_comp,
                            const Values& values)
{
  const auto found = binding.FoundPoints();
  const auto results = binding.Results<Dim>();
  const auto n_valid = static_cast<LO>(found.extent(0));
  const auto up = UpAdjacencies<Dim>(mesh);
  const auto verts_of = VerticesOfEntities<Dim>(mesh);
  if (order == 0) {
    pcms::detail::ForEach(
      "OmegaHLagrangePointEvaluator::EvaluateOrder0", n_valid,
      KOKKOS_LAMBDA(const LO k) {
        const LO elem = ElementOfResult<Dim>(up, results(found(k)));
        for (int c = 0; c < n_comp; ++c) {
          values(found(k), c) = dof_data(elem, c);
        }
      });
  } else {
    pcms::detail::ForEach(
      "OmegaHLagrangePointEvaluator::EvaluateOrder1", n_valid,
      KOKKOS_LAMBDA(const LO k) {
        LO vertices[Dim + 1];
        Real weights[Dim + 1];
        const int nv = EntityVertexWeights<Dim>(verts_of, results(found(k)),
                                                vertices, weights);
        for (int c = 0; c < n_comp; ++c) {
          T val = T{};
          for (int v = 0; v < nv; ++v) {
            val += static_cast<T>(weights[v]) * dof_data(vertices[v], c);
          }
          values(found(k), c) = val;
        }
      });
  }
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
      binding_(std::move(binding))
  {
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
    Omega_h::Mesh& mesh = const_cast<Omega_h::Mesh&>(layout_->GetMesh());
    const int order = layout_->GetOrder();
    const int n_comp = layout_->GetNumComponents();
    if (mesh.dim() == 2) {
      detail::EvaluateOmegaHLagrange<2, T>(mesh, order, *binding_, dof_data,
                                           n_comp, values);
    } else {
      detail::EvaluateOmegaHLagrange<3, T>(mesh, order, *binding_, dof_data,
                                           n_comp, values);
    }
  }

private:
  std::shared_ptr<const OmegaHLagrangeLayout> layout_;
  std::shared_ptr<const pcms::detail::OmegaHMeshBinding> binding_;
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
