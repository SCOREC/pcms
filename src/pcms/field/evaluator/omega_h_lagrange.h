#ifndef PCMS_OMEGA_H_LAGRANGE_EVALUATOR_FACTORY_H
#define PCMS_OMEGA_H_LAGRANGE_EVALUATOR_FACTORY_H

#include "pcms/field/evaluator/omega_h_geometry.hpp"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/field_data.h"
#include "pcms/field/data/omega_h.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/types.h"

#include <Omega_h_mesh.hpp>
#include <Kokkos_Core.hpp>
#include <memory>
#include <utility>

namespace pcms
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
} // namespace detail

// OmegaHLagrangePointEvaluator<T> implements PointEvaluator<T> for simplex
// meshes backed by Omega_h, using the element IDs and barycentric coordinates
// of the bound source.
//
// Output shape: [num_query_points][num_components].
//
// The values are read from the field's Omega_h tag
// (OmegaHFieldData::GetTagName). Its owned entries are always current; its
// ghost entries, which elements on the partition boundary read, are correct
// only after the field's last write was followed by Synchronize.
template <typename T,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class OmegaHLagrangePointEvaluator : public PointEvaluator<T, LayoutPolicy>
{
public:
  /// @param layout Omega_h Lagrange layout of the evaluated fields
  /// @param binding elements and barycentric coordinates of the query points
  OmegaHLagrangePointEvaluator(
    std::shared_ptr<const OmegaHLagrangeLayout> layout,
    std::shared_ptr<const OmegaHMeshBinding> binding)
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
    const auto found = binding_->FoundPoints();
    const auto elem_ids = binding_->ElementIds();
    const auto bary = binding_->Barycentric();
    LO n_valid = static_cast<LO>(found.extent(0));
    int n_comp = layout_->GetNumComponents();

    if (layout_->GetOrder() == 0) {
      Kokkos::parallel_for(
        "OmegaHLagrangePointEvaluator::EvaluateOrder0",
        Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_valid),
        KOKKOS_LAMBDA(LO k) {
          for (int c = 0; c < n_comp; ++c) {
            values(found(k), c) = dof_data(elem_ids(k), c);
          }
        });

    } else {
      // Order-1: barycentric interpolation over element vertices
      Omega_h::Mesh& mesh = const_cast<Omega_h::Mesh&>(layout_->GetMesh());
      int mesh_dim = mesh.dim();
      int nvpe = mesh_dim + 1;
      auto elem_verts = mesh.ask_elem_verts();

      Kokkos::parallel_for(
        "OmegaHLagrangePointEvaluator::EvaluateOrder1",
        Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_valid),
        KOKKOS_LAMBDA(LO k) {
          const LO elem = elem_ids(k);
          for (int c = 0; c < n_comp; ++c) {
            T val = T{};
            for (int v = 0; v < nvpe; ++v) {
              LO vert = elem_verts[elem * nvpe + v];
              val += static_cast<T>(bary(k, v)) * dof_data(vert, c);
            }
            values(found(k), c) = val;
          }
        });
    }
  }

private:
  std::shared_ptr<const OmegaHLagrangeLayout> layout_;
  std::shared_ptr<const OmegaHMeshBinding> binding_;
};

// ---------------------------------------------------------------------------
// EvaluatorFactory
// ---------------------------------------------------------------------------

/// @brief Point evaluation of Lagrange fields on an Omega_h simplex mesh.
/// @tparam T field value type
template <typename T>
class OmegaHLagrangeEvaluatorFactory
  : public TypedPointEvaluatorFactory<T, OmegaHMeshBinding>
{
public:
  /// @param geometry geometry of the layout's mesh
  /// @param layout Omega_h Lagrange layout of the evaluated fields
  OmegaHLagrangeEvaluatorFactory(
    std::shared_ptr<const OmegaHGeometry> geometry,
    std::shared_ptr<const OmegaHLagrangeLayout> layout)
    : TypedPointEvaluatorFactory<T, OmegaHMeshBinding>(std::move(geometry)),
      layout_(std::move(layout))
  {
  }

protected:
  std::unique_ptr<PointEvaluator<T>> CreateEvaluator(
    std::shared_ptr<const OmegaHMeshBinding> bound,
    DifferentialOperator /*op*/) const override
  {
    return std::make_unique<OmegaHLagrangePointEvaluator<T>>(layout_,
                                                             std::move(bound));
  }

private:
  std::shared_ptr<const OmegaHLagrangeLayout> layout_;
};

} // namespace pcms

#endif // PCMS_OMEGA_H_LAGRANGE_EVALUATOR_FACTORY_H
