#ifndef PCMS_BACKENDS_OMEGA_H_MESHFIELDS_EVALUATOR_HPP
#define PCMS_BACKENDS_OMEGA_H_MESHFIELDS_EVALUATOR_HPP

#include "pcms/backends/omega_h_meshfields/evaluator_backend.hpp"
#include "pcms/backends/omega_h/geometry.hpp"
#include "pcms/backends/omega_h_meshfields/field_data.hpp"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/utility/detail/profile.hpp"
#include "pcms/utility/arrays.hpp"
#include <memory>
#include <tuple>

namespace pcms::backend::omega_h_meshfields
{

namespace detail
{

/// @brief Groups found points by their element (CSR), for MeshFields.
///
/// A point's element is the lowest element adjacent to the entity of its search
/// record; its local coordinates are the first Dim barycentric coordinates in
/// that element, the reduced coordinates of the record put in the order of the
/// element's vertices (zero on the vertices the entity does not have).
/// @tparam Dim mesh dimension
/// @param found query point of each found point
/// @param results search record of each query point
/// @return offsets into the rows of each element's points (nelems + 1
///         entries), the local coordinates of each row, and the query point of
///         each row
template <int Dim>
std::tuple<Kokkos::View<LO*>, Kokkos::View<Real**>, Kokkos::View<LO*>>
GroupFoundPointsByElement(
  Omega_h::Mesh& mesh, Kokkos::View<const LO*, DeviceMemorySpace> found,
  Kokkos::View<const typename omega_h::PointLocalizationSearch<Dim>::Result*,
               DeviceMemorySpace>
    results)
{
  const auto num_found = static_cast<LO>(found.extent(0));
  const LO nelems = mesh.nelems();
  const auto up = omega_h::detail::UpAdjacencies<Dim>(mesh);
  const auto verts_of = omega_h::detail::VerticesOfEntities<Dim>(mesh);
  using Range = Kokkos::RangePolicy<DeviceMemorySpace::execution_space>;

  Kokkos::View<LO*> elem_counts("elem_counts", nelems);
  pcms::detail::ForEach(
    "CountPerElement", num_found, KOKKOS_LAMBDA(const LO k) {
      const LO elem =
        omega_h::detail::ElementOfResult<Dim>(up, results(found(k)));
      Kokkos::atomic_add(&elem_counts(elem), 1);
    });
  Kokkos::View<LO*> offsets("offsets", nelems + 1);
  LO total = 0;
  Kokkos::parallel_scan(
    "ComputeOffsets", Range(0, nelems),
    KOKKOS_LAMBDA(const LO e, LO& partial, const bool is_final) {
      if (is_final) {
        offsets(e) = partial;
      }
      partial += elem_counts(e);
    },
    total);
  Kokkos::deep_copy(Kokkos::subview(offsets, nelems), total);

  Kokkos::View<Real**> local_coords("local_coords", num_found, Dim);
  Kokkos::View<LO*> indices("indices", num_found);
  pcms::detail::ForEach(
    "FillLocalCoordinates", num_found, KOKKOS_LAMBDA(const LO k) {
      const auto& result = results(found(k));
      const LO elem = omega_h::detail::ElementOfResult<Dim>(up, result);
      const LO count =
        Kokkos::atomic_fetch_add(&elem_counts(elem), static_cast<LO>(-1));
      const LO index = offsets(elem) + count - 1;
      LO vertices[Dim + 1];
      Real weights[Dim + 1];
      const int nv = omega_h::detail::EntityVertexWeights<Dim>(
        verts_of, result, vertices, weights);
      Real bary[Dim + 1] = {};
      for (int j = 0; j < nv; ++j) {
        for (int v = 0; v <= Dim; ++v) {
          if (verts_of[Dim][elem * (Dim + 1) + v] == vertices[j]) {
            bary[v] = weights[j];
          }
        }
      }
      for (int j = 0; j < Dim; ++j) {
        local_coords(index, j) = bary[j];
      }
      indices(index) = found(k);
    });
  return {offsets, local_coords, indices};
}

} // namespace detail

/// @brief Evaluates MeshFields fields at the bound points with the field's
///        MeshFields backend.
template <typename T,
          typename LayoutPolicy =
            pcms::detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class MeshFieldsPointEvaluator : public PointEvaluator<T, LayoutPolicy>
{
public:
  /// @brief Groups the found points by element, with their first `dim`
  ///        barycentric coordinates as the coordinates MeshFields evaluates
  ///        at.
  /// @param layout MeshFields layout of the evaluated fields
  /// @param binding the search records of the query points
  MeshFieldsPointEvaluator(
    std::shared_ptr<const MeshFieldsAdapterLayout> layout,
    std::shared_ptr<const omega_h::OmegaHMeshBinding> binding)
    : PointEvaluator<T, LayoutPolicy>(binding, DifferentialOperator::Value,
                                      layout),
      layout_(std::move(layout))
  {
    auto& mesh = layout_->GetMesh();
    if (mesh.dim() == 2) {
      std::tie(offsets_, local_coords_, indices_) =
        detail::GroupFoundPointsByElement<2>(mesh, binding->FoundPoints(),
                                             binding->Results<2>());
    } else {
      std::tie(offsets_, local_coords_, indices_) =
        detail::GroupFoundPointsByElement<3>(mesh, binding->FoundPoints(),
                                             binding->Results<3>());
    }
  }

  // Public because CUDA forbids extended lambdas in non-public member
  // functions. User should call Evaluate
  // TODO: add passkey or other mechanism to avoid user calling this
  void EvaluateImpl(
    const Field<T>& field,
    ValueView<T, DeviceMemorySpace, LayoutPolicy> out) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto values = out.GetValues();
    auto const* mesh_field_data =
      dynamic_cast<const MeshFieldsFieldData<T>*>(&field.GetData());
    if (!mesh_field_data) {
      throw pcms_error(
        "MeshFieldsPointEvaluator::Evaluate: incompatible FieldData type");
    }

    auto eval_results =
      mesh_field_data->GetMeshFieldBackend()->evaluate(local_coords_, offsets_);

    int ncomp = layout_->GetNumComponents();
    int num_eval = static_cast<int>(eval_results.extent(0));
    const auto indices = indices_;
    pcms::detail::ForEach(
      "CopyEvalResultsToValues", num_eval, KOKKOS_LAMBDA(const LO i) {
        for (int c = 0; c < ncomp; ++c) {
          values(indices(i), c) = eval_results(i, c);
        }
      });
  }

private:
  std::shared_ptr<const MeshFieldsAdapterLayout> layout_;
  // the found points grouped by element: offsets into local_coords_ of each
  // element's points (number of elements + 1 entries), the points'
  // coordinates in their element, and the query point of each row
  Kokkos::View<LO*> offsets_;
  Kokkos::View<Real**> local_coords_;
  Kokkos::View<LO*> indices_;
};

/// @brief Point evaluation of MeshFields fields on a 2D simplex mesh.
/// @tparam T field value type
template <typename T>
class MeshFieldsEvaluatorFactory
  : public TypedPointEvaluatorFactory<T, omega_h::OmegaHMeshBinding>
{
public:
  /// @param geometry geometry of the layout's mesh
  /// @param layout MeshFields layout of the evaluated fields
  MeshFieldsEvaluatorFactory(
    std::shared_ptr<const omega_h::OmegaHGeometry> geometry,
    std::shared_ptr<const MeshFieldsAdapterLayout> layout)
    : TypedPointEvaluatorFactory<T, omega_h::OmegaHMeshBinding>(
        std::move(geometry)),
      layout_(std::move(layout))
  {
    if (layout_->GetMesh().dim() == 3) {
      throw pcms_error("MeshFieldsEvaluatorFactory does not support 3D meshes");
    }
  }

protected:
  std::unique_ptr<PointEvaluator<T>> CreateEvaluator(
    std::shared_ptr<const omega_h::OmegaHMeshBinding> bound,
    DifferentialOperator /*op*/) const override
  {
    return std::make_unique<MeshFieldsPointEvaluator<T>>(layout_,
                                                         std::move(bound));
  }

private:
  std::shared_ptr<const MeshFieldsAdapterLayout> layout_;
};

} // namespace pcms::backend::omega_h_meshfields

#endif // PCMS_BACKENDS_OMEGA_H_MESHFIELDS_EVALUATOR_HPP
