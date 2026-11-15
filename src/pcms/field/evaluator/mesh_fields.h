#ifndef PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_EVALUATOR_FACTORY_H
#define PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_EVALUATOR_FACTORY_H

#include "pcms/field/evaluator/mesh_fields_backend.h"
#include "pcms/field/evaluator/omega_h_geometry.hpp"
#include "pcms/field/data/mesh_fields.h"
#include "pcms/field/point_evaluator_factory.hpp"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/field_data.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/arrays.h"
#include <memory>
#include <tuple>

namespace pcms
{

namespace detail
{

/// @brief Groups found points by their element (CSR), for MeshFields.
/// @param found query point of each found point
/// @param elem_ids element of each found point
/// @param bary barycentric coordinates of each found point ([point][vertex])
/// @param dim mesh dimension; the first `dim` barycentric coordinates are the
///        local coordinates
/// @param nelems number of elements
/// @return offsets into the rows of each element's points (nelems + 1
///         entries), the local coordinates of each row, and the query point of
///         each row
inline std::tuple<Kokkos::View<LO*>, Kokkos::View<Real**>, Kokkos::View<LO*>>
GroupFoundPointsByElement(Kokkos::View<const LO*, DeviceMemorySpace> found,
                          Kokkos::View<const LO*, DeviceMemorySpace> elem_ids,
                          Kokkos::View<const Real**, DeviceMemorySpace> bary,
                          int dim, LO nelems)
{
  const auto num_found = static_cast<LO>(found.extent(0));
  using Range = Kokkos::RangePolicy<DeviceMemorySpace::execution_space>;

  Kokkos::View<LO*> elem_counts("elem_counts", nelems);
  detail::ForEach(
    "CountPerElement", num_found, KOKKOS_LAMBDA(const LO k) {
      Kokkos::atomic_add(&elem_counts(elem_ids(k)), 1);
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

  Kokkos::View<Real**> local_coords("local_coords", num_found, dim);
  Kokkos::View<LO*> indices("indices", num_found);
  detail::ForEach(
    "FillLocalCoordinates", num_found, KOKKOS_LAMBDA(const LO k) {
      const LO elem = elem_ids(k);
      const LO count =
        Kokkos::atomic_fetch_add(&elem_counts(elem), static_cast<LO>(-1));
      const LO index = offsets(elem) + count - 1;
      for (int j = 0; j < dim; ++j) {
        local_coords(index, j) = bary(k, j);
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
            detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class MeshFieldsPointEvaluator : public PointEvaluator<T, LayoutPolicy>
{
public:
  /// @brief Groups the found points by element, with their first `dim`
  ///        barycentric coordinates as the coordinates MeshFields evaluates
  ///        at.
  /// @param layout MeshFields layout of the evaluated fields
  /// @param binding elements and barycentric coordinates of the query points
  MeshFieldsPointEvaluator(
    std::shared_ptr<const MeshFieldsAdapterLayout> layout,
    std::shared_ptr<const OmegaHMeshBinding> binding)
    : PointEvaluator<T, LayoutPolicy>(binding, DifferentialOperator::Value,
                                      layout),
      layout_(std::move(layout))
  {
    auto& mesh = layout_->GetMesh();
    std::tie(offsets_, local_coords_, indices_) =
      detail::GroupFoundPointsByElement(
        binding->FoundPoints(), binding->ElementIds(), binding->Barycentric(),
        mesh.dim(), mesh.nelems());
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
    detail::ForEach(
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
  : public TypedPointEvaluatorFactory<T, OmegaHMeshBinding>
{
public:
  /// @param geometry geometry of the layout's mesh
  /// @param layout MeshFields layout of the evaluated fields
  MeshFieldsEvaluatorFactory(
    std::shared_ptr<const OmegaHGeometry> geometry,
    std::shared_ptr<const MeshFieldsAdapterLayout> layout)
    : TypedPointEvaluatorFactory<T, OmegaHMeshBinding>(std::move(geometry)),
      layout_(std::move(layout))
  {
    if (layout_->GetMesh().dim() == 3) {
      throw pcms_error("MeshFieldsEvaluatorFactory does not support 3D meshes");
    }
  }

protected:
  std::unique_ptr<PointEvaluator<T>> CreateEvaluator(
    std::shared_ptr<const OmegaHMeshBinding> bound,
    DifferentialOperator /*op*/) const override
  {
    return std::make_unique<MeshFieldsPointEvaluator<T>>(layout_,
                                                         std::move(bound));
  }

private:
  std::shared_ptr<const MeshFieldsAdapterLayout> layout_;
};

} // namespace pcms

#endif // PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_EVALUATOR_FACTORY_H
