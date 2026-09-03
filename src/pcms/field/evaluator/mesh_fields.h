#ifndef PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_EVALUATOR_FACTORY_H
#define PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_EVALUATOR_FACTORY_H

#include "pcms/field/evaluator/mesh_fields_backend.h"
#include "pcms/field/data/mesh_fields.h"
#include "pcms/field/field_evaluator_factory.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/field_data.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/arrays.h"
#include <memory>

namespace pcms
{

// MeshFieldsPointEvaluator<T> implements PointEvaluator<T> for
// MeshFields-backed simplex meshes. Each Evaluate call loads DOF data from the
// FieldData argument into the MeshFieldBackend's shape_field_ via SetData, then
// calls evaluate().
template <typename T,
          typename LayoutPolicy =
            detail::default_layout_for_memory_space_t<DeviceMemorySpace>>
class MeshFieldsPointEvaluator : public PointEvaluator<T, LayoutPolicy>
{
public:
  MeshFieldsPointEvaluator(
    std::shared_ptr<const MeshFieldsAdapterLayout> layout,
    MeshFieldsAdapter2LocalizationHint hint, Real fill_value)
    : layout_(std::move(layout)),
      hint_(std::move(hint)),
      fill_value_(fill_value)
  {
  }

  void Evaluate(
    const Field<T>& field,
    Rank2View<T, DeviceMemorySpace, LayoutPolicy> values) const override
  {
    PCMS_FUNCTION_TIMER;
    PCMS_ALWAYS_ASSERT(values.extent(0) ==
                       hint_.coordinates_.extent(0) + hint_.num_missing_);
    PCMS_ALWAYS_ASSERT(values.extent(1) ==
                       static_cast<size_t>(layout_->GetNumComponents()));
    auto const* mesh_field_data =
      dynamic_cast<const MeshFieldsFieldData<T>*>(&field.GetData());
    if (!mesh_field_data) {
      throw pcms_error(
        "MeshFieldsPointEvaluator::Evaluate: incompatible FieldData type");
    }

    auto eval_results = mesh_field_data->GetMeshFieldBackend()->evaluate(
      hint_.coordinates_d_, hint_.offsets_d_);

    int ncomp = layout_->GetNumComponents();
    int num_eval = static_cast<int>(eval_results.extent(0));
    Kokkos::parallel_for(
      "CopyEvalResultsToValues",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_eval, ncomp}),
      KOKKOS_CLASS_LAMBDA(LO i, int c) {
        values(hint_.indices_d_(i), c) = eval_results(i, c);
      });

    if (hint_.num_missing_ > 0 && hint_.mode_ == OutOfBoundsMode::FILL) {
      T fill_val = static_cast<T>(fill_value_);
      int num_missing = static_cast<int>(hint_.num_missing_);
      Kokkos::parallel_for(
        "FillMissingValues",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {num_missing, ncomp}),
        KOKKOS_CLASS_LAMBDA(LO i, int c) {
          values(hint_.missing_indices_d_(i), c) = fill_val;
        });
    }
  }

private:
  std::shared_ptr<const MeshFieldsAdapterLayout> layout_;
  MeshFieldsAdapter2LocalizationHint hint_;
  Real fill_value_;
};

// MeshFieldsEvaluatorFactory<T> implements FieldEvaluatorFactory<T> for
// MeshFields-backed simplex meshes. It owns the spatial search structure and
// the MeshFieldBackend (which holds the internal shape field).
template <typename T>
class MeshFieldsEvaluatorFactory : public FieldEvaluatorFactory<T>
{
public:
  explicit MeshFieldsEvaluatorFactory(
    std::shared_ptr<const MeshFieldsAdapterLayout> layout)
    : layout_(std::move(layout)),
      mesh_(layout_->GetMesh()),
      search_(mesh_, 10, 10)
  {
    if (mesh_.dim() == 3) {
      throw pcms_error("MeshFieldsEvaluatorFactory does not support 3D meshes");
    }
  }

  const FieldLayout& GetLayout() const override { return *layout_; }

  CoordinateSystem GetCoordinateSystem() const override
  {
    return layout_->GetDOFHolderCoordinates().GetCoordinateSystem();
  }

  bool HasDOFHolderCoordinates() const override { return true; }

  bool SupportsNearestBoundary() const override { return false; }

  std::unique_ptr<PointEvaluator<T>> CreatePointEvaluator(
    const EvaluationRequest& request) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto coords = request.coords;
    const auto policy = request.policy;
    if (coords.GetCoordinateSystem() != GetCoordinateSystem()) {
      throw pcms_error(
        "MeshFieldsEvaluatorFactory: coordinate system mismatch");
    }
    if (policy.mode == OutOfBoundsMode::NEAREST_BOUNDARY) {
      throw pcms_error(
        "MeshFieldsEvaluatorFactory: NearestBoundary is not supported");
    }

    auto coordinates = coords.GetValues();
    Kokkos::View<Real* [2], DeviceMemorySpace> coords_search(
      "coords_search", coordinates.extent(0));
    Kokkos::parallel_for(
      "copy_coords",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
        0, coordinates.extent(0)),
      KOKKOS_LAMBDA(const int i) {
        coords_search(i, 0) = coordinates(i, 0);
        coords_search(i, 1) = coordinates(i, 1);
      });
    auto results = search_(coords_search);
    auto owning_ids = search_.GetOwningElementIds(results);

    MeshFieldsAdapter2LocalizationHint hint(mesh_, results, coords_search,
                                            owning_ids, policy.mode);

    return std::make_unique<MeshFieldsPointEvaluator<T>>(
      layout_, std::move(hint), policy.fill_value);
  }

  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override
  {
    return layout_->GetDOFHolderCoordinates();
  }

private:
  std::shared_ptr<const MeshFieldsAdapterLayout> layout_;
  Omega_h::Mesh& mesh_;
  mutable GridPointSearch2D search_;
};

} // namespace pcms

#endif // PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_EVALUATOR_FACTORY_H
