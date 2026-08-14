#ifndef PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_REDUCED_QUINTIC_EVALUATOR_FACTORY_H
#define PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_REDUCED_QUINTIC_EVALUATOR_FACTORY_H

#include "pcms/field/evaluator/mesh_fields_backend.h"
#include "pcms/field/evaluator/mesh_fields.h"
#include "pcms/field/data/mesh_fields_reduced_quintic.h"
#include "pcms/field/field_evaluator_factory.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/field/point_evaluator.h"
#include "pcms/field/field_data.h"
#include "pcms/utility/profile.h"
#include "pcms/utility/arrays.h"
#include <memory>

namespace pcms
{

// MeshFieldsReducedQuinticEvaluatorFactory<T> implements
// FieldEvaluatorFactory<T> for MeshFields-backed reduced quintic elements on
// simplex meshes. It owns the spatial search structure and a
// MeshFieldBackend<T> which internally creates the shape field with 6 DOFs
// per vertex for the reduced quintic element.
//
// The reduced quintic element uses 6 DOFs per vertex:
//   DOF 0: function value
//   DOF 1: df/dx
//   DOF 2: df/dy
//   DOF 3: d2f/dx2
//   DOF 4: d2f/dxy
//   DOF 5: d2f/dy2
// This provides C1 continuity on triangle meshes.
//
// Note: The PointEvaluator is shared with the standard Lagrange backend
// (MeshFieldsPointEvaluator) since the evaluation strategy is already
// encapsulated in the MeshFieldBackend via the policy-based design.
template <typename T>
class MeshFieldsReducedQuinticEvaluatorFactory
  : public FieldEvaluatorFactory<T>
{
public:
  explicit MeshFieldsReducedQuinticEvaluatorFactory(
    std::shared_ptr<const MeshFieldsAdapterLayout> layout)
    : layout_(std::move(layout)),
      mesh_(layout_->GetMesh()),
      search_(mesh_, 10, 10)
  {
    if (mesh_.dim() == 3) {
      throw pcms_error(
        "MeshFieldsReducedQuinticEvaluatorFactory does not support 3D meshes");
    }
    if (layout_->GetNumComponents() != 1) {
      throw pcms_error(
        "MeshFieldsReducedQuinticEvaluatorFactory only supports "
        "single-component fields");
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

    auto coordinates = coords.GetCoordinates();
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

    MeshFieldsAdapter2LocalizationHint hint(mesh_, results, policy.mode);

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

#endif // PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_REDUCED_QUINTIC_EVALUATOR_FACTORY_H