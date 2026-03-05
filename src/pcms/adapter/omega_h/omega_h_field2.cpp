#include "omega_h_field2.h"
#include "omega_h_field.h"
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include <memory>

namespace pcms
{
template <int Dim, int Order>
class MeshFieldBackendImpl : public MeshFieldBackend
{
public:
  MeshFieldBackendImpl(Omega_h::Mesh& mesh)
    : mesh_(mesh),
      mesh_field_(mesh),
      shape_field_(mesh_field_.template CreateLagrangeField<Real, Order, 1>())
  {
  }

  Kokkos::View<Real* [1]> evaluate(Kokkos::View<Real**> localCoords,
                                   Kokkos::View<LO*> offsets) const override
  {
    auto self = const_cast<MeshFieldBackendImpl<Dim, Order>*>(this);
    return self->mesh_field_.triangleLocalPointEval(localCoords, offsets,
                                                    shape_field_);
  }

  void SetData(Rank1View<const Real, HostMemorySpace> data, size_t num_nodes,
               size_t num_components, int dim) override
  {
    size_t stride = num_nodes * num_components;
    auto topo = static_cast<MeshField::Mesh_Topology>(dim);
    Kokkos::View<Real*, DefaultExecutionSpace::memory_space> data_d(
      "data_d", data.size());
    Kokkos::deep_copy(data_d, Kokkos::View<const Real*, HostMemorySpace>(
                                data.data_handle(), data.size()));
    Kokkos::parallel_for(
      mesh_.nents(dim), KOKKOS_CLASS_LAMBDA(size_t ent) {
        for (size_t n = 0; n < num_nodes; ++n) {
          for (size_t c = 0; c < num_components; ++c) {
            shape_field_(ent, n, c, topo) =
              data_d[ent * stride + n * num_components + c];
          }
        }
      });
  }

  void GetData(Rank1View<Real, HostMemorySpace> data, size_t num_nodes,
               size_t num_components, int dim) const override
  {
    size_t stride = num_nodes * num_components;
    auto topo = static_cast<MeshField::Mesh_Topology>(dim);
    Kokkos::View<Real*, DefaultExecutionSpace::memory_space> data_d(
      "data_d", data.size());
    Kokkos::parallel_for(
      mesh_.nents(dim), KOKKOS_CLASS_LAMBDA(size_t ent) {
        for (size_t n = 0; n < num_nodes; ++n) {
          for (size_t c = 0; c < num_components; ++c) {
            data_d[ent * stride + n * num_components + c] =
              shape_field_(ent, n, c, topo);
          }
        }
      });
    Kokkos::deep_copy(
      Kokkos::View<Real*, HostMemorySpace>(data.data_handle(), data.size()),
      data_d);
  }

private:
  Omega_h::Mesh& mesh_;
  MeshField::OmegahMeshField<DefaultExecutionSpace, Dim> mesh_field_;
  using ShapeField =
    decltype(mesh_field_.template CreateLagrangeField<Real, Order, 1>());
  ShapeField shape_field_;
};

struct ComputeOffsetsFunctor
{
  Kokkos::View<LO*, HostMemorySpace> offsets_;
  Kokkos::View<LO*, HostMemorySpace> elem_counts_;

  ComputeOffsetsFunctor(Kokkos::View<LO*, HostMemorySpace> offsets,
                        Kokkos::View<LO*, HostMemorySpace> elem_counts)
    : offsets_(offsets), elem_counts_(elem_counts)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i, LO& partial, bool is_final) const
  {
    if (is_final) {
      offsets_(i) = partial;
    }
    partial += elem_counts_(i);
  }
};

struct CountPointsPerElementFunctor
{
  Kokkos::View<LO*> elem_counts_;
  Kokkos::View<GridPointSearch::Result*> search_results_;

  CountPointsPerElementFunctor(
    Kokkos::View<LO*> elem_counts,
    Kokkos::View<GridPointSearch::Result*> search_results)
    : elem_counts_(elem_counts), search_results_(search_results)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    auto [dim, elem_idx, coord] = search_results_(i);
    Kokkos::atomic_add(&elem_counts_(elem_idx), 1);
  }
};

struct FillCoordinatesAndIndicesFunctor
{
  Omega_h::Mesh& mesh_;
  Kokkos::View<LO*> elem_counts_;
  Kokkos::View<LO*> offsets_;
  Kokkos::View<Real**> coordinates_;
  Kokkos::View<LO*> indices_;
  Kokkos::View<GridPointSearch::Result*> search_results_;
  Omega_h::Int dim_;

  FillCoordinatesAndIndicesFunctor(
    Omega_h::Mesh& mesh, Kokkos::View<LO*> elem_counts,
    Kokkos::View<LO*> offsets, Kokkos::View<Real**> coordinates,
    Kokkos::View<LO*> indices,
    Kokkos::View<GridPointSearch::Result*> search_results)
    : mesh_(mesh),
      elem_counts_(elem_counts),
      offsets_(offsets),
      coordinates_(coordinates),
      indices_(indices),
      search_results_(search_results),
      dim_(mesh.dim())
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    auto [dim, elem_idx, coord] = search_results_(i);
    // disable the host assertion macro for device code
    // currently don't handle case where point is on a boundary
    // PCMS_ALWAYS_ASSERT(static_cast<int>(dim) == mesh_.dim());
    // element should be inside the domain (positive)
    // PCMS_ALWAYS_ASSERT(elem_idx >= 0 && elem_idx < mesh_.nelems());
    LO count = Kokkos::atomic_sub_fetch(&elem_counts_(elem_idx), 1);
    LO index = offsets_(elem_idx) + count - 1;
    for (int j = 0; j < (dim_ + 1); ++j) {
      coordinates_(index, j) = coord[j];
    }
    indices_(index) = i;
  }
};

struct OmegaHField2LocalizationHint
{
  OmegaHField2LocalizationHint(
    Omega_h::Mesh& mesh,
    Kokkos::View<GridPointSearch::Result*, HostMemorySpace> search_results,
    OutOfBoundsMode mode)
    : mode_(mode), num_valid_(0), num_missing_(0)
  {
    // First pass: count valid and invalid points
    std::vector<size_t> valid_point_indices;
    std::vector<size_t> missing_point_indices;

    if (mode_ == OutOfBoundsMode::ERROR) {
      // Error mode - throw error immediately if any point is out of bounds
      for (size_t i = 0; i < search_results.size(); ++i) {
        auto [dim, elem_idx, coord] = search_results(i);
        bool is_missing =
          (static_cast<int>(dim) != mesh.dim()) || (elem_idx < 0);
        PCMS_ALWAYS_ASSERT(!is_missing && "Points found outside mesh domain");
        valid_point_indices.push_back(i);
      }
    } else {
      // Other modes - collect valid and missing points separately
      for (size_t i = 0; i < search_results.size(); ++i) {
        auto [dim, elem_idx, coord] = search_results(i);
        bool is_missing =
          (static_cast<int>(dim) != mesh.dim()) || (elem_idx < 0);
        if (is_missing) {
          missing_point_indices.push_back(i);
        } else {
          valid_point_indices.push_back(i);
        }
      }
    }

    num_valid_ = valid_point_indices.size();
    num_missing_ = missing_point_indices.size();

    // Handle missing points based on mode
    if (num_missing_ > 0 && mode_ == OutOfBoundsMode::NEAREST_BOUNDARY) {
      PCMS_ALWAYS_ASSERT(false && "NEAREST_BOUNDARY mode not implemented yet");
    }

    // Allocate arrays for valid points only
    offsets_ = Kokkos::View<LO*, HostMemorySpace>("offsets", mesh.nelems() + 1);
    coordinates_ = Kokkos::View<Real**, HostMemorySpace>(
      "coordinates", num_valid_, mesh.dim() + 1);
    indices_ = Kokkos::View<LO*, HostMemorySpace>("indices", num_valid_);

    // Store missing point indices
    if (num_missing_ > 0) {
      missing_indices_ =
        Kokkos::View<LO*, HostMemorySpace>("missing_indices", num_missing_);
      for (size_t i = 0; i < num_missing_; ++i) {
        missing_indices_(i) = static_cast<LO>(missing_point_indices[i]);
      }
    }

    // Count points per element (valid points only)
    Kokkos::View<LO*, HostMemorySpace> elem_counts("elem_counts",
                                                   mesh.nelems());
    for (size_t i = 0; i < num_valid_; ++i) {
      auto [dim, elem_idx, coord] = search_results(valid_point_indices[i]);
      elem_counts[elem_idx] += 1;
    }

    // Compute offsets
    LO total;
    ComputeOffsetsFunctor functor(offsets_, elem_counts);
    Kokkos::parallel_scan(
      "ComputeOffsets",
      Kokkos::RangePolicy<HostMemorySpace::execution_space>(0, mesh.nelems()),
      functor, total);
    offsets_(mesh.nelems()) = total;

    // Fill coordinates and indices for valid points
    for (size_t i = 0; i < num_valid_; ++i) {
      size_t orig_idx = valid_point_indices[i];
      auto [dim, elem_idx, coord] = search_results(orig_idx);
      elem_counts(elem_idx) -= 1;
      LO index = offsets_(elem_idx) + elem_counts(elem_idx);
      for (int j = 0; j < (mesh.dim() + 1); ++j) {
        coordinates_(index, j) = coord[j];
      }
      indices_(index) = static_cast<LO>(orig_idx);
    }
  }

  OutOfBoundsMode mode_;
  size_t num_valid_;
  size_t num_missing_;

  // offsets is the number of points in each element
  Kokkos::View<LO*, HostMemorySpace> offsets_;
  // coordinates are the parametric coordinates of each point
  Kokkos::View<Real**, HostMemorySpace> coordinates_;
  // indices are the index of the original point (for valid points)
  Kokkos::View<LO*, HostMemorySpace> indices_;
  // indices of points not found in mesh
  Kokkos::View<LO*, HostMemorySpace> missing_indices_;
};

/*
 * Field
 */
OmegaHField2::OmegaHField2(const OmegaHFieldLayout& layout)
  : layout_(layout),
    mesh_(layout.GetMesh()),
    search_(mesh_, 10, 10),
    dof_holder_data_("dof_holder_data", static_cast<size_t>(layout.OwnedSize()))
{
  auto nodes_per_dim = layout.GetNodesPerDim();
  if (nodes_per_dim[2] == 0 && nodes_per_dim[3] == 0) {
    if (nodes_per_dim[0] == 1 && nodes_per_dim[1] == 0) {
      switch (mesh_.dim()) {
        case 1:
          mesh_field_ = std::make_unique<MeshFieldBackendImpl<1, 1>>(mesh_);
          break;
        case 2:
          mesh_field_ = std::make_unique<MeshFieldBackendImpl<2, 1>>(mesh_);
          break;
        default: break; // backend is null
      }
    } else if (nodes_per_dim[0] == 1 && nodes_per_dim[1] == 1) {
      switch (mesh_.dim()) {
        case 2:
          mesh_field_ = std::make_unique<MeshFieldBackendImpl<2, 2>>(mesh_);
          break;
        case 3:
          mesh_field_ = std::make_unique<MeshFieldBackendImpl<3, 2>>(mesh_);
          break;
        default: break; // backend is null
      }
    }
  }
}

Rank1View<const Real, HostMemorySpace> OmegaHField2::GetDOFHolderData() const
{
  PCMS_FUNCTION_TIMER;
  auto nodes_per_dim = layout_.GetNodesPerDim();
  auto num_components = layout_.GetNumComponents();
  size_t offset = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    if (nodes_per_dim[i]) {
      size_t len = static_cast<size_t>(mesh_.nents(i)) *
                   static_cast<size_t>(nodes_per_dim[i]) *
                   static_cast<size_t>(num_components);
      Rank1View<Real, HostMemorySpace> subspan{
        std::data(dof_holder_data_) + offset, len};
      mesh_field_->GetData(subspan, nodes_per_dim[i], num_components, i);
      offset += len;
    }
  }

  return make_const_array_view(dof_holder_data_);
};

void OmegaHField2::SetDOFHolderData(Rank1View<const Real, HostMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;

  auto nodes_per_dim = layout_.GetNodesPerDim();
  auto num_components = layout_.GetNumComponents();
  PCMS_ALWAYS_ASSERT(static_cast<LO>(data.size()) ==
                     layout_.GetNumOwnedDofHolder() * num_components);
  size_t offset = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    if (nodes_per_dim[i]) {
      size_t len = static_cast<size_t>(mesh_.nents(i)) *
                   static_cast<size_t>(nodes_per_dim[i]) *
                   static_cast<size_t>(num_components);
      Rank1View<const Real, HostMemorySpace> subspan{
        data.data_handle() + offset, len};
      mesh_field_->SetData(subspan, nodes_per_dim[i], num_components, i);
      offset += len;
    }
  }
}

LocalizationHint OmegaHField2::GetLocalizationHint(
  CoordinateView<HostMemorySpace> coordinate_view) const
{
  PCMS_FUNCTION_TIMER;
  // TODO decide if we want to implicitly perform the coordinate transformations
  // when possible
  if (coordinate_view.GetCoordinateSystem() !=
      layout_.GetDOFHolderCoordinates().GetCoordinateSystem()) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Coordinate system mismatch");
  }

  auto coordinates = coordinate_view.GetCoordinates();
  Kokkos::View<Real* [2]> coords("coords", coordinates.size() / 2);
  auto coordinates_host = Kokkos::View<const Real**, HostMemorySpace>(
    coordinates.data_handle(), coordinates.extent(0), coordinates.extent(1));
  deep_copy_mismatch_layouts(coords, coordinates_host);
  auto results = search_(coords);
  Kokkos::View<GridPointSearch::Result*, HostMemorySpace> results_h(
    "results_h", results.size());
  Kokkos::deep_copy(results_h, results);
  auto hint = std::make_shared<OmegaHField2LocalizationHint>(
    mesh_, results_h, out_of_bounds_mode_);

  return LocalizationHint{hint};
}

void OmegaHField2::Evaluate(LocalizationHint location,
                            FieldDataView<Real, HostMemorySpace> results) const
{
  PCMS_FUNCTION_TIMER;
  // TODO decide if we want to implicitly perform the coordinate transformations
  // when possible
  if (results.GetCoordinateSystem() !=
      layout_.GetDOFHolderCoordinates().GetCoordinateSystem()) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Coordinate system mismatch");
  }

  OmegaHField2LocalizationHint hint =
    *reinterpret_cast<OmegaHField2LocalizationHint*>(location.data.get());

  Kokkos::View<Real**> coordinates_d(
    "coordinates_d", hint.coordinates_.extent(0), hint.coordinates_.extent(1));
  deep_copy_mismatch_layouts(coordinates_d, hint.coordinates_);
  Kokkos::View<LO*> offsets_d("offsets_d", hint.offsets_.extent(0));
  Kokkos::deep_copy(offsets_d, hint.offsets_);
  auto eval_results = mesh_field_->evaluate(coordinates_d, offsets_d);
  Kokkos::View<Real**, HostMemorySpace> eval_results_h(
    "eval_results_h", eval_results.extent(0), eval_results.extent(1));
  deep_copy_mismatch_layouts(eval_results_h, eval_results);
  Rank1View<Real, HostMemorySpace> values = results.GetValues();

  // Copy results for valid points
  Kokkos::parallel_for(
    "CopyEvalResultsToValues",
    Kokkos::RangePolicy<HostMemorySpace::execution_space>(
      0, eval_results_h.extent(0)),
    KOKKOS_LAMBDA(LO i) { values[hint.indices_(i)] = eval_results_h(i, 0); });

  // Handle missing points based on mode
  if (hint.num_missing_ > 0 && hint.mode_ == OutOfBoundsMode::FILL) {
    auto fill_val = fill_value_;
    Kokkos::parallel_for(
      "FillMissingValues",
      Kokkos::RangePolicy<HostMemorySpace::execution_space>(0,
                                                            hint.num_missing_),
      KOKKOS_LAMBDA(LO i) { values[hint.missing_indices_(i)] = fill_val; });
  }
}

void OmegaHField2::EvaluateGradient(
  FieldDataView<Real, HostMemorySpace> /* unused */)
{
  // TODO when moved to PCMS throw PCMS exception
  throw std::runtime_error("Not implemented");
}

const FieldLayout& OmegaHField2::GetLayout() const
{
  return layout_;
}

bool OmegaHField2::CanEvaluateGradient()
{
  // TODO compute the gradient field using element shape functions
  return false;
}

int OmegaHField2::Serialize(
  Rank1View<Real, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const
{
  PCMS_FUNCTION_TIMER;
  // host copy of filtered field data array
  const auto array_h = GetDOFHolderData();
  if (buffer.size() > 0) {
    auto owned = layout_.GetOwned();
    for (size_t i = 0; i < array_h.size(); i++) {
      if (owned[i])
        buffer[permutation[i]] = array_h[i];
    }
  }
  return array_h.size();
}

void OmegaHField2::Deserialize(
  Rank1View<const Real, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
{
  PCMS_FUNCTION_TIMER;
  Omega_h::HostWrite<Real> sorted_buffer(permutation.size());
  auto owned = layout_.GetOwned();
  for (LO i = 0; i < sorted_buffer.size(); ++i) {
    if (owned[i])
      sorted_buffer[i] = buffer[permutation[i]];
  }

  SetDOFHolderData(pcms::make_const_array_view(sorted_buffer));
}

} // namespace pcms
