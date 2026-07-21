#ifndef PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_BACKEND_H
#define PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_BACKEND_H

// Backend types shared between MeshFieldsAdapter2 and
// MeshFieldsEvaluatorFactory. Extracted here to avoid a circular include.

#include <Kokkos_Core.hpp>
#include <MeshField.hpp>
#include <memory>

#include "pcms/field/layout/mesh_fields.h"
#include "pcms/utility/types.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/arrays.h"
#include "pcms/localization/point_search.h"
#include "pcms/field/field.h" // OutOfBoundsMode

namespace pcms
{

// ---------------------------------------------------------------------------
// Abstract backend interface
// ---------------------------------------------------------------------------
template <typename T>
class MeshFieldBackend
{
public:
  virtual ~MeshFieldBackend() = default;
  virtual Kokkos::View<T* [1]> evaluate(Kokkos::View<T**> localCoords,
                                        Kokkos::View<LO*> offsets) const = 0;
  virtual void SetData(Rank1View<const T, DeviceMemorySpace> data,
                       size_t num_nodes, size_t num_components, int dim) = 0;
  virtual void GetData(Rank1View<T, DeviceMemorySpace> data, size_t num_nodes,
                       size_t num_components, int dim) const = 0;
};

// ---------------------------------------------------------------------------
// Concrete backend implementation
// ---------------------------------------------------------------------------
template <typename T, int Dim, int Order>
class MeshFieldBackendImpl : public MeshFieldBackend<T>
{
public:
  MeshFieldBackendImpl(Omega_h::Mesh& mesh)
    : mesh_(mesh),
      mesh_field_(mesh),
      shape_field_(mesh_field_.template CreateLagrangeField<T, Order, 1>())
  {
  }

  Kokkos::View<T* [1]> evaluate(Kokkos::View<T**> localCoords,
                                Kokkos::View<LO*> offsets) const override
  {
    auto self = const_cast<MeshFieldBackendImpl<T, Dim, Order>*>(this);
    return self->mesh_field_.triangleLocalPointEval(localCoords, offsets,
                                                    shape_field_);
  }

  void SetData(Rank1View<const T, DeviceMemorySpace> data, size_t num_nodes,
               size_t num_components, int dim) override
  {
    size_t stride = num_nodes * num_components;
    auto topo = static_cast<MeshField::Mesh_Topology>(dim);
    Kokkos::parallel_for(
      mesh_.nents(dim), KOKKOS_CLASS_LAMBDA(size_t ent) {
        for (size_t n = 0; n < num_nodes; ++n) {
          for (size_t c = 0; c < num_components; ++c) {
            shape_field_(ent, n, c, topo) =
              data[ent * stride + n * num_components + c];
          }
        }
      });
  }

  void GetData(Rank1View<T, DeviceMemorySpace> data, size_t num_nodes,
               size_t num_components, int dim) const override
  {
    size_t stride = num_nodes * num_components;
    auto topo = static_cast<MeshField::Mesh_Topology>(dim);
    Kokkos::parallel_for(
      mesh_.nents(dim), KOKKOS_CLASS_LAMBDA(size_t ent) {
        for (size_t n = 0; n < num_nodes; ++n) {
          for (size_t c = 0; c < num_components; ++c) {
            data[ent * stride + n * num_components + c] =
              shape_field_(ent, n, c, topo);
          }
        }
      });
  }

private:
  Omega_h::Mesh& mesh_;
  MeshField::OmegahMeshField<DefaultExecutionSpace, Dim> mesh_field_;
  using ShapeField =
    decltype(mesh_field_.template CreateLagrangeField<T, Order, 1>());
  ShapeField shape_field_;
};

// ---------------------------------------------------------------------------
// Factory function: create a MeshFieldBackend from a layout
// ---------------------------------------------------------------------------
template <typename T>
std::shared_ptr<MeshFieldBackend<T>> MakeMeshFieldBackend(
  const MeshFieldsAdapterLayout& layout)
{
  if constexpr (!std::is_same_v<T, MeshField::Real4> &&
                !std::is_same_v<T, MeshField::Real8>) {
    throw pcms_error(
      "MeshFieldBackend only supports the MeshFields scalar types enabled in "
      "this build");
  }

  Omega_h::Mesh& mesh = layout.GetMesh();
  if (mesh.dim() == 3) {
    throw pcms_error("MeshFieldBackend does not support 3D meshes");
  }
  auto nodes_per_dim = layout.GetNodesPerDim();
  if (nodes_per_dim[0] == 1 && nodes_per_dim[1] == 0 && nodes_per_dim[2] == 0 &&
      nodes_per_dim[3] == 0) {
    switch (mesh.dim()) {
      case 1: return std::make_shared<MeshFieldBackendImpl<T, 1, 1>>(mesh);
      case 2: return std::make_shared<MeshFieldBackendImpl<T, 2, 1>>(mesh);
      default: break;
    }
  } else if (nodes_per_dim[0] == 1 && nodes_per_dim[1] == 1 &&
             nodes_per_dim[2] == 0 && nodes_per_dim[3] == 0) {
    switch (mesh.dim()) {
      case 2: return std::make_shared<MeshFieldBackendImpl<T, 2, 2>>(mesh);
      case 3: return std::make_shared<MeshFieldBackendImpl<T, 3, 2>>(mesh);
      default: break;
    }
  }
  return nullptr;
}

// ---------------------------------------------------------------------------
// Helper functors used in MeshFieldsAdapter2LocalizationHint
// ---------------------------------------------------------------------------
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

// Device-side version
struct ComputeOffsetsDeviceFunctor
{
  Kokkos::View<LO*> offsets_;
  Kokkos::View<LO*> elem_counts_;

  ComputeOffsetsDeviceFunctor(Kokkos::View<LO*> offsets,
                              Kokkos::View<LO*> elem_counts)
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
  Kokkos::View<GridPointSearch2D::Result*> search_results_;
  Kokkos::View<LO*> owning_elem_ids_;

  CountPointsPerElementFunctor(
    Kokkos::View<LO*> elem_counts,
    Kokkos::View<GridPointSearch2D::Result*> search_results,
    Kokkos::View<LO*> owning_elem_ids)
    : elem_counts_(elem_counts),
      search_results_(search_results),
      owning_elem_ids_(owning_elem_ids)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    const auto owner_idx = owning_elem_ids_(i);
    Kokkos::atomic_add(&elem_counts_(owner_idx), 1);
  }
};

struct FillCoordinatesAndIndicesFunctor
{
  Omega_h::Mesh& mesh_;
  Kokkos::View<LO*> elem_counts_;
  Kokkos::View<LO*> offsets_;
  Kokkos::View<Real**> coordinates_;
  Kokkos::View<LO*> indices_;
  Kokkos::View<GridPointSearch2D::Result*> search_results_;
  Kokkos::View<LO*> owning_elem_ids_;
  Omega_h::Int dim_;

  FillCoordinatesAndIndicesFunctor(
    Omega_h::Mesh& mesh, Kokkos::View<LO*> elem_counts,
    Kokkos::View<LO*> offsets, Kokkos::View<Real**> coordinates,
    Kokkos::View<LO*> indices,
    Kokkos::View<GridPointSearch2D::Result*> search_results,
    Kokkos::View<LO*> owning_elem_ids)
    : mesh_(mesh),
      elem_counts_(elem_counts),
      offsets_(offsets),
      coordinates_(coordinates),
      indices_(indices),
      search_results_(search_results),
      owning_elem_ids_(owning_elem_ids),
      dim_(mesh.dim())
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    auto [dim, elem_idx, coord] = search_results_(i);
    (void)dim;
    (void)elem_idx;
    const auto owner_idx = owning_elem_ids_(i);
    LO count = Kokkos::atomic_sub_fetch(&elem_counts_(owner_idx), 1);
    LO index = offsets_(owner_idx) + count - 1;
    for (int j = 0; j < (dim_ + 1); ++j) {
      coordinates_(index, j) = coord[j];
    }
    indices_(index) = i;
  }
};

// Device-side functors for processing search results
struct FilterValidPointsFunctor
{
  Kokkos::View<GridPointSearch2D::Result*> search_results_;
  Kokkos::View<LO*> valid_flags_;
  Kokkos::View<LO*> owning_elem_ids_;
  Omega_h::Int mesh_dim_;
  OutOfBoundsMode mode_;

  FilterValidPointsFunctor(Kokkos::View<GridPointSearch2D::Result*> results,
                           Kokkos::View<LO*> flags,
                           Kokkos::View<LO*> owning_elem_ids,
                           Omega_h::Int mesh_dim, OutOfBoundsMode mode)
    : search_results_(results),
      valid_flags_(flags),
      owning_elem_ids_(owning_elem_ids),
      mesh_dim_(mesh_dim),
      mode_(mode)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i, LO& num_valid, LO& num_missing) const
  {
    auto [dim, elem_idx, coord] = search_results_(i);
    (void)dim;
    (void)coord;
    const auto owner_idx = owning_elem_ids_(i);
    bool is_missing = (elem_idx < 0) || (owner_idx < 0);

    if (mode_ == OutOfBoundsMode::ERROR && is_missing) {
      Kokkos::abort("Points found outside mesh domain");
    }

    valid_flags_(i) = !is_missing ? 1 : 0;
    if (!is_missing) {
      num_valid++;
    } else {
      num_missing++;
    }
  }
};

struct CompactIndicesFunctor
{
  Kokkos::View<LO*> flags_;
  Kokkos::View<LO*> scan_;
  Kokkos::View<LO*> valid_indices_;
  Kokkos::View<LO*> missing_indices_;
  LO num_valid_;

  CompactIndicesFunctor(Kokkos::View<LO*> flags, Kokkos::View<LO*> scan,
                        Kokkos::View<LO*> valid_indices,
                        Kokkos::View<LO*> missing_indices, LO num_valid)
    : flags_(flags),
      scan_(scan),
      valid_indices_(valid_indices),
      missing_indices_(missing_indices),
      num_valid_(num_valid)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    if (flags_(i) == 1) {
      // Valid point
      valid_indices_(scan_(i)) = i;
    } else {
      // Missing point - use offset from end of valid region
      LO missing_offset = i - scan_(i);
      missing_indices_(missing_offset) = i;
    }
  }
};

struct CountPerElementFunctor
{
  Kokkos::View<GridPointSearch2D::Result*> search_results_;
  Kokkos::View<LO*> valid_indices_;
  Kokkos::View<LO*> owning_elem_ids_;
  Kokkos::View<LO*> elem_counts_;

  CountPerElementFunctor(Kokkos::View<GridPointSearch2D::Result*> results,
                         Kokkos::View<LO*> valid_indices,
                         Kokkos::View<LO*> owning_elem_ids,
                         Kokkos::View<LO*> elem_counts)
    : search_results_(results),
      valid_indices_(valid_indices),
      owning_elem_ids_(owning_elem_ids),
      elem_counts_(elem_counts)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    LO orig_idx = valid_indices_(i);
    const auto owner_idx = owning_elem_ids_(orig_idx);
    Kokkos::atomic_add(&elem_counts_(owner_idx), 1);
  }
};

struct FillCoordinatesDeviceFunctor
{
  Kokkos::View<GridPointSearch2D::Result*> search_results_;
  Kokkos::View<LO*> valid_indices_;
  Kokkos::View<LO*> owning_elem_ids_;
  Kokkos::View<LO*> elem_counts_;
  Kokkos::View<LO*> offsets_;
  Kokkos::View<Real**> coordinates_;
  Kokkos::View<LO*> indices_;
  Omega_h::LOs tris2verts_;
  Omega_h::Reals mesh_coords_;
  Kokkos::View<const Real**, DeviceMemorySpace> global_coords_;
  Omega_h::Mesh& mesh_;
  Omega_h::Int dim_;

  FillCoordinatesDeviceFunctor(
    Kokkos::View<GridPointSearch2D::Result*> results,
    Kokkos::View<LO*> valid_indices, Kokkos::View<LO*> owning_elem_ids,
    Kokkos::View<LO*> elem_counts, Kokkos::View<LO*> offsets,
    Kokkos::View<Real**> coordinates, Kokkos::View<LO*> indices,
    Omega_h::LOs tris2verts, Omega_h::Reals mesh_coords,
    Kokkos::View<const Real**, DeviceMemorySpace> global_coords,
    Omega_h::Mesh& mesh, Omega_h::Int dim)
    : search_results_(results),
      valid_indices_(valid_indices),
      owning_elem_ids_(owning_elem_ids),
      elem_counts_(elem_counts),
      offsets_(offsets),
      coordinates_(coordinates),
      indices_(indices),
      tris2verts_(tris2verts),
      mesh_coords_(mesh_coords),
      global_coords_(global_coords),
      mesh_(mesh),
      dim_(dim)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO i) const
  {
    LO orig_idx = valid_indices_(i);
    const auto owner_idx = owning_elem_ids_(orig_idx);
    LO count =
      Kokkos::atomic_fetch_add(&elem_counts_(owner_idx), static_cast<LO>(-1));
    LO index = offsets_(owner_idx) + count - 1;
    const auto elem_tri2verts =
      Omega_h::gather_verts<3>(tris2verts_, owner_idx);
    const auto vertex_coords =
      Omega_h::gather_vectors<3, 2>(mesh_coords_, elem_tri2verts);
    const Omega_h::Vector<2> point{global_coords_(orig_idx, 0),
                                   global_coords_(orig_idx, 1)};
    const auto local =
      Omega_h::barycentric_from_global<2, 2>(point, vertex_coords);
    for (int j = 0; j < (dim_ + 1); ++j)
      coordinates_(index, j) = local[j];
    indices_(index) = orig_idx;
  }
};

struct ExclusiveScanFunctor
{
  Kokkos::View<LO*> valid_flags_;
  Kokkos::View<LO*> scan_result_;

  ExclusiveScanFunctor(Kokkos::View<LO*> flags, Kokkos::View<LO*> scan)
    : valid_flags_(flags), scan_result_(scan)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const LO i, LO& update, const bool final) const
  {
    if (final) {
      scan_result_(i) = update;
    }
    update += valid_flags_(i);
  }
};

struct SetFinalOffsetFunctor
{
  Kokkos::View<LO*> offsets_;
  LO nelems_;
  LO total_;

  SetFinalOffsetFunctor(Kokkos::View<LO*> offsets, LO nelems, LO total)
    : offsets_(offsets), nelems_(nelems), total_(total)
  {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(LO) const { offsets_(nelems_) = total_; }
};

// ---------------------------------------------------------------------------
// Localization hint for MeshFields
// ---------------------------------------------------------------------------
struct MeshFieldsAdapter2LocalizationHint
{
  // Host-side constructor (legacy, for compatibility)
  MeshFieldsAdapter2LocalizationHint(
    Omega_h::Mesh& mesh,
    Kokkos::View<GridPointSearch2D::Result*, HostMemorySpace> search_results,
    Kokkos::View<const Real**, HostMemorySpace> global_coords,
    Kokkos::View<LO*, HostMemorySpace> owning_elem_ids, OutOfBoundsMode mode)
    : mode_(mode), num_valid_(0), num_missing_(0)
  {
    std::vector<size_t> valid_point_indices;
    std::vector<size_t> missing_point_indices;

    if (mode_ == OutOfBoundsMode::ERROR) {
      for (size_t i = 0; i < search_results.size(); ++i) {
        auto [dim, elem_idx, coord] = search_results(i);
        const auto owner_idx = owning_elem_ids(i);
        (void)dim;
        (void)coord;
        bool is_missing = (elem_idx < 0) || (owner_idx < 0);
        PCMS_ALWAYS_ASSERT(!is_missing && "Points found outside mesh domain");
        valid_point_indices.push_back(i);
      }
    } else {
      for (size_t i = 0; i < search_results.size(); ++i) {
        auto [dim, elem_idx, coord] = search_results(i);
        const auto owner_idx = owning_elem_ids(i);
        (void)dim;
        (void)coord;
        bool is_missing = (elem_idx < 0) || (owner_idx < 0);
        if (is_missing) {
          missing_point_indices.push_back(i);
        } else {
          valid_point_indices.push_back(i);
        }
      }
    }

    num_valid_ = valid_point_indices.size();
    num_missing_ = missing_point_indices.size();

    if (num_missing_ > 0 && mode_ == OutOfBoundsMode::NEAREST_BOUNDARY) {
      PCMS_ALWAYS_ASSERT(false && "NEAREST_BOUNDARY mode not implemented yet");
    }

    offsets_ = Kokkos::View<LO*, HostMemorySpace>("offsets", mesh.nelems() + 1);
    coordinates_ = Kokkos::View<Real**, HostMemorySpace>(
      "coordinates", num_valid_, mesh.dim() + 1);
    indices_ = Kokkos::View<LO*, HostMemorySpace>("indices", num_valid_);

    if (num_missing_ > 0) {
      missing_indices_ =
        Kokkos::View<LO*, HostMemorySpace>("missing_indices", num_missing_);
      for (size_t i = 0; i < num_missing_; ++i) {
        missing_indices_(i) = static_cast<LO>(missing_point_indices[i]);
      }
    }

    Kokkos::View<LO*, HostMemorySpace> elem_counts("elem_counts",
                                                   mesh.nelems());
    Kokkos::deep_copy(elem_counts, 0);
    for (size_t i = 0; i < num_valid_; ++i) {
      const auto owner_idx = owning_elem_ids(valid_point_indices[i]);
      elem_counts[owner_idx] += 1;
    }

    LO total;
    ComputeOffsetsFunctor functor(offsets_, elem_counts);
    Kokkos::parallel_scan(
      "ComputeOffsets",
      Kokkos::RangePolicy<HostMemorySpace::execution_space>(0, mesh.nelems()),
      functor, total);
    offsets_(mesh.nelems()) = total;

    const auto tris2verts = mesh.ask_elem_verts();
    const auto mesh_coords = mesh.coords();
    for (size_t i = 0; i < num_valid_; ++i) {
      size_t orig_idx = valid_point_indices[i];
      const auto owner_idx = owning_elem_ids(orig_idx);
      elem_counts(owner_idx) -= 1;
      LO index = offsets_(owner_idx) + elem_counts(owner_idx);
      const auto elem_tri2verts =
        Omega_h::gather_verts<3>(tris2verts, owner_idx);
      const auto vertex_coords =
        Omega_h::gather_vectors<3, 2>(mesh_coords, elem_tri2verts);
      const Omega_h::Vector<2> point{global_coords(orig_idx, 0),
                                     global_coords(orig_idx, 1)};
      const auto local =
        Omega_h::barycentric_from_global<2, 2>(point, vertex_coords);
      for (int j = 0; j < (mesh.dim() + 1); ++j)
        coordinates_(index, j) = local[j];
      indices_(index) = static_cast<LO>(orig_idx);
    }

    // Copy all data to device once during construction
    offsets_d_ = Kokkos::View<LO*>("offsets_d", offsets_.extent(0));
    Kokkos::deep_copy(offsets_d_, offsets_);

    coordinates_d_ = Kokkos::View<Real**>(
      "coordinates_d", coordinates_.extent(0), coordinates_.extent(1));
    DeepCopyMismatchLayouts(coordinates_d_, coordinates_);

    indices_d_ = Kokkos::View<LO*>("indices_d", indices_.extent(0));
    Kokkos::deep_copy(indices_d_, indices_);

    if (num_missing_ > 0) {
      missing_indices_d_ =
        Kokkos::View<LO*>("missing_indices_d", missing_indices_.extent(0));
      Kokkos::deep_copy(missing_indices_d_, missing_indices_);
    }
  }

  // Device-side constructor
  MeshFieldsAdapter2LocalizationHint(
    Omega_h::Mesh& mesh,
    Kokkos::View<GridPointSearch2D::Result*> search_results_d,
    Kokkos::View<const Real**, DeviceMemorySpace> global_coords,
    Kokkos::View<LO*, DeviceMemorySpace> owning_elem_ids, OutOfBoundsMode mode)
    : mode_(mode), num_valid_(0), num_missing_(0)
  {
    const LO n_points = search_results_d.size();

    if (mode_ == OutOfBoundsMode::NEAREST_BOUNDARY) {
      PCMS_ALWAYS_ASSERT(false && "NEAREST_BOUNDARY mode not implemented yet");
    }

    // Step 1: Filter valid/missing points on device
    Kokkos::View<LO*> valid_flags("valid_flags", n_points);
    FilterValidPointsFunctor filter_functor(search_results_d, valid_flags,
                                            owning_elem_ids, mesh.dim(), mode_);
    LO num_valid = 0, num_missing = 0;
    Kokkos::parallel_reduce(
      "FilterValidPoints",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_points),
      filter_functor, num_valid, num_missing);

    num_valid_ = num_valid;
    num_missing_ = num_missing;

    Kokkos::View<LO*> scan_result("scan_result", n_points);
    ExclusiveScanFunctor scan_functor(valid_flags, scan_result);
    Kokkos::parallel_scan(
      "ExclusiveScanFlags",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_points),
      scan_functor);

    Kokkos::View<LO*> valid_indices_d("valid_indices_d", num_valid_);
    Kokkos::View<LO*> missing_indices_tmp("missing_indices_tmp", num_missing_);

    CompactIndicesFunctor compact_functor(valid_flags, scan_result,
                                          valid_indices_d, missing_indices_tmp,
                                          num_valid_);
    Kokkos::parallel_for(
      "CompactIndices",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, n_points),
      compact_functor);

    Kokkos::View<LO*> elem_counts_d("elem_counts_d", mesh.nelems());
    Kokkos::deep_copy(elem_counts_d, 0);

    CountPerElementFunctor count_functor(search_results_d, valid_indices_d,
                                         owning_elem_ids, elem_counts_d);
    Kokkos::parallel_for(
      "CountPerElement",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_valid_),
      count_functor);

    offsets_d_ = Kokkos::View<LO*>("offsets_d", mesh.nelems() + 1);
    ComputeOffsetsDeviceFunctor offsets_functor(offsets_d_, elem_counts_d);
    LO total;
    Kokkos::parallel_scan(
      "ComputeOffsets",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, mesh.nelems()),
      offsets_functor, total);

    // Set final offset
    SetFinalOffsetFunctor set_final_functor(offsets_d_, mesh.nelems(), total);
    Kokkos::parallel_for(
      "SetFinalOffset",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, 1),
      set_final_functor);

    // Step 6: Fill coordinates and indices on device
    coordinates_d_ =
      Kokkos::View<Real**>("coordinates_d", num_valid_, mesh.dim() + 1);
    indices_d_ = Kokkos::View<LO*>("indices_d", num_valid_);
    const auto tris2verts = mesh.ask_elem_verts();
    const auto mesh_coords = mesh.coords();

    FillCoordinatesDeviceFunctor fill_functor(
      search_results_d, valid_indices_d, owning_elem_ids, elem_counts_d,
      offsets_d_, coordinates_d_, indices_d_, tris2verts, mesh_coords,
      global_coords, mesh, mesh.dim());
    Kokkos::parallel_for(
      "FillCoordinatesAndIndices",
      Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, num_valid_),
      fill_functor);

    // Step 7: Handle missing indices
    if (num_missing_ > 0) {
      missing_indices_d_ = missing_indices_tmp;
    }

    // Create host mirrors for compatibility (lazy copy - only if needed)
    offsets_ = Kokkos::create_mirror_view(offsets_d_);
    coordinates_ = Kokkos::View<Real**, HostMemorySpace>(
      "coordinates_", num_valid_, mesh.dim() + 1);
    DeepCopyMismatchLayouts(coordinates_, coordinates_d_);
    indices_ = Kokkos::create_mirror_view(indices_d_);
    if (num_missing_ > 0) {
      missing_indices_ = Kokkos::create_mirror_view(missing_indices_d_);
    }
  }

  OutOfBoundsMode mode_;
  size_t num_valid_;
  size_t num_missing_;

  // Host views (for construction and debugging)
  Kokkos::View<LO*, HostMemorySpace> offsets_;
  Kokkos::View<Real**, HostMemorySpace> coordinates_;
  Kokkos::View<LO*, HostMemorySpace> indices_;
  Kokkos::View<LO*, HostMemorySpace> missing_indices_;

  // Device views (primary data for evaluation - avoid copies)
  Kokkos::View<LO*> offsets_d_;
  Kokkos::View<Real**> coordinates_d_;
  Kokkos::View<LO*> indices_d_;
  Kokkos::View<LO*> missing_indices_d_;
};

} // namespace pcms

#endif // PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_BACKEND_H
