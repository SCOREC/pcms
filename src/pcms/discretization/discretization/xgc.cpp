#include "pcms/discretization/discretization/xgc.hpp"

namespace pcms
{

struct ClassifyVertsFunctor
{
  pcms::DimID geom;
  Kokkos::View<LO*, DeviceMemorySpace> verts;
  Kokkos::View<ClassificationDimension*, DeviceMemorySpace> class_dims_;
  Kokkos::View<ClassificationId*, DeviceMemorySpace> class_ids_;

  ClassifyVertsFunctor(
    DimID geom, Kokkos::View<LO*, DeviceMemorySpace> verts,
    Kokkos::View<ClassificationDimension*, DeviceMemorySpace> class_dims,
    Kokkos::View<ClassificationId*, DeviceMemorySpace> class_ids)
    : geom(geom), verts(verts), class_dims_(class_dims), class_ids_(class_ids)
  {
  }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    LO vert = verts(i);
    if (vert >= 0 && vert < class_dims_.extent(0)) {
      class_dims_(vert) = geom.dim;
      class_ids_(vert) = geom.id;
    }
  };
};

XGCDiscretization::XGCDiscretization(
  const ReverseClassificationVertex& reverse_classification, LO num_plane_nodes)
  : reverse_classification_(&reverse_classification),
    num_plane_nodes_(num_plane_nodes),
    class_dims_("xgc_discretization_class_dims", num_plane_nodes),
    class_ids_("xgc_discretization_class_ids", num_plane_nodes)
{
  Kokkos::deep_copy(class_dims_, static_cast<ClassificationDimension>(-1));
  Kokkos::deep_copy(class_ids_, static_cast<ClassificationId>(-1));

  for (const auto& [geom, verts] : reverse_classification) {
    auto verts_host =
      Kokkos::View<LO*, HostMemorySpace>("verts_host", verts.size());
    int idx = 0;
    for (LO vert : verts)
      verts_host(idx++) = vert;
    auto verts_device =
      Kokkos::View<LO*, DeviceMemorySpace>("verts_device", verts.size());
    Kokkos::deep_copy(verts_device, verts_host);
    Kokkos::parallel_for(
      "ClassifyVerts", Kokkos::RangePolicy<>(0, verts.size()),
      ClassifyVertsFunctor(geom, verts_device, class_dims_, class_ids_));
    Kokkos::fence(); // Wait for kernel to complete before verts_device is
                     // destroyed, better would be to optimize the
                     // reverse_classification data structure to avoid this copy
                     // and synchronization, but this is simpler for now
  }
}

bool XGCDiscretization::SameEntities(const Discretization& other) const noexcept
{
  auto p = dynamic_cast<const XGCDiscretization*>(&other);
  return p != nullptr &&
         p->reverse_classification_ == reverse_classification_ &&
         p->num_plane_nodes_ == num_plane_nodes_;
}

int XGCDiscretization::GetDimension() const
{
  return 2;
}

LO XGCDiscretization::GetNumEntities(int entity_dim) const
{
  return entity_dim == Vertex ? num_plane_nodes_ : 0;
}

Rank1View<const ClassificationDimension, DeviceMemorySpace>
XGCDiscretization::GetEntityClassificationDimensions(int entity_dim) const
{
  return entity_dim == Vertex
           ? make_const_array_view(class_dims_)
           : Rank1View<const ClassificationDimension, DeviceMemorySpace>(
               nullptr, 0);
}

Rank1View<const ClassificationId, DeviceMemorySpace>
XGCDiscretization::GetEntityClassificationIds(int entity_dim) const
{
  return entity_dim == Vertex
           ? make_const_array_view(class_ids_)
           : Rank1View<const ClassificationId, DeviceMemorySpace>(nullptr, 0);
}

} // namespace pcms
