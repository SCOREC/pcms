#include "pcms/field/layout/mesh_fields.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/inclusive_scan.h"
#include "pcms/utility/profile.h"
#include <Omega_h_for.hpp>
#include <MeshField.hpp>
#include <memory>

namespace pcms
{
/*
 * Field Layout
 */

template <typename T>
Omega_h::Write<Omega_h::GO> GetGidsHelper(LO total_ents,
                                          std::array<int, 4> nodes_per_dim,
                                          Omega_h::Mesh& mesh,
                                          const std::string& global_id_name)
{
  PCMS_FUNCTION_TIMER;

  Omega_h::Write<Omega_h::GO> owned_gids(total_ents);
  LO offset = 0;
  for (int i = 0; i <= mesh.dim(); ++i) {
    if (nodes_per_dim[i]) {
      auto dim_gids = mesh.get_array<T>(i, global_id_name);
      Omega_h::parallel_for(
        dim_gids.size(),
        OMEGA_H_LAMBDA(int i) { owned_gids[i + offset] = dim_gids[i]; });
      offset += dim_gids.size();
    }
  }

  PCMS_ALWAYS_ASSERT(offset == total_ents);

  return owned_gids;
}

namespace
{
// Maps a meshFields topology to the dimension.
KOKKOS_INLINE_FUNCTION int TopologyToDim(MeshField::Mesh_Topology topo)
{
  switch (topo) {
    case MeshField::Vertex: return 0;
    case MeshField::Edge: return 1;
    case MeshField::Triangle: return 2;
    case MeshField::Tetrahedron: return 3;
    default: return -1;
  }
}

// Gets the appropriate MeshField element
template <int Dim, int Order>
auto GetMeshFieldElement(Omega_h::Mesh& mesh)
{
  if constexpr (Dim == 2) {
    return MeshField::Omegah::getTriangleElement<Order>(mesh);
  } else {
    return MeshField::Omegah::getTetrahedronElement<Order>(mesh);
  }
}

template <int Dim, int Order>
void BuildDofHolderCoordsFromMeshFieldImpl(
  Omega_h::Mesh& mesh, Kokkos::View<Real**> holder_coords,
  const std::array<int, 4>& nodes_per_dim)
{
  const auto elem = GetMeshFieldElement<Dim, Order>(mesh);

  using ShapeT = std::decay_t<decltype(elem.shp)>;
  constexpr size_t numNodes = ShapeT::numNodes;
  constexpr size_t meshDim = ShapeT::meshEntDim;
  constexpr MeshField::Mesh_Topology elemTopo =
    (Dim == 2) ? MeshField::Triangle : MeshField::Tetrahedron;
  const auto map = elem.map;

  // Compute the starting offset of DOF holders for each entity dimension.
  Kokkos::Array<LO, 4> dof_holder_offset;
  {
    LO cur_offset = 0;
    for (int d = 0; d < 4; ++d) {
      dof_holder_offset[d] = cur_offset;
      if (d <= mesh.dim() && nodes_per_dim[d] != 0)
        cur_offset += mesh.nents(d);
    }
  }

  const LO nelems = mesh.nelems();
  const auto vtx_coords = mesh.coords();
  const auto param_coords = elem.shp.getNodeParametricCoords();

  Kokkos::parallel_for(
    "MeshFieldsDofHolderCoords",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, nelems),
    KOKKOS_LAMBDA(const LO e) {
      constexpr size_t NVerts = Dim + 1;

      Real Xv[NVerts * Dim];
      for (size_t b = 0; b < NVerts; ++b) {
        const auto hv = map(static_cast<MeshField::LO>(b), 0, e, elemTopo);
        const LO v = hv.entity;
        for (size_t d = 0; d < Dim; ++d)
          Xv[b * Dim + d] = vtx_coords[static_cast<size_t>(v) * Dim + d];
      }

      for (size_t n = 0; n < numNodes; ++n) {
        const auto h = map(static_cast<MeshField::LO>(n), 0, e, elemTopo);

        // Compute the barycentric coordinates L[b] at the node's parametric
        // coordinates.
        Real L[NVerts];
        L[0] = 1;
        for (size_t d = 0; d < Dim; ++d) {
          const Real xi_d = param_coords[n * Dim + d];
          L[d + 1] = xi_d;
          L[0] -= xi_d;
        }

        Real X[Dim] = {0};
        for (size_t b = 0; b < NVerts; ++b)
          for (size_t d = 0; d < Dim; ++d)
            X[d] += L[b] * Xv[b * Dim + d];

        const int dim_of_holder = TopologyToDim(h.topo);
        const LO row = dof_holder_offset[dim_of_holder] + h.entity;
        for (size_t d = 0; d < Dim; ++d)
          holder_coords(row, d) = X[d];
      }
    });
}

void BuildDofHolderCoordsFromMeshField(Omega_h::Mesh& mesh,
                                       Kokkos::View<Real**> holder_coords,
                                       const std::array<int, 4>& nodes_per_dim)
{
  int dim = mesh.dim();
  int order = 0;
  for (int i = 0; i <= dim; ++i) {
    if (nodes_per_dim[i] == 1)
      ++order;
    else if (nodes_per_dim[i] != 0) {
      std::cerr << "Unsupported" << std::endl;
      std::abort();
    }
  }

  if (dim == 2 && order == 1)
    BuildDofHolderCoordsFromMeshFieldImpl<2, 1>(mesh, holder_coords,
                                                nodes_per_dim);
  else if (dim == 2 && order == 2)
    BuildDofHolderCoordsFromMeshFieldImpl<2, 2>(mesh, holder_coords,
                                                nodes_per_dim);
  else if (dim == 3 && order == 1)
    BuildDofHolderCoordsFromMeshFieldImpl<3, 1>(mesh, holder_coords,
                                                nodes_per_dim);
  else if (dim == 3 && order == 2)
    BuildDofHolderCoordsFromMeshFieldImpl<3, 2>(mesh, holder_coords,
                                                nodes_per_dim);
  else {
    std::cerr << "Unsupported element/order combination" << std::endl;
    std::abort();
  }
}
} // namespace

struct CopyClassInfoFunctor
{
  Omega_h::Write<Omega_h::ClassId> class_ids_;
  Omega_h::Write<Omega_h::I8> class_dims_;
  Kokkos::View<bool*> owned_;
  Omega_h::Read<Omega_h::ClassId> ids_;
  Omega_h::Read<Omega_h::I8> dims_;
  Omega_h::Read<Omega_h::I8> owned_data_;
  size_t offset_;

  CopyClassInfoFunctor(Omega_h::Write<Omega_h::ClassId> class_ids,
                       Omega_h::Write<Omega_h::I8> class_dims,
                       Kokkos::View<bool*> owned,
                       Omega_h::Read<Omega_h::ClassId> ids,
                       Omega_h::Read<Omega_h::I8> dims,
                       Omega_h::Read<Omega_h::I8> owned_data, size_t offset)
    : class_ids_(class_ids),
      class_dims_(class_dims),
      owned_(owned),
      ids_(ids),
      dims_(dims),
      owned_data_(owned_data),
      offset_(offset)
  {
  }

  OMEGA_H_DEVICE
  void operator()(LO i) const
  {
    class_ids_[offset_ + i] = ids_[i];
    class_dims_[offset_ + i] = dims_[i];
    owned_[offset_ + i] = owned_data_[i];
  }
};

MeshFieldsAdapterLayout::MeshFieldsAdapterLayout(
  Omega_h::Mesh& mesh, std::array<int, 4> nodes_per_dim, int num_components,
  CoordinateSystem coordinate_system, std::string global_id_name)
  : mesh_(mesh),
    global_id_name_(global_id_name),
    num_components_(num_components),
    coordinate_system_(coordinate_system),
    nodes_per_dim_(nodes_per_dim),
    dof_holder_coords_("", GetNumOwnedDofHolder(), mesh_.dim()),
    class_ids_(GetNumEnts()),
    class_dims_(class_ids_.size()),
    owned_("", class_dims_.size()),
    owned_host_("", class_dims_.size())
{
  PCMS_FUNCTION_TIMER;
  LO total_ents = GetNumEnts();

  auto tag = mesh_.get_tagbase(0, global_id_name_);
  if (Omega_h::is<GO>(tag)) {
    gids_ = GetGidsHelper<GO>(total_ents, nodes_per_dim, mesh, global_id_name);
  } else if (Omega_h::is<LO>(tag)) {
    gids_ = GetGidsHelper<LO>(total_ents, nodes_per_dim, mesh, global_id_name);
  } else {
    std::cerr << "Weird tag type for global arrays.\n";
    std::abort();
  }

  BuildDofHolderCoordsFromMeshField(mesh_, dof_holder_coords_, nodes_per_dim);

  size_t offset = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    if (nodes_per_dim_[i]) {
      auto ids = mesh_.get_array<Omega_h::ClassId>(i, "class_id");
      auto dims = mesh_.get_array<Omega_h::I8>(i, "class_dim");
      auto owned = mesh_.owned(i);
      PCMS_ALWAYS_ASSERT(ids.size() == dims.size() &&
                         dims.size() == mesh_.nents(i));

      CopyClassInfoFunctor functor(class_ids_, class_dims_, owned_, ids, dims,
                                   owned, offset);
      Omega_h::parallel_for(mesh_.nents(i), functor);
      offset += mesh.nents(i);
    }
  }
  gids_host_ = Omega_h::HostWrite<Omega_h::GO>(gids_);

  int n = class_ids_.size();
  classification_dims_host_ =
    Kokkos::View<LO*, HostMemorySpace>("classification_dims", n);
  classification_ids_host_ =
    Kokkos::View<LO*, HostMemorySpace>("classification_ids", n);
  auto class_dims_h = Omega_h::HostRead<Omega_h::I8>(class_dims_);
  auto class_ids_h = Omega_h::HostRead<Omega_h::ClassId>(class_ids_);
  for (int i = 0; i < n; ++i) {
    classification_dims_host_(i) =
      static_cast<LO>(static_cast<unsigned char>(class_dims_h[i]));
    classification_ids_host_(i) = static_cast<LO>(class_ids_h[i]);
  }
  discretization_ = std::make_shared<OmegaHDiscretization>(mesh_);
}

std::shared_ptr<const Discretization>
MeshFieldsAdapterLayout::GetDiscretization() const noexcept
{
  return discretization_;
}

int MeshFieldsAdapterLayout::GetNumComponents() const
{
  return num_components_;
}

LO MeshFieldsAdapterLayout::GetNumOwnedDofHolder() const
{
  LO count = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    count += mesh_.nents(i) * nodes_per_dim_[i];
  }
  return count;
}

GO MeshFieldsAdapterLayout::GetNumGlobalDofHolder() const
{
  LO count = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    count += mesh_.nglobal_ents(i) * nodes_per_dim_[i];
  }
  return count;
}

std::array<int, 4> MeshFieldsAdapterLayout::GetNodesPerDim() const
{
  return nodes_per_dim_;
}

Rank1View<const bool, HostMemorySpace> MeshFieldsAdapterLayout::GetOwnedHost()
  const
{
  Kokkos::deep_copy(owned_host_, owned_);
  return make_const_array_view(owned_host_);
}

GlobalIDView<HostMemorySpace> MeshFieldsAdapterLayout::GetGidsHost() const
{
  return GlobalIDView<HostMemorySpace>(gids_host_.data(), gids_host_.size());
}

CoordinateView<DeviceMemorySpace>
MeshFieldsAdapterLayout::GetDOFHolderCoordinates() const
{
  auto coords_view = MakeConstRank2View(dof_holder_coords_);
  return CoordinateView<DeviceMemorySpace>{coordinate_system_, coords_view};
}

bool MeshFieldsAdapterLayout::IsDistributed() const
{
  return true;
}

Omega_h::Read<Omega_h::ClassId> MeshFieldsAdapterLayout::GetClassIDs() const
{
  PCMS_FUNCTION_TIMER;
  return Omega_h::Read(class_ids_);
}

Omega_h::Read<Omega_h::I8> MeshFieldsAdapterLayout::GetClassDims() const
{
  PCMS_FUNCTION_TIMER;
  return Omega_h::Read(class_dims_);
}

size_t MeshFieldsAdapterLayout::GetNumEnts() const
{
  size_t n = 0;
  for (int i = 0; i <= mesh_.dim(); ++i) {
    if (nodes_per_dim_[i])
      n += mesh_.nents(i);
  }
  return n;
}

Omega_h::Mesh& MeshFieldsAdapterLayout::GetMesh() const
{
  return mesh_;
}

EntOffsetsArray MeshFieldsAdapterLayout::GetEntOffsets() const
{
  EntOffsetsArray offsets{};
  size_t offset = 0;
  for (size_t i = 0; i < offsets.size(); ++i) {
    offsets[i] = offset;
    if (i <= static_cast<size_t>(mesh_.dim()) && nodes_per_dim_[i])
      offset += mesh_.nents(i);
  }
  offsets[offsets.size() - 1] = offset;
  return offsets;
}

int MeshFieldsAdapterLayout::GetDimension() const
{
  return mesh_.dim();
}

Rank1View<const LO, HostMemorySpace>
MeshFieldsAdapterLayout::GetDOFHolderClassificationDimensionsHost() const
{
  return make_const_array_view(classification_dims_host_);
}

Rank1View<const LO, HostMemorySpace>
MeshFieldsAdapterLayout::GetDOFHolderClassificationIdsHost() const
{
  return make_const_array_view(classification_ids_host_);
}

} // namespace pcms
