#include "pcms/field/layout/mfem.h"

#include "pcms/discretization/discretization/point_cloud.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/profile.h"

#include <memory>

namespace pcms
{

namespace
{

std::shared_ptr<const Discretization> MakeVertexDiscretization(
  int dim, Kokkos::View<const Real**, HostMemorySpace> coords_host)
{
  return std::make_shared<PointCloudDiscretization>(dim, coords_host);
}

} // namespace

MFEMLayout::MFEMLayout(mfem::ParMesh& pmesh,
                       mfem::ParFiniteElementSpace& pfes,
                       CoordinateSystem coordinate_system)
  : pmesh_(pmesh),
    pfes_(pfes),
    dim_(pmesh.SpaceDimension()),
    coordinate_system_(coordinate_system)
{
  PCMS_FUNCTION_TIMER;

  AssertVertexScalarSpace();

  const int nv = pmesh_.GetNV();

  // Host coordinates, one row per vertex. GetVertex returns a pointer to the
  // contiguous spatial coordinates of the vertex, independent of MFEM's
  // internal storage ordering.
  Kokkos::View<Real**, HostMemorySpace> coords_host("mfem_coords_host", nv,
                                                    dim_);
  for (int v = 0; v < nv; ++v) {
    const double* x = pmesh_.GetVertex(v);
    for (int d = 0; d < dim_; ++d) {
      coords_host(v, d) = static_cast<Real>(x[d]);
    }
  }
  coords_ = Kokkos::View<Real**, DeviceMemorySpace>("mfem_coords", nv, dim_);
  Kokkos::deep_copy(coords_, coords_host);

  // Global vertex ids.
  mfem::Array<HYPRE_BigInt> vertex_gids;
  pmesh_.GetGlobalVertexIndices(vertex_gids);
  gids_host_ = Kokkos::View<GO*, HostMemorySpace>("mfem_gids", nv);

  // Ownership: a vertex DOF holder is owned by this rank iff it maps to a
  // local true DOF. Shared vertices are owned by exactly one rank, so every
  // global DOF holder has a single owner across the communicator.
  owned_host_ = Kokkos::View<bool*, HostMemorySpace>("mfem_owned", nv);

  // Placeholder classification: vertex entity dim (0) and local-index ids.
  // RCB partitioning ignores these; they exist to satisfy the FieldLayout and
  // OverlapMask contracts.
  classification_dims_host_ =
    Kokkos::View<LO*, HostMemorySpace>("mfem_class_dims", nv);
  classification_ids_host_ =
    Kokkos::View<LO*, HostMemorySpace>("mfem_class_ids", nv);

  mfem::Array<int> vdofs;
  for (int v = 0; v < nv; ++v) {
    pfes_.GetVertexDofs(v, vdofs);
    PCMS_ALWAYS_ASSERT(vdofs.Size() == 1);
    const int dof = vdofs[0];

    gids_host_(v) = static_cast<GO>(vertex_gids[v]);
    owned_host_(v) = (pfes_.GetLocalTDofNumber(dof) >= 0);
    classification_dims_host_(v) = 0;
    classification_ids_host_(v) = v;
  }

  discretization_ = MakeVertexDiscretization(dim_, coords_host);
}

void MFEMLayout::AssertVertexScalarSpace() const
{
  // Order-1 H1, single scalar component: exactly one DOF per vertex.
  PCMS_ALWAYS_ASSERT(pfes_.GetVDim() == 1);
  PCMS_ALWAYS_ASSERT(pfes_.GetNDofs() == pmesh_.GetNV());
}

std::shared_ptr<const Discretization> MFEMLayout::GetDiscretization()
  const noexcept
{
  return discretization_;
}

int MFEMLayout::GetNumComponents() const
{
  return 1;
}

LO MFEMLayout::GetNumOwnedDofHolder() const
{
  // Number of local DOF holders (all vertices on this rank). The owned mask
  // distinguishes the single-owner subset used for communication.
  return static_cast<LO>(pmesh_.GetNV());
}

GO MFEMLayout::GetNumGlobalDofHolder() const
{
  // Each globally unique vertex corresponds to one true DOF.
  return static_cast<GO>(pfes_.GlobalTrueVSize());
}

Rank1View<const bool, HostMemorySpace> MFEMLayout::GetOwnedHost() const
{
  return make_const_array_view(owned_host_);
}

GlobalIDView<HostMemorySpace> MFEMLayout::GetGidsHost() const
{
  return GlobalIDView<HostMemorySpace>(gids_host_.data(), gids_host_.size());
}

CoordinateView<DeviceMemorySpace> MFEMLayout::GetDOFHolderCoordinates() const
{
  return CoordinateView<DeviceMemorySpace>{coordinate_system_,
                                           MakeConstRank2View(coords_)};
}

bool MFEMLayout::IsDistributed() const
{
  return true;
}

EntOffsetsArray MFEMLayout::GetEntOffsets() const
{
  EntOffsetsArray offsets{};
  offsets.fill(0);
  const auto n = static_cast<size_t>(GetNumOwnedDofHolder());
  // Vertex DOF holders occupy entity dimension 0.
  for (int i = 1; i < ent_offsets_len; ++i) {
    offsets[i] = n;
  }
  return offsets;
}

int MFEMLayout::GetDimension() const
{
  return dim_;
}

Rank1View<const LO, HostMemorySpace>
MFEMLayout::GetDOFHolderClassificationDimensionsHost() const
{
  return make_const_array_view(classification_dims_host_);
}

Rank1View<const LO, HostMemorySpace>
MFEMLayout::GetDOFHolderClassificationIdsHost() const
{
  return make_const_array_view(classification_ids_host_);
}

Kokkos::View<bool*, HostMemorySpace> MFEMLayout::OverlapMaskFromAttribute(
  mfem::ParMesh& pmesh, int attribute)
{
  PCMS_FUNCTION_TIMER;

  const int nv = pmesh.GetNV();
  Kokkos::View<bool*, HostMemorySpace> overlap("mfem_overlap_mask", nv);
  Kokkos::deep_copy(overlap, false);

  mfem::Array<int> verts;
  for (int e = 0; e < pmesh.GetNE(); ++e) {
    if (pmesh.GetAttribute(e) != attribute) {
      continue;
    }
    pmesh.GetElementVertices(e, verts);
    for (int j = 0; j < verts.Size(); ++j) {
      overlap(verts[j]) = true;
    }
  }

  return overlap;
}

} // namespace pcms
