#include "pcms/field/data/mfem.h"

#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/profile.h"

namespace pcms
{

MFEMVertexFieldData::MFEMVertexFieldData(mfem::ParFiniteElementSpace& pfes,
                                         mfem::ParGridFunction& gf,
                                         FieldMetadata metadata)
  : pfes_(pfes),
    gf_(gf),
    metadata_(metadata),
    host_data_("mfem_field_data_host", pfes.GetNDofs()),
    device_data_("mfem_field_data_device", pfes.GetNDofs())
{
  PCMS_ALWAYS_ASSERT(pfes_.GetVDim() == 1);
}

const FieldMetadata& MFEMVertexFieldData::GetMetadata() const
{
  return metadata_;
}

Rank1View<const Real, HostMemorySpace>
MFEMVertexFieldData::GetDOFHolderDataHost() const
{
  PCMS_FUNCTION_TIMER;

  // Pull owner values into a true-DOF vector via the parallel restriction so
  // every shared vertex reports its owner's value.
  mfem::Vector true_values(pfes_.GetTrueVSize());
  gf_.GetTrueDofs(true_values);

  mfem::Array<int> vdofs;
  const int nv = pfes_.GetNDofs();
  for (int v = 0; v < nv; ++v) {
    pfes_.GetVertexDofs(v, vdofs);
    PCMS_ALWAYS_ASSERT(vdofs.Size() == 1);
    const int lt = pfes_.GetLocalTDofNumber(vdofs[0]);
    host_data_(v) = (lt >= 0) ? static_cast<Real>(true_values[lt]) : Real{0};
  }

  return make_const_array_view(host_data_);
}

void MFEMVertexFieldData::SetDOFHolderDataHost(
  Rank1View<const Real, HostMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;
  PCMS_ALWAYS_ASSERT(static_cast<int>(data.size()) == pfes_.GetNDofs());

  // Scatter received owner values into a true-DOF vector, then distribute to
  // all local DOFs (including shared, non-owned vertices) via the parallel
  // prolongation.
  mfem::Vector true_values(pfes_.GetTrueVSize());
  mfem::Array<int> vdofs;
  const int nv = pfes_.GetNDofs();
  for (int v = 0; v < nv; ++v) {
    pfes_.GetVertexDofs(v, vdofs);
    PCMS_ALWAYS_ASSERT(vdofs.Size() == 1);
    const int lt = pfes_.GetLocalTDofNumber(vdofs[0]);
    if (lt >= 0) {
      true_values[lt] = static_cast<double>(data[v]);
    }
  }

  gf_.SetFromTrueDofs(true_values);
}

Rank1View<const Real, DeviceMemorySpace>
MFEMVertexFieldData::GetDOFHolderData() const
{
  GetDOFHolderDataHost();
  Kokkos::deep_copy(device_data_, host_data_);
  return make_const_array_view(device_data_);
}

void MFEMVertexFieldData::SetDOFHolderData(
  Rank1View<const Real, DeviceMemorySpace> data)
{
  PCMS_ALWAYS_ASSERT(static_cast<int>(data.size()) == pfes_.GetNDofs());
  CopyDeviceRank1ViewToHostView(host_data_, data);
  SetDOFHolderDataHost(make_const_array_view(host_data_));
}

} // namespace pcms
