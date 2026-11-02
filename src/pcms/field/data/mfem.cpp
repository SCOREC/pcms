#include "pcms/field/data/mfem.hpp"

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
    host_data_("mfem_field_data_host", pfes.GetNDofs(), 1),
    device_data_("mfem_field_data_device", pfes.GetNDofs(), 1)
{
  PCMS_ALWAYS_ASSERT(pfes_.GetVDim() == 1);
}

const FieldMetadata& MFEMVertexFieldData::GetMetadata() const
{
  return metadata_;
}

Rank2View<const Real, HostMemorySpace>
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
    host_data_(v, 0) = (lt >= 0) ? static_cast<Real>(true_values[lt]) : Real{0};
  }

  return MakeConstRank2View(host_data_);
}

void MFEMVertexFieldData::SetDOFHolderDataHost(
  Rank2View<const Real, HostMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;
  const int nv = pfes_.GetNDofs();
  PCMS_ALWAYS_ASSERT(static_cast<int>(data.extent(0)) == nv);
  PCMS_ALWAYS_ASSERT(data.extent(1) == 1);

  // Scatter received owner values into a true-DOF vector, then distribute to
  // all local DOFs (including shared, non-owned vertices) via the parallel
  // prolongation.
  mfem::Vector true_values(pfes_.GetTrueVSize());
  mfem::Array<int> vdofs;
  for (int v = 0; v < nv; ++v) {
    pfes_.GetVertexDofs(v, vdofs);
    PCMS_ALWAYS_ASSERT(vdofs.Size() == 1);
    const int lt = pfes_.GetLocalTDofNumber(vdofs[0]);
    if (lt >= 0) {
      true_values[lt] = static_cast<double>(data(v, 0));
    }
  }

  gf_.SetFromTrueDofs(true_values);
}

Rank2View<const Real, DeviceMemorySpace> MFEMVertexFieldData::GetDOFHolderData()
  const
{
  GetDOFHolderDataHost();
  DeepCopyMismatchLayouts(device_data_, host_data_);
  return MakeConstRank2View(device_data_);
}

void MFEMVertexFieldData::SetDOFHolderData(
  Rank2View<const Real, DeviceMemorySpace> data)
{
  CopyDeviceRank2ViewToDeviceView(device_data_, data);
  DeepCopyMismatchLayouts(host_data_, device_data_);
  SetDOFHolderDataHost(MakeConstRank2View(host_data_));
}

} // namespace pcms
