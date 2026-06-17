#ifndef PCMS_FIELD_DATA_MFEM_H
#define PCMS_FIELD_DATA_MFEM_H

#include "pcms/field/field_data.h"
#include "pcms/field/field_metadata.h"
#include "pcms/utility/arrays.h"

#include <mfem.hpp>

namespace pcms
{

// FieldData backend for an MFEM order-1 H1 (vertex) scalar field. The
// coefficient store is the live mfem::ParGridFunction: Get/Set operate on it
// directly so the coupler and the MFEM solver share state.
//
// DOF-holder ordering is the local vertex ordering, matching MFEMLayout. To
// stay consistent across process boundaries, Get/Set round-trip through the
// FE space true DOFs: Get reads owner values via the parallel restriction and
// Set distributes received owner values via the parallel prolongation, filling
// shared (non-owned) vertices on this rank.
class MFEMVertexFieldData : public FieldData<Real>
{
public:
  MFEMVertexFieldData(mfem::ParFiniteElementSpace& pfes,
                      mfem::ParGridFunction& gf, FieldMetadata metadata = {});

  const FieldMetadata& GetMetadata() const override;

  Rank1View<const Real, HostMemorySpace> GetDOFHolderDataHost() const override;
  void SetDOFHolderDataHost(
    Rank1View<const Real, HostMemorySpace> data) override;

  Rank1View<const Real, DeviceMemorySpace> GetDOFHolderData() const override;
  void SetDOFHolderData(Rank1View<const Real, DeviceMemorySpace> data) override;

private:
  mfem::ParFiniteElementSpace& pfes_;
  mfem::ParGridFunction& gf_;
  FieldMetadata metadata_;
  mutable Kokkos::View<Real*, HostMemorySpace> host_data_;
  mutable Kokkos::View<Real*, DeviceMemorySpace> device_data_;
};

} // namespace pcms

#endif // PCMS_FIELD_DATA_MFEM_H
