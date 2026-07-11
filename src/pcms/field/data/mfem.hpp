#ifndef PCMS_FIELD_DATA_MFEM_HPP
#define PCMS_FIELD_DATA_MFEM_HPP

#include "pcms/field/field_data.h"
#include "pcms/field/field_metadata.h"
#include "pcms/utility/arrays.h"

#include <mfem.hpp>

namespace pcms
{

// FieldData backend for an MFEM order-1 H1 (vertex) scalar field.
class MFEMVertexFieldData : public FieldData<Real>
{
public:
  MFEMVertexFieldData(mfem::ParFiniteElementSpace& pfes,
                      mfem::ParGridFunction& gf, FieldMetadata metadata = {});

  const FieldMetadata& GetMetadata() const override;

  Rank2View<const Real, HostMemorySpace> GetDOFHolderDataHost() const override;
  void SetDOFHolderDataHost(
    Rank2View<const Real, HostMemorySpace> data) override;

  Rank2View<const Real, DeviceMemorySpace> GetDOFHolderData() const override;
  void SetDOFHolderData(Rank2View<const Real, DeviceMemorySpace> data) override;

private:
  mfem::ParFiniteElementSpace& pfes_;
  mfem::ParGridFunction& gf_;
  FieldMetadata metadata_;
  mutable Kokkos::View<Real**, HostMemorySpace> host_data_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> device_data_;
};

} // namespace pcms

#endif // PCMS_FIELD_DATA_MFEM_HPP
