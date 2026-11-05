#ifndef PCMS_BACKENDS_MFEM_FIELD_DATA_HPP
#define PCMS_BACKENDS_MFEM_FIELD_DATA_HPP

#include "pcms/field/field_data.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.hpp"

#include <mfem.hpp>

namespace pcms::backend::mfem
{

// FieldData backend for an MFEM order-1 H1 (vertex) scalar field.
class MFEMVertexFieldData
  : public FieldData<Real>
  , public DOFAccess<Real>
{
public:
  MFEMVertexFieldData(::mfem::ParFiniteElementSpace& pfes,
                      ::mfem::ParGridFunction& gf, ValueBasis basis = {});

  FieldValueType GetValueType() const override;
  const ValueBasis& GetValueBasis() const override;

  const DOFAccess<Real>* GetDOFAccess() const noexcept override { return this; }
  DOFAccess<Real>* GetDOFAccess() noexcept override { return this; }

  Rank2View<const Real, DeviceMemorySpace> GetOwnedValues() const override;
  void SetOwnedValues(Rank2View<const Real, DeviceMemorySpace> data) override;

private:
  Rank2View<const Real, HostMemorySpace> GatherHost() const;
  void ScatterHost(Rank2View<const Real, HostMemorySpace> data);

  ::mfem::ParFiniteElementSpace& pfes_;
  ::mfem::ParGridFunction& gf_;
  ValueBasis basis_;
  mutable Kokkos::View<Real**, HostMemorySpace> host_data_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> device_data_;
};

} // namespace pcms::backend::mfem

#endif // PCMS_BACKENDS_MFEM_FIELD_DATA_HPP
