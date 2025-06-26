#ifndef PCMS_OMEGA_H_FIELD2_H
#define PCMS_OMEGA_H_FIELD2_H

#include <Kokkos_Core.hpp>
#include <MeshField.hpp>
#include <memory>

#include "pcms/adapter/omega_h/omega_h_field_layout.h"
#include "pcms/types.h"
#include "pcms/field.h"
#include "pcms/coordinate_system.h"
#include "pcms/point_search.h"

namespace pcms
{
// TODO template over possible OmegaHField2Types
class OmegaHField2 : public FieldT<Real>
{
public:
  OmegaHField2(std::string name, CoordinateSystem coordinate_system,
               const OmegaHFieldLayout& layout, Omega_h::Mesh& mesh);

  const std::string& GetName() const override;

  CoordinateSystem GetCoordinateSystem() const override;

  LocalizationHint GetLocalizationHint(
    CoordinateView<HostMemorySpace> coordinate_view) const override;

  void Evaluate(LocalizationHint location, FieldDataView<double, HostMemorySpace> results) const override;

  void EvaluateGradient(
    FieldDataView<double, HostMemorySpace> results) override;

  const FieldLayout& GetLayout() const override;

  bool CanEvaluateGradient() override;

  int Serialize(Rank1View<double, pcms::HostMemorySpace> buffer,
                Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
    const override;

  void Deserialize(Rank1View<const double, pcms::HostMemorySpace> buffer,
                   Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
    override;

  FieldDataView<const Real, HostMemorySpace> GetDOFHolderData() const override;
  void SetDOFHolderData(FieldDataView<const Real, HostMemorySpace> data) override;
  CoordinateView<HostMemorySpace> GetDOFHolderCoordinates() const;

  ~OmegaHField2() noexcept = default;

private:
  std::string name_;
  CoordinateSystem coordinate_system_;
  const OmegaHFieldLayout& layout_;
  Omega_h::Mesh& mesh_;
  MeshField::OmegahMeshField<DefaultExecutionSpace, 2> mesh_field_;
  decltype(mesh_field_.template CreateLagrangeField<Real, 1>()) shape_field_;
  GridPointSearch search_;
  Kokkos::View<Real*> dof_holder_data_;
};

} // namespace pcms

#endif // PCMS_OMEGA_H_FIELD2_H
