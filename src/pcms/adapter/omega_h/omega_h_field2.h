#ifndef PCMS_OMEGA_H_FIELD2_H
#define PCMS_OMEGA_H_FIELD2_H

#include <Kokkos_Core.hpp>
#include <vector>

#include "pcms/types.h"
#include "pcms/field.h"
#include "pcms/coordinate_system.h"

namespace pcms
{

// TODO template over possible OmegaHField2Types
class OmegaHField2 : public FieldT<Real>
{
public:
  OmegaHField2(std::string name, CoordinateSystem coordinate_system,
               const OmegaHFieldLayout& layout, Omega_h::Mesh& mesh);

  std::string& GetName() const override;

  mesh_entity_type GetEntityType() const override;

  CoordinateSystem GetCoordinateSystem() const override;

  Kokkos::View<const Real*> GetNodalData() const override;

  void SetNodalData(Kokkos::View<const Real*> data) override;

  void SetEvaluationCoordinates(LocalizationHint hint) override;

  LocalizationHint GetLocalizationHint(
    CoordinateView<HostMemorySpace> coordinate_view) override;

  void Evaluate(FieldDataView<double, HostMemorySpace> results) override;

  void EvaluateGradient(
    FieldDataView<double, HostMemorySpace> results) override;

  const FieldLayout& GetLayout() override;

  bool CanEvaluateGradient() override;

  int Serialize(Rank1View<double, pcms::HostMemorySpace> buffer,
                Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
    const override;

  void Deserialize(Rank1View<const double, pcms::HostMemorySpace> buffer,
                   Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
    const override;

  ~OmegaHField2() noexcept = default;

private:
  std::string name_;
  CoordinateSystem coordinate_system_;
  const OmegaHFieldLayout& layout_;
  Omega_h::Mesh& mesh_;
  Kokkos::View<GridPointSearch::Result*> *hint_ = nullptr;
};

} // namespace pcms

#endif // PCMS_OMEGA_H_FIELD2_H
