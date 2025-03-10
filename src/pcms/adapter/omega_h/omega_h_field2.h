#ifndef PCMS_OMEGA_H_FIELD2_H
#define PCMS_OMEGA_H_FIELD2_H

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

  void SetEvaluationCoordinates(CoordinateView<HostMemorySpace> coordinates,
                                LocalizationHint hint = {}) override;

  LocalizationHint GetLocalizationHint(
    CoordinateView<HostMemorySpace> coordinates) override;

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
  CoordinateSystem coordinate_system_;
  const OmegaHFieldLayout& layout_;
  Omega_h::Mesh& mesh_;
};

} // namespace pcms

#endif // PCMS_OMEGA_H_FIELD2_H
