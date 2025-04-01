#ifndef PCMS_OMEGA_H_FIELD2_H
#define PCMS_OMEGA_H_FIELD2_H

#include <Kokkos_Core.hpp>
#include <memory>

#include "pcms/adapter/omega_h/omega_h_field_layout.h"
#include "pcms/types.h"
#include "pcms/field.h"
#include "pcms/coordinate_system.h"
#include "pcms/point_search.h"

namespace pcms
{

class OmegaHFIeld2LocalizationHint
{
public:
  OmegaHFIeld2LocalizationHint();

  OmegaHFIeld2LocalizationHint(
    Kokkos::View<GridPointSearch::Result*> search_results);

  Kokkos::View<GridPointSearch::Result*> GetResults();

private:
  Kokkos::View<GridPointSearch::Result*> search_results_;
};

// TODO template over possible OmegaHField2Types
class OmegaHField2 : public FieldT<Real>
{
public:
  OmegaHField2(std::string name, CoordinateSystem coordinate_system,
               const OmegaHFieldLayout& layout, Omega_h::Mesh& mesh);

  const std::string& GetName() const override;

  CoordinateSystem GetCoordinateSystem() const override;

  Rank1View<const Real, pcms::HostMemorySpace> GetNodalData() const override;

  Rank1View<const Real, pcms::HostMemorySpace> GetNodalCoordinates()
    const override;

  void SetNodalData(Rank1View<const Real, pcms::HostMemorySpace> data) override;

  LocalizationHint GetLocalizationHint(
    CoordinateView<HostMemorySpace> coordinate_view) override;

  void Evaluate(LocalizationHint location, FieldDataView<double, HostMemorySpace> results) override;

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

  void SetDOFHolderData(FieldDataView<Real, HostMemorySpace> data);
  CoordinateView<HostMemorySpace> GetDOFHolderCoordinates();

  ~OmegaHField2() noexcept = default;

private:
  std::string name_;
  CoordinateSystem coordinate_system_;
  const OmegaHFieldLayout& layout_;
  Omega_h::Mesh& mesh_;
  GridPointSearch search_;
};

} // namespace pcms

#endif // PCMS_OMEGA_H_FIELD2_H
