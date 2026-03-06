#ifndef PCMS_UNIFORM_GRID_FIELD_H
#define PCMS_UNIFORM_GRID_FIELD_H

#include "pcms/adapter/uniform_grid/uniform_grid_field_layout.h"
#include "pcms/utility/types.h"
#include "pcms/field.h"
#include "pcms/coordinate_system.h"
#include "pcms/uniform_grid.h"
#include <memory>

namespace pcms
{
template <unsigned Dim = 2>
class UniformGridField : public FieldT<Real>
{
public:
  UniformGridField(const UniformGridFieldLayout<Dim>& layout);

  LocalizationHint GetLocalizationHint(
    CoordinateView<HostMemorySpace> coordinate_view) const override;

  void Evaluate(LocalizationHint location,
                FieldDataView<Real, HostMemorySpace> results) const override;

  void EvaluateGradient(FieldDataView<Real, HostMemorySpace> results) override;

  const FieldLayout& GetLayout() const override;

  bool CanEvaluateGradient() override;

  int Serialize(Rank1View<Real, pcms::HostMemorySpace> buffer,
                Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
    const override;

  void Deserialize(
    Rank1View<const Real, pcms::HostMemorySpace> buffer,
    Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) override;

  Rank1View<const Real, HostMemorySpace> GetDOFHolderData() const override;
  void SetDOFHolderData(Rank1View<const Real, HostMemorySpace> data) override;

  View<Dim, Real, HostMemorySpace> to_mdspan();
  View<Dim, const Real, HostMemorySpace> to_mdspan() const;

  ~UniformGridField() noexcept = default;

private:
  const UniformGridFieldLayout<Dim>& layout_;
  UniformGrid<Dim>& grid_;
  Kokkos::View<Real*, HostMemorySpace> dof_holder_data_;
};

using UniformGridField2D = UniformGridField<2>;

} // namespace pcms

#endif // PCMS_UNIFORM_GRID_FIELD_H
