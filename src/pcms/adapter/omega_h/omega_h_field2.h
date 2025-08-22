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
class MeshFieldBackend
{
public:
  virtual ~MeshFieldBackend() = default;
  virtual Kokkos::View<Real* [1]> evaluate(Kokkos::View<Real**> localCoords,
                                           Kokkos::View<LO*> offsets) const = 0;
  virtual void SetData(Rank1View<const Real, HostMemorySpace> data,
                       size_t num_nodes, size_t num_components, int dim) = 0;
  virtual void GetData(Rank1View<Real, HostMemorySpace> data, size_t num_nodes,
                       size_t num_components, int dim) const = 0;
};

// TODO template over possible OmegaHField2Types
class OmegaHField2 : public FieldT<Real>
{
public:
  OmegaHField2(std::string name, const OmegaHFieldLayout& layout, Omega_h::Mesh& mesh);

  const std::string& GetName() const override;


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

  ~OmegaHField2() noexcept = default;

private:
  std::string name_;
  const OmegaHFieldLayout& layout_;
  Omega_h::Mesh& mesh_;
  std::unique_ptr<MeshFieldBackend> mesh_field_;
  GridPointSearch search_;
  Kokkos::View<Real*> dof_holder_data_;
};

} // namespace pcms

#endif // PCMS_OMEGA_H_FIELD2_H
