#ifndef PCMS_ADAPTER_OMEGA_H_LAGRANGE_LAYOUT_H
#define PCMS_ADAPTER_OMEGA_H_LAGRANGE_LAYOUT_H

#include <Omega_h_mesh.hpp>
#include "pcms/utility/arrays.h"
#include "pcms/discretization/discretization/omega_h.hpp"
#include "pcms/field/field_layout.h"
#include "pcms/field/field.h"

namespace pcms
{

// Layout for native Omega_h Lagrange fields.
//   order 0: one DOF holder per element (centroid coordinates)
//   order 1: one DOF holder per vertex (barycentric interpolation)
// Throws std::invalid_argument for any other order.
class OmegaHLagrangeLayout : public FieldLayout
{
public:
  OmegaHLagrangeLayout(
    Omega_h::Mesh& mesh, int order, int num_components,
    std::shared_ptr<const CoordinateSystem> coordinate_system,
    std::string global_id_name = "global");
  OmegaHLagrangeLayout(
    Omega_h::Mesh& mesh, int order, int num_components,
    std::shared_ptr<const CoordinateSystem> coordinate_system,
    Omega_h::Read<Omega_h::I8> owned_mask,
    std::string global_id_name = "global");

  std::shared_ptr<const Discretization> GetDiscretization()
    const noexcept override;

  int GetNumComponents() const override;
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override;
  GlobalIDView<HostMemorySpace> GetGidsHost() const override;
  GlobalIDView<DeviceMemorySpace> GetGids() const;
  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override;

  [[nodiscard]] bool IsDistributed() const override;

  EntOffsetsArray GetEntOffsets() const override;

  int GetDimension() const override;

  Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override;

  Rank1View<const LO, HostMemorySpace> GetDOFHolderClassificationIdsHost()
    const override;

  int GetOrder() const;
  Omega_h::Mesh& GetMesh() const;

private:
  Omega_h::Mesh& mesh_;
  int order_;
  int num_components_;
  std::string global_id_name_;

  Omega_h::Write<Omega_h::GO> gids_;
  Omega_h::HostWrite<Omega_h::GO> gids_host_;
  Kokkos::View<Real**, DeviceMemorySpace> coords_2d_;
  Omega_h::Read<Real> coords_; // device coordinates
  Kokkos::View<bool*, DeviceMemorySpace> owned_;
  Kokkos::View<bool*, HostMemorySpace> owned_host_;
  Omega_h::Read<Omega_h::ClassId> class_ids_;
  Omega_h::Read<Omega_h::I8> class_dims_;
  Kokkos::View<LO*, HostMemorySpace> classification_dims_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_ids_host_;
  std::shared_ptr<const Discretization> discretization_;
};

} // namespace pcms
#endif // PCMS_ADAPTER_OMEGA_H_LAGRANGE_LAYOUT_H
