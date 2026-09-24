#ifndef PCMS_FIELD_LAYOUT_OMEGA_H_ENTITY_H
#define PCMS_FIELD_LAYOUT_OMEGA_H_ENTITY_H

#include <Omega_h_mesh.hpp>

#include "pcms/discretization/discretization/omega_h.hpp"
#include "pcms/field/coordinate_system.h"
#include "pcms/field/field_layout.h"

namespace pcms
{

// Layout for fields with one DOF holder on each entity of a single Omega_h
// mesh dimension. Unlike OmegaHLagrangeLayout, this is not limited to
// Lagrange orders and directly represents "one value per chosen entity"
// layouts such as face-centroid reconstruction sites.
class OmegaHEntityLayout : public FieldLayout
{
public:
  OmegaHEntityLayout(Omega_h::Mesh& mesh, int entity_dim, int num_components,
                     CoordinateSystem coordinate_system,
                     std::string global_id_name = "global");

  std::shared_ptr<const Discretization> GetDiscretization()
    const noexcept override;

  int GetNumComponents() const override;
  LO GetNumLocalDofHolder() const override;
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override;
  GlobalIDView<HostMemorySpace> GetGidsHost() const override;
  GlobalIDView<DeviceMemorySpace> GetGids() const override;
  GlobalIDView<HostMemorySpace> GetOwnedGidsHost() const override;
  GlobalIDView<DeviceMemorySpace> GetOwnedGids() const override;
  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override;
  CoordinateView<DeviceMemorySpace> GetOwnedDOFHolderCoordinates()
    const override;
  Kokkos::View<const LO*, HostMemorySpace> GetOwnedToLocalHost() const override;
  Kokkos::View<const LO*, DeviceMemorySpace> GetOwnedToLocal() const override;

  [[nodiscard]] bool IsDistributed() const override;
  EntOffsetsArray GetEntOffsets() const override;
  int GetDimension() const override;
  int GetDOFHolderEntityDim() const override;

  Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override;

  Rank1View<const LO, HostMemorySpace> GetDOFHolderClassificationIdsHost()
    const override;

private:
  int dimension_;
  int entity_dim_;
  int num_components_;
  GO num_global_dof_holder_;
  CoordinateSystem coordinate_system_;

  Omega_h::Write<Omega_h::GO> gids_;
  Omega_h::HostWrite<Omega_h::GO> gids_host_;
  Omega_h::Read<Real> coords_; // device coordinates (1D flattened)
  Kokkos::View<Real**, DeviceMemorySpace> coords_2d_; // device coordinates (2D)
  Omega_h::Read<Omega_h::ClassId> class_ids_;
  Omega_h::Read<Omega_h::I8> class_dims_;
  Kokkos::View<bool*, DeviceMemorySpace> owned_;
  Kokkos::View<bool*, HostMemorySpace> owned_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_dims_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_ids_host_;
  LO num_owned_ = 0;
  Kokkos::View<GO*, HostMemorySpace> owned_gids_host_;
  Kokkos::View<Real**, DeviceMemorySpace> owned_coords_2d_;
  Kokkos::View<LO*, HostMemorySpace> owned_to_local_host_;
  Kokkos::View<GO*, DeviceMemorySpace> owned_gids_;
  Kokkos::View<LO*, DeviceMemorySpace> owned_to_local_;

  void BuildOwnedViews();

  std::shared_ptr<const Discretization> discretization_;
};

} // namespace pcms

#endif // PCMS_FIELD_LAYOUT_OMEGA_H_ENTITY_H
