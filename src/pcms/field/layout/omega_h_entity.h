#ifndef PCMS_FIELD_LAYOUT_OMEGA_H_ENTITY_H
#define PCMS_FIELD_LAYOUT_OMEGA_H_ENTITY_H

#include <Omega_h_mesh.hpp>

#include "pcms/discretization/discretization/omega_h.hpp"
#include "pcms/field/field_layout.h"

namespace pcms
{

class OmegaHEntityLayout : public FieldLayout
{
public:
  OmegaHEntityLayout(Omega_h::Mesh& mesh, int entity_dim, int num_components,
                     std::shared_ptr<const CoordinateSystem> coordinate_system,
                     std::string global_id_name = "global");

  std::shared_ptr<const Discretization> GetDiscretization()
    const noexcept override;

  int GetNumComponents() const override;
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override;
  GlobalIDView<HostMemorySpace> GetGidsHost() const override;
  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override;

  [[nodiscard]] bool IsDistributed() const override;
  EntOffsetsArray GetEntOffsets() const override;
  int GetDimension() const override;

  Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override;

  Rank1View<const LO, HostMemorySpace> GetDOFHolderClassificationIdsHost()
    const override;

private:
  int dimension_;
  int entity_dim_;
  int num_components_;
  GO num_global_dof_holder_;

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
  std::shared_ptr<const Discretization> discretization_;
};

} // namespace pcms

#endif // PCMS_FIELD_LAYOUT_OMEGA_H_ENTITY_H
