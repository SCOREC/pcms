#ifndef PCMS_FIELD_LAYOUT_OMEGA_H_ENTITY_H
#define PCMS_FIELD_LAYOUT_OMEGA_H_ENTITY_H

#include <Omega_h_mesh.hpp>

#include "pcms/field/evaluator/omega_h_geometry.hpp"
#include "pcms/field/field_layout.hpp"

namespace pcms
{

class OmegaHEntityLayout : public FieldLayout
{
public:
  /// The DOF holders are the entities of dimension `entity_dim` of the
  /// geometry's mesh. The layout is on `geometry`, and takes its coordinate
  /// system and its global id tag name.
  /// @throws pcms_error if `geometry` is null
  OmegaHEntityLayout(std::shared_ptr<const OmegaHGeometry> geometry,
                     int entity_dim, int num_components);

  int GetNumComponents() const override;
  LO GetNumDofHolder() const override;
  [[nodiscard]] bool HasGhosts() const override;
  // Number of DOF holders, owned and ghost: the size of the field data
  // this layout's backends store.
  LO GetNumDofHolderWithGhosts() const;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override;
  GlobalIDView<HostMemorySpace> GetGidsHost() const override;
  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override;
  EntOffsetsArray GetEntOffsets() const override;
  Kokkos::View<const LO*, HostMemorySpace> GetOwnedIndicesHost() const override;
  Kokkos::View<const LO*, DeviceMemorySpace> GetOwnedIndices() const override;

  [[nodiscard]] bool IsDistributed() const override;
  int GetDimension() const override;
  int GetDOFHolderEntityDim() const override;

  Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override;
  Rank1View<const LO, HostMemorySpace> GetDOFHolderClassificationIdsHost()
    const override;

private:
  Omega_h::Mesh& mesh_;
  int dimension_;
  int entity_dim_;
  int num_components_;
  GO num_global_dof_holder_;

  Omega_h::Write<Omega_h::GO> gids_;
  Omega_h::Read<Real> coords_; // device coordinates (1D flattened)
  Kokkos::View<Real**, DeviceMemorySpace> coords_2d_; // device coordinates (2D)
  Omega_h::Read<Omega_h::ClassId> class_ids_;
  Omega_h::Read<Omega_h::I8> class_dims_;
  Kokkos::View<bool*, DeviceMemorySpace> owned_;
  Kokkos::View<bool*, HostMemorySpace> owned_host_;
  LO num_owned_ = 0;
  Kokkos::View<GO*, HostMemorySpace> owned_gids_host_;
  Kokkos::View<Real**, DeviceMemorySpace> owned_coords_2d_;
  Kokkos::View<LO*, HostMemorySpace> owned_indices_host_;
  Kokkos::View<LO*, DeviceMemorySpace> owned_indices_;

  Kokkos::View<LO*, HostMemorySpace> owned_classification_dims_host_;
  Kokkos::View<LO*, HostMemorySpace> owned_classification_ids_host_;

  // Entity offsets over the DOF holders with ghosts.
  EntOffsetsArray EntOffsetsWithGhosts() const;
  void BuildOwnedViews();
  // Requires the owned views and the classification arrays with ghosts.
  void BuildOwnedClassification();
};

} // namespace pcms

#endif // PCMS_FIELD_LAYOUT_OMEGA_H_ENTITY_H
