#ifndef PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_ADAPTER_LAYOUT_H
#define PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_ADAPTER_LAYOUT_H

#include <Omega_h_mesh.hpp>

#include "pcms/utility/arrays.hpp"
#include "pcms/backends/omega_h/geometry.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/field.hpp"

#include <array>

namespace pcms
{
class MeshFieldsAdapterLayout : public FieldLayout
{
public:
  /// The layout is on `geometry`, whose mesh it lays out, and takes its
  /// coordinate system and its global id tag name; it does not apply an owned
  /// mask.
  /// @throws pcms_error if `geometry` is null
  /// @param name the layout's name (FieldLayout::GetName)
  MeshFieldsAdapterLayout(std::shared_ptr<const OmegaHGeometry> geometry,
                          std::array<int, 4> nodes_per_dim, int num_components,
                          std::string name = "");

  int GetNumComponents() const override;
  // nodes for standard lagrange FEM
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

  // returns true if the field layout is distributed: on several ranks it may
  // hold ghost DOF holders in addition to the owned ones (see HasGhosts)
  [[nodiscard]] bool IsDistributed() const override;

  int GetDimension() const override;
  int GetDOFHolderEntityDim() const override;

  Rank1View<const LO, HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override;
  Rank1View<const LO, HostMemorySpace> GetDOFHolderClassificationIdsHost()
    const override;

  std::array<int, 4> GetNodesPerDim() const;
  size_t GetNumEnts() const;
  Omega_h::Mesh& GetMesh() const;

private:
  Omega_h::Read<Omega_h::ClassId> GetClassIDs() const;
  Omega_h::Read<Omega_h::I8> GetClassDims() const;

  Omega_h::Mesh& mesh_;
  Omega_h::Write<Omega_h::GO> gids_;
  int num_components_;
  std::array<int, 4> nodes_per_dim_;
  Kokkos::View<Real**> dof_holder_coords_;
  Omega_h::Write<Omega_h::ClassId> class_ids_;
  Omega_h::Write<Omega_h::I8> class_dims_;
  Omega_h::HostWrite<Omega_h::ClassId> class_ids_host_;
  Omega_h::HostWrite<Omega_h::I8> class_dims_host_;
  Kokkos::View<bool*> owned_;
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
#endif // PCMS_ADAPTER_MESHFIELDS_MESH_FIELDS_ADAPTER_LAYOUT_H
