#ifndef PCMS_BACKENDS_OMEGA_H_MESHFIELDS_LAYOUT_HPP
#define PCMS_BACKENDS_OMEGA_H_MESHFIELDS_LAYOUT_HPP

#include <Omega_h_mesh.hpp>

#include "pcms/utility/arrays.hpp"
#include "pcms/backends/omega_h/discretization.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/field.hpp"

#include <array>

namespace pcms::backend::omega_h_meshfields
{
class MeshFieldsAdapterLayout : public FieldLayout
{
public:
  MeshFieldsAdapterLayout(
    Omega_h::Mesh& mesh, std::array<int, 4> nodes_per_dim, int num_components,
    std::shared_ptr<const CoordinateSystem> coordinate_system,
    std::string global_id_name = "global");

  std::shared_ptr<const Discretization> GetDiscretization()
    const noexcept override;

  int GetNumComponents() const override;
  // nodes for standard lagrange FEM
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override;
  GlobalIDView<HostMemorySpace> GetGidsHost() const override;
  CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates() const override;

  // returns true if the field layout is distributed
  // if the field layout is distributed, the owned and global dofs are the same
  [[nodiscard]] bool IsDistributed() const override;

  EntOffsetsArray GetEntOffsets() const override;

  int GetDimension() const override;

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
  Omega_h::HostWrite<Omega_h::GO> gids_host_;
  std::string global_id_name_;
  int num_components_;
  std::array<int, 4> nodes_per_dim_;
  Kokkos::View<Real**> dof_holder_coords_;
  Omega_h::Write<Omega_h::ClassId> class_ids_;
  Omega_h::Write<Omega_h::I8> class_dims_;
  Omega_h::HostWrite<Omega_h::ClassId> class_ids_host_;
  Omega_h::HostWrite<Omega_h::I8> class_dims_host_;
  Kokkos::View<bool*> owned_;
  Kokkos::View<bool*, HostMemorySpace> owned_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_dims_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_ids_host_;
  std::shared_ptr<const Discretization> discretization_;
};

} // namespace pcms::backend::omega_h_meshfields
#endif // PCMS_BACKENDS_OMEGA_H_MESHFIELDS_LAYOUT_HPP
