#ifndef PCMS_BACKENDS_OMEGA_H_LAGRANGE_LAYOUT_HPP
#define PCMS_BACKENDS_OMEGA_H_LAGRANGE_LAYOUT_HPP

#include <Omega_h_mesh.hpp>
#include <memory>
#include "pcms/utility/arrays.hpp"
#include "pcms/backends/omega_h/geometry.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/field.hpp"

namespace pcms::backend::omega_h
{

// Layout for native Omega_h Lagrange fields.
//   order 0: one DOF holder per element (centroid coordinates)
//   order 1: one DOF holder per vertex (barycentric interpolation)
// Throws std::invalid_argument for any other order.
class OmegaHLagrangeLayout : public FieldLayout
{
public:
  /// The layout is on `geometry`, whose mesh it lays out, and takes its
  /// coordinate system, its global id tag name and its owned mask. The mask is
  /// per vertex: it applies to the DOF holders of an order-1 layout and not to
  /// those of an order-0 layout.
  /// @throws pcms_error if `geometry` is null
  /// @param name the layout's name (FieldLayout::GetName)
  /// @throws std::invalid_argument if `order` is neither 0 nor 1
  OmegaHLagrangeLayout(std::shared_ptr<const OmegaHGeometry> geometry,
                       int order, int num_components, std::string name = "");

  int GetNumComponents() const override;
  LO GetNumDofHolder() const override;
  [[nodiscard]] bool HasGhosts() const override;
  // Number of DOF holders, owned and ghost: the size of the field data
  // this layout's backends store.
  LO GetNumDofHolderWithGhosts() const;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwnedHost() const override;
  GlobalIDView<HostMemorySpace> GetGidsHost() const override;
  GlobalIDView<DeviceMemorySpace> GetGids() const;
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

  int GetOrder() const;
  Omega_h::Mesh& GetMesh() const;

private:
  Omega_h::Mesh& mesh_;
  int order_;
  int num_components_;

  Omega_h::Write<Omega_h::GO> gids_;
  Kokkos::View<Real**, DeviceMemorySpace> coords_2d_;
  Omega_h::Read<Real> coords_; // device coordinates
  Kokkos::View<bool*, DeviceMemorySpace> owned_;
  Kokkos::View<bool*, HostMemorySpace> owned_host_;
  Omega_h::Read<Omega_h::ClassId> class_ids_;
  Omega_h::Read<Omega_h::I8> class_dims_;
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
  // Requires the owned views and the classification arrays.
  void BuildOwnedClassification();
};

} // namespace pcms::backend::omega_h
#endif // PCMS_BACKENDS_OMEGA_H_LAGRANGE_LAYOUT_HPP
