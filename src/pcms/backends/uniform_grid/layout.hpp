#ifndef PCMS_BACKENDS_UNIFORM_GRID_LAYOUT_HPP
#define PCMS_BACKENDS_UNIFORM_GRID_LAYOUT_HPP

#include "pcms/utility/arrays.hpp"
#include "pcms/backends/uniform_grid/geometry.hpp"
#include "pcms/field/field_layout.hpp"
#include "pcms/field/field.hpp"
#include "pcms/utility/uniform_grid.hpp"

#include <array>

namespace pcms::backend::uniform_grid
{
template <unsigned Dim = 2>
class UniformGridFieldLayout : public FieldLayout
{
public:
  /// The layout of order 0 (one DOF holder per cell) or 1 (one per vertex) on
  /// the grid of `geometry`, which it is on; it takes the geometry's
  /// coordinate system.
  /// @throws pcms_error if `geometry` is null
  UniformGridFieldLayout(
    std::shared_ptr<const UniformGridGeometry<Dim>> geometry,
    int num_components, int order = 1, std::string name = "");

  int GetNumComponents() const override;
  LO GetNumDofHolder() const override;
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

  const UniformGrid<Dim>& GetGrid() const;
  LO GetNumCells() const;
  LO GetNumVertices() const;
  int GetOrder() const;

private:
  LO CountDofHolders() const;

  UniformGrid<Dim> grid_;
  int num_components_;
  int order_;
  Kokkos::View<GO*, DeviceMemorySpace> gids_;
  Kokkos::View<bool*, DeviceMemorySpace> owned_;
  Kokkos::View<Real**, DeviceMemorySpace> dof_holder_coords_;
  Kokkos::View<GO*, HostMemorySpace> gids_host_;
  Kokkos::View<bool*, HostMemorySpace> owned_host_;
  Kokkos::View<LO*, DeviceMemorySpace> classification_dims_;
  Kokkos::View<LO*, DeviceMemorySpace> classification_ids_;
  Kokkos::View<LO*, HostMemorySpace> classification_dims_host_;
  Kokkos::View<LO*, HostMemorySpace> classification_ids_host_;
};

using UniformGridFieldLayout2D = UniformGridFieldLayout<2>;

} // namespace pcms::backend::uniform_grid
#endif // PCMS_BACKENDS_UNIFORM_GRID_LAYOUT_HPP
