#ifndef PCMS_UNIFORM_GRID_FIELD_LAYOUT_H
#define PCMS_UNIFORM_GRID_FIELD_LAYOUT_H

#include "pcms/utility/arrays.hpp"
#include "pcms/backends/uniform_grid/discretization.hpp"
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
  UniformGridFieldLayout(
    UniformGrid<Dim> grid, int num_components,
    std::shared_ptr<const CoordinateSystem> coordinate_system, int order = 1);

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

  const UniformGrid<Dim>& GetGrid() const;
  LO GetNumCells() const;
  LO GetNumVertices() const;
  int GetOrder() const;

private:
  LO GetNumDofHolders() const;

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
  std::shared_ptr<const Discretization> discretization_;
};

using UniformGridFieldLayout2D = UniformGridFieldLayout<2>;

} // namespace pcms::backend::uniform_grid
#endif // PCMS_UNIFORM_GRID_FIELD_LAYOUT_H
