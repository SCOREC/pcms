#ifndef PCMS_UNIFORM_GRID_FIELD_LAYOUT_H
#define PCMS_UNIFORM_GRID_FIELD_LAYOUT_H

#include "pcms/utility/arrays.h"
#include "pcms/field_layout.h"
#include "pcms/coordinate_system.h"
#include "pcms/field.h"
#include "pcms/uniform_grid.h"

#include <array>

namespace pcms
{
template <unsigned Dim = 2>
class UniformGridFieldLayout : public FieldLayout
{
public:
  UniformGridFieldLayout(UniformGrid<Dim>& grid, int num_components,
                         CoordinateSystem coordinate_system);

  std::unique_ptr<FieldT<Real>> CreateField() const override;

  int GetNumComponents() const override;
  LO GetNumOwnedDofHolder() const override;
  GO GetNumGlobalDofHolder() const override;

  Rank1View<const bool, HostMemorySpace> GetOwned() const override;
  GlobalIDView<HostMemorySpace> GetGids() const override;
  CoordinateView<HostMemorySpace> GetDOFHolderCoordinates() const override;

  bool IsDistributed() override;

  EntOffsetsArray GetEntOffsets() const override;

  ReversePartitionMap2 GetReversePartitionMap(
    const redev::Partition& partition) const override;

  UniformGrid<Dim>& GetGrid() const;
  LO GetNumCells() const;
  LO GetNumVertices() const;

private:
  UniformGrid<Dim>& grid_;
  int num_components_;
  CoordinateSystem coordinate_system_;
  Kokkos::View<GO*, HostMemorySpace> gids_;
  Kokkos::View<Real**, HostMemorySpace> dof_holder_coords_;
  Kokkos::View<bool*, HostMemorySpace> owned_;
};

using UniformGridFieldLayout2D = UniformGridFieldLayout<2>;

} // namespace pcms
#endif // PCMS_UNIFORM_GRID_FIELD_LAYOUT_H
