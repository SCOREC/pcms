#ifndef PCMS_DISCRETIZATION_DISCRETIZATION_UNIFORM_GRID_H
#define PCMS_DISCRETIZATION_DISCRETIZATION_UNIFORM_GRID_H

#include "pcms/discretization/discretization.h"
#include "pcms/utility/uniform_grid.h"

#include <Kokkos_Core.hpp>

namespace pcms
{

template <unsigned Dim>
class UniformGridDiscretization : public Discretization
{
public:
  explicit UniformGridDiscretization(UniformGrid<Dim> grid);

  bool SameEntities(const Discretization& other) const noexcept override;

  int GetDimension() const override;

  LO GetNumEntities(int entity_dim) const override;

  Rank1View<const ClassificationDimension, DeviceMemorySpace>
  GetEntityClassificationDimensions(int entity_dim) const override;

  Rank1View<const ClassificationId, DeviceMemorySpace>
  GetEntityClassificationIds(int entity_dim) const override;

private:
  static constexpr int CellEntityDim = static_cast<int>(Dim);

  UniformGrid<Dim> grid_;
  Kokkos::View<ClassificationDimension*, DeviceMemorySpace> vertex_class_dims_;
  Kokkos::View<ClassificationId*, DeviceMemorySpace> vertex_class_ids_;
  Kokkos::View<ClassificationDimension*, DeviceMemorySpace> cell_class_dims_;
  Kokkos::View<ClassificationId*, DeviceMemorySpace> cell_class_ids_;
};

template <unsigned Dim>
UniformGridDiscretization<Dim>::UniformGridDiscretization(UniformGrid<Dim> grid)
  : grid_(std::move(grid)),
    vertex_class_dims_("uniform_grid_vertex_class_dims",
                       [&]() {
                         LO n = 1;
                         for (unsigned d = 0; d < Dim; ++d)
                           n *= (grid_.divisions[d] + 1);
                         return n;
                       }()),
    vertex_class_ids_("uniform_grid_vertex_class_ids",
                      vertex_class_dims_.extent(0)),
    cell_class_dims_("uniform_grid_cell_class_dims", grid_.GetNumCells()),
    cell_class_ids_("uniform_grid_cell_class_ids", cell_class_dims_.extent(0))
{
  Kokkos::deep_copy(vertex_class_dims_,
                    static_cast<ClassificationDimension>(Vertex));
  Kokkos::deep_copy(vertex_class_ids_, static_cast<ClassificationId>(0));
  Kokkos::deep_copy(cell_class_dims_,
                    static_cast<ClassificationDimension>(CellEntityDim));
  Kokkos::deep_copy(cell_class_ids_, static_cast<ClassificationId>(0));
}

template <unsigned Dim>
bool UniformGridDiscretization<Dim>::SameEntities(
  const Discretization& other) const noexcept
{
  auto p = dynamic_cast<const UniformGridDiscretization*>(&other);
  if (p == nullptr)
    return false;
  for (unsigned d = 0; d < Dim; ++d) {
    if (grid_.bot_left[d] != p->grid_.bot_left[d] ||
        grid_.edge_length[d] != p->grid_.edge_length[d] ||
        grid_.divisions[d] != p->grid_.divisions[d]) {
      return false;
    }
  }
  return true;
}

template <unsigned Dim>
int UniformGridDiscretization<Dim>::GetDimension() const
{
  return static_cast<int>(Dim);
}

template <unsigned Dim>
LO UniformGridDiscretization<Dim>::GetNumEntities(int entity_dim) const
{
  if (entity_dim == Vertex)
    return static_cast<LO>(vertex_class_dims_.extent(0));
  if (entity_dim == CellEntityDim)
    return static_cast<LO>(cell_class_dims_.extent(0));
  return 0;
}

template <unsigned Dim>
Rank1View<const ClassificationDimension, DeviceMemorySpace>
UniformGridDiscretization<Dim>::GetEntityClassificationDimensions(
  int entity_dim) const
{
  if (entity_dim == Vertex)
    return make_const_array_view(vertex_class_dims_);
  if (entity_dim == CellEntityDim)
    return make_const_array_view(cell_class_dims_);
  return Rank1View<const ClassificationDimension, DeviceMemorySpace>(nullptr,
                                                                     0);
}

template <unsigned Dim>
Rank1View<const ClassificationId, DeviceMemorySpace>
UniformGridDiscretization<Dim>::GetEntityClassificationIds(int entity_dim) const
{
  if (entity_dim == Vertex)
    return make_const_array_view(vertex_class_ids_);
  if (entity_dim == CellEntityDim)
    return make_const_array_view(cell_class_ids_);
  return Rank1View<const ClassificationId, DeviceMemorySpace>(nullptr, 0);
}

} // namespace pcms

#endif // PCMS_DISCRETIZATION_DISCRETIZATION_UNIFORM_GRID_H
