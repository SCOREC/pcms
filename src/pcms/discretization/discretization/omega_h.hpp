#ifndef PCMS_DISCRETIZATION_DISCRETIZATION_OMEGA_H_H
#define PCMS_DISCRETIZATION_DISCRETIZATION_OMEGA_H_H

#include "pcms/discretization/discretization.h"

#include <Omega_h_array.hpp>
#include <Omega_h_mesh.hpp>

#include <array>

namespace pcms
{

class OmegaHDiscretization : public Discretization
{
public:
  explicit OmegaHDiscretization(Omega_h::Mesh& mesh);

  bool SameEntities(const Discretization& other) const noexcept override;

  int GetDimension() const override;

  LO GetNumEntities(int entity_dim) const override;

  Rank1View<const ClassificationDimension, DeviceMemorySpace>
  GetEntityClassificationDimensions(int entity_dim) const override;

  Rank1View<const ClassificationId, DeviceMemorySpace>
  GetEntityClassificationIds(int entity_dim) const override;

  // The Omega_h mesh this discretization is defined on.
  Omega_h::Mesh& GetMesh() const noexcept { return mesh_; }

private:
  static constexpr int max_entity_dim_ = Region;

  Omega_h::Mesh& mesh_;
  std::array<Omega_h::Read<Omega_h::I8>, max_entity_dim_ + 1>
    class_dims_by_dim_;
  std::array<Omega_h::Read<Omega_h::ClassId>, max_entity_dim_ + 1>
    class_ids_by_dim_;
};

} // namespace pcms

#endif // PCMS_DISCRETIZATION_DISCRETIZATION_OMEGA_H_H
