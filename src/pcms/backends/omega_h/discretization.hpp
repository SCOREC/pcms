#ifndef PCMS_BACKENDS_OMEGA_H_DISCRETIZATION_HPP
#define PCMS_BACKENDS_OMEGA_H_DISCRETIZATION_HPP

#include "pcms/discretization/discretization.hpp"

#include <Omega_h_array.hpp>
#include <Omega_h_mesh.hpp>

#include <array>
#include <optional>
#include <string>

namespace pcms::backend::omega_h
{

class OmegaHDiscretization : public Discretization
{
public:
  explicit OmegaHDiscretization(Omega_h::Mesh& mesh);

  // global_id_name: tag holding global vertex/element ids; owned_mask:
  // optional per-vertex ownership used instead of the mesh's partition.
  OmegaHDiscretization(Omega_h::Mesh& mesh, std::string global_id_name,
                       std::optional<Omega_h::Read<Omega_h::I8>> owned_mask =
                         std::nullopt);

  bool SameEntities(const Discretization& other) const noexcept override;

  int GetDimension() const override;

  LO GetNumEntities(int entity_dim) const override;

  Rank1View<const ClassificationDimension, DeviceMemorySpace>
  GetEntityClassificationDimensions(int entity_dim) const override;

  Rank1View<const ClassificationId, DeviceMemorySpace>
  GetEntityClassificationIds(int entity_dim) const override;

  // The Omega_h mesh this discretization is defined on.
  Omega_h::Mesh& GetMesh() const noexcept { return mesh_; }

  const std::string& GetGlobalIdName() const noexcept { return global_id_name_; }

  const std::optional<Omega_h::Read<Omega_h::I8>>& GetOwnedMask()
    const noexcept
  {
    return owned_mask_;
  }

private:
  static constexpr int max_entity_dim_ = Region;

  Omega_h::Mesh& mesh_;
  std::string global_id_name_ = "global";
  std::optional<Omega_h::Read<Omega_h::I8>> owned_mask_;
  std::array<Omega_h::Read<Omega_h::I8>, max_entity_dim_ + 1>
    class_dims_by_dim_;
  std::array<Omega_h::Read<Omega_h::ClassId>, max_entity_dim_ + 1>
    class_ids_by_dim_;
};

} // namespace pcms::backend::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_DISCRETIZATION_HPP
