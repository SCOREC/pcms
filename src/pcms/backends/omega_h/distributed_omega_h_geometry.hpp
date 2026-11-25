#ifndef PCMS_BACKENDS_OMEGA_H_DISTRIBUTED_OMEGA_H_GEOMETRY_HPP
#define PCMS_BACKENDS_OMEGA_H_DISTRIBUTED_OMEGA_H_GEOMETRY_HPP

#include "pcms/backends/omega_h/distributed_point_search.hpp"
#include "pcms/backends/omega_h/omega_h_geometry.hpp"
#include "pcms/geometry/distributed_geometry.hpp"
#include "pcms/utility/types.hpp"

#include <Omega_h_mesh.hpp>

#include <memory>
#include <variant>
#include <vector>

namespace pcms::omega_h
{

/**
 * \brief The DistributedGeometry of an Omega_h mesh partitioned over ranks:
 * the owner of a query point is the rank that owns the element containing it
 * (DistributedPointSearch). Created by CreateDistributedGeometry.
 */
class DistributedOmegaHGeometry final : public DistributedGeometry
{
public:
  using SearchVariant =
    std::variant<std::unique_ptr<DistributedPointSearch<2>>,
                 std::unique_ptr<DistributedPointSearch<3>>>;

  /**
   * Collective: publishes this rank's part to the routing table.
   * \param local the geometry of this rank's part of the mesh
   * \param grid_divisions cells per axis of the routing grid; 0 derives it
   *        from the global element count. Must be the same on every rank.
   */
  explicit DistributedOmegaHGeometry(
    std::shared_ptr<const OmegaHGeometry> local, LO grid_divisions = 0);

  [[nodiscard]] Omega_h::Mesh& GetMesh() const noexcept { return mesh_; }

  /** The search the owners are located with. */
  [[nodiscard]] const SearchVariant& GetSearch() const noexcept
  {
    return search_;
  }

  /**
   * Starts the search of each point from the element its hint gives
   * (ElementGuesses); the entities are the elements found.
   */
  [[nodiscard]] std::vector<PointOwner> LocateOwners(
    Kokkos::View<const Real**, DeviceMemorySpace> points,
    const OwnerHints& owners) const override;

  [[nodiscard]] int OwnerEntityDim() const noexcept override
  {
    return mesh_.dim();
  }

private:
  Omega_h::Mesh& mesh_;
  SearchVariant search_;
};

} // namespace pcms::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_DISTRIBUTED_OMEGA_H_GEOMETRY_HPP
