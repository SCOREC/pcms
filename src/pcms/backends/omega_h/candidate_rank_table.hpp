#ifndef PCMS_BACKENDS_OMEGA_H_CANDIDATE_RANK_TABLE_HPP
#define PCMS_BACKENDS_OMEGA_H_CANDIDATE_RANK_TABLE_HPP

#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/localization/candidate_rank_table.hpp"
#include "pcms/utility/bounding_box.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/detail/uniform_grid.hpp"
#include "pcms/utility/uniform_grid.hpp"

#include <Omega_h_bbox.hpp>
#include <Omega_h_mesh.hpp>

#include <optional>
#include <stdexcept>

namespace pcms::omega_h
{

/**
 * Cells per axis giving about one element per cell of the whole
 * decomposition, the same on every rank (`nelems()` counts ghosts, so it
 * differs between ranks). Collective.
 */
template <int Dim>
[[nodiscard]] LO DivisionsPerAxisForGlobalMesh(Omega_h::Mesh& mesh,
                                               LO min_divisions = 10)
{
  return pcms::detail::DivisionsPerAxisForMesh<Dim>(
    static_cast<LO>(mesh.nglobal_ents(mesh.dim())), min_divisions);
}

/**
 * The shared routing grid for `mesh`, identical on every rank: the globally
 * reduced bounding box, about one element per cell of the global element
 * count. Collective.
 */
template <int Dim>
[[nodiscard]] UniformGrid<Dim> SharedRoutingGrid(Omega_h::Mesh& mesh,
                                                 LO min_divisions = 10)
{
  return CreateUniformGridFromMesh<Dim>(
    mesh, DivisionsPerAxisForGlobalMesh<Dim>(mesh, min_divisions));
}

/**
 * The routing grid for `mesh` with `divisions` cells per axis, or the default
 * resolution (`SharedRoutingGrid`) when `divisions` is 0. Collective.
 *
 * `divisions` must be the same on every rank, or cell ids would name different
 * regions of space on different ranks.
 */
template <int Dim>
[[nodiscard]] UniformGrid<Dim> RoutingGrid(Omega_h::Mesh& mesh, LO divisions)
{
  if (divisions > 0) {
    return CreateUniformGridFromMesh<Dim>(mesh, divisions);
  }
  return SharedRoutingGrid<Dim>(mesh);
}

/**
 * The bounding box of this rank's part of `mesh`, or none when it holds no
 * vertices: a rank with no local coordinates owns no elements either, so there
 * is nothing for a query to find there.
 */
template <int Dim>
[[nodiscard]] std::optional<AABBox<Dim>> LocalBounds(Omega_h::Mesh& mesh)
{
  if (mesh.nverts() == 0) {
    return std::nullopt;
  }
  const auto bbox = Omega_h::find_bounding_box<Dim>(mesh.coords());
  AABBox<Dim> box{};
  for (int d = 0; d < Dim; ++d) {
    box.center[d] = (bbox.min[d] + bbox.max[d]) / 2.0;
    box.half_width[d] = (bbox.max[d] - bbox.min[d]) / 2.0;
  }
  return box;
}

/**
 * The candidate rank table of the parts of `mesh` on `grid` (the shared
 * routing grid by default). Collective.
 * \throws std::invalid_argument if the mesh dimension is not `Dim`
 */
template <int Dim>
[[nodiscard]] CandidateRankTable<Dim> CreateCandidateRankTable(
  Omega_h::Mesh& mesh, const UniformGrid<Dim>& grid)
{
  if (mesh.dim() != Dim) {
    throw std::invalid_argument(
      "CreateCandidateRankTable: mesh dimension does not match the table");
  }
  return CandidateRankTable<Dim>(mesh.comm()->get_impl(), grid,
                                 LocalBounds<Dim>(mesh));
}

/** \see CreateCandidateRankTable */
template <int Dim>
[[nodiscard]] CandidateRankTable<Dim> CreateCandidateRankTable(
  Omega_h::Mesh& mesh)
{
  return CreateCandidateRankTable<Dim>(mesh, SharedRoutingGrid<Dim>(mesh));
}

} // namespace pcms::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_CANDIDATE_RANK_TABLE_HPP
