// Distributed point search, layer by layer.
//
// 1. The routing table (`CandidateRankTable`): the global search grid is the
//    same on every rank, and a point inside an element a rank owns always
//    lists that rank as a candidate.
// 2. The local pass (`ClassifyLocally`): a rank answers only for elements it
//    owns, never for ghosts.
// 3. The full search (`Locate`): every point is attributed to its owning rank
//    and element, including points that had to be routed, and a warm cache
//    routes nothing.
#include "distributed_test_utils.h"

#include "pcms/localization/candidate_rank_table.h"
#include "pcms/localization/distributed_point_search.h"
#include "pcms/utility/uniform_grid.h"

#include <Omega_h_bbox.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

using distributed_test::require;
using distributed_test::require_same_on_all_ranks;
using distributed_test::require_same_real_on_all_ranks;
using pcms::LO;
using pcms::Real;

namespace
{

// ---- 1. routing table ----

/// Points of a simplex on a barycentric lattice of denominator `denom`.
///
/// Enumerating barycentric weights covers every vertex, every edge and face
/// interior point, and the interior, in one sweep, which is what the routing
/// guarantee has to hold for.
template <int Dim>
std::vector<std::array<Real, Dim>> SimplexSamples(
  const std::array<std::array<Real, Dim>, Dim + 1>& vertices, int denom)
{
  std::vector<std::array<Real, Dim>> samples;
  std::array<int, Dim + 1> weight{};
  while (true) {
    int head_sum = 0;
    for (int v = 0; v < Dim; ++v) {
      head_sum += weight[v];
    }
    if (head_sum <= denom) {
      weight[Dim] = denom - head_sum;
      std::array<Real, Dim> point{};
      for (int v = 0; v <= Dim; ++v) {
        for (int d = 0; d < Dim; ++d) {
          point[d] += static_cast<Real>(weight[v]) * vertices[v][d];
        }
      }
      for (int d = 0; d < Dim; ++d) {
        point[d] /= static_cast<Real>(denom);
      }
      samples.push_back(point);
    }
    // Odometer over the first Dim weights; the last is determined by them.
    int axis = 0;
    for (; axis < Dim; ++axis) {
      if (++weight[axis] <= denom) {
        break;
      }
      weight[axis] = 0;
    }
    if (axis == Dim) {
      break;
    }
  }
  return samples;
}

/// The routing grid must be the shared one, or a cell id means something
/// different on each rank and every lookup below is meaningless.
///
/// Its resolution therefore has to come from the global element count. The
/// local count, `mesh.nelems()`, includes ghost copies and differs per rank,
/// so this only means something on a mesh that really has ghosts.
template <int Dim>
void check_shared_grid(Omega_h::Mesh& mesh,
                       const pcms::CandidateRankTable<Dim>& table)
{
  auto comm = mesh.comm();
  const auto& grid = table.GetGrid();
  const auto expected = pcms::CandidateRankTable<Dim>::SharedRoutingGrid(mesh);

  const auto global_count = pcms::GlobalElementCount(mesh);
  require(comm->allreduce(mesh.nelems(), Omega_h_Op::OMEGA_H_SUM) >
            global_count,
          "no rank holds ghost elements, so the shared-grid check is vacuous");
  require_same_on_all_ranks(comm, global_count,
                            "global element count differs between ranks");

  for (int d = 0; d < Dim; ++d) {
    require_same_on_all_ranks(comm, grid.divisions[d],
                              "grid divisions differ between ranks");
    require_same_real_on_all_ranks(comm, grid.bot_left[d], 1e-12,
                                   "grid bottom-left differs between ranks");
    require_same_real_on_all_ranks(comm, grid.edge_length[d], 1e-12,
                                   "grid extent differs between ranks");
    // ...and the table is on that grid rather than one of its own making.
    require(grid.divisions[d] == expected.divisions[d],
            "table did not use the shared routing grid");
    require(grid.bot_left[d] == expected.bot_left[d],
            "table did not use the shared routing grid");
    require(grid.edge_length[d] == expected.edge_length[d],
            "table did not use the shared routing grid");
  }

  if (comm->rank() == 0) {
    std::printf("routing grid: %lld x %lld cells over [%.3f, %.3f]\n",
                static_cast<long long>(grid.divisions[0]),
                static_cast<long long>(grid.divisions[Dim - 1]),
                grid.bot_left[0], grid.bot_left[0] + grid.edge_length[0]);
  }
}

/// Helper: this rank's element vertices, and which of them it owns.
struct ElementGeom
{
  std::vector<std::array<std::array<Real, 3>, 4>> vertices; // [elem][vert][dim]
  std::vector<char> owned;
};

ElementGeom GatherElementGeom(Omega_h::Mesh& mesh)
{
  const int dim = mesh.dim();
  const int nverts_per_elem = dim + 1;
  const auto n_elems = mesh.nelems();
  const auto elem_verts = Omega_h::HostRead<LO>(mesh.ask_elem_verts());
  const auto coords = Omega_h::HostRead<Real>(mesh.coords());
  const auto owned = Omega_h::HostRead<Omega_h::I8>(mesh.owned(dim));

  ElementGeom geom;
  geom.vertices.resize(static_cast<size_t>(n_elems));
  geom.owned.resize(static_cast<size_t>(n_elems));
  for (LO e = 0; e < n_elems; ++e) {
    geom.owned[static_cast<size_t>(e)] = owned[e] != 0;
    for (int v = 0; v < nverts_per_elem; ++v) {
      const LO vert = elem_verts[e * nverts_per_elem + v];
      for (int d = 0; d < 3; ++d) {
        geom.vertices[static_cast<size_t>(e)][static_cast<size_t>(v)]
                     [static_cast<size_t>(d)] =
          d < dim ? coords[vert * dim + d] : 0.0;
      }
    }
  }
  return geom;
}

/// The guarantee: a point inside an element this rank owns must route back to
/// this rank. Sampled over a barycentric lattice so vertices, edge and face
/// interiors, and the element interior are all covered.
template <int Dim>
void check_owned_points_route_to_owner(
  Omega_h::Mesh& mesh, const pcms::CandidateRankTable<Dim>& table, int denom)
{
  auto comm = mesh.comm();
  const int rank = comm->rank();
  const auto geom = GatherElementGeom(mesh);

  LO n_owned = 0;
  LO n_probes = 0;
  LO n_uncovered = 0;
  std::array<Real, 3> first_uncovered{};
  for (size_t e = 0; e < geom.vertices.size(); ++e) {
    if (!geom.owned[e]) {
      continue;
    }
    ++n_owned;
    std::array<std::array<Real, Dim>, Dim + 1> vertices{};
    for (int v = 0; v <= Dim; ++v) {
      for (int d = 0; d < Dim; ++d) {
        vertices[static_cast<size_t>(v)][static_cast<size_t>(d)] =
          geom.vertices[e][static_cast<size_t>(v)][static_cast<size_t>(d)];
      }
    }
    for (const auto& sample : SimplexSamples<Dim>(vertices, denom)) {
      ++n_probes;
      Omega_h::Vector<Dim> point;
      for (int d = 0; d < Dim; ++d) {
        point[d] = sample[static_cast<size_t>(d)];
      }
      const auto candidates = table.CandidatesForPoint(point);
      if (std::find(candidates.begin(), candidates.end(), rank) ==
          candidates.end()) {
        if (n_uncovered == 0) {
          for (int d = 0; d < 3; ++d) {
            first_uncovered[static_cast<size_t>(d)] =
              sample[static_cast<size_t>(d < Dim ? d : 0)];
          }
        }
        ++n_uncovered;
      }
    }
  }

  const LO total_probes = comm->allreduce(n_probes, Omega_h_Op::OMEGA_H_SUM);
  const LO total_uncovered =
    comm->allreduce(n_uncovered, Omega_h_Op::OMEGA_H_SUM);
  const LO total_owned = comm->allreduce(n_owned, Omega_h_Op::OMEGA_H_SUM);

  if (n_uncovered > 0) {
    std::fprintf(stderr, "[rank %d] first uncovered probe (%.17g, %.17g)\n",
                 rank, first_uncovered[0], first_uncovered[1]);
  }
  if (rank == 0) {
    std::printf(
      "%dD: %lld owned elements, %lld probes, %lld not routed to their owner\n",
      Dim, static_cast<long long>(total_owned),
      static_cast<long long>(total_probes),
      static_cast<long long>(total_uncovered));
  }
  require(total_probes > 0, "no probes were generated");
  require(total_uncovered == 0,
          "a point inside an owned element did not route to its owner");
}

/// The routing decision must be a function of the point alone: every rank holds
/// the same table, so every rank must produce the same candidate list for the
/// same point. A gather that mismatched ranks would show up here.
template <int Dim>
void check_candidate_lists_agree_across_ranks(
  Omega_h::Mesh& mesh, const pcms::CandidateRankTable<Dim>& table)
{
  auto comm = mesh.comm();
  const int nranks = table.NumRanks();
  const auto& grid = table.GetGrid();

  // Cell centers are a deterministic, rank-independent point set because the
  // grid is shared.
  const LO num_cells = grid.GetNumCells();
  const LO stride = std::max<LO>(1, num_cells / 64);
  std::vector<LO> cells;
  for (LO cell = 0; cell < num_cells; cell += stride) {
    cells.push_back(cell);
  }

  const size_t width = cells.size() * static_cast<size_t>(nranks);
  std::vector<int> local(width, -1);
  for (size_t i = 0; i < cells.size(); ++i) {
    const auto bbox = grid.GetCellBBOX(cells[i]);
    Omega_h::Vector<Dim> point;
    for (int d = 0; d < Dim; ++d) {
      point[d] = bbox.center[static_cast<size_t>(d)];
    }
    const auto candidates = table.CandidatesForPoint(point);
    for (size_t j = 0; j < candidates.size(); ++j) {
      local[i * static_cast<size_t>(nranks) + j] = candidates[j];
    }
  }

  std::vector<int> all(width * static_cast<size_t>(nranks), 0);
  MPI_Allgather(local.data(), static_cast<int>(width), MPI_INT, all.data(),
                static_cast<int>(width), MPI_INT, comm->get_impl());
  for (int r = 1; r < nranks; ++r) {
    for (size_t k = 0; k < width; ++k) {
      if (all[k] != all[static_cast<size_t>(r) * width + k]) {
        std::fprintf(stderr,
                     "[rank %d] rank %d computed a different candidate list "
                     "for sampled cell %zu\n",
                     comm->rank(), r, k / static_cast<size_t>(nranks));
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
  }

  if (comm->rank() == 0) {
    std::printf("%dD: %zu sampled cells give identical candidate lists on all "
                "%d ranks\n",
                Dim, cells.size(), nranks);
  }
}

/// The published range must be exactly what the local bounding box implies.
///
/// Recomputed here from first principles -- box corners binned by the grid's
/// own per-axis rule, with one cell of padding -- rather than by calling the
/// table's helper, so this catches a wrong pad, a wrong axis order, or a range
/// built from the wrong box. It also validates the gather: the range a rank
/// sees for rank `r` must equal what rank `r` computed for itself, which is
/// exactly what a neighbour-scoped gather silently gets wrong by returning a
/// short buffer.
template <int Dim>
void check_published_ranges(Omega_h::Mesh& mesh,
                            const pcms::CandidateRankTable<Dim>& table,
                            Real tolerance_reach)
{
  auto comm = mesh.comm();
  const auto& grid = table.GetGrid();

  typename pcms::CandidateRankTable<Dim>::CellRange expected{};
  if (mesh.nverts() == 0) {
    expected = pcms::CandidateRankTable<Dim>::CellRange::Empty();
  } else {
    const auto bbox = Omega_h::find_bounding_box<Dim>(mesh.coords());
    for (int d = 0; d < Dim; ++d) {
      const Real lo_coord = bbox.min[d] - tolerance_reach;
      const Real hi_coord = bbox.max[d] + tolerance_reach;
      expected.lo[d] = std::max<LO>(grid.AxisCellIndex(d, lo_coord) - 1, 0);
      expected.hi[d] = std::min<LO>(grid.AxisCellIndex(d, hi_coord) + 1,
                                    grid.divisions[d] - 1);
    }
  }

  for (int d = 0; d < Dim; ++d) {
    require(table.LocalRange().lo[d] == expected.lo[d],
            "published range's lower bound does not match the local bbox");
    require(table.LocalRange().hi[d] == expected.hi[d],
            "published range's upper bound does not match the local bbox");
  }

  // Independent expectations, gathered and compared against the table: this is
  // the check that a rank-indexed gather went to all ranks rather than only to
  // the mesh's communication-graph neighbours.
  Kokkos::Array<LO, 2 * static_cast<size_t>(Dim)> local{};
  for (int d = 0; d < Dim; ++d) {
    local[static_cast<size_t>(d)] = expected.lo[d];
    local[static_cast<size_t>(Dim) + d] = expected.hi[d];
  }
  const auto all =
    pcms::AllGatherAllRanks(comm->get_impl(), local.data(), 2 * Dim);
  for (int r = 0; r < table.NumRanks(); ++r) {
    for (int d = 0; d < Dim; ++d) {
      require(table.RangeForRank(r).lo[d] ==
                all[static_cast<size_t>(r) * 2 * Dim + d],
              "rank r's gathered range does not match rank r's own range");
      require(table.RangeForRank(r).hi[d] ==
                all[static_cast<size_t>(r) * 2 * Dim + Dim + d],
              "rank r's gathered range does not match rank r's own range");
    }
  }

  if (comm->rank() == 0) {
    std::printf("%dD: published ranges verified against each rank's own bbox "
                "and against a rank-indexed gather\n",
                Dim);
  }
}

/// Widening the published ranges by the containment tolerance may only ever add
/// candidate ranks, and must actually add them once the tolerance exceeds a
/// cell size. Also checks the case the widening exists for: a point just
/// outside an element, which the search may still accept, must route to that
/// element's owner.
template <int Dim>
void check_tolerance_widening(Omega_h::Mesh& mesh)
{
  auto comm = mesh.comm();
  const int rank = comm->rank();

  // Use a grid fine enough that a rank's bare range does not already cover
  // every cell. On a coarse grid with few ranks both ranges can span the whole
  // domain, and then widening cannot add a candidate and the property is
  // unobservable.
  const auto shared = pcms::CandidateRankTable<Dim>::SharedRoutingGrid(mesh);
  const LO fine_divisions = 8 * shared.divisions[0];
  const auto fine_grid =
    pcms::CreateUniformGridFromMesh<Dim>(mesh, fine_divisions);
  const Real cell_size = fine_grid.edge_length[0] / fine_grid.divisions[0];

  // Widen by several cells so the effect is unambiguous, and probe a face at
  // half that distance: more than the single cell of padding a bare range
  // already has, so a bare range must miss the probe while a widened one
  // catches it.
  const Real reach = 4.0 * cell_size;
  const pcms::CandidateRankTable<Dim> strict(mesh, fine_grid, 0.0);
  const pcms::CandidateRankTable<Dim> wide(mesh, fine_grid, reach);

  LO pairs_strict = 0;
  LO pairs_wide = 0;
  LO cells_grown = 0;
  LO subset_violations = 0;
  for (LO cell = 0; cell < strict.GetGrid().GetNumCells(); ++cell) {
    const auto s = strict.CandidatesForCell(cell);
    const auto w = wide.CandidatesForCell(cell);
    pairs_strict += static_cast<LO>(s.size());
    pairs_wide += static_cast<LO>(w.size());
    if (w.size() > s.size()) {
      ++cells_grown;
    }
    // Both lists are ascending, so this is a subset test.
    if (!std::includes(w.begin(), w.end(), s.begin(), s.end())) {
      ++subset_violations;
    }
  }

  // The case the widening exists for: a point the search may still accept
  // because it lies just outside this rank's geometry. Probe just beyond each
  // face of the rank's own bounding box, at half the tolerance, so it is inside
  // the tolerance-widened box and must route back here. Probbing relative to
  // the box rather than to an element makes this deterministic: an interior
  // element's face can sit far inside the box, where the bare range already
  // covers the probe and no widening is needed.
  LO n_face_probes = 0;
  LO n_missed_wide = 0;
  LO n_missed_strict = 0;
  if (mesh.nverts() > 0) {
    const auto bbox = Omega_h::find_bounding_box<Dim>(mesh.coords());
    for (int d = 0; d < Dim; ++d) {
      for (int side = 0; side < 2; ++side) {
        Omega_h::Vector<Dim> probe;
        for (int e = 0; e < Dim; ++e) {
          probe[e] = (bbox.min[e] + bbox.max[e]) / 2.0;
        }
        // Half a tolerance outside a face of the box: within `reach` of the
        // rank's geometry, so the widened range must cover it.
        probe[d] =
          side == 0 ? bbox.min[d] - reach / 2.0 : bbox.max[d] + reach / 2.0;
        ++n_face_probes;
        const auto wide_candidates = wide.CandidatesForPoint(probe);
        if (std::find(wide_candidates.begin(), wide_candidates.end(), rank) ==
            wide_candidates.end()) {
          ++n_missed_wide;
        }
        const auto strict_candidates = strict.CandidatesForPoint(probe);
        if (std::find(strict_candidates.begin(), strict_candidates.end(),
                      rank) == strict_candidates.end()) {
          ++n_missed_strict;
        }
      }
    }
  }

  const LO total_pairs_strict =
    comm->allreduce(pairs_strict, Omega_h_Op::OMEGA_H_SUM);
  const LO total_pairs_wide =
    comm->allreduce(pairs_wide, Omega_h_Op::OMEGA_H_SUM);
  const LO total_grown = comm->allreduce(cells_grown, Omega_h_Op::OMEGA_H_SUM);
  const LO total_violations =
    comm->allreduce(subset_violations, Omega_h_Op::OMEGA_H_SUM);
  const LO total_probes =
    comm->allreduce(n_face_probes, Omega_h_Op::OMEGA_H_SUM);
  const LO total_missed_wide =
    comm->allreduce(n_missed_wide, Omega_h_Op::OMEGA_H_SUM);
  const LO total_missed_strict =
    comm->allreduce(n_missed_strict, Omega_h_Op::OMEGA_H_SUM);

  if (rank == 0) {
    std::printf(
      "%dD tolerance widening (reach %.4g = %.1f cells): %lld -> %lld "
      "candidate "
      "pairs, %lld grown cells; %lld probes outside a face, missed by widened "
      "ranges %lld, by bare ranges %lld\n",
      Dim, static_cast<double>(reach), static_cast<double>(reach / cell_size),
      static_cast<long long>(total_pairs_strict),
      static_cast<long long>(total_pairs_wide),
      static_cast<long long>(total_grown), static_cast<long long>(total_probes),
      static_cast<long long>(total_missed_wide),
      static_cast<long long>(total_missed_strict));
  }

  require(total_violations == 0, "widening removed a candidate rank");
  // The widened range must contain the bare range on every axis, and the
  // widened candidate list must contain the bare one. Both hold whether or not
  // the partition happens to saturate the grid.
  for (int d = 0; d < Dim; ++d) {
    require(wide.LocalRange().lo[d] <= strict.LocalRange().lo[d],
            "widening moved a range's lower bound outward");
    require(wide.LocalRange().hi[d] >= strict.LocalRange().hi[d],
            "widening shrank a range's upper bound");
  }
  require(total_missed_wide == 0,
          "a point just outside an owned element was not routed to its owner");

  /// Widening is unobservable exactly when every non-empty rank's range already
  /// spans the whole grid: there is no cell left to add, so the candidate lists
  /// cannot grow no matter how large the tolerance. A range can span the grid
  /// because the partition interleaves (each rank's coordinate extent covers
  /// the domain) or because the grid is coarse relative to the decomposition.
  bool any_rank_can_grow = false;
  for (int r = 0; r < strict.NumRanks(); ++r) {
    const auto& range = strict.RangeForRank(r);
    if (range.IsEmpty()) {
      continue;
    }
    bool spans = true;
    for (int d = 0; d < Dim; ++d) {
      if (range.lo[d] != 0 || range.hi[d] != fine_grid.divisions[d] - 1) {
        spans = false;
      }
    }
    if (!spans) {
      any_rank_can_grow = true;
    }
  }
  if (rank == 0) {
    std::printf("%dD widening is %s on this partition\n", Dim,
                any_rank_can_grow ? "observable" : "unobservable");
  }
  if (any_rank_can_grow) {
    require(total_pairs_wide > total_pairs_strict,
            "widening by several cells did not add any candidate");
    require(total_grown > 0, "widening did not grow any cell's candidate list");
  }
}

/// A rank with no local geometry can find nothing, so it must publish an empty
/// range and never be a routing candidate.
template <int Dim>
void check_empty_ranks(Omega_h::Mesh& mesh,
                       const pcms::CandidateRankTable<Dim>& table)
{
  auto comm = mesh.comm();
  const int nranks = table.NumRanks();
  const auto nverts_all =
    pcms::AllGatherAllRanks(comm->get_impl(), static_cast<LO>(mesh.nverts()));

  LO n_empty = 0;
  for (int r = 0; r < nranks; ++r) {
    if (nverts_all[r] != 0) {
      continue;
    }
    ++n_empty;
    require(table.RangeForRank(r).IsEmpty(),
            "a rank with no geometry published a non-empty range");
    for (LO cell = 0; cell < table.GetGrid().GetNumCells(); ++cell) {
      const auto candidates = table.CandidatesForCell(cell);
      require(std::find(candidates.begin(), candidates.end(), r) ==
                candidates.end(),
              "an empty rank is a routing candidate");
    }
  }

  // `n_empty` is already global: every rank holds the same gathered view, so it
  // must not be reduced again.
  if (comm->rank() == 0) {
    std::printf("%dD: %lld rank(s) with no local geometry\n", Dim,
                static_cast<long long>(n_empty));
  }
}

// ---- 2. local pass, 3. full search ----

/// A centroid lies strictly inside its own element, so the geometrically right
/// answer is that element. That makes ownership of the answer unambiguous: an
/// owned element must be answerable here, a ghost must not be.
template <int Dim>
void check_local_classification(Omega_h::Mesh& mesh)
{
  using Status = pcms::DistributedPointSearch<Dim>::LocalStatus;

  const auto n_elems = mesh.nelems();
  const auto n_owned = pcms::OwnedElementCount(mesh);
  require(mesh.nghost_layers() > 0, "mesh is not ghosted");

  const auto divisions = pcms::DivisionsPerAxisForGlobalMesh<Dim>(mesh);
  const auto local_search = distributed_test::MakeLocalSearch<Dim>(mesh,
                                                                    divisions);
  const pcms::DistributedPointSearch<Dim> search{mesh, local_search};
  require(search.HasGhostLayer(), "distributed search missed the ghost layer");

  const auto centroids =
    Omega_h::HostRead<Real>(pcms::get_entity_centroids(mesh, Dim));
  Kokkos::View<pcms::Real* [Dim]> points("centroids", n_elems);
  auto points_h = Kokkos::create_mirror_view(points);
  for (LO i = 0; i < n_elems; ++i) {
    for (int d = 0; d < Dim; ++d) {
      points_h(i, d) = centroids[i * Dim + d];
    }
  }
  Kokkos::deep_copy(points, points_h);

  const auto classified = search.ClassifyLocally(points);
  require(classified.size() == static_cast<size_t>(n_elems),
          "one classification per point expected");

  LO n_here = 0;
  LO n_elsewhere = 0;
  LO n_missing = 0;
  LO n_owned_wrong = 0;
  LO n_ghost_wrong = 0;
  for (LO i = 0; i < n_elems; ++i) {
    const auto& c = classified[static_cast<size_t>(i)];
    const bool owned = search.OwnsLocalElement(i);
    switch (c.status) {
      case Status::LOCATED_HERE:
        ++n_here;
        // Only an owned element may be answered here, and it must be the
        // element the point came from.
        if (!owned || c.element_id != i) {
          ++n_owned_wrong;
        }
        break;
      case Status::LOCATED_ELSEWHERE:
        ++n_elsewhere;
        // A ghost answer must name a non-owned element, and element_id must be
        // a valid local element.
        if (owned || c.element_id < 0 || c.element_id >= n_elems) {
          ++n_ghost_wrong;
        }
        break;
      case Status::NOT_FOUND_HERE: ++n_missing; break;
    }
  }

  const LO total_here = mesh.comm()->allreduce(n_here, Omega_h_Op::OMEGA_H_SUM);
  const LO total_elsewhere =
    mesh.comm()->allreduce(n_elsewhere, Omega_h_Op::OMEGA_H_SUM);
  const LO total_missing =
    mesh.comm()->allreduce(n_missing, Omega_h_Op::OMEGA_H_SUM);
  const LO total_owned_wrong =
    mesh.comm()->allreduce(n_owned_wrong, Omega_h_Op::OMEGA_H_SUM);
  const LO total_ghost_wrong =
    mesh.comm()->allreduce(n_ghost_wrong, Omega_h_Op::OMEGA_H_SUM);
  const LO total_owned =
    mesh.comm()->allreduce(n_owned, Omega_h_Op::OMEGA_H_SUM);
  // A rank whose partition is isolated holds no ghosts and simply contributes
  // nothing here; what matters is that the decomposition as a whole exercised
  // the ghost classification. That is checked below via `total_elsewhere`.

  if (mesh.comm()->rank() == 0) {
    std::printf(
      "classification of %lld owned element centroids: %lld answerable here, "
      "%lld answerable elsewhere (ghost), %lld not found, %lld/%lld errors\n",
      static_cast<long long>(total_owned), static_cast<long long>(total_here),
      static_cast<long long>(total_elsewhere),
      static_cast<long long>(total_missing),
      static_cast<long long>(total_owned_wrong) + total_ghost_wrong,
      static_cast<long long>(total_here) + total_elsewhere +
        total_missing);
  }

  require(
    total_owned_wrong == 0,
    "an element not owned by this rank was classified as answerable here");
  require(total_ghost_wrong == 0,
          "a ghost answer did not name a non-owned element");
  require(total_here == total_owned,
          "not every owned element was answerable on its owner");
  require(
    total_missing == 0,
    "a local element centroid was not found at all by the rank holding it");
  require(total_elsewhere > 0,
          "no point was classified as answerable elsewhere, so the ownership "
          "distinction is not being exercised");
}

/// Every rank locates every owned element centroid in the decomposition.
///
/// Each point is therefore answered by the rank that owns it, locally, and by
/// every other rank, through routing. Both paths must agree on the owning rank
/// and on the element's global id -- that is the cross-path consistency the
/// merge rule exists to guarantee.
template <int Dim>
void check_distributed_locate(Omega_h::Mesh& mesh)
{
  using ResultStatus = pcms::DistributedPointSearch<Dim>::ResultStatus;

  const int rank = mesh.comm()->rank();
  require(mesh.has_tag(Dim, "global"),
          "test needs element global ids to check which element was found");

  // Every owned element centroid in the decomposition, plus one point far
  // outside the domain on every axis.
  auto gathered = distributed_test::GatherOwnedCentroids(mesh);
  const int total = gathered.total;
  const auto& global_ids = gathered.global_id;
  const auto& expected_owner = gathered.contributor;
  const LO n_local = static_cast<LO>(
    std::count(expected_owner.begin(), expected_owner.end(), rank));
  for (int d = 0; d < Dim; ++d) {
    gathered.coords.push_back(2.0);
  }
  const LO n_points = static_cast<LO>(total) + 1;
  const auto points = distributed_test::ToPointView<Dim>(
    gathered.coords, "global element centroids");

  const auto divisions = pcms::DivisionsPerAxisForGlobalMesh<Dim>(mesh);
  const auto local_search = distributed_test::MakeLocalSearch<Dim>(mesh,
                                                                    divisions);
  pcms::DistributedPointSearch<Dim> search{mesh, local_search};

  const auto results = search.Locate(points);
  require(results.size() == static_cast<size_t>(n_points),
          "one result per point expected");

  LO n_wrong_owner = 0;
  LO n_wrong_element = 0;
  LO n_not_found = 0;
  LO n_bad_bary = 0;
  for (int i = 0; i < total; ++i) {
    const auto& r = results[static_cast<size_t>(i)];
    if (r.status != ResultStatus::LOCATED) {
      ++n_not_found;
      continue;
    }
    if (r.owner_rank != expected_owner[static_cast<size_t>(i)]) {
      ++n_wrong_owner;
    }
    if (r.global_element_id != global_ids[static_cast<size_t>(i)]) {
      ++n_wrong_element;
    }
    // A centroid is strictly inside its element, so every barycentric weight is
    // positive, and the weights sum to one. This also catches a result
    // assembled from the wrong slots.
    Real weight = 0;
    for (int d = 0; d <= Dim; ++d) {
      weight += r.barycentric[d];
      if (r.barycentric[d] <= 0) {
        ++n_bad_bary;
      }
    }
    if (std::fabs(weight - 1.0) > 1e-8) {
      ++n_bad_bary;
    }
  }
  // A point no rank placed must not name an owner or an element, so a caller
  // cannot mistake it for a location.
  const auto& outside = results[static_cast<size_t>(total)];
  const bool outside_not_found = outside.status == ResultStatus::NOT_FOUND &&
                                 outside.owner_rank == -1 &&
                                 outside.element_id == -1 &&
                                 outside.global_element_id == -1;

  const LO t_not_found =
    mesh.comm()->allreduce(n_not_found, Omega_h_Op::OMEGA_H_SUM);
  const LO t_wrong_owner =
    mesh.comm()->allreduce(n_wrong_owner, Omega_h_Op::OMEGA_H_SUM);
  const LO t_wrong_element =
    mesh.comm()->allreduce(n_wrong_element, Omega_h_Op::OMEGA_H_SUM);
  const LO t_bad_bary =
    mesh.comm()->allreduce(n_bad_bary, Omega_h_Op::OMEGA_H_SUM);
  const LO t_points = mesh.comm()->allreduce(n_points, Omega_h_Op::OMEGA_H_SUM);
  const LO t_routed_away = mesh.comm()->allreduce(
    static_cast<LO>(total) - n_local,
    Omega_h_Op::OMEGA_H_SUM);

  if (rank == 0) {
    std::printf(
      "located %lld gathered points (%lld queries per rank left their rank): "
      "%lld not found, %lld wrong owner, %lld wrong element, %lld bad "
      "barycentric\n",
      static_cast<long long>(t_points), static_cast<long long>(t_routed_away),
      static_cast<long long>(t_not_found),
      static_cast<long long>(t_wrong_owner),
      static_cast<long long>(t_wrong_element),
      static_cast<long long>(t_bad_bary));
  }

  require(outside_not_found,
          "a point outside the domain was reported as located, or its result "
          "named an owner or element");
  require(t_not_found == 0, "a point inside the domain was not located");
  require(t_wrong_owner == 0,
          "a point was attributed to the wrong owning rank");
  require(t_wrong_element == 0, "a point resolved to the wrong element");
  require(t_bad_bary == 0,
          "a located point has non-positive barycentric weight");
  require(t_routed_away > 0,
          "no query left its rank, so routing was not exercised");

  // Determinism and caching: the same batch twice must give the same answer,
  // and a second pass over a warm cache must not re-search anything it
  // resolved.
  std::vector<typename pcms::DistributedPointSearch<Dim>::Result> cache;
  const auto pass_one = search.Locate(points, &cache);
  const auto routed_one = search.LastRoutedPointCount();
  const auto pass_two = search.Locate(points, &cache);
  const auto routed_two = search.LastRoutedPointCount();

  LO n_unresolved_results = 0;
  for (const auto& r : pass_one) {
    if (r.status != ResultStatus::LOCATED) {
      ++n_unresolved_results;
    }
  }

  LO n_unstable = 0;
  for (size_t i = 0; i < results.size(); ++i) {
    if (pass_two[i].status != results[i].status ||
        pass_two[i].owner_rank != results[i].owner_rank ||
        pass_two[i].global_element_id != results[i].global_element_id) {
      ++n_unstable;
    }
  }

  if (rank == 0) {
    std::printf("routed %zu points on the first pass, %zu on the second "
                "(only unresolved points are retried, and points outside the "
                "mesh bounding box are never routed)\n",
                routed_one, routed_two);
  }
  require(routed_one > 0, "the first pass routed nothing");
  // The only point still unresolved on a warm cache is the one outside the
  // grid, and that one needs no message, so a warm pass routes nothing at all.
  require(n_unresolved_results > 0, "the out-of-domain probe was located");
  require(routed_two == 0, "a cached pass still routed points");
  require(routed_two < routed_one, "caching did not reduce the routing work");
  require(mesh.comm()->allreduce(n_unstable, Omega_h_Op::OMEGA_H_SUM) == 0,
          "a cached pass gave a different answer");

  // A batch of points outside the global bounding box needs no communication at
  // all: no element can contain them, so not a single query is routed.
  // One past each end of the first axis, then one past the far corner.
  Kokkos::View<Real* [Dim]> far_points("far outside", 3);
  auto far_h = Kokkos::create_mirror_view(far_points);
  for (int d = 0; d < Dim; ++d) {
    far_h(0, d) = (d == 0) ? 2.0 : 0.5;
    far_h(1, d) = (d == 0) ? -1.0 : 0.5;
    far_h(2, d) = 0.5;
  }
  far_h(2, Dim - 1) = 5.0;
  Kokkos::deep_copy(far_points, far_h);
  const auto far_results = search.Locate(far_points);
  const std::size_t far_routed = search.LastRoutedPointCount();
  LO n_far_located = 0;
  for (const auto& r : far_results) {
    if (r.status != ResultStatus::NOT_FOUND) {
      ++n_far_located;
    }
  }
  require(far_routed == 0, "a point outside the mesh bounding box was routed");
  require(mesh.comm()->allreduce(n_far_located, Omega_h_Op::OMEGA_H_SUM) == 0,
          "a point outside the mesh bounding box was reported as located");
}

/// The global search grid's resolution can be chosen by the caller. Any
/// resolution must give the same answers -- it only changes how many candidate
/// ranks a point is sent to -- and a resolution that differs between ranks must
/// be rejected on every rank, since cell ids would then mean different things.
template <int Dim>
void check_grid_divisions(Omega_h::Mesh& mesh)
{
  using Search = pcms::DistributedPointSearch<Dim>;
  const auto gathered = distributed_test::GatherOwnedCentroids(mesh);
  const auto points = distributed_test::ToPointView<Dim>(gathered.coords);
  const auto local_search = distributed_test::MakeLocalSearch<Dim>(
    mesh, pcms::DivisionsPerAxisForGlobalMesh<Dim>(mesh));

  Search by_default{mesh, local_search};
  const auto expected = by_default.Locate(points);

  LO n_differ = 0;
  for (const LO divisions : {LO{3}, LO{37}}) {
    Search custom{mesh, local_search, 0.0, divisions};
    for (int d = 0; d < Dim; ++d) {
      require(custom.GetRoutingTable().GetGrid().divisions[d] == divisions,
              "the search did not use the requested grid resolution");
    }
    const auto got = custom.Locate(points);
    for (size_t i = 0; i < got.size(); ++i) {
      if (got[i].status != expected[i].status ||
          got[i].owner_rank != expected[i].owner_rank ||
          got[i].global_element_id != expected[i].global_element_id) {
        ++n_differ;
      }
    }
  }
  require(mesh.comm()->allreduce(n_differ, Omega_h_Op::OMEGA_H_SUM) == 0,
          "a custom grid resolution changed where a point was located");

  bool threw = false;
  try {
    Search mismatched{mesh, local_search, 0.0, 10 + mesh.comm()->rank()};
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  require(threw, "a grid resolution that differs between ranks was accepted");

  if (mesh.comm()->rank() == 0) {
    std::printf("%dD: custom grid resolutions give the same answers; a "
                "mismatched one is rejected\n",
                Dim);
  }
}

template <int Dim>
void run_all(Omega_h::CommPtr world)
{
  auto mesh = distributed_test::BuildGhostedBox<Dim>(world);
  if (world->rank() == 0) {
    std::printf("--- %dD, %d ranks ---\n", Dim, world->size());
  }

  const pcms::CandidateRankTable<Dim> table(mesh);
  check_shared_grid(mesh, table);
  check_published_ranges<Dim>(mesh, table, 0.0);
  check_owned_points_route_to_owner(mesh, table, /*denom=*/Dim == 2 ? 4 : 3);
  check_candidate_lists_agree_across_ranks(mesh, table);
  check_tolerance_widening<Dim>(mesh);
  check_empty_ranks(mesh, table);

  check_local_classification<Dim>(mesh);
  check_distributed_locate<Dim>(mesh);
  check_grid_divisions<Dim>(mesh);
}

} // namespace

int main(int argc, char** argv)
{
  try {
    MPI_Init(&argc, &argv);
    {
      Omega_h::Library library(&argc, &argv);
      auto world = library.world();
      require(world->size() >= 2, "test requires at least 2 MPI ranks");

      run_all<2>(world);
      run_all<3>(world);

      if (world->rank() == 0) {
        std::printf("distributed search test passed\n");
      }
    }
    MPI_Finalize();
    return 0;
  } catch (const std::exception& e) {
    distributed_test::AbortAllRanks(e.what());
  } catch (...) {
    distributed_test::AbortAllRanks("unknown exception");
  }
}
