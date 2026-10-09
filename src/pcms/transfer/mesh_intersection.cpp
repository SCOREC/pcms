#include "pcms/transfer/mesh_intersection.hpp"
#include "pcms/field/evaluator/omega_h_lagrange.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/mesh_geometry.h"
#include "pcms/utility/omega_h_array_utils.h"

namespace pcms
{
namespace
{
// Construct the source-mesh containing-element search appropriate to the
// spatial dimension: a background grid sized so cells hold ~1 element.
template <int Dim>
auto MakeGridPointSearch(Omega_h::Mesh& source_mesh)
{
  const auto n = pcms::DivisionsPerAxisForMesh<Dim>(source_mesh.nelems());
  if constexpr (Dim == 3) {
    return pcms::GridPointSearch3D(source_mesh, n, n, n);
  } else {
    return pcms::GridPointSearch2D(source_mesh, n, n);
  }
}

// For each target element, the source element containing its centroid: where
// the BFS over the source dual graph starts.
template <int Dim>
Kokkos::View<const LO*> LocateTargetCentroids(
  Omega_h::Mesh& target_mesh, const PointLocalizationSearch<Dim>& source_search)
{
  const auto flat_centroids = pcms::get_entity_centroids(target_mesh, Dim);
  // Convert layout_right 1D Omega_h array to 2D Kokkos view with correct layout
  auto centroids = ConvertCoordsTo2D(flat_centroids, target_mesh.nelems(), Dim);
  auto results = source_search(centroids);
  return source_search.GetOwningElementIds(results);
}
} // namespace

template <int Dim>
void FindIntersections::adjBasedIntersectSearch(
  const Kokkos::View<const LO*>& start_elements,
  const Omega_h::LOs& tgt2src_offsets,
  Omega_h::Write<Omega_h::LO>& nIntersections,
  Omega_h::Write<Omega_h::LO>& tgt2src_indices, bool is_count_only,
  bool use_prefilter)
{
  // Element entity dimension equals the spatial dimension (FACE for 2D, REGION
  // for 3D); measures are triangle areas (2D) or tet volumes (3D).
  const auto& tgt_coords = target_mesh_.coords();
  const auto& src_coords = source_mesh_.coords();
  const auto& tgt_elems2nodes = target_mesh_.ask_down(Dim, Omega_h::VERT).ab2b;
  const auto& src_elems2nodes = source_mesh_.ask_down(Dim, Omega_h::VERT).ab2b;
  const auto& src_elem_measures = measure_elements_real(&source_mesh_);
  const auto& tgt_elem_measures = measure_elements_real(&target_mesh_);
  const auto& t2t =
    source_mesh_.ask_dual(); // gives connected element neighbors
  const auto& t2tt = t2t.a2ab;
  const auto& tt2t = t2t.ab2b;

  auto nelems_target = target_mesh_.nelems();
  PCMS_ALWAYS_ASSERT(static_cast<Omega_h::LO>(start_elements.extent(0)) ==
                     nelems_target);
  Omega_h::parallel_for(
    nelems_target,
    OMEGA_H_LAMBDA(const Omega_h::LO id) {
      Queue queue;
      Track visited;

      auto current_cell_id = start_elements(id);
      auto current_tgt_elm_measure = tgt_elem_measures[id];

      OMEGA_H_CHECK_PRINTF(current_cell_id >= 0,
                           "ERROR: source cell id not found for given target "
                           "centroid %d\n",
                           id);

      auto tgt_elm_vert_coords =
        get_vert_coords_of_elem<Dim>(tgt_coords, tgt_elems2nodes, id);

      Omega_h::LO start_counter;
      if (!is_count_only) {
        start_counter = tgt2src_offsets[id];
      }

      int count = 0;

      count++;
      visited.push_back(current_cell_id);
      queue.push_back(current_cell_id);

      if (!is_count_only) {
        int idx_count = count - 1;
        tgt2src_indices[start_counter + idx_count] = current_cell_id;
      }

      while (!queue.isEmpty()) {
        Omega_h::LO currentElm = queue.front();
        queue.pop_front();
        auto start = t2tt[currentElm];
        auto end = t2tt[currentElm + 1];

        for (int i = start; i < end; ++i) {
          auto neighborElmId = tt2t[i];

          if (visited.notVisited(neighborElmId)) {
            // If the visited buffer is full, skip this neighbor so the BFS
            // terminates. Without this, an unrecorded neighbor stays "not
            // visited" and is re-queued forever (infinite loop). Mirrors the
            // guard in adj_search.cpp. Raise PCMS_INTERSECTION_TRACK_SIZE if
            // this fires.
            if (!visited.push_back(neighborElmId)) {
              printf("ERROR: visited buffer full "
                     "(PCMS_INTERSECTION_TRACK_SIZE=%d) for "
                     "target %d; some intersections may be missed\n",
                     PCMS_INTERSECTION_TRACK_SIZE, id);
              continue;
            }
            auto elm_vert_coords = get_vert_coords_of_elem<Dim>(
              src_coords, src_elems2nodes, neighborElmId);

            // Most neighbors reached through the dual graph share a face with
            // an element that already overlaps and so contribute no volume.
            // The plane test settles those far more cheaply than a clip and
            // only ever rejects provably degenerate overlaps. It is not a pure
            // speed-up, though: on meshes with non-dyadic coordinates the clip
            // itself accepts many of these shared-face pairs as roundoff
            // slivers (volume ~1e-12 of the element, right at eps), so the
            // filtered map is strictly smaller and the more correct one.
            if (use_prefilter && simplices_have_degenerate_overlap<Dim>(
                                   tgt_elm_vert_coords, elm_vert_coords)) {
              continue;
            }

            r3d::Polytope<Dim> intersection;
            r3d::intersect_simplices(intersection, tgt_elm_vert_coords,
                                     elm_vert_coords);
            // Take the magnitude: r3d::measure is signed by the orientation of
            // the target simplex used to initialize the polytope, which can be
            // negative for tetrahedra. This mirrors the fabs applied in the
            // sub-simplex decomposition and mass assembly; without it a
            // negatively-oriented target element would reject all of its real
            // overlaps and break conservation.
            auto intersected_measure = Kokkos::fabs(r3d::measure(intersection));
            auto current_src_elm_measure = src_elem_measures[neighborElmId];
            auto scale =
              Kokkos::fmax(current_tgt_elm_measure, current_src_elm_measure);
            auto eps = Kokkos::fmax(PCMS_INTERSECTION_ABS_TOL,
                                    PCMS_INTERSECTION_REL_TOL * scale);
            // A valid intersection is a non-degenerate simplex-simplex overlap:
            // at least Dim+1 vertices (a polygon in 2D, a polyhedron in 3D).
            if (intersection.nverts >= Dim + 1 && intersected_measure >= eps) {
              count++;

              OMEGA_H_CHECK_PRINTF(
                count < PCMS_INTERSECTION_QUEUE_SIZE,
                "intersection count for target %d reached the cap %d; raise "
                "PCMS_INTERSECTION_QUEUE_SIZE",
                id, PCMS_INTERSECTION_QUEUE_SIZE);

              queue.push_back(neighborElmId);

              if (!is_count_only) {
                Omega_h::LO idx_count = count - 1;
                tgt2src_indices[start_counter + idx_count] = neighborElmId;

              } // end of tgt2src_indices check

            } // end of intersection with bbox check

          } // end of not visited check

        } // end of loop over adj elements to the current element

      } // end of while loop

      nIntersections[id] = count;
    }, // end of lambda
    "count the number of intersections for each target element");
}

// Explicit instantiations for the supported spatial dimensions.
template void FindIntersections::adjBasedIntersectSearch<2>(
  const Kokkos::View<const LO*>&, const Omega_h::LOs&,
  Omega_h::Write<Omega_h::LO>&, Omega_h::Write<Omega_h::LO>&, bool, bool);
template void FindIntersections::adjBasedIntersectSearch<3>(
  const Kokkos::View<const LO*>&, const Omega_h::LOs&,
  Omega_h::Write<Omega_h::LO>&, Omega_h::Write<Omega_h::LO>&, bool, bool);

namespace
{
template <int Dim>
IntersectionResults intersectTargetsImpl(
  Omega_h::Mesh& source_mesh, Omega_h::Mesh& target_mesh,
  const PointLocalizationSearch<Dim>& source_search, bool use_prefilter)
{
  FindIntersections intersect(source_mesh, target_mesh);

  auto nelems_target = target_mesh.nelems();

  const auto start_elements =
    LocateTargetCentroids<Dim>(target_mesh, source_search);

  Omega_h::Write<Omega_h::LO> nIntersections(
    nelems_target, 0, "number of intersections in each target element");

  Omega_h::Write<Omega_h::LO> tgt2src_indices;

  intersect.adjBasedIntersectSearch<Dim>(start_elements, Omega_h::LOs(),
                                         nIntersections, tgt2src_indices, true,
                                         use_prefilter);

  Kokkos::fence();
  auto tgt2src_offsets = Omega_h::offset_scan(Omega_h::Read(nIntersections),
                                              "offsets for intersections");
  auto ntotal_intersections = tgt2src_offsets.last();

  Kokkos::fence();

  tgt2src_indices = Omega_h::Write<Omega_h::LO>(
    ntotal_intersections, 0,
    "indices of the source elements that intersect the given target element");

  intersect.adjBasedIntersectSearch<Dim>(start_elements, tgt2src_offsets,
                                         nIntersections, tgt2src_indices, false,
                                         use_prefilter);
  return {.tgt2src_offsets = tgt2src_offsets,
          .tgt2src_indices = Omega_h::read(tgt2src_indices)};
}

template <int Dim>
const PointLocalizationSearch<Dim>& RequireSearchDimension(
  const GridPointSearchVariant& source_search)
{
  using SearchT =
    std::conditional_t<Dim == 3, GridPointSearch3D, GridPointSearch2D>;
  const auto* search = std::get_if<SearchT>(&source_search);
  if (search == nullptr) {
    throw pcms_error("intersectTargets: the supplied source search has a "
                     "different spatial dimension than the meshes");
  }
  return *search;
}
} // namespace

IntersectionResults intersectTargets(Omega_h::Mesh& source_mesh,
                                     Omega_h::Mesh& target_mesh,
                                     bool use_prefilter)
{
  OMEGA_H_CHECK(source_mesh.dim() == target_mesh.dim());
  if (source_mesh.dim() == 3) {
    const auto search = MakeGridPointSearch<3>(source_mesh);
    return intersectTargetsImpl<3>(source_mesh, target_mesh, search,
                                   use_prefilter);
  }
  const auto search = MakeGridPointSearch<2>(source_mesh);
  return intersectTargetsImpl<2>(source_mesh, target_mesh, search,
                                 use_prefilter);
}

IntersectionResults intersectTargets(
  Omega_h::Mesh& source_mesh, Omega_h::Mesh& target_mesh,
  const GridPointSearchVariant& source_search, bool use_prefilter)
{
  OMEGA_H_CHECK(source_mesh.dim() == target_mesh.dim());
  if (source_mesh.dim() == 3) {
    return intersectTargetsImpl<3>(source_mesh, target_mesh,
                                   RequireSearchDimension<3>(source_search),
                                   use_prefilter);
  }
  return intersectTargetsImpl<2>(source_mesh, target_mesh,
                                 RequireSearchDimension<2>(source_search),
                                 use_prefilter);
}
namespace
{
std::shared_ptr<const OmegaHDiscretization> RequireOmegaHDiscretization(
  const FunctionSpace& space, const char* role)
{
  auto discretization = std::dynamic_pointer_cast<const OmegaHDiscretization>(
    space.GetDiscretization());
  if (!discretization) {
    throw pcms_error(std::string("OmegaHMeshIntersection: the ") + role +
                     " space is not on an Omega_h discretization");
  }
  return discretization;
}
} // namespace

const GridPointSearchVariant* SourceSearchFromSpace(const FunctionSpace& space)
{
  const auto* lagrange = dynamic_cast<const LagrangeFunctionSpace*>(&space);
  if (lagrange == nullptr) {
    return nullptr;
  }
  const auto* factory =
    dynamic_cast<const OmegaHLagrangeEvaluatorFactory<Real>*>(
      lagrange->GetEvaluatorFactory().get());
  return factory ? &factory->GetSearch() : nullptr;
}

OmegaHMeshIntersection::OmegaHMeshIntersection(
  const FunctionSpace& source_space, const FunctionSpace& target_space,
  bool use_prefilter)
  : OmegaHMeshIntersection(RequireOmegaHDiscretization(source_space, "source"),
                           RequireOmegaHDiscretization(target_space, "target"),
                           SourceSearchFromSpace(source_space), use_prefilter)
{
}

OmegaHMeshIntersection::OmegaHMeshIntersection(
  std::shared_ptr<const OmegaHDiscretization> source,
  std::shared_ptr<const OmegaHDiscretization> target,
  const GridPointSearchVariant* source_search, bool use_prefilter)
  : source_(std::move(source)), target_(std::move(target))
{
  if (!source_ || !target_) {
    throw pcms_error("OmegaHMeshIntersection: discretizations must be set");
  }
  if (source_->GetMesh().dim() != target_->GetMesh().dim()) {
    throw pcms_error(
      "OmegaHMeshIntersection: source and target mesh dimensions differ");
  }
  results_ =
    source_search
      ? intersectTargets(source_->GetMesh(), target_->GetMesh(), *source_search,
                         use_prefilter)
      : intersectTargets(source_->GetMesh(), target_->GetMesh(), use_prefilter);
}

std::shared_ptr<const Discretization>
OmegaHMeshIntersection::GetSourceDiscretization() const noexcept
{
  return source_;
}

std::shared_ptr<const Discretization>
OmegaHMeshIntersection::GetTargetDiscretization() const noexcept
{
  return target_;
}

const IntersectionResults& OmegaHMeshIntersection::GetTargetToSource()
  const noexcept
{
  return results_;
}

Omega_h::Mesh& OmegaHMeshIntersection::GetSourceMesh() const noexcept
{
  return source_->GetMesh();
}

Omega_h::Mesh& OmegaHMeshIntersection::GetTargetMesh() const noexcept
{
  return target_->GetMesh();
}

std::shared_ptr<MeshIntersection> IntersectMeshes(
  const FunctionSpace& source_space, const FunctionSpace& target_space,
  bool use_prefilter)
{
  const bool omega_h = std::dynamic_pointer_cast<const OmegaHDiscretization>(
                         source_space.GetDiscretization()) &&
                       std::dynamic_pointer_cast<const OmegaHDiscretization>(
                         target_space.GetDiscretization());
  if (omega_h) {
    return std::make_shared<OmegaHMeshIntersection>(source_space, target_space,
                                                    use_prefilter);
  }
  throw pcms_error("IntersectMeshes: no mesh intersection implementation for "
                   "this pair of discretizations");
}
} // namespace pcms
