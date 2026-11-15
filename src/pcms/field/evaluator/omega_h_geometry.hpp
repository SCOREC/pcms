#ifndef PCMS_FIELD_EVALUATOR_OMEGA_H_GEOMETRY_HPP
#define PCMS_FIELD_EVALUATOR_OMEGA_H_GEOMETRY_HPP

#include "pcms/geometry/geometry.hpp"
#include "pcms/localization/point_search.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"
#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>
#include <memory>
#include <mutex>
#include <string>
#include <utility>

namespace pcms
{

namespace detail
{

/// @brief The element to start the search of each point from, given owner
///        hints: the hinted element for an element hint, an element around the
///        hinted vertex for a vertex hint (the first of the vertex's adjacent
///        elements); empty for no hint, or one of another dimension, and a
///        -1 for a vertex without elements.
/// @param owners the hints, as handed to Bind or Neighbourhood (already
///        checked to name the mesh's geometry and one entity per point), or
///        none
/// @param caller name of the calling function, for the error message
/// @throws pcms_error if a used id is not an entity of the mesh
Kokkos::View<const LO*, DeviceMemorySpace> ElementGuesses(
  Omega_h::Mesh& mesh, const OwnerHints& owners, const std::string& caller);

/// The element a search result lifts to, and the point's barycentric
/// coordinates in it ([vertex], in the element's vertex order). An element
/// (face in 2D, region in 3D) result is its own element, with the barycentric
/// coordinates its reduced coordinates give. A result on a lower-dimensional
/// entity (shared by several elements) lifts to the lowest adjacent element,
/// with the barycentric coordinates the reduced coordinates give on the
/// vertices of the entity and zero on the others. A point outside the mesh
/// lifts through its nearest entity; NO_INTERSECT lifts to element -1.
/// @param entities2elems per entity dimension below Dim, the adjacency from
///        the entities to the elements around them
/// @param entities2verts per entity dimension from 1 below Dim, the
///        adjacency from the entities to their vertices
/// @param elem_verts the vertices of each element
/// @param result the search's result
template <int Dim>
KOKKOS_INLINE_FUNCTION Kokkos::pair<LO, Omega_h::Vector<Dim + 1>> LiftToElement(
  const Omega_h::Few<Omega_h::Adj, Dim>& entities2elems,
  const Omega_h::Few<Omega_h::Adj, Dim>& entities2verts,
  const Omega_h::LOs& elem_verts,
  const typename PointLocalizationSearch<Dim>::Result& result)
{
  const int entity_dim = static_cast<int>(result.dimensionality);
  const LO entity = result.EntityId();
  auto bary = Omega_h::zero_vector<Dim + 1>();
  if (entity < 0 || entity_dim > Dim) {
    return {-1, bary};
  }
  // the barycentric coordinates of the point in the entity, in the order of
  // the entity's vertices
  Real weights[Dim + 1];
  weights[0] = 1.0;
  for (int j = 0; j < entity_dim; ++j) {
    weights[j + 1] = result.parametric[j];
    weights[0] -= result.parametric[j];
  }
  if (entity_dim == Dim) {
    for (int k = 0; k <= Dim; ++k) {
      bary[k] = weights[k];
    }
    return {entity, bary};
  }
  const auto& up = entities2elems[entity_dim];
  const LO begin = up.a2ab[entity];
  const LO end = up.a2ab[entity + 1];
  if (begin >= end) {
    return {-1, bary};
  }
  LO elem = up.ab2b[begin];
  for (LO i = begin + 1; i < end; ++i) {
    elem = Kokkos::min(elem, up.ab2b[i]);
  }
  for (int j = 0; j <= entity_dim; ++j) {
    const LO vertex =
      entity_dim == 0
        ? entity
        : entities2verts[entity_dim].ab2b[entity * (entity_dim + 1) + j];
    for (int k = 0; k <= Dim; ++k) {
      if (elem_verts[elem * (Dim + 1) + k] == vertex) {
        bary[k] = weights[j];
      }
    }
  }
  return {elem, bary};
}

/// The element and barycentric coordinates ([point][vertex]) every result of a
/// search lifts to, see LiftToElement. Without `with_bary` only the elements
/// are written (the coordinates view is empty).
template <int Dim>
std::pair<Kokkos::View<LO*, DeviceMemorySpace>,
          Kokkos::View<Real**, DeviceMemorySpace>>
LiftResults(Omega_h::Mesh& mesh,
            Kokkos::View<const typename PointLocalizationSearch<Dim>::Result*,
                         DeviceMemorySpace>
              results,
            bool with_bary = true)
{
  const auto n = static_cast<LO>(results.extent(0));
  Kokkos::View<LO*, DeviceMemorySpace> elem_ids(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "omega_h_lift_elem_ids"),
    n);
  Kokkos::View<Real**, DeviceMemorySpace> bary(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "omega_h_lift_bary"),
    with_bary ? n : 0, Dim + 1);
  Omega_h::Few<Omega_h::Adj, Dim> entities2elems;
  Omega_h::Few<Omega_h::Adj, Dim> entities2verts;
  for (int d = 0; d < Dim; ++d) {
    entities2elems[d] = mesh.ask_up(d, Dim);
    if (d > 0) {
      entities2verts[d] = mesh.ask_down(d, Omega_h::VERT);
    }
  }
  const auto elem_verts = mesh.ask_elem_verts();
  detail::ForEach(
    "LiftResults", n, KOKKOS_LAMBDA(const LO i) {
      const auto lifted = LiftToElement<Dim>(entities2elems, entities2verts,
                                             elem_verts, results(i));
      elem_ids(i) = lifted.first;
      if (with_bary) {
        for (int k = 0; k <= Dim; ++k) {
          bary(i, k) = lifted.second[k];
        }
      }
    });
  return {elem_ids, bary};
}

/// The element every result of a search lifts to (LiftResults without the
/// barycentric coordinates), e.g. to seed a neighbourhood search.
template <int Dim>
Kokkos::View<LO*, DeviceMemorySpace> LiftElements(
  Omega_h::Mesh& mesh,
  Kokkos::View<const typename PointLocalizationSearch<Dim>::Result*,
               DeviceMemorySpace>
    results)
{
  return LiftResults<Dim>(mesh, results, /*with_bary=*/false).first;
}

/// The element and barycentric coordinates of the found points, gathered from
/// per-point lifts: row k is the lift of point rows(k). The results go to
/// `elem_ids` and `bary`, reallocated only when too small (ReserveRows).
/// (A free function: CUDA does not allow device lambdas inside the generic
/// lambda that dispatches on the search.)
template <int Dim>
void GatherFoundLifts(Kokkos::View<const LO*, DeviceMemorySpace> rows,
                      Kokkos::View<const LO*, DeviceMemorySpace> lifted_elems,
                      Kokkos::View<const Real**, DeviceMemorySpace> lifted_bary,
                      Kokkos::View<LO*, DeviceMemorySpace>& elem_ids,
                      Kokkos::View<Real**, DeviceMemorySpace>& bary)
{
  const auto n = static_cast<LO>(rows.extent(0));
  ReserveRows(elem_ids, "omega_h_elem_ids", n);
  ReserveRows(bary, "omega_h_bary", n, Dim + 1);
  const auto elem_out = elem_ids;
  const auto bary_out = bary;
  detail::ForEach(
    "GatherFoundLifts", n, KOKKOS_LAMBDA(const LO k) {
      elem_out(k) = lifted_elems(rows(k));
      for (int v = 0; v <= Dim; ++v) {
        bary_out(k, v) = lifted_bary(rows(k), v);
      }
    });
}

} // namespace detail

/// Bound source of an Omega_h simplex mesh: element and barycentric
/// coordinates of each found point, in FoundPoints() order.
///
/// The search reports the entity a point lies on and its position on it; the
/// binding's policy picks the element for the evaluators, which need one:
/// shared entities (a point on an edge or vertex) take the lowest adjacent
/// element, hinted or not.
class OmegaHMeshBinding final : public BoundSource
{
public:
  /// Locates `request.coords` with `search`, starting from the elements the
  /// request's owner hints give (ElementGuesses); a wrong hint costs a grid
  /// search for that point, and a hinted point can be missing. `search` and
  /// `mesh` are the geometry's; a rebind searches with them again.
  OmegaHMeshBinding(std::shared_ptr<const Geometry> geometry,
                    const BindRequest& request,
                    const GridPointSearchVariant& search, Omega_h::Mesh& mesh);

  [[nodiscard]] Kokkos::View<const LO*, DeviceMemorySpace> ElementIds()
    const noexcept
  {
    return Kokkos::subview(elem_ids_, Kokkos::make_pair(LO{0}, NumFound()));
  }

  /// [found point][vertex]
  [[nodiscard]] Kokkos::View<const Real**, DeviceMemorySpace> Barycentric()
    const noexcept
  {
    return Kokkos::subview(bary_, Kokkos::make_pair(LO{0}, NumFound()),
                           Kokkos::ALL);
  }

protected:
  void RebindImpl() override;

  [[nodiscard]] bool UnsharedImpl() const noexcept override
  {
    return elem_ids_.use_count() <= 1 && bary_.use_count() <= 1;
  }

private:
  [[nodiscard]] LO NumFound() const noexcept
  {
    return static_cast<LO>(FoundPoints().extent(0));
  }

  /// Locates the points with the search, from the guessed elements (empty
  /// for none), sets the found points and lifts their results to elements
  /// (detail::LiftResults) into elem_ids_ and bary_.
  void Locate(Kokkos::View<const LO*, DeviceMemorySpace> guesses);

  const GridPointSearchVariant& search_;
  Omega_h::Mesh& mesh_;
  Kokkos::View<LO*, DeviceMemorySpace> elem_ids_;
  Kokkos::View<Real**, DeviceMemorySpace> bary_;
};

/// @brief Options of an Omega_h geometry.
struct OmegaHGeometryOptions
{
  /// divisions and tolerance of the containing-element search
  SearchOptions search;
  /// the tag holding global entity ids
  std::string global_id_name = "global";
  /// per-vertex ownership override: a vertex the mask is 0 for is not owned,
  /// whatever the mesh says; empty means the mesh's own ownership. Layouts
  /// whose DOF holders are the mesh vertices apply it; layouts of other
  /// entities do not.
  Omega_h::Read<Omega_h::I8> owned_mask;
};

/// Geometry of a 2D or 3D Omega_h simplex mesh with
/// straight-sided elements. The containing-element search over the mesh is
/// built by the first Bind or Neighbourhood, which costs about as much as a
/// bind: create one geometry per mesh and share it, so that every layout and
/// space on the mesh is on the same entities. Create it with
/// CreateOmegaHGeometry.
class OmegaHGeometry final : public Geometry
{
public:
  /// Own it by std::shared_ptr, as CreateOmegaHGeometry does: Bind uses
  /// shared_from_this().
  /// @throws pcms_error if the options' owned mask is not empty and does not
  ///         have one entry per mesh vertex
  OmegaHGeometry(Omega_h::Mesh& mesh,
                 std::shared_ptr<const CoordinateSystem> system,
                 OmegaHGeometryOptions options = {});

  [[nodiscard]] const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const noexcept override
  {
    return system_;
  }

  [[nodiscard]] Omega_h::Mesh& GetMesh() const noexcept { return mesh_; }

  /// @return the name of the tag that holds the global entity ids the layouts
  ///         on this geometry use
  [[nodiscard]] const std::string& GetGlobalIdName() const noexcept
  {
    return global_id_name_;
  }

  /// @return the per-vertex ownership override, empty if the mesh's own
  ///         ownership applies
  [[nodiscard]] Omega_h::Read<Omega_h::I8> GetOwnedMask() const noexcept
  {
    return owned_mask_;
  }

  /// One source element per query point, to seed adjacency searches. Points
  /// must be in the mesh's native coordinate system. `guesses`, empty or an
  /// element per point (-1 for none), makes the search start there, see
  /// GridPointSearch2D::operator()(points, guesses, results); they are
  /// verified and do not change the result. An outside point uses
  /// the nearest entity reported by the grid search, lifted to an element;
  /// -1 means the search cell has no candidate. Unlike Bind, this preserves
  /// outside-point seeds: an outside target centroid can still belong to an
  /// element that overlaps the mesh.
  [[nodiscard]] Kokkos::View<const LO*, DeviceMemorySpace> FindElementSeeds(
    Kokkos::View<const Real**, DeviceMemorySpace> points,
    Kokkos::View<const LO*, DeviceMemorySpace> guesses = {}) const;

  /// Supports among the mesh vertices (`site_dim` 0) or the centroids of its
  /// elements (`site_dim` the mesh dimension), by the mesh-adjacency search
  /// from the element Bind's search finds for each point, which starts from the
  /// cell or vertex hint of this geometry when there is one (verified, see
  /// OwnerHints: a bad hint is corrected). A point outside the mesh starts from
  /// the element nearest it, so it gets supports near the boundary; a point
  /// whose search cell holds no element gets none.
  /// @throws pcms_error if `site_dim` is neither 0 nor the mesh dimension
  [[nodiscard]] SupportResults Neighbourhood(
    CoordinateView<DeviceMemorySpace> coords, int site_dim,
    const NeighbourhoodOptions& options,
    const OwnerHints& owners = {}) const override;

  /// The per-element constant Jacobian of the affine simplex map, computed
  /// for all elements on the first call.
  [[nodiscard]] DifferentialView Differential(
    const BoundSource& bound) const override;

protected:
  [[nodiscard]] std::shared_ptr<BoundSource> BindImpl(
    const BindRequest& request) const override;

private:
  /// The containing-element search, built on first use.
  const GridPointSearchVariant& Search() const;

  Omega_h::Mesh& mesh_;
  std::string global_id_name_;
  Omega_h::Read<Omega_h::I8> owned_mask_;
  std::shared_ptr<const CoordinateSystem> system_;
  SearchOptions search_options_;
  // The search and Differential are built once, lazily; the geometry is const
  // and shared, hence the once_flags. These are the mutable members and they
  // are caches.
  mutable std::once_flag search_once_;
  mutable std::unique_ptr<GridPointSearchVariant> search_;
  mutable std::once_flag differential_once_;
  mutable DifferentialView differential_;
};

/// The geometry of an Omega_h simplex mesh (2D or 3D). Create one per mesh and
/// share it: the layouts and spaces on one geometry object share its entities.
/// @param mesh the mesh; must outlive the geometry
/// @param system the mesh's coordinate system, or null for Cartesian
/// @param options search, global id tag and ownership
std::shared_ptr<const OmegaHGeometry> CreateOmegaHGeometry(
  Omega_h::Mesh& mesh, std::shared_ptr<const CoordinateSystem> system,
  OmegaHGeometryOptions options = {});

/// The mesh of an Omega_h geometry.
/// @throws pcms_error if `geometry` is not an Omega_h geometry
Omega_h::Mesh& GetMesh(const Geometry& geometry);

/// The grid search over a mesh, with options.divisions divisions per axis when
/// nonzero and DivisionsPerAxisForMesh otherwise.
/// @throws std::invalid_argument unless the mesh is 2D or 3D
GridPointSearchVariant MakeGridSearch(Omega_h::Mesh& mesh,
                                      const SearchOptions& options);

} // namespace pcms

#endif // PCMS_FIELD_EVALUATOR_OMEGA_H_GEOMETRY_HPP
