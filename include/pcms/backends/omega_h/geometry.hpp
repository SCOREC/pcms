#ifndef PCMS_BACKENDS_OMEGA_H_GEOMETRY_HPP
#define PCMS_BACKENDS_OMEGA_H_GEOMETRY_HPP

#include "pcms/geometry/geometry.hpp"
#include "pcms/backends/omega_h/point_search.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <variant>

namespace pcms::backend::omega_h
{

namespace detail
{

/// @brief The element to start the search of each point from, given owner
///        hints: the hinted element for an element hint, else the lowest
///        element adjacent to the hinted entity (LowestAdjacentElement), -1
///        for an entity without elements; empty for no hint.
/// @param owners the hints, as handed to Bind or Neighbourhood (already
///        checked to name the mesh's geometry and one entity per point), or
///        none
/// @param caller name of the calling function, for the error message
/// @throws pcms_error if a used id is not an entity of the mesh
Kokkos::View<const LO*, DeviceMemorySpace> ElementGuesses(
  Omega_h::Mesh& mesh, const OwnerHints& owners, const std::string& caller);

/// @return per entity dimension below Dim, the adjacency from the entities to
///         the elements around them
template <int Dim>
Omega_h::Few<Omega_h::Adj, Dim> UpAdjacencies(Omega_h::Mesh& mesh)
{
  Omega_h::Few<Omega_h::Adj, Dim> up;
  for (int d = 0; d < Dim; ++d) {
    up[d] = mesh.ask_up(d, Dim);
  }
  return up;
}

/// @return per entity dimension from 1 to Dim, the vertices of each entity
///         (entry 0 is unused: a vertex is its own vertex)
template <int Dim>
Omega_h::Few<Omega_h::LOs, Dim + 1> VerticesOfEntities(Omega_h::Mesh& mesh)
{
  Omega_h::Few<Omega_h::LOs, Dim + 1> verts;
  for (int d = 1; d <= Dim; ++d) {
    verts[d] = mesh.ask_verts_of(d);
  }
  return verts;
}

/// The lowest element adjacent to an entity: an element is its own, -1 for an
/// entity of another dimension than 0..Dim, a negative id or an entity without
/// elements.
/// @param up per entity dimension below Dim, the adjacency from the entities
///        to the elements around them (UpAdjacencies)
template <int Dim>
KOKKOS_INLINE_FUNCTION LO LowestAdjacentElement(
  const Omega_h::Few<Omega_h::Adj, Dim>& up, int entity_dim, LO entity)
{
  if (entity < 0 || entity_dim < 0 || entity_dim > Dim) {
    return -1;
  }
  if (entity_dim == Dim) {
    return entity;
  }
  const auto& adj = up[entity_dim];
  const LO begin = adj.a2ab[entity];
  const LO end = adj.a2ab[entity + 1];
  if (begin >= end) {
    return -1;
  }
  LO elem = adj.ab2b[begin];
  for (LO i = begin + 1; i < end; ++i) {
    elem = Kokkos::min(elem, adj.ab2b[i]);
  }
  return elem;
}

/// The element the record of a search lowest-adjoins: the lowest element
/// adjacent to the entity it names (for a point outside the mesh, its nearest
/// entity); -1 for NO_INTERSECT.
template <int Dim>
KOKKOS_INLINE_FUNCTION LO
ElementOfResult(const Omega_h::Few<Omega_h::Adj, Dim>& up,
                const typename PointLocalizationSearch<Dim>::Result& result)
{
  return LowestAdjacentElement(up, static_cast<int>(result.dimensionality),
                               result.EntityId());
}

/// The vertices of the entity a search record names and the barycentric
/// coordinates of the point on them, in the entity's vertex order: with xi the
/// reduced coordinates, w_0 = 1 - sum(xi) and w_j = xi_j. The entity is the
/// record's, also for a point outside the mesh (its nearest entity).
/// @param verts_of VerticesOfEntities
/// @param vertices out: the Dim + 1 vertices at most
/// @param weights out: the weights, as many as vertices
/// @return the number of vertices: the entity's dimension plus 1, 0 for
///         NO_INTERSECT
template <int Dim>
KOKKOS_INLINE_FUNCTION int EntityVertexWeights(
  const Omega_h::Few<Omega_h::LOs, Dim + 1>& verts_of,
  const typename PointLocalizationSearch<Dim>::Result& result, LO* vertices,
  Real* weights)
{
  const int d = static_cast<int>(result.dimensionality);
  const LO entity = result.EntityId();
  if (entity < 0 || d > Dim) {
    return 0;
  }
  weights[0] = 1.0;
  for (int j = 0; j < d; ++j) {
    weights[j + 1] = result.parametric[j];
    weights[0] -= result.parametric[j];
  }
  for (int j = 0; j <= d; ++j) {
    vertices[j] = d == 0 ? entity : verts_of[d][entity * (d + 1) + j];
  }
  return d + 1;
}

/// Whether each search record is in the mesh (not outside it and not
/// NO_INTERSECT), as the points the bound source evaluates.
template <int Dim>
void FoundFlags(
  Kokkos::View<const typename PointLocalizationSearch<Dim>::Result*,
               DeviceMemorySpace>
    results,
  Kokkos::View<bool*, DeviceMemorySpace> found)
{
  using Result = typename PointLocalizationSearch<Dim>::Result;
  pcms::detail::ForEach(
    "FoundFlags", static_cast<LO>(results.extent(0)),
    KOKKOS_LAMBDA(const LO i) {
      found(i) =
        !results(i).IsOutside() &&
        results(i).dimensionality != Result::Dimensionality::NO_INTERSECT;
    });
}

/// Writes to `elements(i)` the element of record i (ElementOfResult), for the
/// first `num_records` records and -1 beyond them, as the guesses of a search
/// of as many points as `elements` has.
template <int Dim>
void ElementsOfResults(
  Omega_h::Mesh& mesh,
  Kokkos::View<const typename PointLocalizationSearch<Dim>::Result*,
               DeviceMemorySpace>
    results,
  LO num_records, Kokkos::View<LO*, DeviceMemorySpace> elements)
{
  const auto up = UpAdjacencies<Dim>(mesh);
  pcms::detail::ForEach(
    "ElementsOfResults", static_cast<LO>(elements.extent(0)),
    KOKKOS_LAMBDA(const LO i) {
      elements(i) = i < num_records ? ElementOfResult<Dim>(up, results(i)) : -1;
    });
}

/// Writes to `elements(i)` the lowest element adjacent to the entity of
/// dimension `entity_dim` and id `entities(i)` (LowestAdjacentElement).
template <int Dim>
void ElementsOfEntities(Omega_h::Mesh& mesh, int entity_dim,
                        Kokkos::View<const LO*, DeviceMemorySpace> entities,
                        Kokkos::View<LO*, DeviceMemorySpace> elements)
{
  const auto up = UpAdjacencies<Dim>(mesh);
  pcms::detail::ForEach(
    "ElementsOfEntities", static_cast<LO>(elements.extent(0)),
    KOKKOS_LAMBDA(const LO i) {
      elements(i) = LowestAdjacentElement<Dim>(up, entity_dim, entities(i));
    });
}

/// The element every result of a search lowest-adjoins (ElementOfResult), e.g.
/// to seed a neighbourhood search.
template <int Dim>
Kokkos::View<LO*, DeviceMemorySpace> LiftElements(
  Omega_h::Mesh& mesh,
  Kokkos::View<const typename PointLocalizationSearch<Dim>::Result*,
               DeviceMemorySpace>
    results)
{
  Kokkos::View<LO*, DeviceMemorySpace> elements(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "omega_h_lift_elem_ids"),
    results.extent(0));
  ElementsOfResults<Dim>(mesh, results, static_cast<LO>(results.extent(0)),
                         elements);
  return elements;
}

} // namespace detail

/// Bound source of an Omega_h simplex mesh: the search's record of each query
/// point, the entity it lies on and its position on it (Results()). The binding
/// picks no element: an evaluator applies its own rule to the entity, e.g.
/// interpolating on the entity's own vertices, or taking the lowest adjacent
/// element of a shared entity.
///
/// A rebind starts from the element of each point's record in the last locate,
/// a guess that the search verifies like any hint, so the records and found
/// points are those of a new bind. The records are kept across rebinds and
/// reallocated only when the points outgrow them (after which the next rebind
/// has no guesses); the first bind has the guesses its hints give, else none.
///
/// Implementer API: not covered by the public-API stability promise.
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

  /// [query point]: the search's record of each point, the entity it lies on
  /// and its reduced coordinates on it (PointLocalizationSearch<Dim>::Result;
  /// `Dim` is the mesh dimension). The rows of missing points hold the
  /// search's outside or NO_INTERSECT record; read them only through
  /// FoundPoints(). Evaluators apply their own element rule to it (see
  /// detail::LowestAdjacentElement, detail::EntityVertexWeights).
  template <int Dim>
  [[nodiscard]] Kokkos::View<
    const typename PointLocalizationSearch<Dim>::Result*, DeviceMemorySpace>
  Results() const
  {
    using Results = Kokkos::View<typename PointLocalizationSearch<Dim>::Result*,
                                 DeviceMemorySpace>;
    return Kokkos::subview(std::get<Results>(point_results_),
                           Kokkos::make_pair(LO{0}, NumPoints()));
  }

protected:
  void RebindImpl() override;

  [[nodiscard]] bool UnsharedImpl() const noexcept override
  {
    return std::visit([](const auto& r) { return r.use_count() <= 1; },
                      point_results_);
  }

private:
  /// Locates the points with the search, starting from `hints`, an element
  /// per point or empty, else from the element of each point's record of the
  /// last locate (the lowest adjacent element of its entity); records the
  /// results in point_results_ and sets the found points.
  void Locate(Kokkos::View<const LO*, DeviceMemorySpace> hints);

  const GridPointSearchVariant& search_;
  Omega_h::Mesh& mesh_;
  // per query point, the search's record, kept across rebinds; the first
  // num_records_ are the last locate's
  std::variant<Kokkos::View<GridPointSearch2D::Result*, DeviceMemorySpace>,
               Kokkos::View<GridPointSearch3D::Result*, DeviceMemorySpace>>
    point_results_;
  Kokkos::View<LO*, DeviceMemorySpace> point_guesses_;
  Kokkos::View<bool*, DeviceMemorySpace> point_found_;
  LO num_records_ = 0;
};

/// Options of an Omega_h geometry.
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
/// CreateGeometry.
class OmegaHGeometry final : public Geometry
{
public:
  /// Own it by std::shared_ptr, as CreateGeometry
  /// (pcms/backends/omega_h/omega_h.hpp) does: Bind uses shared_from_this().
  /// @throws pcms_error if the options' owned mask is not empty and does not
  ///         have one entry per mesh vertex
  OmegaHGeometry(Omega_h::Mesh& mesh,
                 std::shared_ptr<const CoordinateSystem> system,
                 OmegaHGeometryOptions options = {});

  OmegaHGeometry(std::shared_ptr<Omega_h::Mesh> mesh,
                 std::shared_ptr<Omega_h::Library> library,
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
  ///
  /// Implementer API: not covered by the public-API stability promise.
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
  ///
  /// Implementer API: not covered by the public-API stability promise.
  [[nodiscard]] DifferentialView Differential(
    const BoundSource& bound) const override;

protected:
  /// Implementer API: not covered by the public-API stability promise.
  [[nodiscard]] std::shared_ptr<BoundSource> BindImpl(
    const BindRequest& request) const override;

private:
  /// The containing-element search, built on first use.
  const GridPointSearchVariant& Search() const;

  // The library must outlive the mesh and the cached search data.
  std::shared_ptr<Omega_h::Library> library_owner_;
  std::shared_ptr<Omega_h::Mesh> mesh_owner_;
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

/// The grid search over a mesh, with options.divisions divisions per axis when
/// nonzero and DivisionsPerAxisForMesh otherwise.
/// @throws std::invalid_argument unless the mesh is 2D or 3D
[[nodiscard]] GridPointSearchVariant MakeGridSearch(
  Omega_h::Mesh& mesh, const SearchOptions& options);

} // namespace pcms::backend::omega_h

#endif // PCMS_BACKENDS_OMEGA_H_GEOMETRY_HPP
