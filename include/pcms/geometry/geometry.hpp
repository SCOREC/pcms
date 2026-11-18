#ifndef PCMS_GEOMETRY_GEOMETRY_HPP
#define PCMS_GEOMETRY_GEOMETRY_HPP

#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/coordinate_view.hpp"
#include "pcms/field/out_of_bounds_policy.hpp"
#include "pcms/geometry/bound_source.hpp"
#include "pcms/geometry/supports.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"
#include <Kokkos_Core.hpp>
#include <memory>
#include <string>

namespace pcms
{

class CoordinateMap; // pcms/field/coordinate_map.hpp
class FieldLayout;   // pcms/field/field_layout.hpp

/// @brief Options of the containing-cell search a geometry builds.
struct SearchOptions
{
  /// grid divisions per axis; 0 derives them from the element count
  int divisions = 0;
  Real tolerance = 1e-12;
};

/// @brief Options of Geometry::Neighbourhood.
struct NeighbourhoodOptions
{
  Real radius = 0.5;
  int min_supports = 10;
  /// upper bound of the supports of the mesh-adjacency search; 0 means the
  /// kind's own bound, 3 * min_supports, which is the only one supported yet
  /// (the kinds throw otherwise)
  int max_supports = 0;
  bool adapt_radius = true;
};

/// @brief Jacobian of the geometry map at bound points.
///
/// jacobian(r, i, j) = d x_j / d xi_i, determinant(r) = det. When
/// per_element is true, r is an element id; otherwise r indexes the bound
/// source's FoundPoints().
struct DifferentialView
{
  Kokkos::View<const Real***, DeviceMemorySpace> jacobian;
  Kokkos::View<const Real*, DeviceMemorySpace> determinant;
  bool per_element = false;
  [[nodiscard]] bool Empty() const noexcept { return jacobian.extent(0) == 0; }
};

/// @brief The coordinate function of a set of entities (a mesh, a grid, a
///        point cloud) and its inverse.
///
/// A Geometry has a native coordinate system, the search that locates
/// physical points among its entities and, when the kind has one, the
/// Jacobian of the map from parametric to physical coordinates. Function
/// spaces are created on a geometry, and their layouts carry it
/// (FieldLayout::GetGeometry); Bind produces bound sources every space on the
/// geometry accepts.
///
/// Lifetime: bound sources and the spaces created on a geometry hold it by
/// std::shared_ptr, so it outlives them; what the geometry refers to (e.g. a
/// mesh) must outlive the geometry. A bound source is immutable and
/// stays valid after later binds on the same geometry; evaluators hold their
/// bound source. Bind uses shared_from_this(), so a geometry must be owned by
/// a std::shared_ptr (create it with the kind's factory:
/// CreateUniformGridGeometry, CreatePointCloudGeometry or
/// backend::omega_h::CreateGeometry). Spaces and layouts are on the same
/// entities when they are on the same Geometry object: OwnerHints name entities
/// by the geometry's identity. A kind's factory (CreateGeometry, ...)
/// returns the most specific type it knows, as a pointer to const since a
/// geometry is immutable once built; more generally, a factory returns a
/// pointer to const for a type with no non-const public member (geometries and
/// function spaces), the object being constructed non-const.
class Geometry : public std::enable_shared_from_this<Geometry>
{
public:
  virtual ~Geometry() noexcept = default;
  Geometry(const Geometry&) = delete;
  Geometry& operator=(const Geometry&) = delete;

  /// @return the native coordinate system: Bind takes points in it
  [[nodiscard]] virtual const std::shared_ptr<const CoordinateSystem>&
  GetCoordinateSystem() const noexcept = 0;

  /// @return whether Bind accepts OutOfBoundsMode::NEAREST_BOUNDARY
  [[nodiscard]] virtual bool SupportsNearestBoundary() const noexcept
  {
    return false;
  }

  /// @brief Localizes physical query points.
  /// @param coords query points, in GetCoordinateSystem() or, with `map`, in
  ///        its source coordinate system; may be released after
  /// @param policy what to do with points outside the geometry
  /// @param map map from the query coordinate system to GetCoordinateSystem(),
  ///        or null; the points are located at their images, and evaluators
  ///        created from the result re-express their values through it
  /// @return bound source for the evaluator sources on this geometry; its
  ///         LocatedPointValues() are the query points (or their images)
  /// @throws pcms_error if the query coordinate system is not
  ///         GetCoordinateSystem() (without a map) or the map does not take it
  ///         to GetCoordinateSystem(), or the policy is NEAREST_BOUNDARY and
  ///         SupportsNearestBoundary() is false, or the policy is ERROR and a
  ///         point is outside
  [[nodiscard]] std::shared_ptr<const BoundSource> Bind(
    CoordinateView<DeviceMemorySpace> coords, OutOfBoundsPolicy policy = {},
    std::shared_ptr<const CoordinateMap> map = nullptr) const;

  /// @brief Bind for points with claimed owners, e.g. points generated in the
  ///        elements of this geometry.
  /// @param coords query points, as for the coordinates overload
  /// @param owners the claimed owner of each query point, to start from
  ///        instead of searching for it (a cell hint on a mesh; vertex and cell
  ///        hints for the supports of a reconstruction geometry); ignored if
  ///        it names another geometry, is not one entity per point, or is
  ///        given with a map (it describes the unmapped points). The kind
  ///        verifies the hint, see OwnerHints: a wrong one costs a search for
  ///        that point, and hinted points that lie outside the geometry are
  ///        missing, so the out-of-bounds policy applies to them as to any
  ///        other.
  /// @throws pcms_error as the coordinates overload, or if a used hint names
  ///         an entity the geometry does not have
  /// @param policy what to do with points outside the geometry
  /// @param map map from the query coordinate system to GetCoordinateSystem(),
  ///        or null, as for the coordinates overload
  [[nodiscard]] std::shared_ptr<const BoundSource> Bind(
    CoordinateView<DeviceMemorySpace> coords, const OwnerHints& owners,
    OutOfBoundsPolicy policy = {},
    std::shared_ptr<const CoordinateMap> map = nullptr) const;

  /// @return whether `owners` is a hint Bind uses for `num_points` points: one
  ///         entity per point, naming this geometry (OwnerHints::Names). A
  ///         geometry built on the entities of another (a point cloud over a
  ///         mesh's entities) also takes hints that name that one.
  [[nodiscard]] virtual bool TakesHints(const OwnerHints& owners,
                                        LO num_points) const noexcept
  {
    return owners.Names(*this, num_points);
  }

  /// @brief Localizes new query points in place of `previous`, a bound source
  ///        of this geometry: same policy and coordinate map.
  ///
  /// When nothing else holds `previous` (its shared_ptr and every Kokkos view
  /// it owns are uniquely held) its storage is reused, so a loop that
  /// releases its evaluators before each call reallocates none of the bound
  /// source's arrays in steady state when the point count does not grow (the
  /// kind's search may still allocate temporaries); otherwise this is an
  /// ordinary Bind. Binds with a coordinate map are not recycled in this
  /// version: they are bound anew with the map. An mdspan taken from a bound
  /// source is valid while a reference to that source is held, as always.
  /// @param previous bound source of this geometry; consumed
  /// @param coords query points, in the coordinate system Bind took them in
  /// @return the bound source of `coords`: `previous`, reused, or a new one
  /// @throws pcms_error if `previous` is null or belongs to another geometry,
  ///         or as Bind
  [[nodiscard]] std::shared_ptr<const BoundSource> Rebind(
    std::shared_ptr<const BoundSource> previous,
    CoordinateView<DeviceMemorySpace> coords) const;

  /// @brief Support sets of query points among this geometry's sites, for
  ///        reconstruction spaces.
  /// @throws pcms_error unless the kind overrides it, or if the coordinate
  ///         system does not have the identity metric (the distances are
  ///         Euclidean), or `options.max_supports` is neither 0 nor
  ///         3 * `options.min_supports` (other bounds are not honoured yet)
  /// @param site_dim the dimension of the entities whose coordinates are the
  ///        sites: 0 for a mesh's vertices, its dimension for its cell
  ///        centroids; the kinds that do not search among the entities of
  ///        their geometry (a point cloud, whose sites are its site layout's)
  ///        ignore it
  /// @param owners the claimed owner of each query point, as for Bind: the
  ///        kinds with an adjacency search start their search for the point
  ///        there (a vertex hint starts from an element around the vertex);
  ///        verified, so a wrong hint gives the supports of the unhinted
  ///        search
  [[nodiscard]] virtual SupportResults Neighbourhood(
    CoordinateView<DeviceMemorySpace> coords, int site_dim,
    const NeighbourhoodOptions& options, const OwnerHints& owners = {}) const;

  /// @brief Jacobian at the bound points.
  /// @param bound a bound source of this geometry
  /// @return an empty view unless the kind overrides it
  /// @throws pcms_error from the kinds that override it if `bound` belongs to
  ///         another geometry
  [[nodiscard]] virtual DifferentialView Differential(
    const BoundSource& bound) const;

protected:
  Geometry() = default;

  /// @param request points to locate, already checked to be in
  ///        GetCoordinateSystem() (mapped there if the bind had a map), and the
  ///        query points the bound source keeps
  /// @return the bound source of the kind
  [[nodiscard]] virtual std::shared_ptr<BoundSource> BindImpl(
    const BindRequest& request) const = 0;

private:
  // creates the passkey of BoundSource::Rebind
  struct RebindKey : BoundSource::RebindKey
  {
    RebindKey() = default;
  };
};

namespace detail
{

/// @brief The checks every Neighbourhood makes of its arguments.
/// @param geometry geometry whose Neighbourhood is called
/// @param coords the query points
/// @param options the search options
/// @param caller name of the calling function, for the message
/// @throws pcms_error if `coords` is not in the geometry's coordinate system,
///         that system does not have the identity metric, or
///         `options.max_supports` is neither 0 nor 3 * `options.min_supports`
void CheckNeighbourhoodArguments(
  const Geometry& geometry, const CoordinateView<DeviceMemorySpace>& coords,
  const NeighbourhoodOptions& options, const std::string& caller);

/// @brief Checks that every id of a hint is one of `num_entities` entities.
/// @param ids entity ids, e.g. OwnerHints::entities
/// @param num_entities number of entities the ids index
/// @param caller name of the checking function, for the message
/// @param entity what the ids index, e.g. "element", for the message
/// @throws pcms_error naming the range of the ids if one is outside
///         [0, num_entities)
void CheckEntityIds(Kokkos::View<const LO*, DeviceMemorySpace> ids,
                    LO num_entities, const std::string& caller,
                    const std::string& entity);

} // namespace detail

} // namespace pcms

#endif // PCMS_GEOMETRY_GEOMETRY_HPP
