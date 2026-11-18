#include "pcms/backends/omega_h/geometry.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/backends/omega_h/adj_search.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/omega_h_array_utils.hpp"
#include "pcms/utility/detail/profile.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/utility/uniform_grid.hpp"
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>

namespace pcms
{
using utility::ConvertCoordsTo2D;
using utility::flatten_to_omega_h_reals;
} // namespace pcms

namespace pcms::backend::omega_h
{
// pcms and Omega_h exchange data without conversion
static_assert(std::is_same_v<Real, Omega_h::Real>,
              "pcms and Omega_h real types must match");
static_assert(std::is_same_v<LO, Omega_h::LO>,
              "pcms and Omega_h LO types must match");
static_assert(std::is_same_v<GO, Omega_h::GO>,
              "pcms and Omega_h GO types must match");

} // namespace pcms::backend::omega_h

namespace pcms::detail
{

OmegaHMeshBinding::OmegaHMeshBinding(
  std::shared_ptr<const Geometry> geometry, const BindRequest& request,
  const backend::omega_h::GridPointSearchVariant& search, Omega_h::Mesh& mesh)
  : BoundSource(std::move(geometry), request), search_(search), mesh_(mesh)
{
  Locate(backend::omega_h::detail::ElementGuesses(mesh, request.owners,
                                                  "OmegaHMeshBinding"));
}

void OmegaHMeshBinding::RebindImpl()
{
  Locate({});
}

void OmegaHMeshBinding::Locate(Kokkos::View<const LO*, DeviceMemorySpace> hints)
{
  std::visit(
    [&](const auto& s) {
      constexpr int Dim = std::decay_t<decltype(s)>::DIM;
      using Results = Kokkos::View<typename std::decay_t<decltype(s)>::Result*,
                                   DeviceMemorySpace>;
      const LO n = NumPoints();
      const auto first = Kokkos::make_pair(LO{0}, n);
      if (!std::holds_alternative<Results>(point_results_)) {
        point_results_ = Results{};
      }
      auto& results = std::get<Results>(point_results_);
      if (results.extent(0) < static_cast<std::size_t>(n)) {
        num_records_ = 0; // the records do not survive a reallocation
      }
      pcms::detail::ReserveRows(results, "omega_h_point_results", n);
      pcms::detail::ReserveRows(point_guesses_, "omega_h_point_guesses", n);
      pcms::detail::ReserveRows(point_found_, "omega_h_point_found", n);
      const auto point_results = Kokkos::subview(results, first);
      const auto point_guesses = Kokkos::subview(point_guesses_, first);
      Kokkos::View<const LO*, DeviceMemorySpace> guesses;
      if (hints.extent(0) > 0) {
        Kokkos::deep_copy(point_guesses, hints);
        guesses = point_guesses;
      } else if (num_records_ > 0) {
        // the elements of the last records; the search overwrites the records
        backend::omega_h::detail::ElementsOfResults<Dim>(
          mesh_, point_results, num_records_, point_guesses);
        guesses = point_guesses;
      }
      s(LocatedPointValues(), guesses, point_results);
      const auto found = Kokkos::subview(point_found_, first);
      backend::omega_h::detail::FoundFlags<Dim>(point_results, found);
      num_records_ = n;
      SetFound(found);
    },
    search_);
}

} // namespace pcms::detail

namespace pcms::backend::omega_h
{

namespace detail
{

Kokkos::View<const LO*, DeviceMemorySpace> ElementGuesses(
  Omega_h::Mesh& mesh, const OwnerHints& owners, const std::string& caller)
{
  const int dim = mesh.dim();
  const int entity_dim = owners.entity_dim;
  if (owners.entities.extent(0) == 0 || entity_dim < 0 || entity_dim > dim) {
    return {};
  }
  static const char* const names[] = {"vertex", "edge", "face", "region"};
  pcms::detail::CheckEntityIds(owners.entities, mesh.nents(entity_dim), caller,
                               entity_dim == dim ? "element"
                                                 : names[entity_dim]);
  if (entity_dim == dim) {
    return owners.entities;
  }
  Kokkos::View<LO*, DeviceMemorySpace> guesses(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "omega_h_element_guesses"),
    owners.entities.extent(0));
  if (dim == 2) {
    ElementsOfEntities<2>(mesh, entity_dim, owners.entities, guesses);
  } else {
    ElementsOfEntities<3>(mesh, entity_dim, owners.entities, guesses);
  }
  return guesses;
}

} // namespace detail

GridPointSearchVariant MakeGridSearch(Omega_h::Mesh& mesh,
                                      const SearchOptions& options)
{
  if (mesh.dim() == 2) {
    const LO n = options.divisions > 0
                   ? options.divisions
                   : pcms::detail::DivisionsPerAxisForMesh<2>(mesh.nelems());
    return GridPointSearch2D(mesh, n, n, options.tolerance);
  }
  if (mesh.dim() == 3) {
    const LO n = options.divisions > 0
                   ? options.divisions
                   : pcms::detail::DivisionsPerAxisForMesh<3>(mesh.nelems());
    return GridPointSearch3D(mesh, n, n, n, options.tolerance);
  }
  throw std::invalid_argument(
    "MakeGridSearch: only 2D and 3D meshes are supported");
}

OmegaHGeometry::OmegaHGeometry(Omega_h::Mesh& mesh,
                               std::shared_ptr<const CoordinateSystem> system,
                               OmegaHGeometryOptions options)
  : mesh_(mesh),
    global_id_name_(std::move(options.global_id_name)),
    owned_mask_(std::move(options.owned_mask)),
    system_(
      pcms::detail::ResolveCoordinateSystem(std::move(system), mesh.dim())),
    search_options_(options.search)
{
  if (owned_mask_.size() != 0 && owned_mask_.size() != mesh_.nverts()) {
    throw pcms_error("OmegaHGeometry: the owned mask must have one entry per "
                     "mesh vertex");
  }
}

namespace
{
Omega_h::Mesh& RequireOwnedMesh(
  const std::shared_ptr<Omega_h::Mesh>& mesh,
  const std::shared_ptr<Omega_h::Library>& library)
{
  if (!mesh || !library || mesh->library() != library.get()) {
    throw pcms_error("OmegaHGeometry: expected a mesh and its library");
  }
  return *mesh;
}
} // namespace

OmegaHGeometry::OmegaHGeometry(std::shared_ptr<Omega_h::Mesh> mesh,
                               std::shared_ptr<Omega_h::Library> library,
                               std::shared_ptr<const CoordinateSystem> system,
                               OmegaHGeometryOptions options)
  : OmegaHGeometry(RequireOwnedMesh(mesh, library), std::move(system),
                   std::move(options))
{
  library_owner_ = std::move(library);
  mesh_owner_ = std::move(mesh);
}

const GridPointSearchVariant& OmegaHGeometry::Search() const
{
  std::call_once(search_once_, [this] {
    search_ = std::make_unique<GridPointSearchVariant>(
      MakeGridSearch(mesh_, search_options_));
  });
  return *search_;
}

std::shared_ptr<BoundSource> OmegaHGeometry::BindImpl(
  const BindRequest& request) const
{
  PCMS_FUNCTION_TIMER;
  return std::make_shared<pcms::detail::OmegaHMeshBinding>(
    shared_from_this(), request, Search(), mesh_);
}

Kokkos::View<const LO*, DeviceMemorySpace> OmegaHGeometry::FindElementSeeds(
  Kokkos::View<const Real**, DeviceMemorySpace> points,
  Kokkos::View<const LO*, DeviceMemorySpace> guesses) const
{
  if (points.extent(1) != static_cast<std::size_t>(mesh_.dim())) {
    throw pcms_error("OmegaHGeometry::FindElementSeeds: point dimension "
                     "differs from the mesh dimension");
  }
  return std::visit(
    [&](const auto& s) -> Kokkos::View<const LO*, DeviceMemorySpace> {
      constexpr int Dim = std::decay_t<decltype(s)>::DIM;
      Kokkos::View<typename std::decay_t<decltype(s)>::Result*,
                   DeviceMemorySpace>
        results(Kokkos::view_alloc(Kokkos::WithoutInitializing,
                                   "point search result"),
                points.extent(0));
      s(points, guesses, results);
      return detail::LiftElements<Dim>(
        mesh_, Kokkos::View<const typename std::decay_t<decltype(s)>::Result*,
                            DeviceMemorySpace>(results));
    },
    Search());
}

SupportResults OmegaHGeometry::Neighbourhood(
  CoordinateView<DeviceMemorySpace> coords, int site_dim,
  const NeighbourhoodOptions& options, const OwnerHints& owners) const
{
  pcms::detail::CheckNeighbourhoodArguments(*this, coords, options,
                                            "OmegaHGeometry::Neighbourhood");
  if (site_dim != 0 && site_dim != mesh_.dim()) {
    throw pcms_error("OmegaHGeometry::Neighbourhood: the sites must be the "
                     "mesh vertices (site_dim 0) or the element centroids "
                     "(site_dim " +
                     std::to_string(mesh_.dim()) + "), got " +
                     std::to_string(site_dim));
  }
  const auto flat_points = flatten_to_omega_h_reals(coords.GetValues());
  const LO n = static_cast<LO>(coords.GetValues().extent(0));
  const int dim = mesh_.dim();
  // the element Bind's search finds for each point, starting from the hint
  // (for a point outside the mesh, the element nearest it; -1, and so no
  // supports, if its grid cell has none)
  const auto points = ConvertCoordsTo2D(flat_points, n, dim);
  const auto seeds = FindElementSeeds(
    points, detail::ElementGuesses(
              mesh_, owners.Names(*this, n) ? owners : OwnerHints{},
              "OmegaHGeometry::Neighbourhood"));
  return searchNeighbors(mesh_, site_dim, flat_points, dim, seeds,
                         options.radius * options.radius, options.min_supports,
                         3 * options.min_supports, options.adapt_radius);
}

namespace
{

/// Determinant of a Dim x Dim matrix (Dim 2 or 3).
/// (A function of its own: CUDA does not allow a device lambda to capture a
/// variable for the first time inside an if constexpr.)
template <int Dim>
KOKKOS_INLINE_FUNCTION Real SimplexDeterminant(const Real (&j)[Dim][Dim])
{
  if constexpr (Dim == 2) {
    return j[0][0] * j[1][1] - j[0][1] * j[1][0];
  } else {
    return j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1]) -
           j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0]) +
           j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
  }
}

/// Per-element Jacobian rows J_i = x_{i+1} - x_0 of the affine simplex map
/// and its determinant.
template <int Dim>
DifferentialView SimplexDifferential(Omega_h::Mesh& mesh)
{
  const LO nelems = mesh.nelems();
  Kokkos::View<Real***, DeviceMemorySpace> jacobian("omega_h_jacobian", nelems,
                                                    Dim, Dim);
  Kokkos::View<Real*, DeviceMemorySpace> determinant("omega_h_jacobian_det",
                                                     nelems);
  const auto elem_verts = mesh.ask_elem_verts();
  const auto coords = mesh.coords();
  Kokkos::parallel_for(
    "OmegaHGeometry::Differential",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(0, nelems),
    KOKKOS_LAMBDA(const LO e) {
      const LO v0 = elem_verts[e * (Dim + 1)];
      Real j[Dim][Dim];
      for (int i = 0; i < Dim; ++i) {
        const LO vi = elem_verts[e * (Dim + 1) + i + 1];
        for (int d = 0; d < Dim; ++d) {
          j[i][d] = coords[vi * Dim + d] - coords[v0 * Dim + d];
          jacobian(e, i, d) = j[i][d];
        }
      }
      determinant(e) = SimplexDeterminant<Dim>(j);
    });
  return DifferentialView{jacobian, determinant, /*per_element=*/true};
}

} // namespace

DifferentialView OmegaHGeometry::Differential(const BoundSource& bound) const
{
  if (bound.GetGeometry().get() != this) {
    throw pcms_error(
      "Differential: the bound source belongs to another geometry");
  }
  std::call_once(differential_once_, [this] {
    differential_ = mesh_.dim() == 2 ? SimplexDifferential<2>(mesh_)
                                     : SimplexDifferential<3>(mesh_);
  });
  return differential_;
}

} // namespace pcms::backend::omega_h
