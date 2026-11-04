#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_shape.hpp>
#include <pcms/configuration.hpp>
#include <pcms/field/coordinate_systems/cartesian.hpp>
#include <pcms/field/coordinate_systems/cylindrical.hpp>
#include <pcms/field/evaluator/omega_h_geometry.hpp>
#include <pcms/field/evaluator/point_cloud_geometry.hpp>
#include <pcms/field/evaluator/uniform_grid_geometry.hpp>
#include <pcms/field/function_space/lagrange.hpp>
#include <pcms/field/function_space/polynomial_reconstruction.hpp>
#include <pcms/field/layout/omega_h_entity.hpp>
#include <pcms/field/layout/omega_h_lagrange.hpp>
#include <pcms/field/layout/empty.hpp>
#include <pcms/field/layout/point_cloud.hpp>
#include <pcms/geometry/geometry.hpp>
#include <pcms/localization/adj_search.hpp>
#include <pcms/localization/mls_support_helpers.hpp>
#include <pcms/utility/mesh_geometry.hpp>
#include <pcms/utility/omega_h_array_utils.hpp>
#include "field_test_utils.hpp"
#include <algorithm>
#include <memory>
#include <random>
#include <string>
#include <variant>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using pcms::DeviceMemorySpace;
using pcms::FunctionSpace;
using pcms::HostMemorySpace;
using pcms::LagrangeBackend;
using pcms::LagrangeOptions;
using pcms::LO;
using pcms::Real;

namespace
{

std::shared_ptr<const pcms::Geometry> MakeGeometry(
  Omega_h::Mesh& mesh, pcms::SearchOptions options = {})
{
  return pcms::CreateOmegaHGeometry(
    mesh, pcms::csys::Cartesian::Create(2),
    pcms::OmegaHGeometryOptions{.search = options});
}

/// The element the grid search of `mesh` finds for each point, lifted from the
/// entity it reports: where the adjacency search of the point starts.
Kokkos::View<LO*> StartElements(Omega_h::Mesh& mesh,
                                Kokkos::View<const Real**> points)
{
  return std::visit(
    [&](const auto& s) {
      constexpr int Dim = std::decay_t<decltype(s)>::DIM;
      const auto results = s(points);
      return pcms::detail::LiftElements<Dim>(
        mesh, Kokkos::View<const typename std::decay_t<decltype(s)>::Result*,
                           pcms::DeviceMemorySpace>(results));
    },
    pcms::MakeGridSearch(mesh, {}));
}

std::shared_ptr<const FunctionSpace> MakeSpace(
  const std::shared_ptr<const pcms::Geometry>& geometry, int order,
  int num_components)
{
  return CreateLagrangeSpace(
    geometry, order, num_components,
    LagrangeOptions{.backend = LagrangeBackend::OmegaH});
}

/// Checks the orientation convention of the Jacobian of element `e`: row i is
/// x_{i+1} - x_0 for the element's vertices x_0, ..., x_Dim.
template <int Dim>
void CheckJacobianRows(Omega_h::Mesh& mesh,
                       const pcms::DifferentialView& differential, LO e)
{
  const auto jacobian = Kokkos::create_mirror_view_and_copy(
    HostMemorySpace(), differential.jacobian);
  const auto elem_verts = Omega_h::HostRead<Omega_h::LO>(mesh.ask_elem_verts());
  const auto coords = Omega_h::HostRead<Omega_h::Real>(mesh.coords());
  const LO v0 = elem_verts[(Dim + 1) * e];
  for (int i = 0; i < Dim; ++i) {
    const LO vi = elem_verts[(Dim + 1) * e + i + 1];
    for (int d = 0; d < Dim; ++d) {
      CAPTURE(e, i, d);
      REQUIRE_THAT(jacobian(e, i, d),
                   WithinAbs(coords[Dim * vi + d] - coords[Dim * v0 + d],
                             pcms::test::ExactTol));
    }
  }
}

std::shared_ptr<const pcms::Geometry> MakeMeshGeometry(Omega_h::Mesh& mesh)
{
  return pcms::CreateOmegaHGeometry(mesh,
                                    pcms::csys::Cartesian::Create(mesh.dim()));
}

/// What a bind of an Omega_h geometry found, on the host.
struct BindSummary
{
  std::vector<LO> found;
  std::vector<LO> missing;
  std::vector<pcms::test::HostRecord> records; // of the found points
};

template <typename View>
std::vector<typename View::non_const_value_type> ToVector(const View& view)
{
  const auto host =
    Kokkos::create_mirror_view_and_copy(HostMemorySpace(), view);
  return {host.data(), host.data() + host.size()};
}

BindSummary Summarize(const std::shared_ptr<const pcms::BoundSource>& bound)
{
  return {ToVector(bound->FoundPoints()), ToVector(bound->MissingPoints()),
          pcms::test::FoundRecords(*bound)};
}

void CheckSameBind(const BindSummary& actual, const BindSummary& expected)
{
  REQUIRE(actual.found == expected.found);
  REQUIRE(actual.missing == expected.missing);
  pcms::test::CheckSameRecords(actual.records, expected.records);
}

Kokkos::View<LO*, DeviceMemorySpace> ToDevice(const std::vector<LO>& ids)
{
  Kokkos::View<LO*, DeviceMemorySpace> device("ids", ids.size());
  auto host = Kokkos::create_mirror_view(device);
  for (size_t i = 0; i < ids.size(); ++i) {
    host(i) = ids[i];
  }
  Kokkos::deep_copy(device, host);
  return device;
}

/// `n` random points of [-0.2, 1.2]^dim (some outside the unit box mesh),
/// interleaved.
std::vector<Real> RandomPoints(int dim, int n, unsigned seed)
{
  std::mt19937 gen(seed);
  std::uniform_real_distribution<Real> uni(-0.2, 1.2);
  std::vector<Real> pts(static_cast<size_t>(dim) * n);
  for (auto& x : pts) {
    x = uni(gen);
  }
  return pts;
}

KOKKOS_INLINE_FUNCTION Real Quadratic(Real x, Real y)
{
  return 1.0 + x - 2.0 * y + 3.0 * x * x - x * y + 0.5 * y * y;
}

} // namespace

TEST_CASE("One bind serves every space on a geometry")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto geometry = MakeGeometry(mesh);
  auto scalars = MakeSpace(geometry, 1, 1);
  auto vectors = MakeSpace(geometry, 1, 3);
  auto cells = MakeSpace(geometry, 0, 1);
  REQUIRE(scalars->GetGeometry() == geometry);
  REQUIRE(vectors->GetGeometry() == geometry);
  REQUIRE(cells->GetGeometry() == geometry);

  auto f = scalars->CreateFunction<Real>("f");
  pcms::test::SetField(
    f, OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  const auto pts = pcms::test::StandardEvalCoords2D();
  const auto bound = pcms::test::Bind(*geometry, pts);
  REQUIRE(bound->GetGeometry() == geometry);
  pcms::test::CheckEvaluation(
    *scalars->CreatePointEvaluator<Real>(bound), f, pts,
    OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  // an evaluator still takes only fields of the space that created it
  REQUIRE_THROWS_WITH(
    pcms::test::EvaluateToHost(*cells->CreatePointEvaluator<Real>(bound), f,
                               f.GetData().GetValueBasis()),
    ContainsSubstring("not a field of the function space"));

  // the three components of the vector space see the same points
  auto g = vectors->CreateFunction<Real>("g");
  pcms::test::SetFieldComponents(g, [](Real x, Real y, Real* out) {
    for (int c = 0; c < 3; ++c) {
      out[c] = pcms::test::linear_f(x, y) + c;
    }
  });
  const auto values =
    pcms::test::EvaluateToHost(*vectors->CreatePointEvaluator<Real>(bound), g,
                               g.GetData().GetValueBasis(), 3);
  for (size_t i = 0; i < pts.size() / 2; ++i) {
    for (int c = 0; c < 3; ++c) {
      CAPTURE(i, c);
      REQUIRE_THAT(
        values(i, c),
        WithinAbs(pcms::test::linear_f(pts[2 * i], pts[2 * i + 1]) + c,
                  pcms::test::EvalTol));
    }
  }
}

TEST_CASE("Spaces on a geometry share it")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 4);
  auto geometry = MakeGeometry(mesh);
  auto scalars = MakeSpace(geometry, 1, 1);
  auto cells = MakeSpace(geometry, 0, 1);
  REQUIRE(scalars->GetLayout()->GetGeometry() == geometry);
  REQUIRE(cells->GetLayout()->GetGeometry() == geometry);
  REQUIRE(scalars->GetGeometry() == geometry);
#ifdef PCMS_ENABLE_MESHFIELDS
  auto p2 = CreateLagrangeSpace(
    geometry, 2, 1, LagrangeOptions{.backend = LagrangeBackend::MeshFields});
  REQUIRE(p2->GetLayout()->GetGeometry() == geometry);
#endif
  // a space on another geometry of the same mesh has the same entities, but
  // is on another geometry: owner hints name geometries by identity
  auto other = MakeSpace(MakeGeometry(mesh), 1, 1);
  REQUIRE(other->GetLayout()->GetGeometry() != geometry);
  // the factory of a kind returns the kind's type, which a layout takes
  // directly
  const pcms::OmegaHEntityLayout direct(
    pcms::CreateOmegaHGeometry(mesh, pcms::csys::Cartesian::Create(2)),
    pcms::Vertex, 1);
  REQUIRE(direct.GetGeometry() != nullptr);
  // a layout needs its geometry
  REQUIRE_THROWS_WITH(pcms::OmegaHLagrangeLayout(nullptr, 1, 1),
                      ContainsSubstring("the geometry must not be null"));
  REQUIRE_THROWS_WITH(pcms::OmegaHEntityLayout(nullptr, pcms::Vertex, 1),
                      ContainsSubstring("the geometry must not be null"));
}

namespace
{

int num_grid_allocations = 0;

void CountGridAllocation(const Kokkos::Profiling::SpaceHandle,
                         const char* label, const void*, const uint64_t)
{
  if (std::string(label) == "uniform grid") {
    ++num_grid_allocations;
  }
}

} // namespace

// the variance signatures of the common value types are constants
static_assert(pcms::values::Scalar.Rank() == 0);
static_assert(pcms::values::Vector.Rank() == 1);
static_assert(pcms::values::Tensor.Rank() == 2);
static_assert(pcms::values::Vector == pcms::values::Vector);
static_assert(!(pcms::values::Vector == pcms::values::Covector));
static_assert(pcms::values::Scalar.WithCovariantIndexAppended() ==
              pcms::values::Covector);

TEST_CASE("A layout is named by the options of its space")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  auto geometry = MakeGeometry(mesh);
  REQUIRE(MakeSpace(geometry, 1, 1)->GetLayout()->GetName().empty());
  const auto named = CreateLagrangeSpace(
    geometry, 1, 1,
    LagrangeOptions{.backend = LagrangeBackend::OmegaH, .layout_name = "gids"});
  REQUIRE(named->GetLayout()->GetName() == "gids");
#ifdef PCMS_ENABLE_MESHFIELDS
  REQUIRE(
    CreateLagrangeSpace(geometry, 1, 1,
                        LagrangeOptions{.backend = LagrangeBackend::MeshFields,
                                        .layout_name = "mf"})
      ->GetLayout()
      ->GetName() == "mf");
#endif
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {2, 2};
  REQUIRE(CreateLagrangeSpace(pcms::CreateUniformGridGeometry(
                                grid, pcms::csys::Cartesian::Create(2)),
                              1, 1, LagrangeOptions{.layout_name = "grid"})
            ->GetLayout()
            ->GetName() == "grid");
  // a layout built directly takes its name as a constructor argument
  REQUIRE(pcms::OmegaHLagrangeLayout(pcms::test::MakeOmegaHGeometry(mesh), 1, 1,
                                     "direct")
            .GetName() == "direct");
  REQUIRE(pcms::EmptyFieldLayout(pcms::csys::Cartesian::Create(2), 1, "empty")
            .GetName() == "empty");
}

TEST_CASE("The search of a geometry is built when first needed")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 4);
  num_grid_allocations = 0;
  Kokkos::Tools::Experimental::set_allocate_data_callback(CountGridAllocation);
  // a geometry, and layouts, spaces and functions on it, need no search: one
  // made for a coupler only costs nothing
  const auto geometry = pcms::test::MakeOmegaHGeometry(mesh);
  auto space = MakeSpace(geometry, 1, 1);
  auto f = space->CreateFunction<Real>("f");
  const int before_bind = num_grid_allocations;
  const auto pts = pcms::test::StandardEvalCoords2D();
  const auto bound = pcms::test::Bind(*geometry, pts);
  const int after_bind = num_grid_allocations;
  // later binds and neighbourhoods reuse it
  const auto again = pcms::test::Bind(*geometry, pts);
  const auto coords = pcms::test::CreateDeviceCoordinateView(
    pts, geometry->GetCoordinateSystem());
  const auto supports = geometry->Neighbourhood(
    coords.coordinate_view, 0, pcms::NeighbourhoodOptions{.min_supports = 3});
  const int after_all = num_grid_allocations;
  Kokkos::Tools::Experimental::set_allocate_data_callback(nullptr);
  REQUIRE(before_bind == 0);
  REQUIRE(after_bind > 0);
  REQUIRE(after_all == after_bind);
}

TEST_CASE("Correct owner hints give the search's bind")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 4);
  auto geometry = MakeGeometry(mesh);
  auto space = MakeSpace(geometry, 1, 1);
  auto f = space->CreateFunction<Real>("f");
  pcms::test::SetField(
    f, OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });

  // the element centroids, with their elements as owners
  const LO nelems = mesh.nelems();
  const auto centroids =
    pcms::ConvertCoordsTo2D(pcms::get_entity_centroids(mesh, 2), nelems, 2);
  const auto coords =
    pcms::test::MakeCoords(centroids, geometry->GetCoordinateSystem());
  Kokkos::View<LO*, DeviceMemorySpace> elements("elements", nelems);
  pcms::iota_view<LO>(elements);
  const auto searched = geometry->Bind(coords);
  const auto hinted =
    geometry->Bind(coords, pcms::OwnerHints{geometry, 2, elements});
  REQUIRE(hinted->FoundPoints().extent(0) == static_cast<size_t>(nelems));
  // inside the mesh, the hint gives what the search finds
  const auto searched_ids = pcms::test::FoundRecords(*searched);
  const auto hinted_ids = pcms::test::FoundRecords(*hinted);
  const auto searched_values =
    pcms::test::EvaluateToHost(*space->CreatePointEvaluator<Real>(searched), f,
                               f.GetData().GetValueBasis());
  const auto hinted_values = pcms::test::EvaluateToHost(
    *space->CreatePointEvaluator<Real>(hinted), f, f.GetData().GetValueBasis());
  for (LO e = 0; e < nelems; ++e) {
    CAPTURE(e);
    REQUIRE(searched_ids[e].dim == 2);
    REQUIRE(searched_ids[e].entity == e);
    REQUIRE(hinted_ids[e].dim == 2);
    REQUIRE(hinted_ids[e].entity == e);
    REQUIRE_THAT(hinted_values(e, 0),
                 WithinAbs(searched_values(e, 0), pcms::test::ExactTol));
  }

  // ids must be entities of the geometry
  Kokkos::View<LO*, DeviceMemorySpace> beyond("beyond", nelems);
  Kokkos::deep_copy(beyond, elements);
  Kokkos::deep_copy(Kokkos::subview(beyond, nelems - 1), nelems);
  REQUIRE_THROWS_WITH(
    geometry->Bind(coords, pcms::OwnerHints{geometry, 2, beyond}),
    ContainsSubstring("element ids span [0, " + std::to_string(nelems) +
                      "], outside [0, " + std::to_string(nelems) + ")"));

  // a hint for another geometry is ignored: the points are searched
  Kokkos::View<LO*, DeviceMemorySpace> wrong("wrong", nelems);
  const auto ignored =
    geometry->Bind(coords, pcms::OwnerHints{MakeGeometry(mesh), 2, wrong});
  const auto ignored_ids = pcms::test::FoundRecords(*ignored);
  for (LO e = 0; e < nelems; ++e) {
    REQUIRE(ignored_ids[e].entity == e);
  }
}

TEST_CASE("Differential gives the element measure")
{
  auto lib = Omega_h::Library{};
  SECTION("triangles")
  {
    auto mesh = pcms::test::BuildSquareMesh(lib, 5);
    auto geometry = MakeGeometry(mesh);
    const auto pts = pcms::test::StandardEvalCoords2D();
    const auto differential =
      geometry->Differential(*pcms::test::Bind(*geometry, pts));
    REQUIRE(differential.per_element);
    REQUIRE(differential.jacobian.extent(0) ==
            static_cast<size_t>(mesh.nelems()));
    const auto det = Kokkos::create_mirror_view_and_copy(
      HostMemorySpace(), differential.determinant);
    const auto measures =
      Omega_h::HostRead<Omega_h::Real>(Omega_h::measure_elements_real(&mesh));
    for (LO e = 0; e < mesh.nelems(); ++e) {
      CAPTURE(e);
      REQUIRE_THAT(det(e) / 2.0, WithinAbs(measures[e], pcms::test::ExactTol));
    }
    CheckJacobianRows<2>(mesh, differential, mesh.nelems() / 2);
    // the bound source must be one of this geometry
    auto other = MakeGeometry(mesh);
    REQUIRE_THROWS_WITH(geometry->Differential(*pcms::test::Bind(*other, pts)),
                        ContainsSubstring("belongs to another geometry"));
  }
  SECTION("tetrahedra")
  {
    auto mesh = pcms::test::BuildUnitCube(lib, 2);
    auto geometry =
      pcms::CreateOmegaHGeometry(mesh, pcms::csys::Cartesian::Create(3));
    const auto differential =
      geometry->Differential(*pcms::test::Bind(*geometry, {0.5, 0.5, 0.5}));
    const auto det = Kokkos::create_mirror_view_and_copy(
      HostMemorySpace(), differential.determinant);
    const auto measures =
      Omega_h::HostRead<Omega_h::Real>(Omega_h::measure_elements_real(&mesh));
    for (LO e = 0; e < mesh.nelems(); ++e) {
      CAPTURE(e);
      REQUIRE_THAT(det(e) / 6.0, WithinAbs(measures[e], pcms::test::ExactTol));
    }
    CheckJacobianRows<3>(mesh, differential, mesh.nelems() - 1);
  }
}

TEST_CASE("Search options are honoured")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto coarse = MakeGeometry(mesh, pcms::SearchOptions{3, 1e-12});
  auto derived = MakeGeometry(mesh);
  auto coarse_space = MakeSpace(coarse, 1, 1);
  auto derived_space = MakeSpace(derived, 1, 1);
  auto f = coarse_space->CreateFunction<Real>("f");
  pcms::test::SetField(
    f, OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  auto g = derived_space->CreateFunction<Real>("g");
  pcms::test::SetField(
    g, OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });

  // inside and outside points
  std::vector<Real> pts = pcms::test::StandardEvalCoords2D();
  for (Real x : pcms::test::StandardOutsideCoords2D()) {
    pts.push_back(x);
  }
  const pcms::OutOfBoundsPolicy fill{pcms::OutOfBoundsMode::FILL, -1.0};
  const auto coarse_bound = pcms::test::Bind(*coarse, pts, fill);
  const auto derived_bound = pcms::test::Bind(*derived, pts, fill);
  const auto coarse_found = Kokkos::create_mirror_view_and_copy(
    HostMemorySpace(), coarse_bound->FoundPoints());
  const auto derived_found = Kokkos::create_mirror_view_and_copy(
    HostMemorySpace(), derived_bound->FoundPoints());
  REQUIRE(coarse_found.extent(0) == 5);
  REQUIRE(coarse_found.extent(0) == derived_found.extent(0));
  for (size_t k = 0; k < coarse_found.extent(0); ++k) {
    REQUIRE(coarse_found(k) == derived_found(k));
  }
  const auto a = pcms::test::EvaluateToHost(
    *coarse_space->CreatePointEvaluator<Real>(coarse_bound), f,
    f.GetData().GetValueBasis());
  const auto b = pcms::test::EvaluateToHost(
    *derived_space->CreatePointEvaluator<Real>(derived_bound), g,
    g.GetData().GetValueBasis());
  for (size_t i = 0; i < a.extent(0); ++i) {
    CAPTURE(i);
    REQUIRE_THAT(a(i, 0), WithinAbs(b(i, 0), pcms::test::ExactTol));
  }

  // the tolerance decides whether a point just outside the boundary is found
  const std::vector<Real> near = {0.5, 0.5, 1.0 + 1e-8, 0.5};
  const auto loose = MakeGeometry(mesh, pcms::SearchOptions{.tolerance = 1e-6});
  const auto tight =
    MakeGeometry(mesh, pcms::SearchOptions{.tolerance = 1e-12});
  REQUIRE(pcms::test::Bind(*loose, near, fill)->FoundPoints().extent(0) == 2);
  const auto tight_bound = pcms::test::Bind(*tight, near, fill);
  REQUIRE(tight_bound->FoundPoints().extent(0) == 1);
  const auto tight_missing = Kokkos::create_mirror_view_and_copy(
    HostMemorySpace(), tight_bound->MissingPoints());
  REQUIRE(tight_missing.extent(0) == 1);
  REQUIRE(tight_missing(0) == 1);
}

TEST_CASE("Neighbourhood returns vertex supports")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  auto geometry = MakeGeometry(mesh);
  auto coords = pcms::test::CreateDeviceCoordinateView(
    pcms::test::StandardEvalCoords2D(), geometry->GetCoordinateSystem());
  pcms::NeighbourhoodOptions options{.radius = 0.3, .min_supports = 6};
  const pcms::SupportResults actual =
    geometry->Neighbourhood(coords.coordinate_view, 0, options);
  Omega_h::Real radius_sq = 0.09;
  const pcms::SupportResults expected = pcms::searchNeighbors(
    mesh, 0, pcms::flatten_to_omega_h_reals(coords.coordinate_view.GetValues()),
    2, StartElements(mesh, coords.view), radius_sq, 6, 18);
  pcms::test::CheckSupportResultsEquivalent(actual, expected);
  // a point outside the mesh starts from the element nearest it
  auto outside = pcms::test::CreateDeviceCoordinateView(
    {1.05, 0.5}, geometry->GetCoordinateSystem());
  REQUIRE(geometry->Neighbourhood(outside.coordinate_view, 0, options)
            .supports_idx.size() >= 6);
  // the options matter: the defaults give other supports
  const pcms::SupportResults defaults = geometry->Neighbourhood(
    coords.coordinate_view, 0, pcms::NeighbourhoodOptions{});
  REQUIRE(defaults.supports_idx.size() != actual.supports_idx.size());

  // other upper bounds than 3 * min_supports are not honoured yet, so they
  // are refused
  options.max_supports = 12;
  REQUIRE_THROWS_WITH(
    geometry->Neighbourhood(coords.coordinate_view, 0, options),
    ContainsSubstring(
      "max_supports must be 0 or 3 * min_supports (18), got 12"));
  // the query points must be in the geometry's coordinate system
  auto cylindrical = pcms::test::CreateDeviceCoordinateView(
    pcms::test::StandardEvalCoords2D(), pcms::csys::CylindricalRZ::Create());
  REQUIRE_THROWS_WITH(
    geometry->Neighbourhood(cylindrical.coordinate_view, 0, {}),
    ContainsSubstring("is not the geometry's coordinate system"));
}

TEST_CASE("PointCloudGeometry::Neighbourhood takes the search of its sites")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 8);
  const auto system = pcms::csys::Cartesian::Create(2);
  const auto mesh_geometry = pcms::test::MakeOmegaHGeometry(mesh);
  auto vertices =
    std::make_shared<pcms::OmegaHEntityLayout>(mesh_geometry, pcms::Vertex, 1);
  auto faces =
    std::make_shared<pcms::OmegaHEntityLayout>(mesh_geometry, pcms::Face, 1);
  const auto centroids = pcms::get_entity_centroids(mesh, pcms::Face);
  auto cloud = std::make_shared<pcms::PointCloudLayout>(
    2, pcms::ConvertCoordsTo2D(centroids, mesh.nfaces(), 2), system);
  const auto query = vertices->GetDOFHolderCoordinates();
  const auto query_points =
    pcms::ConvertCoordsTo2D(mesh.coords(), mesh.nverts(), 2);
  pcms::NeighbourhoodOptions options{.radius = 0.25, .min_supports = 8};
  Omega_h::Real radius_sq = options.radius * options.radius;

  // a bare point cloud has no geometry to search: the distance search
  pcms::test::CheckSupportResultsEquivalent(
    pcms::CreatePointCloudGeometry(cloud)->Neighbourhood(query, 0, options),
    pcms::BuildPointCloudSupports(centroids, mesh.coords(), 2, options.radius,
                                  8, true));
  // sites that are entities of a mesh geometry's mesh, vertices or element
  // centroids: the adjacency search of the mesh from the element of each
  // query point (of the entity dimension the sites are)
  const auto start = StartElements(mesh, query_points);
  const auto vertex_sites = pcms::CreatePointCloudGeometry(vertices);
  const auto centroid_sites = pcms::CreatePointCloudGeometry(faces);
  const auto vertex_supports =
    pcms::searchNeighbors(mesh, 0, mesh.coords(), 2, start, radius_sq, 8, 24);
  const auto centroid_supports =
    pcms::searchNeighbors(mesh, 2, mesh.coords(), 2, start, radius_sq, 8, 24);
  pcms::test::CheckSupportResultsEquivalent(
    vertex_sites->Neighbourhood(query, 0, options), vertex_supports);
  pcms::test::CheckSupportResultsEquivalent(
    centroid_sites->Neighbourhood(query, 0, options), centroid_supports);

  // the DOF holders of a layout on the mesh geometry are hints for the
  // cloud: they name the geometry its sites are entities of; the cloud takes
  // them, and a hint naming the cloud itself, and not those of another
  // geometry
  const auto vertex_owners = pcms::DOFHolderOwners(*vertices);
  const auto n = static_cast<LO>(mesh.nverts());
  REQUIRE(vertex_owners.entity_dim == pcms::Vertex);
  REQUIRE(centroid_sites->TakesHints(vertex_owners, n));
  REQUIRE(vertex_sites->TakesHints(vertex_owners, n));
  REQUIRE(centroid_sites->TakesHints(
    pcms::OwnerHints{centroid_sites, 0, vertex_owners.entities}, n));
  REQUIRE_FALSE(centroid_sites->TakesHints(
    pcms::OwnerHints{MakeGeometry(mesh), 0, vertex_owners.entities}, n));
  REQUIRE_FALSE(MakeGeometry(mesh)->TakesHints(vertex_owners, n));
  pcms::test::CheckSupportResultsEquivalent(
    centroid_sites->Neighbourhood(query, 0, options, vertex_owners),
    centroid_supports);
  // the owners are verified: a vertex far from the point, for the geometry or
  // for the cloud, is corrected by the search and gives the supports of the
  // unhinted search
  auto origin = pcms::test::CreateDeviceCoordinateView({0.0, 0.0}, system);
  Kokkos::View<LO*, DeviceMemorySpace> far("far", 1);
  Kokkos::deep_copy(far, mesh.nverts() - 1);
  const pcms::NeighbourhoodOptions fixed{
    .radius = 0.25, .min_supports = 8, .adapt_radius = false};
  const auto unhinted =
    centroid_sites->Neighbourhood(origin.coordinate_view, 0, fixed);
  REQUIRE(unhinted.supports_idx.size() > 0);
  pcms::test::CheckSupportResultsEquivalent(
    centroid_sites->Neighbourhood(
      origin.coordinate_view, 0, fixed,
      pcms::OwnerHints{mesh_geometry, pcms::Vertex, far}),
    unhinted);
  pcms::test::CheckSupportResultsEquivalent(
    centroid_sites->Neighbourhood(origin.coordinate_view, 0, fixed,
                                  pcms::OwnerHints{centroid_sites, 0, far}),
    unhinted);

  // a bind of the cloud through hints of the mesh geometry's layout is the
  // bind without them
  const auto hinted_bind = centroid_sites->Bind(query, vertex_owners);
  pcms::test::CheckSupportResultsEquivalent(
    dynamic_cast<const pcms::MLSBinding&>(*hinted_bind).Supports(),
    dynamic_cast<const pcms::MLSBinding&>(*centroid_sites->Bind(query))
      .Supports());
}

TEST_CASE("OmegaHGeometry::Neighbourhood searches among vertices or cells")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 6);
  auto geometry = MakeGeometry(mesh);
  const auto points = pcms::test::StandardEvalCoords2D();
  const auto coords = pcms::test::CreateDeviceCoordinateView(
    points, geometry->GetCoordinateSystem());
  const pcms::NeighbourhoodOptions options{.radius = 0.3, .min_supports = 6};
  const auto start = StartElements(mesh, coords.view);
  const auto flat =
    pcms::flatten_to_omega_h_reals(coords.coordinate_view.GetValues());
  // the sites are the element centroids: what the adjacency search finds
  // among them from the same elements
  pcms::test::CheckSupportResultsEquivalent(
    geometry->Neighbourhood(coords.coordinate_view, 2, options),
    pcms::searchNeighbors(mesh, 2, flat, 2, start, 0.09, 6, 18));
  pcms::test::CheckSupportResultsEquivalent(
    geometry->Neighbourhood(coords.coordinate_view, 0, options),
    pcms::searchNeighbors(mesh, 0, flat, 2, start, 0.09, 6, 18));
  // other entities are not sites of the adjacency search
  REQUIRE_THROWS_WITH(
    geometry->Neighbourhood(coords.coordinate_view, 1, options),
    ContainsSubstring("the sites must be the mesh vertices"));
}

TEST_CASE("Points the adjacency search cannot start from are missing")
{
  auto lib = Omega_h::Library{};
  // two triangles far apart: the search cells between them hold no element
  Omega_h::Mesh mesh(&lib);
  Omega_h::build_from_elems_and_coords(
    &mesh, OMEGA_H_SIMPLEX, 2, Omega_h::LOs({0, 1, 2, 3, 4, 5}),
    Omega_h::Reals({0, 0, 1, 0, 0, 1, 9, 9, 10, 9, 9, 10}));
  pcms::test::AddDefaultClassification(mesh);
  const auto system = pcms::csys::Cartesian::Create(2);
  auto vertices = std::make_shared<pcms::OmegaHEntityLayout>(
    pcms::test::MakeOmegaHGeometry(mesh), pcms::Vertex, 1);
  auto space = pcms::CreatePolynomialReconstructionSpace(
    pcms::CreatePointCloudGeometry(
      vertices, pcms::NeighbourhoodOptions{.radius = 0.5, .min_supports = 3}),
    pcms::MLSOptions{.degree = 1});
  auto f = space->CreateFunction<Real>();
  pcms::test::SetField(
    f, OMEGA_H_LAMBDA(Real x, Real y) { return pcms::test::linear_f(x, y); });
  const std::vector<Real> pts{0.25, 0.25, 5.5, 5.5};
  auto coords = pcms::test::CreateDeviceCoordinateView(pts, system);
  const auto bound = space->GetGeometry()->Bind(
    coords.coordinate_view,
    pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, -1.0});
  REQUIRE(bound->MissingPoints().extent(0) == 1);
  const auto values = pcms::test::EvaluateToHost(
    *space->CreatePointEvaluator<Real>(bound), f, f.GetData().GetValueBasis());
  REQUIRE_THAT(values(0, 0), WithinAbs(pcms::test::linear_f(0.25, 0.25),
                                       pcms::test::EvalTol));
  REQUIRE(values(1, 0) == -1.0);
}

TEST_CASE("An empty bind evaluates to nothing")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 2);
  auto geometry = MakeGeometry(mesh);
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {1.0, 1.0};
  grid.divisions = {2, 2};
  std::vector<std::shared_ptr<const FunctionSpace>> spaces{
    MakeSpace(geometry, 1, 1),
    CreateLagrangeSpace(
      pcms::CreateUniformGridGeometry(grid, pcms::csys::Cartesian::Create(2)),
      1, 1)};
#ifdef PCMS_ENABLE_MESHFIELDS
  spaces.push_back(CreateLagrangeSpace(
    geometry, 1, 1, LagrangeOptions{.backend = LagrangeBackend::MeshFields}));
#endif
  const Kokkos::View<Real**, DeviceMemorySpace> none("none", 0, 2);
  for (const auto& space : spaces) {
    auto f = space->CreateFunction<Real>("f");
    const auto bound = space->GetGeometry()->Bind(
      pcms::test::MakeCoords(none, space->GetCoordinateSystem()));
    REQUIRE(
      pcms::test::EvaluateToHost(*space->CreatePointEvaluator<Real>(bound), f,
                                 f.GetData().GetValueBasis())
        .extent(0) == 0);
  }
}

TEST_CASE("DOFHolderOwners gives the entity of each DOF holder")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 3);
  auto geometry = MakeGeometry(mesh);
  const auto entities = [&](const pcms::FunctionSpace& space, int dim) {
    const auto owners = pcms::DOFHolderOwners(*space.GetLayout());
    REQUIRE(owners.geometry == geometry);
    REQUIRE(owners.entity_dim == dim);
    const auto ids =
      Kokkos::create_mirror_view_and_copy(HostMemorySpace(), owners.entities);
    REQUIRE(ids.extent(0) == static_cast<size_t>(mesh.nents(dim)));
    for (LO i = 0; i < mesh.nents(dim); ++i) {
      REQUIRE(ids(i) == i);
    }
  };
  entities(*MakeSpace(geometry, 1, 1), pcms::Vertex);
  entities(*MakeSpace(geometry, 0, 1), pcms::Face);
#ifdef PCMS_ENABLE_MESHFIELDS
  entities(
    *CreateLagrangeSpace(
      geometry, 1, 1, LagrangeOptions{.backend = LagrangeBackend::MeshFields}),
    pcms::Vertex);
#endif
  // with ghosts, DOF holders are not the entities
  Omega_h::Write<Omega_h::I8> mask(mesh.nverts(), 1);
  Kokkos::deep_copy(Kokkos::View<Omega_h::I8*>(mask.data(), 1), 0);
  const auto ghosted = CreateLagrangeSpace(
    pcms::CreateOmegaHGeometry(
      mesh, pcms::csys::Cartesian::Create(2),
      pcms::OmegaHGeometryOptions{.owned_mask = Omega_h::Read(mask)}),
    1, 1, LagrangeOptions{.backend = LagrangeBackend::OmegaH});
  REQUIRE(ghosted->GetLayout()->HasGhosts());
  REQUIRE(pcms::DOFHolderOwners(*ghosted->GetLayout()).entities.extent(0) == 0);
  // without a geometry there are no owners to name: a bare point cloud
  // (whose DOF holders are not entities of a mesh) and an empty layout
  const pcms::PointCloudLayout cloud(2, Kokkos::View<Real**>("sites", 3, 2),
                                     pcms::csys::Cartesian::Create(2));
  REQUIRE(cloud.GetGeometry() == nullptr);
  REQUIRE(pcms::DOFHolderOwners(cloud).entities.extent(0) == 0);
  const pcms::EmptyFieldLayout none(pcms::csys::Cartesian::Create(2));
  REQUIRE(none.GetGeometry() == nullptr);
  REQUIRE(pcms::DOFHolderOwners(none).entities.extent(0) == 0);
  // the hint names the geometry of the layout, so a layout on another
  // geometry of the same mesh names that one
  const auto other = MakeSpace(MakeGeometry(mesh), 1, 1);
  REQUIRE(
    pcms::DOFHolderOwners(*other->GetLayout()).geometry !=
    pcms::DOFHolderOwners(*MakeSpace(geometry, 1, 1)->GetLayout()).geometry);
}

TEST_CASE("Hinted neighbourhoods equal searched ones inside the mesh")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 6);
  auto geometry = MakeGeometry(mesh);
  const pcms::NeighbourhoodOptions options{.radius = 0.3, .min_supports = 6};
  for (const int dim : {pcms::Vertex, pcms::Face}) {
    CAPTURE(dim);
    const auto points = pcms::ConvertCoordsTo2D(
      pcms::get_entity_centroids(mesh, dim), mesh.nents(dim), 2);
    const auto coords =
      pcms::test::MakeCoords(points, geometry->GetCoordinateSystem());
    Kokkos::View<LO*, DeviceMemorySpace> ids("ids", mesh.nents(dim));
    pcms::iota_view<LO>(ids);
    pcms::test::CheckSupportResultsEquivalent(
      geometry->Neighbourhood(coords, 0, options,
                              pcms::OwnerHints{geometry, dim, ids}),
      geometry->Neighbourhood(coords, 0, options));
  }
}

TEST_CASE("Binding with owners keeps the policy")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 4);
  auto geometry = MakeGeometry(mesh);
  auto space = MakeSpace(geometry, 1, 1);
  auto f = space->CreateFunction<Real>("f");
  // the second point is outside; without a usable hint it is searched
  const std::vector<Real> pts{0.5, 0.5, 2.0, 2.0};
  auto coords = pcms::test::CreateDeviceCoordinateView(
    pts, geometry->GetCoordinateSystem());
  const auto bound =
    geometry->Bind(coords.coordinate_view, pcms::OwnerHints{},
                   pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, -1.0});
  REQUIRE(bound->Policy().fill_value == -1.0);
  REQUIRE(bound->MissingPoints().extent(0) == 1);
  const auto values = pcms::test::EvaluateToHost(
    *space->CreatePointEvaluator<Real>(bound), f, f.GetData().GetValueBasis());
  REQUIRE(values(1, 0) == -1.0);
}

TEST_CASE("Stale owner hints bind like the search", "[owner_hints]")
{
  auto lib = Omega_h::Library{};
  for (const int dim : {2, 3}) {
    CAPTURE(dim);
    auto mesh = dim == 2 ? pcms::test::BuildSquareMesh(lib, 6)
                         : pcms::test::BuildUnitCube(lib, 3);
    auto geometry = MakeMeshGeometry(mesh);
    const int n = 200;
    const auto pts = RandomPoints(dim, n, 17);
    const auto coords = pcms::test::CreateDeviceCoordinateView(
      pts, geometry->GetCoordinateSystem(), dim);
    const pcms::OutOfBoundsPolicy fill{pcms::OutOfBoundsMode::FILL, -1.0};
    const auto cold = Summarize(geometry->Bind(coords.coordinate_view, fill));
    // some of the random points are outside, some inside
    REQUIRE(!cold.found.empty());
    REQUIRE(!cold.missing.empty());

    const auto bind_with = [&](int entity_dim, const std::vector<LO>& ids) {
      return Summarize(geometry->Bind(
        coords.coordinate_view,
        pcms::OwnerHints{geometry, entity_dim, ToDevice(ids)}, fill));
    };
    // the elements the search finds, 0 for the points outside
    std::vector<LO> right(n, 0);
    for (size_t k = 0; k < cold.found.size(); ++k) {
      // inside the mesh, away from shared entities: the record is an element
      REQUIRE(cold.records[k].dim == dim);
      right[cold.found[k]] = cold.records[k].entity;
    }
    {
      INFO("the right elements");
      CheckSameBind(bind_with(dim, right), cold);
    }
    {
      INFO("the elements shuffled");
      auto stale = right;
      std::shuffle(stale.begin(), stale.end(), std::mt19937(3));
      CheckSameBind(bind_with(dim, stale), cold);
    }
    {
      INFO("one element for every point");
      CheckSameBind(bind_with(dim, std::vector<LO>(n, mesh.nelems() - 1)),
                    cold);
    }
    {
      INFO("vertices, right or wrong");
      std::mt19937 gen(9);
      std::uniform_int_distribution<LO> any(0, mesh.nverts() - 1);
      std::vector<LO> vertices(n);
      for (auto& v : vertices) {
        v = any(gen);
      }
      CheckSameBind(bind_with(pcms::Vertex, vertices), cold);
    }
    for (int hint_dim = 1; hint_dim < dim; ++hint_dim) {
      INFO("entities of dimension " << hint_dim << ", right or wrong");
      std::mt19937 gen(13);
      std::uniform_int_distribution<LO> any(0, mesh.nents(hint_dim) - 1);
      std::vector<LO> entities(n);
      for (auto& e : entities) {
        e = any(gen);
      }
      CheckSameBind(bind_with(hint_dim, entities), cold);
    }
  }
}

TEST_CASE("Owner hints of points outside the mesh follow the policy",
          "[owner_hints]")
{
  auto lib = Omega_h::Library{};
  for (const int dim : {2, 3}) {
    CAPTURE(dim);
    auto mesh = dim == 2 ? pcms::test::BuildSquareMesh(lib, 4)
                         : pcms::test::BuildUnitCube(lib, 2);
    auto geometry = MakeMeshGeometry(mesh);
    // an inside point and a point beyond the mesh, both hinted with element 0
    std::vector<Real> pts(static_cast<size_t>(2 * dim), 0.0);
    for (int d = 0; d < dim; ++d) {
      pts[d] = 0.45;
      pts[dim + d] = 2.5;
    }
    const auto coords = pcms::test::CreateDeviceCoordinateView(
      pts, geometry->GetCoordinateSystem(), dim);
    const pcms::OwnerHints hint{geometry, dim, ToDevice({0, 0})};

    const pcms::OutOfBoundsPolicy fill{pcms::OutOfBoundsMode::FILL, -1.0};
    const auto filled = geometry->Bind(coords.coordinate_view, hint, fill);
    REQUIRE(ToVector(filled->FoundPoints()) == std::vector<LO>{0});
    REQUIRE(ToVector(filled->MissingPoints()) == std::vector<LO>{1});
    CheckSameBind(Summarize(filled),
                  Summarize(geometry->Bind(coords.coordinate_view, fill)));

    REQUIRE_THROWS_WITH(geometry->Bind(coords.coordinate_view, hint),
                        ContainsSubstring("outside the source's domain"));
    REQUIRE_THROWS_WITH(geometry->Bind(coords.coordinate_view),
                        ContainsSubstring("outside the source's domain"));
  }
}

TEST_CASE("DOF holder hints of an order-1 layout bind like the search",
          "[owner_hints]")
{
  auto lib = Omega_h::Library{};
  for (const int dim : {2, 3}) {
    CAPTURE(dim);
    auto mesh = dim == 2 ? pcms::test::BuildSquareMesh(lib, 5)
                         : pcms::test::BuildUnitCube(lib, 3);
    auto geometry = MakeMeshGeometry(mesh);
    auto space = CreateLagrangeSpace(
      geometry, 1, 1, LagrangeOptions{.backend = LagrangeBackend::OmegaH});
    const auto owners = pcms::DOFHolderOwners(*space->GetLayout());
    REQUIRE(owners.entity_dim == pcms::Vertex);
    const auto dofs = space->GetLayout()->GetDOFHolderCoordinates();
    // the DOF holders are the mesh vertices, on the boundary of the domain:
    // the points sit on the tolerance band of the elements they are hinted in
    CheckSameBind(Summarize(geometry->Bind(dofs, owners)),
                  Summarize(geometry->Bind(dofs)));
    pcms::test::CheckSupportResultsEquivalent(
      geometry->Neighbourhood(dofs, 0, {.radius = 0.4, .min_supports = 6},
                              owners),
      geometry->Neighbourhood(dofs, 0, {.radius = 0.4, .min_supports = 6}));
  }
}

TEST_CASE("Stale owner hints give the neighbourhood of the search",
          "[owner_hints]")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 6);
  auto geometry = MakeGeometry(mesh);
  const pcms::NeighbourhoodOptions options{.radius = 0.3, .min_supports = 6};
  const int n = 100;
  const auto pts = RandomPoints(2, n, 5);
  const auto coords = pcms::test::CreateDeviceCoordinateView(
    pts, geometry->GetCoordinateSystem());
  const auto cold = geometry->Neighbourhood(coords.coordinate_view, 0, options);
  std::mt19937 gen(2);
  for (const int dim : {pcms::Vertex, pcms::Face}) {
    CAPTURE(dim);
    std::uniform_int_distribution<LO> any(0, mesh.nents(dim) - 1);
    std::vector<LO> ids(n);
    for (auto& id : ids) {
      id = any(gen);
    }
    pcms::test::CheckSupportResultsEquivalent(
      geometry->Neighbourhood(coords.coordinate_view, 0, options,
                              pcms::OwnerHints{geometry, dim, ToDevice(ids)}),
      cold);
  }
}

namespace
{

/// The centroids of all entities of `mesh` (vertices, edge midpoints, face and
/// element centres): points on every kind of shared entity, interleaved.
std::vector<Real> EntityCentroids(Omega_h::Mesh& mesh)
{
  std::vector<Real> pts;
  for (int d = 0; d <= mesh.dim(); ++d) {
    const auto c = Omega_h::HostRead<Real>(pcms::get_entity_centroids(mesh, d));
    pts.insert(pts.end(), c.data(), c.data() + c.size());
  }
  return pts;
}

/// P1 and P0 values at `pts`, from a cold bind and, after the first
/// evaluation is done with, from a rebind started from the last records.
struct ColdAndWarm
{
  Kokkos::View<Real**, HostMemorySpace> cold;
  Kokkos::View<Real**, HostMemorySpace> warm;
};

ColdAndWarm EvaluateColdAndWarm(
  const std::shared_ptr<const pcms::Geometry>& geometry,
  const std::shared_ptr<const FunctionSpace>& space, const pcms::Field<Real>& f,
  const std::vector<Real>& pts)
{
  const int dim = geometry->GetCoordinateSystem()->Dimension();
  const auto coords = pcms::test::CreateDeviceCoordinateView(
    pts, geometry->GetCoordinateSystem(), dim);
  std::shared_ptr<const pcms::BoundSource> bound =
    geometry->Bind(coords.coordinate_view);
  ColdAndWarm out;
  out.cold = pcms::test::EvaluateToHost(
    *space->CreatePointEvaluator<Real>(bound), f, f.GetData().GetValueBasis());
  bound = geometry->Rebind(std::move(bound), coords.coordinate_view);
  out.warm = pcms::test::EvaluateToHost(
    *space->CreatePointEvaluator<Real>(bound), f, f.GetData().GetValueBasis());
  return out;
}

} // namespace

TEST_CASE("Omega_h evaluators use the entity the search finds",
          "[omega_h_entities]")
{
  auto lib = Omega_h::Library{};
  for (const int dim : {2, 3}) {
    CAPTURE(dim);
    auto mesh = dim == 2 ? pcms::test::BuildSquareMesh(lib, 4)
                         : pcms::test::BuildUnitCube(lib, 2);
    auto geometry = MakeMeshGeometry(mesh);
    const auto pts = EntityCentroids(mesh);
    const auto n = static_cast<LO>(pts.size() / dim);

    // P1 interpolates on the entity itself: a linear field is reproduced
    // exactly at vertices, edge midpoints, face and element centres, and a cold
    // and a warm bind agree
    auto p1 = MakeSpace(geometry, 1, 1);
    auto f1 = p1->CreateFunction<Real>("f1");
    const auto linear = [dim](const std::vector<Real>& x, LO i) {
      const auto p = x.data() + static_cast<size_t>(dim) * i;
      return 1.0 + 2.0 * p[0] - 3.0 * p[1] + (dim == 3 ? 0.5 * p[2] : 0.0);
    };
    if (dim == 2) {
      pcms::test::SetField(
        f1, OMEGA_H_LAMBDA(Real x, Real y) { return 1.0 + 2.0 * x - 3.0 * y; });
    } else {
      pcms::test::SetField(
        f1, OMEGA_H_LAMBDA(Real x, Real y, Real z) {
          return 1.0 + 2.0 * x - 3.0 * y + 0.5 * z;
        });
    }
    const auto values = EvaluateColdAndWarm(geometry, p1, f1, pts);
    for (LO i = 0; i < n; ++i) {
      CAPTURE(i);
      REQUIRE_THAT(values.cold(i, 0),
                   WithinAbs(linear(pts, i), pcms::test::EvalTol));
      REQUIRE_THAT(values.warm(i, 0),
                   WithinAbs(values.cold(i, 0), pcms::test::ExactTol));
    }

    // P0 takes the value of the lowest element adjacent to a shared entity
    auto p0 = MakeSpace(geometry, 0, 1);
    auto f0 = p0->CreateFunction<Real>("f0");
    std::vector<Real> element_values(mesh.nelems());
    for (LO e = 0; e < mesh.nelems(); ++e) {
      element_values[e] = 10.0 + e;
    }
    f0.SetDOFHolderDataUncheckedHost(
      pcms::Rank2View<const Real, HostMemorySpace>(element_values.data(),
                                                   mesh.nelems(), 1));
    const auto p0_values = EvaluateColdAndWarm(geometry, p0, f0, pts);
    // the entities of the centroids, in the order of EntityCentroids
    LO i = 0;
    for (int d = 0; d <= dim; ++d) {
      const auto up = d < dim ? mesh.ask_up(d, dim) : Omega_h::Adj{};
      for (LO entity = 0; entity < mesh.nents(d); ++entity, ++i) {
        CAPTURE(d, entity);
        LO lowest = d == dim ? entity : mesh.nelems();
        if (d < dim) {
          const auto a2ab = Omega_h::HostRead<LO>(up.a2ab);
          const auto ab2b = Omega_h::HostRead<LO>(up.ab2b);
          for (LO a = a2ab[entity]; a < a2ab[entity + 1]; ++a) {
            lowest = std::min(lowest, ab2b[a]);
          }
        }
        REQUIRE(p0_values.cold(i, 0) == 10.0 + lowest);
        REQUIRE(p0_values.warm(i, 0) == p0_values.cold(i, 0));
      }
    }
  }
}

#ifdef PCMS_ENABLE_MESHFIELDS
TEST_CASE("MeshFields values at shared entities equal the Omega_h P1 values",
          "[omega_h_entities]")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 4);
  auto geometry = MakeMeshGeometry(mesh);
  const auto pts = EntityCentroids(mesh);
  const auto n = static_cast<LO>(pts.size() / 2);
  auto omega_h = MakeSpace(geometry, 1, 1);
  auto mesh_fields = CreateLagrangeSpace(
    geometry, 1, 1, LagrangeOptions{.backend = LagrangeBackend::MeshFields});
  auto f_omega_h = omega_h->CreateFunction<Real>("f_omega_h");
  auto f_mesh_fields = mesh_fields->CreateFunction<Real>("f_mesh_fields");
  const auto quadratic_ish = OMEGA_H_LAMBDA(Real x, Real y)
  {
    return 1.0 + 2.0 * x - 3.0 * y + 0.25 * x * y;
  };
  pcms::test::SetField(f_omega_h, quadratic_ish);
  pcms::test::SetField(f_mesh_fields, quadratic_ish);
  const auto coords = pcms::test::CreateDeviceCoordinateView(
    pts, geometry->GetCoordinateSystem());
  const auto bound = geometry->Bind(coords.coordinate_view);
  const auto a =
    pcms::test::EvaluateToHost(*omega_h->CreatePointEvaluator<Real>(bound),
                               f_omega_h, f_omega_h.GetData().GetValueBasis());
  const auto b = pcms::test::EvaluateToHost(
    *mesh_fields->CreatePointEvaluator<Real>(bound), f_mesh_fields,
    f_mesh_fields.GetData().GetValueBasis());
  for (LO i = 0; i < n; ++i) {
    CAPTURE(i);
    REQUIRE_THAT(b(i, 0), WithinAbs(a(i, 0), pcms::test::EvalTol));
  }
}
#endif

#ifdef PCMS_ENABLE_MESHFIELDS
TEST_CASE("MeshFields evaluates P2 fields through the geometry's bind")
{
  // P2 reproduces a quadratic exactly, so any change in the element
  // coordinates MeshFields receives shows up here.
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildSquareMesh(lib, 6);
  auto geometry = MakeGeometry(mesh);
  auto p2 = CreateLagrangeSpace(
    geometry, 2, 1, LagrangeOptions{.backend = LagrangeBackend::MeshFields});
  auto f2 = p2->CreateFunction<Real>("f2");
  pcms::test::SetField(
    f2, OMEGA_H_LAMBDA(Real x, Real y) { return Quadratic(x, y); });

  const auto pts = pcms::test::StandardEvalCoords2D();
  const auto bound = pcms::test::Bind(*geometry, pts);
  pcms::test::CheckEvaluation(
    *p2->CreatePointEvaluator<Real>(bound), f2, pts,
    OMEGA_H_LAMBDA(Real x, Real y) { return Quadratic(x, y); },
    pcms::test::ExactTol);
}
#endif
