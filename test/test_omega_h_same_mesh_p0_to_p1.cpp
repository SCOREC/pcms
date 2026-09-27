#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>

#include <pcms/transfer/mass_matrix_type.hpp>
#include <pcms/transfer/omega_h_conservative_projection.hpp>
#include <pcms/transfer/omega_h_form_integrator_utils.hpp>
#include <pcms/utility/arrays.h>
#include "field_test_utils.h"

#include <array>
#include <cmath>
#include <random>
#include <utility>
#include <vector>

// Same-mesh P0->P1 projection has closed-form answers, which makes it a sharp
// test of the intersection quadrature: every target element is clipped against
// itself, the pieces must tile it exactly, and the lumped projection must equal
// the element-nodal average to roundoff. Both dimensions run through the same
// decomposition code except for the boundary walk, so a geometric slip in the
// shared path shows up in both.

namespace
{

// Tets taken verbatim from the PUMI-Tally wistell-D mesh. Clipping each one
// against itself with r3d gives signed distances of ~1e-13 on every face plane,
// and r3d then interpolates new vertices at O(1) fractions along the edges
// between two on-plane vertices. The resulting polytope has the right volume
// but a folded face graph whose triangles come out in both orientations, so a
// star decomposition that takes |volume| per piece over-counts (by 28.6%,
// 1.7% and 9.0% respectively). Any element sharing a face plane with a source
// element can hit this; the same-mesh case hits it for every element.
const std::vector<double> kFoldingTetCoords = {
  30.198030285537243,
  -936.1389388516545,
  483.16848456836186,
  20.132020190358162,
  -979.0909772444453,
  485.85648258043625,
  0.0,
  -966.3369691371918,
  483.16848456836186,
  15.099015142768621,
  -951.2379539944232,
  468.06946942559324,

  -452.97045428305864,
  392.57439371198416,
  241.58424228406392,
  -457.87167295326776,
  423.1168047352564,
  234.21738729666106,
  -482.6122426477242,
  412.69480663896707,
  261.65815839972544,
  -460.4863470953233,
  422.81883403751397,
  272.4806339440611,

  -620.6927499466186,
  1017.0028114430277,
  82.31454644833995,
  -634.1586359962821,
  1056.9310599938033,
  90.59409085637769,
  -637.7105881441516,
  1022.8656314424136,
  102.40466297536494,
  -603.9377425704466,
  1022.2805800406807,
  87.57395423730232,
};

template <int Dim>
constexpr int kVertsPerElem = Dim + 1;

// Signed measure of the simplex whose Dim+1 vertices are packed in `c`.
template <int Dim>
double SignedSimplexMeasure(const double* c)
{
  if constexpr (Dim == 2) {
    const double bx = c[2] - c[0], by = c[3] - c[1];
    const double dx = c[4] - c[0], dy = c[5] - c[1];
    return 0.5 * (bx * dy - by * dx);
  } else {
    const double b[3] = {c[3] - c[0], c[4] - c[1], c[5] - c[2]};
    const double d[3] = {c[6] - c[0], c[7] - c[1], c[8] - c[2]};
    const double e[3] = {c[9] - c[0], c[10] - c[1], c[11] - c[2]};
    return (b[0] * (d[1] * e[2] - d[2] * e[1]) -
            b[1] * (d[0] * e[2] - d[2] * e[0]) +
            b[2] * (d[0] * e[1] - d[1] * e[0])) /
           6.0;
  }
}

// Random, positively oriented, non-sliver simplices at the coordinate scale of
// a tokamak/stellarator mesh (~1e3), so shared-plane roundoff is
// representative.
template <int Dim>
std::vector<double> RandomSimplexCoords(int n, unsigned seed)
{
  constexpr int nv = kVertsPerElem<Dim>;
  constexpr double min_measure = (Dim == 3) ? 1e3 : 1e2;
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> base(-1900.0, 1900.0);
  std::uniform_real_distribution<double> offset(-50.0, 50.0);
  std::vector<double> coords;
  coords.reserve(Dim * nv * n);
  while (static_cast<int>(coords.size()) < Dim * nv * n) {
    double c[Dim * nv];
    for (int d = 0; d < Dim; ++d) {
      c[d] = base(gen);
    }
    for (int v = 1; v < nv; ++v) {
      for (int d = 0; d < Dim; ++d) {
        c[Dim * v + d] = c[d] + offset(gen);
      }
    }
    const double measure = SignedSimplexMeasure<Dim>(c);
    if (std::fabs(measure) < min_measure) {
      continue;
    }
    if (measure < 0.0) {
      for (int d = 0; d < Dim; ++d) {
        std::swap(c[Dim * (nv - 2) + d], c[Dim * (nv - 1) + d]);
      }
    }
    coords.insert(coords.end(), c, c + Dim * nv);
  }
  return coords;
}

// One mesh holding mutually disjoint simplices.
template <int Dim>
Omega_h::Mesh BuildDisjointSimplices(Omega_h::Library& lib,
                                     const std::vector<double>& coords)
{
  constexpr int nv = kVertsPerElem<Dim>;
  const int nelems = static_cast<int>(coords.size()) / (Dim * nv);
  Omega_h::HostWrite<Omega_h::LO> ev2v(nv * nelems);
  for (int i = 0; i < nv * nelems; ++i) {
    ev2v[i] = i;
  }
  Omega_h::HostWrite<Omega_h::Real> x(static_cast<Omega_h::LO>(coords.size()));
  for (std::size_t i = 0; i < coords.size(); ++i) {
    x[static_cast<Omega_h::LO>(i)] = coords[i];
  }
  Omega_h::Mesh mesh(&lib);
  Omega_h::build_from_elems_and_coords(&mesh, OMEGA_H_SIMPLEX, Dim,
                                       Omega_h::LOs(ev2v.write()),
                                       Omega_h::Reals(x.write()));
  pcms::test::AddDefaultClassification(mesh);
  return mesh;
}

// Unit box tessellation moved to non-dyadic coordinates of magnitude ~1e3. On
// the unit box every on-plane distance is exactly zero and r3d never folds;
// this is what a real mesh looks like to the clipper.
template <int Dim>
Omega_h::Mesh BuildScaledBox(Omega_h::Library& lib, int n)
{
  auto mesh = (Dim == 3) ? Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0,
                                              1.0, 1.0, n, n, n)
                         : Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0,
                                              1.0, 0.0, n, n, 0);
  const auto coords = Omega_h::HostRead<Omega_h::Real>(mesh.coords());
  Omega_h::HostWrite<Omega_h::Real> scaled(coords.size());
  constexpr double scale = 937.31;
  constexpr std::array<double, 3> shift = {-1894.78, 211.7, -598.52};
  for (Omega_h::LO v = 0; v < mesh.nverts(); ++v) {
    for (int d = 0; d < Dim; ++d) {
      scaled[Dim * v + d] = coords[Dim * v + d] * scale + shift[d];
    }
  }
  mesh.set_coords(Omega_h::Reals(scaled.write()));
  return mesh;
}

// Sum of the sub-simplex measures handed to the integrand for target element
// e when the source mesh is the target mesh itself and the intersection map is
// the identity. Exactly |e| for a correct decomposition.
template <int Dim>
Omega_h::HostRead<Omega_h::Real> SelfDecompositionMeasures(Omega_h::Mesh& mesh)
{
  const int nelems = mesh.nelems();
  const auto coords = mesh.coords();
  const auto e2n = mesh.ask_down(Dim, Omega_h::VERT).ab2b;
  Omega_h::Write<Omega_h::LO> offsets(nelems + 1);
  Omega_h::Write<Omega_h::LO> indices(nelems);
  Omega_h::parallel_for(nelems + 1, OMEGA_H_LAMBDA(int i) { offsets[i] = i; });
  Omega_h::parallel_for(nelems, OMEGA_H_LAMBDA(int i) { indices[i] = i; });
  const pcms::IntersectionResults self{Omega_h::LOs(offsets),
                                       Omega_h::LOs(indices)};

  Omega_h::Write<Omega_h::Real> sums(nelems, 0.0);
  Omega_h::parallel_for(
    nelems, OMEGA_H_LAMBDA(int e) {
      Omega_h::Real acc = 0.0;
      pcms::detail::ForEachIntersectionSubsimplex<Dim>(
        e, self, coords, coords, e2n, e2n,
        [&](const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>&, int,
            Omega_h::Real measure) { acc += measure; });
      sums[e] = acc;
    });
  return Omega_h::HostRead<Omega_h::Real>(sums);
}

template <int Dim>
void RequireSelfDecompositionExact(Omega_h::Mesh& mesh)
{
  const auto sums = SelfDecompositionMeasures<Dim>(mesh);
  const auto sizes =
    Omega_h::HostRead<Omega_h::Real>(Omega_h::measure_elements_real(&mesh));
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    REQUIRE(sums[e] == Catch::Approx(std::fabs(sizes[e])).epsilon(1e-12));
  }
}

// The lumped P1 mass matrix is diagonal with M_ii = sum_{e in i} |e|/(Dim+1)
// and the exact P0 load vector is b_i = sum_{e in i} |e|/(Dim+1) f_e, so the
// lumped projection of a P0 field onto P1 on the same mesh is the
// measure-weighted average of the incident element values. Nothing here is
// approximate; the operator must reproduce it to roundoff.
template <int Dim>
std::vector<double> ElementNodalAverage(Omega_h::Mesh& mesh,
                                        const std::vector<double>& elem_values)
{
  constexpr int nv = kVertsPerElem<Dim>;
  const auto sizes =
    Omega_h::HostRead<Omega_h::Real>(Omega_h::measure_elements_real(&mesh));
  const auto e2v = Omega_h::HostRead<Omega_h::LO>(mesh.ask_elem_verts());
  std::vector<double> numerator(mesh.nverts(), 0.0);
  std::vector<double> denominator(mesh.nverts(), 0.0);
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    const double w = std::fabs(sizes[e]) / nv;
    for (int k = 0; k < nv; ++k) {
      const auto v = e2v[nv * e + k];
      numerator[v] += w * elem_values[e];
      denominator[v] += w;
    }
  }
  for (Omega_h::LO v = 0; v < mesh.nverts(); ++v) {
    numerator[v] /= denominator[v];
  }
  return numerator;
}

template <typename View>
std::vector<double> ToVector(const View& values)
{
  std::vector<double> out(values.size());
  for (std::size_t i = 0; i < out.size(); ++i) {
    out[i] = values[i];
  }
  return out;
}

// Zeroth and first moments (integrals of 1 and of each coordinate times the
// field) of a P0 element field and of a P1 vertex field, exact on simplices.
template <int Dim>
std::array<double, Dim + 1> P0Moments(Omega_h::Mesh& mesh,
                                      const std::vector<double>& f)
{
  constexpr int nv = kVertsPerElem<Dim>;
  const auto sizes =
    Omega_h::HostRead<Omega_h::Real>(Omega_h::measure_elements_real(&mesh));
  const auto e2v = Omega_h::HostRead<Omega_h::LO>(mesh.ask_elem_verts());
  const auto x = Omega_h::HostRead<Omega_h::Real>(mesh.coords());
  std::array<double, Dim + 1> m{};
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    const double measure = std::fabs(sizes[e]);
    m[0] += measure * f[e];
    for (int d = 0; d < Dim; ++d) {
      double centroid = 0.0;
      for (int k = 0; k < nv; ++k) {
        centroid += x[Dim * e2v[nv * e + k] + d];
      }
      m[1 + d] += measure * f[e] * centroid / nv;
    }
  }
  return m;
}

template <int Dim>
std::array<double, Dim + 1> P1Moments(Omega_h::Mesh& mesh,
                                      const std::vector<double>& u)
{
  constexpr int nv = kVertsPerElem<Dim>;
  // int_T phi_i phi_j = |T| (1 + delta_ij) / ((Dim+1)(Dim+2)), so
  // int_T x u_h = |T| / ((Dim+1)(Dim+2)) * ((sum_i u_i)(sum_j x_j) +
  //                                          sum_i u_i x_i).
  constexpr double mass_factor = 1.0 / ((Dim + 1) * (Dim + 2));
  const auto sizes =
    Omega_h::HostRead<Omega_h::Real>(Omega_h::measure_elements_real(&mesh));
  const auto e2v = Omega_h::HostRead<Omega_h::LO>(mesh.ask_elem_verts());
  const auto x = Omega_h::HostRead<Omega_h::Real>(mesh.coords());
  std::array<double, Dim + 1> m{};
  for (Omega_h::LO e = 0; e < mesh.nelems(); ++e) {
    const double measure = std::fabs(sizes[e]);
    double sum_u = 0.0;
    for (int k = 0; k < nv; ++k) {
      sum_u += u[e2v[nv * e + k]];
    }
    m[0] += measure * sum_u / nv;
    for (int d = 0; d < Dim; ++d) {
      double sum_x = 0.0;
      double sum_ux = 0.0;
      for (int k = 0; k < nv; ++k) {
        const double xk = x[Dim * e2v[nv * e + k] + d];
        sum_x += xk;
        sum_ux += u[e2v[nv * e + k]] * xk;
      }
      m[1 + d] += measure * mass_factor * (sum_u * sum_x + sum_ux);
    }
  }
  return m;
}

// A P0 field that varies from element to element at the ~10% level, so a
// misattributed or double-counted piece of any element is visible. The arity
// of operator() selects the dimension in pcms::test::SetField.
template <int Dim>
struct SourceFunction;

template <>
struct SourceFunction<2>
{
  KOKKOS_INLINE_FUNCTION pcms::Real operator()(pcms::Real x, pcms::Real y) const
  {
    return 10.0 + Kokkos::sin(0.011 * x) + 0.5 * Kokkos::cos(0.017 * y);
  }
};

template <>
struct SourceFunction<3>
{
  KOKKOS_INLINE_FUNCTION pcms::Real operator()(pcms::Real x, pcms::Real y,
                                               pcms::Real z) const
  {
    return 10.0 + Kokkos::sin(0.011 * x) + 0.5 * Kokkos::cos(0.017 * y) +
           0.001 * z;
  }
};

template <int Dim>
std::vector<double> ProjectSameMesh(Omega_h::Mesh& mesh,
                                    pcms::MassMatrixType mass_type,
                                    std::vector<double>& source_values)
{
  auto source_space = pcms::test::MakeP0Space(mesh);
  auto target_space = pcms::test::MakeP1Space(mesh);
  auto source = source_space->CreateFunction<pcms::Real>();
  auto target = target_space->CreateFunction<pcms::Real>();
  pcms::test::SetField(source, SourceFunction<Dim>{});

  pcms::OmegaHConservativeProjection projection(*source_space, *target_space,
                                                mass_type);
  projection.Apply(source, target);

  source_values =
    ToVector(pcms::FlattenToRank1View(source.GetDOFHolderDataHost()));
  return ToVector(pcms::FlattenToRank1View(target.GetDOFHolderDataHost()));
}

template <int Dim>
void RequireLumpedEqualsNodalAverage(Omega_h::Mesh& mesh)
{
  std::vector<double> f;
  const auto u = ProjectSameMesh<Dim>(mesh, pcms::MassMatrixType::Lumped, f);
  REQUIRE(static_cast<Omega_h::LO>(u.size()) == mesh.nverts());
  const auto expected = ElementNodalAverage<Dim>(mesh, f);
  for (Omega_h::LO v = 0; v < mesh.nverts(); ++v) {
    REQUIRE(u[v] == Catch::Approx(expected[v]).epsilon(1e-12));
  }
  REQUIRE(P1Moments<Dim>(mesh, u)[0] ==
          Catch::Approx(P0Moments<Dim>(mesh, f)[0]).epsilon(1e-12));
}

template <int Dim>
void RequireConsistentPreservesMoments(Omega_h::Mesh& mesh)
{
  std::vector<double> f;
  const auto u =
    ProjectSameMesh<Dim>(mesh, pcms::MassMatrixType::Consistent, f);
  const auto m0 = P0Moments<Dim>(mesh, f);
  const auto m1 = P1Moments<Dim>(mesh, u);
  for (int k = 0; k <= Dim; ++k) {
    REQUIRE(m1[k] == Catch::Approx(m0[k]).epsilon(1e-9));
  }
}

// The polytope r3d produced when clipping tet 38125 above against itself in
// an optimized build, captured verbatim (positions to 1e-6, graph exact). The
// four original vertices appear three or four times each and three new
// vertices sit at 1/2 and 2/3 along edges; two of the faces the walk emits are
// mirror images that must cancel. Whether a given build reproduces this fold
// from the coordinates depends on the compiler's floating-point contraction,
// so the polytope is pinned here rather than regenerated.
r3d::Polytope<3> FoldedSelfClipPolytope()
{
  constexpr int nverts = 14;
  constexpr double pos[nverts][3] = {
    {0.000000, -966.336969, 483.168485},  {20.132020, -979.090977, 485.856483},
    {7.549508, -958.787462, 475.618977},  {0.000000, -966.336969, 483.168485},
    {30.198030, -936.138939, 483.168485}, {15.099015, -951.237954, 468.069469},
    {22.648523, -943.688446, 475.618977}, {30.198030, -936.138939, 483.168485},
    {0.000000, -966.336969, 483.168485},  {20.132020, -979.090977, 485.856483},
    {20.132020, -979.090977, 485.856483}, {13.421347, -974.839641, 484.960483},
    {20.132020, -979.090977, 485.856483}, {30.198030, -936.138939, 483.168485},
  };
  constexpr int pnbrs[nverts][3] = {
    {3, 2, 8},  {4, 9, 10}, {0, 5, 11},  {0, 4, 6},   {1, 7, 3},
    {2, 6, 12}, {3, 7, 5},  {4, 13, 6},  {0, 11, 9},  {1, 8, 10},
    {1, 9, 13}, {2, 12, 8}, {5, 13, 11}, {7, 10, 12},
  };
  r3d::Polytope<3> poly;
  poly.nverts = nverts;
  for (int v = 0; v < nverts; ++v) {
    for (int d = 0; d < 3; ++d) {
      poly.verts[v].pos[d] = pos[v][d];
      poly.verts[v].pnbrs[d] = pnbrs[v][d];
    }
  }
  return poly;
}

// The 2D analogue, built by hand: the triangle A B C with the cycle detouring
// from B out to D (on edge BC) and straight back, A -> B -> D -> B -> C -> A.
// This is the pattern r3d leaves after clipping on a plane through B: a
// spliced edge vertex reached twice. The two spur edges B->D and D->B lift to
// mirrored triangles whose signed areas cancel, so the cycle still encloses
// the area of ABC; summing |area| per piece adds twice the spur triangle.
r3d::Polytope<2> FoldedPolygon()
{
  constexpr int nverts = 5;
  constexpr double pos[nverts][2] = {
    {0.0, 0.0}, {4.0, 0.0}, {2.0, 1.5}, {4.0, 0.0}, {0.0, 3.0}};
  r3d::Polytope<2> poly;
  poly.nverts = nverts;
  for (int v = 0; v < nverts; ++v) {
    poly.verts[v].pos[0] = pos[v][0];
    poly.verts[v].pos[1] = pos[v][1];
    poly.verts[v].pnbrs[0] = (v + 1) % nverts;
    poly.verts[v].pnbrs[1] = (v + nverts - 1) % nverts;
  }
  return poly;
}

template <int Dim>
double StarDecompositionMeasure(const r3d::Polytope<Dim>& poly)
{
  const double eps = PCMS_INTERSECTION_ABS_TOL +
                     PCMS_INTERSECTION_REL_TOL * std::fabs(r3d::measure(poly));
  double sum = 0.0;
  pcms::detail::ForEachPolytopeStarSimplex<Dim>(
    poly, eps,
    [&](const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>&,
        Omega_h::Real measure) { sum += measure; });
  return sum;
}

// A clean simplex star-decomposes from its first vertex into exactly one
// piece, itself: the Dim facets incident to the apex have zero measure and are
// dropped. One piece means one quadrature rule per same-mesh element, which is
// the cost floor for the RHS build.
template <int Dim>
void RequireCleanSimplexDecomposition(const double* coords)
{
  r3d::Few<r3d::Vector<Dim>, Dim + 1> verts;
  for (int v = 0; v < Dim + 1; ++v) {
    for (int d = 0; d < Dim; ++d) {
      verts[v][d] = coords[Dim * v + d];
    }
  }
  r3d::Polytope<Dim> poly;
  r3d::init(poly, verts);
  const double measure = std::fabs(r3d::measure(poly));
  int pieces = 0;
  pcms::detail::ForEachPolytopeStarSimplex<Dim>(
    poly, 0.0,
    [&](const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>&,
        Omega_h::Real piece) {
      REQUIRE(piece == Catch::Approx(measure).epsilon(1e-12));
      ++pieces;
    });
  REQUIRE(pieces == 1);
}

} // namespace

TEST_CASE("star decomposition of a folded r3d polytope matches its measure",
          "[transfer][mesh_intersection]")
{
  SECTION("2D, clean triangle: one piece, itself")
  {
    constexpr double tri[6] = {0.0, 0.0, 4.0, 0.0, 0.0, 3.0};
    RequireCleanSimplexDecomposition<2>(tri);
  }
  SECTION("3D, clean tetrahedron: one piece, itself")
  {
    RequireCleanSimplexDecomposition<3>(kFoldingTetCoords.data());
  }
  SECTION("2D, folded polygon")
  {
    const auto poly = FoldedPolygon();
    const double area = std::fabs(r3d::measure(poly));
    REQUIRE(area == Catch::Approx(6.0));
    REQUIRE(StarDecompositionMeasure<2>(poly) ==
            Catch::Approx(area).epsilon(1e-12));
  }
  SECTION("3D, folded self-clip polytope")
  {
    const auto poly = FoldedSelfClipPolytope();
    const double volume = std::fabs(r3d::measure(poly));
    // Volume of the original tet (positions were rounded to 1e-6): r3d's own
    // signed integration gets it right despite the fold.
    REQUIRE(volume == Catch::Approx(2499.121742).epsilon(1e-6));
    REQUIRE(StarDecompositionMeasure<3>(poly) ==
            Catch::Approx(volume).epsilon(1e-12));
  }
}

TEST_CASE("self-clip star decomposition sums to the element measure",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  SECTION("3D, tets whose r3d self-clip is known to fold")
  {
    auto mesh = BuildDisjointSimplices<3>(lib, kFoldingTetCoords);
    REQUIRE(mesh.nelems() == 3);
    RequireSelfDecompositionExact<3>(mesh);
  }
  SECTION("2D, random triangles at realistic coordinate scale")
  {
    auto mesh =
      BuildDisjointSimplices<2>(lib, RandomSimplexCoords<2>(400, 20260903u));
    RequireSelfDecompositionExact<2>(mesh);
  }
  SECTION("3D, random tets at realistic coordinate scale")
  {
    auto mesh =
      BuildDisjointSimplices<3>(lib, RandomSimplexCoords<3>(400, 20260902u));
    RequireSelfDecompositionExact<3>(mesh);
  }
  SECTION("2D, scaled box mesh")
  {
    auto mesh = BuildScaledBox<2>(lib, 12);
    RequireSelfDecompositionExact<2>(mesh);
  }
  SECTION("3D, scaled box mesh")
  {
    auto mesh = BuildScaledBox<3>(lib, 5);
    RequireSelfDecompositionExact<3>(mesh);
  }
}

TEST_CASE("same-mesh P0->P1 lumped projection equals the element-nodal average",
          "[transfer][mesh_intersection]")
{
  Omega_h::Library lib;

  SECTION("3D, tets whose r3d self-clip is known to fold")
  {
    auto mesh = BuildDisjointSimplices<3>(lib, kFoldingTetCoords);
    RequireLumpedEqualsNodalAverage<3>(mesh);
  }
  SECTION("2D, scaled box mesh")
  {
    auto mesh = BuildScaledBox<2>(lib, 12);
    RequireLumpedEqualsNodalAverage<2>(mesh);
  }
  SECTION("3D, scaled box mesh")
  {
    auto mesh = BuildScaledBox<3>(lib, 5);
    RequireLumpedEqualsNodalAverage<3>(mesh);
  }
}

TEST_CASE("same-mesh P0->P1 consistent projection preserves zeroth and first "
          "moments",
          "[transfer][mesh_intersection]")
{
  // Galerkin orthogonality with v = 1 and each coordinate in P1: the
  // consistent L2 projection matches every moment the target space represents.
  Omega_h::Library lib;

  SECTION("2D, scaled box mesh")
  {
    auto mesh = BuildScaledBox<2>(lib, 12);
    RequireConsistentPreservesMoments<2>(mesh);
  }
  SECTION("3D, scaled box mesh")
  {
    auto mesh = BuildScaledBox<3>(lib, 5);
    RequireConsistentPreservesMoments<3>(mesh);
  }
}
