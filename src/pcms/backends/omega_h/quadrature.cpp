#include "pcms/backends/omega_h/omega_h_geometry.hpp"
#include "pcms/backends/omega_h/simplex_quadrature.hpp"
#include "pcms/configuration.hpp"
#include "pcms/geometry/detail/geometry.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/random.hpp"
#include <Kokkos_Random.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_map.hpp>
#include <Omega_h_shape.hpp>
#include <cstdint>
#include <string>
#include <variant>

namespace pcms::omega_h
{

namespace
{

using Generator = Kokkos::Random_XorShift64<DefaultExecutionSpace>;

// The global id of each element, from the geometry's global id tag.
Omega_h::GOs ElementGlobalIds(Omega_h::Mesh& mesh, const std::string& name)
{
  const int dim = mesh.dim();
  if (Omega_h::is<Omega_h::LO>(mesh.get_tagbase(dim, name))) {
    const auto ids = mesh.get_array<Omega_h::LO>(dim, name);
    Omega_h::Write<Omega_h::GO> gids(ids.size());
    Omega_h::parallel_for(
      ids.size(), OMEGA_H_LAMBDA(int i) { gids[i] = ids[i]; });
    return gids;
  }
  return mesh.get_array<Omega_h::GO>(dim, name);
}

// Normalized exponential draws are uniform on the simplex (Dirichlet(1, ...)).
template <int Dim>
KOKKOS_INLINE_FUNCTION Omega_h::Vector<Dim + 1> UniformBarycentric(
  Generator& gen)
{
  Omega_h::Vector<Dim + 1> b;
  Real sum = 0;
  for (int k = 0; k <= Dim; ++k) {
    b[k] = -Kokkos::log(1.0 - gen.drand());
    sum += b[k];
  }
  for (int k = 0; k <= Dim; ++k) {
    b[k] /= sum;
  }
  return b;
}

// The reference rule of `order` (order > 0) or `samples` uniform random
// samples (samples > 0) on each owned element.
template <int Dim>
Quadrature CreateSimplexQuadrature(const OmegaHGeometry& geometry, int order,
                                   int samples, std::uint64_t seed)
{
  static_assert(Dim == 2 || Dim == 3, "Omega_h meshes are 2D or 3D");
  // the measure of the reference simplex, 1 / Dim!
  constexpr Real ReferenceMeasure = [] {
    Real factorial = 1;
    for (int k = 2; k <= Dim; ++k) {
      factorial *= k;
    }
    return 1 / factorial;
  }();
  const bool random = samples > 0;
  int per_cell = samples;
  Kokkos::View<Real * [Dim + 1], DeviceMemorySpace> bary;
  Kokkos::View<Real*, DeviceMemorySpace> rule_weights;
  Omega_h::GOs gids;
  if (random) {
    gids = ElementGlobalIds(geometry.GetMesh(), geometry.GetGlobalIdName());
  } else {
#ifdef PCMS_ENABLE_OMEGA_H_MESHFIELDS
    const detail::IntegrationData<Dim> rule(order);
    bary = rule.bary_coords;
    rule_weights = rule.weights;
    per_cell = rule.size();
#else
    throw pcms_error("OmegaHGeometry::CreateQuadrature: the reference rules "
                     "need PCMS_ENABLE_OMEGA_H_MESHFIELDS");
#endif
  }
  Omega_h::Mesh& mesh = geometry.GetMesh();
  // a ghosted mesh also holds copies of its neighbours' elements
  const auto owned = Omega_h::collect_marked(mesh.owned(Dim));
  const LO num_owned = owned.size();
  const auto num_points = static_cast<std::size_t>(num_owned) * per_cell;
  Kokkos::View<LO*, DeviceMemorySpace> elements("quadrature_elements",
                                                num_points);
  Kokkos::View<Real**, DeviceMemorySpace> coords("quadrature_coords",
                                                 num_points, Dim);
  Kokkos::View<Real*, DeviceMemorySpace> weights("quadrature_weights",
                                                 num_points);
  const auto elem_verts = mesh.ask_elem_verts();
  const auto x = mesh.coords();
  const auto determinant = geometry.CellJacobian().determinant;
  Kokkos::parallel_for(
    "OmegaHGeometry::CreateQuadrature",
    Kokkos::RangePolicy<DefaultExecutionSpace>(0, num_owned),
    KOKKOS_LAMBDA(const LO k) {
      const LO e = owned[k];
      const auto verts = Omega_h::gather_verts<Dim + 1>(elem_verts, e);
      const Real measure = Kokkos::fabs(determinant(e));
      const auto fill = [&](int q, const Omega_h::Vector<Dim + 1>& b, Real w) {
        const LO i = k * per_cell + q;
        elements(i) = e;
        weights(i) = w * measure;
        for (int d = 0; d < Dim; ++d) {
          coords(i, d) = 0;
          for (int k = 0; k <= Dim; ++k) {
            coords(i, d) += b[k] * x[verts[k] * Dim + d];
          }
        }
      };
      if (random) {
        auto gen = pcms::detail::EntityGenerator(seed, gids[e]);
        for (int q = 0; q < per_cell; ++q) {
          fill(q, UniformBarycentric<Dim>(gen), ReferenceMeasure / per_cell);
        }
      } else {
        for (int q = 0; q < per_cell; ++q) {
          Omega_h::Vector<Dim + 1> b;
          for (int k = 0; k <= Dim; ++k) {
            b[k] = bary(q, k);
          }
          fill(q, b, rule_weights(q));
        }
      }
    });
  return pcms::detail::WeightedQuadrature(
    geometry,
    CoordinateView<DeviceMemorySpace>(geometry.GetCoordinateSystem(),
                                      MakeConstRank2View(coords)),
    weights, OwnerHints{geometry.shared_from_this(), Dim, elements});
}

} // namespace

Quadrature OmegaHGeometry::CreateQuadrature(
  const pcms::detail::QuadratureRule& rule) const
{
  int order = 0;
  int samples = 0;
  std::uint64_t seed = 0;
  if (const auto* gauss = std::get_if<pcms::detail::GaussRule>(&rule)) {
    if (gauss->order < 1) {
      throw pcms_error("OmegaHGeometry::CreateQuadrature: the order must be "
                       "positive");
    }
    order = gauss->order;
  } else if (const auto* mc =
               std::get_if<pcms::detail::MonteCarloRule>(&rule)) {
    samples = mc->samples_per_cell;
    seed = mc->seed;
  } else {
    return Geometry::CreateQuadrature(rule);
  }
  return detail::VisitMeshDimension(mesh_, [&]<int Dim>() {
    return CreateSimplexQuadrature<Dim>(*this, order, samples, seed);
  });
}

} // namespace pcms::omega_h
