#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

#include "distributed_test_utils.hpp"
#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/field/function_spaces.hpp"
#include "pcms/geometry/quadrature.hpp"
#include "pcms/utility/arrays.hpp"
#ifdef PCMS_ENABLE_CONSERVATIVE_PROJECTION
#include "pcms/transfer/omega_h_intersection_quadrature.hpp"
#endif

#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

using pcms::GO;
using pcms::LO;
using pcms::Real;

namespace
{

using Host2D = Kokkos::View<Real**, pcms::HostMemorySpace>;
using Device2D = Kokkos::View<Real**, pcms::DeviceMemorySpace>;

int failures = 0;

int Rank()
{
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

void Check(bool condition, const std::string& what)
{
  if (!condition) {
    ++failures;
    std::fprintf(stderr, "[rank %d] check failed: %s\n", Rank(), what.c_str());
  }
}

Real SumOverRanks(Real local)
{
  Real global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

long long SumOverRanks(long long local)
{
  long long global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

// Linear in the coordinates, a different function per component.
Real Linear(const Real* x, int dim, int c)
{
  Real value = 1.0 + c;
  for (int d = 0; d < dim; ++d) {
    value += (d + 1.0 + 0.5 * c) * x[d];
  }
  return value;
}

// A DOF value that depends only on the holder's global id, so that every rank
// knows the value of every holder; computed exactly, so that every call site
// gives the same value however the compiler contracts it.
Real DofOfGid(GO gid, int c)
{
  return static_cast<Real>((gid * 37 + GO{11} * c) % 101) / 128.0 - 0.4;
}

std::vector<Real> Centroids(Omega_h::Mesh& mesh, int entity_dim)
{
  const int dim = mesh.dim();
  const auto coords = Omega_h::HostRead<Real>(mesh.coords());
  const auto verts = Omega_h::HostRead<LO>(mesh.ask_verts_of(entity_dim));
  const int per = entity_dim + 1;
  const LO n = mesh.nents(entity_dim);
  std::vector<Real> points(static_cast<size_t>(n) * dim, 0.0);
  for (LO e = 0; e < n; ++e) {
    for (int k = 0; k < per; ++k) {
      for (int d = 0; d < dim; ++d) {
        points[e * dim + d] += coords[verts[e * per + k] * dim + d] / per;
      }
    }
  }
  return points;
}

// The local points: element centroids, and for order 1 also the vertices and
// edge midpoints, which include the points on the partition boundary.
std::vector<Real> LocalPoints(Omega_h::Mesh& mesh, int order)
{
  auto points = Centroids(mesh, mesh.dim());
  if (order == 1) {
    for (int entity_dim : {0, 1}) {
      const auto coords = Omega_h::HostRead<Real>(mesh.coords());
      const auto more =
        entity_dim == 0
          ? std::vector<Real>(coords.data(), coords.data() + coords.size())
          : Centroids(mesh, entity_dim);
      points.insert(points.end(), more.begin(), more.end());
    }
  }
  return points;
}

struct Points
{
  Device2D storage;
  std::shared_ptr<const pcms::BoundPoints> bound;
};

Points Bind(const pcms::Geometry& geometry, const std::vector<Real>& points,
            int dim)
{
  const LO n = static_cast<LO>(points.size() / dim);
  Device2D storage("points", n, dim);
  auto host = Kokkos::create_mirror_view(storage);
  for (LO i = 0; i < n; ++i) {
    for (int d = 0; d < dim; ++d) {
      host(i, d) = points[i * dim + d];
    }
  }
  Kokkos::deep_copy(storage, host);
  auto bound = geometry.Bind(pcms::CoordinateView<pcms::DeviceMemorySpace>(
    geometry.GetCoordinateSystem(), pcms::MakeConstRank2View(storage)));
  return {storage, bound};
}

Device2D ToDevice(const Host2D& host)
{
  Device2D device("device", host.extent(0), host.extent(1));
  auto mirror = Kokkos::create_mirror_view(device);
  for (size_t i = 0; i < host.extent(0); ++i) {
    for (size_t j = 0; j < host.extent(1); ++j) {
      mirror(i, j) = host(i, j);
    }
  }
  Kokkos::deep_copy(device, mirror);
  return device;
}

Host2D ToHost(const Device2D& device)
{
  auto mirror =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, device);
  Host2D host("host", device.extent(0), device.extent(1));
  for (size_t i = 0; i < host.extent(0); ++i) {
    for (size_t j = 0; j < host.extent(1); ++j) {
      host(i, j) = mirror(i, j);
    }
  }
  return host;
}

void SetOwned(pcms::Field<Real>& field, const Host2D& values)
{
  field.SetDOFHolderDataHost(pcms::ValueView<const Real, pcms::HostMemorySpace>(
    field.GetData().GetValueBasis(), pcms::MakeConstRank2View(values)));
  field.Synchronize();
}

Host2D Evaluate(const pcms::PointEvaluator<Real>& evaluator,
                const pcms::Field<Real>& field, int nc)
{
  Device2D values("values", evaluator.NumPoints(), nc);
  evaluator.Evaluate(field, distributed_test::TagLike(field, values));
  return ToHost(values);
}

Host2D OwnedHolderCoordinates(const pcms::FieldLayout& layout, int dim)
{
  const auto coords = layout.GetDOFHolderCoordinates().GetValues();
  Device2D device("holder_coordinates", layout.GetNumDofHolder(), dim);
  Kokkos::parallel_for(
    "OwnedHolderCoordinates",
    Kokkos::MDRangePolicy<pcms::DefaultExecutionSpace, Kokkos::Rank<2>>(
      {0, 0},
      {static_cast<int64_t>(device.extent(0)), static_cast<int64_t>(dim)}),
    KOKKOS_LAMBDA(const int64_t i, const int64_t d) {
      device(i, d) = coords(i, d);
    });
  return ToHost(device);
}

std::string Case(bool ghosted, int dim, int order, int nc)
{
  return std::string(ghosted ? "ghosted" : "element-based") + " dim " +
         std::to_string(dim) + " order " + std::to_string(order) +
         " components " + std::to_string(nc);
}

// Evaluate reproduces a field exact in the space at every local point.
void CheckEvaluate(const pcms::FunctionSpace& space,
                   const std::vector<Real>& points, const Points& bound,
                   int dim, const std::string& name)
{
  const auto& layout = *space.GetLayout();
  const int nc = layout.GetNumComponents();
  const auto holders = OwnedHolderCoordinates(layout, dim);
  Host2D owned("owned", layout.GetNumDofHolder(), nc);
  for (LO h = 0; h < layout.GetNumDofHolder(); ++h) {
    for (int c = 0; c < nc; ++c) {
      owned(h, c) = Linear(&holders(h, 0), dim, c);
    }
  }
  auto field = space.CreateFunction<Real>();
  SetOwned(field, owned);
  const auto evaluator = space.CreatePointEvaluator<Real>(bound.bound);
  const auto values = Evaluate(*evaluator, field, nc);
  int wrong = 0;
  for (LO p = 0; p < evaluator->NumPoints(); ++p) {
    for (int c = 0; c < nc; ++c) {
      const Real expected =
        Linear(&points[static_cast<std::size_t>(p) * dim], dim, c);
      if (std::abs(values(p, c) - expected) > 1e-12) {
        if (wrong == 0) {
          std::fprintf(stderr, "[rank %d] %s: point %d comp %d: %g != %g\n",
                       Rank(), name.c_str(), p, c, values(p, c), expected);
        }
        ++wrong;
      }
    }
  }
  Check(wrong == 0, name + ": Evaluate of a field exact in the space");
}

// The adjoint identity summed over the ranks, the overwrite of a second
// EvaluateTranspose, and the matrix against Evaluate.
void CheckTransposeAndMatrix(const pcms::FunctionSpace& space,
                             const Points& bound, const std::string& name)
{
  const auto& layout = *space.GetLayout();
  const int nc = layout.GetNumComponents();
  const LO n_owned = layout.GetNumDofHolder();
  const auto gids = layout.GetGidsHost();
  Host2D x("x", n_owned, nc);
  for (LO h = 0; h < n_owned; ++h) {
    for (int c = 0; c < nc; ++c) {
      x(h, c) = DofOfGid(gids(h), c);
    }
  }
  auto field = space.CreateFunction<Real>();
  SetOwned(field, x);
  const auto evaluator = space.CreatePointEvaluator<Real>(bound.bound);
  const LO n_points = evaluator->NumPoints();
  const auto e_x = Evaluate(*evaluator, field, nc);

  Host2D y("y", n_points, nc);
  for (LO p = 0; p < n_points; ++p) {
    for (int c = 0; c < nc; ++c) {
      y(p, c) = std::cos(0.37 * (p + 1) * (Rank() + 1) + 1.1 * c);
    }
  }
  const auto y_device = ToDevice(y);
  auto functional = space.CreateFunction<Real>();
  const pcms::ValueView<const Real, pcms::DeviceMemorySpace> y_view(
    functional.GetData().GetValueBasis(), pcms::MakeConstRank2View(y_device));
  evaluator->EvaluateTranspose(y_view, functional);
  Host2D once("once", n_owned, nc);
  {
    const auto values = functional.GetDOFHolderDataHost().GetValues();
    for (LO h = 0; h < n_owned; ++h) {
      for (int c = 0; c < nc; ++c) {
        once(h, c) = values(h, c);
      }
    }
  }
  Real lhs = 0, rhs = 0, scale = 0;
  for (LO p = 0; p < n_points; ++p) {
    for (int c = 0; c < nc; ++c) {
      lhs += e_x(p, c) * y(p, c);
      scale += std::abs(e_x(p, c) * y(p, c));
    }
  }
  for (LO h = 0; h < n_owned; ++h) {
    for (int c = 0; c < nc; ++c) {
      rhs += x(h, c) * once(h, c);
    }
  }
  lhs = SumOverRanks(lhs);
  rhs = SumOverRanks(rhs);
  scale = SumOverRanks(scale);
  if (std::abs(lhs - rhs) > 1e-12 * std::max(1.0, scale)) {
    std::fprintf(stderr, "[rank %d] %s: <Ex,y> = %.17g, <x,E^T y> = %.17g\n",
                 Rank(), name.c_str(), lhs, rhs);
  }
  Check(std::abs(lhs - rhs) <= 1e-12 * std::max(1.0, scale),
        name + ": adjoint identity over the ranks");

  evaluator->EvaluateTranspose(y_view, functional);
  {
    const auto again = functional.GetDOFHolderDataHost().GetValues();
    int differ = 0;
    for (LO h = 0; h < n_owned; ++h) {
      for (int c = 0; c < nc; ++c) {
        differ += std::abs(again(h, c) - once(h, c)) > 1e-14;
      }
    }
    Check(differ == 0, name + ": a second EvaluateTranspose overwrites");
  }

  pcms::detail::EvaluateTranspose(*evaluator, y_view, functional,
                                  pcms::Real(1));
  {
    const auto twice = functional.GetDOFHolderDataHost().GetValues();
    int differ = 0;
    for (LO h = 0; h < n_owned; ++h) {
      for (int c = 0; c < nc; ++c) {
        differ += std::abs(twice(h, c) - 2 * once(h, c)) > 1e-14;
      }
    }
    Check(differ == 0, name + ": beta = 1 adds E^T y to the functional");
  }

  // E x from the CSR rows, with x looked up by global id
  const auto matrix = evaluator->GetMatrix();
  const auto offsets = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace{}, matrix.row_offsets);
  const auto columns = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace{}, matrix.columns);
  const auto entries =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, matrix.values);
  Check(offsets.extent(0) == static_cast<size_t>(n_points) + 1,
        name + ": one CSR row per point");
  int wrong = 0;
  for (LO p = 0; p + 1 < static_cast<LO>(offsets.extent(0)); ++p) {
    for (int c = 0; c < nc; ++c) {
      Real sum = 0;
      for (LO k = offsets(p); k < offsets(p + 1); ++k) {
        sum += entries(k) * DofOfGid(columns(k), c);
      }
      if (std::abs(sum - e_x(p, c)) > 1e-13) {
        if (wrong < 3) {
          std::fprintf(stderr,
                       "[rank %d] %s: point %d of %d comp %d: rows %.17g, "
                       "Evaluate %.17g\n",
                       Rank(), name.c_str(), p, n_points, c, sum, e_x(p, c));
          for (LO k = offsets(p); k < offsets(p + 1); ++k) {
            std::fprintf(stderr, "    k %d gid %lld weight %.17g value %.17g\n",
                         k, static_cast<long long>(columns(k)), entries(k),
                         DofOfGid(columns(k), c));
          }
        }
        ++wrong;
      }
    }
  }
  if (wrong > 0) {
    std::fprintf(stderr, "[rank %d] %s: %d of %d values differ\n", Rank(),
                 name.c_str(), wrong, n_points * nc);
  }
  Check(wrong == 0, name + ": the CSR rows times x equal Evaluate");
}

// A quadrature sums to the domain measure over the ranks and, given its
// points per element (> 0), has each element on exactly one rank.
void CheckQuadrature(const pcms::Quadrature& quadrature, Omega_h::Mesh& mesh,
                     int points_per_element, const std::string& name)
{
  const auto weights = Kokkos::create_mirror_view_and_copy(
    pcms::HostMemorySpace{}, quadrature.weights);
  Real sum = 0;
  for (size_t i = 0; i < weights.extent(0); ++i) {
    sum += weights(i);
  }
  sum = SumOverRanks(sum);
  const long long n_points =
    SumOverRanks(static_cast<long long>(weights.extent(0)));
  const long long expected_points =
    static_cast<long long>(mesh.nglobal_ents(mesh.dim())) * points_per_element;
  const bool count_ok = points_per_element <= 0 || n_points == expected_points;
  if (std::abs(sum - 1.0) > 1e-12 || !count_ok) {
    std::fprintf(stderr,
                 "[rank %d] %s: measure %.17g (expected 1), %lld points "
                 "(expected %lld)\n",
                 Rank(), name.c_str(), sum, n_points, expected_points);
  }
  Check(std::abs(sum - 1.0) <= 1e-12,
        name + ": the weights sum to the measure");
  Check(count_ok, name + ": each element is integrated on one rank");
}

// The points per element of the Gauss rule of `order` on a serial mesh.
int PointsPerElement(Omega_h::Library& library, int dim, int order)
{
  auto mesh =
    Omega_h::build_box(library.self(), OMEGA_H_SIMPLEX, 1.0, 1.0,
                       dim == 3 ? 1.0 : 0.0, 1, 1, dim == 3 ? 1 : 0, false);
  auto geometry =
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
  return static_cast<int>(
    pcms::CreateGaussQuadrature(geometry, order).weights.extent(0) /
    mesh.nelems());
}

void RunCases(Omega_h::Library& library, int dim, bool ghosted)
{
  auto world = library.world();
  const int n = (dim == 2 ? 4 : 2) * world->size();
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, dim == 3 ? 1.0 : 0.0,
                       n, n, dim == 3 ? n : 0, false);
  if (ghosted) {
    mesh.set_parting(OMEGA_H_GHOSTED, 1, false);
  } else {
    mesh.set_parting(OMEGA_H_ELEM_BASED, 0, false);
  }
  auto geometry =
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
  for (int order : {0, 1}) {
    const auto points = LocalPoints(mesh, order);
    const auto bound = Bind(*geometry, points, dim);
    for (int nc : {1, 2}) {
      const auto name = Case(ghosted, dim, order, nc);
      const auto space = pcms::CreateLagrangeSpace(
        geometry, order, nc,
        pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH});
      CheckEvaluate(*space, points, bound, dim, name);
      CheckTransposeAndMatrix(*space, bound, name);
      if (nc == 1) {
        CheckQuadrature(space->GetQuadrature(), mesh,
                        PointsPerElement(library, dim, std::max(1, 2 * order)),
                        name + ": GetQuadrature");
      }
    }
  }
  const auto gauss = pcms::CreateGaussQuadrature(geometry, 2);
  CheckQuadrature(gauss, mesh, PointsPerElement(library, dim, 2),
                  Case(ghosted, dim, 1, 1) + ": CreateGaussQuadrature");
  CheckQuadrature(pcms::CreateMonteCarloQuadrature(geometry, 3), mesh, 3,
                  Case(ghosted, dim, 1, 1) + ": CreateMonteCarloQuadrature");
#ifdef PCMS_ENABLE_CONSERVATIVE_PROJECTION
  // the intersection of the mesh with itself; a 3D element may be cut into
  // several sub-simplices
  const auto omega_h_geometry =
    std::dynamic_pointer_cast<const pcms::omega_h::OmegaHGeometry>(geometry);
  CheckQuadrature(pcms::omega_h::CreateIntersectionQuadrature(
                    omega_h_geometry, omega_h_geometry, 2),
                  mesh, 0,
                  Case(ghosted, dim, 1, 1) + ": CreateIntersectionQuadrature");
#endif
}

} // namespace

int main(int argc, char** argv)
{
  try {
    MPI_Init(&argc, &argv);
    {
      Omega_h::Library library(&argc, &argv);
      auto world = library.world();
      distributed_test::require(world->size() >= 2,
                                "test requires at least 2 MPI ranks");
      for (int dim : {2, 3}) {
        for (bool ghosted : {false, true}) {
          RunCases(library, dim, ghosted);
        }
      }
    }
    const long long total = SumOverRanks(static_cast<long long>(failures));
    if (Rank() == 0) {
      std::printf("distributed evaluation test: %lld failed checks\n", total);
    }
    MPI_Finalize();
    return total == 0 ? 0 : 1;
  } catch (const std::exception& e) {
    distributed_test::AbortAllRanks(e.what());
  } catch (...) {
    distributed_test::AbortAllRanks("unknown exception");
  }
}
