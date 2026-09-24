#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>

#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/types.h"

#include <cmath>
#include <cstdio>
#include <type_traits>
#include <vector>

using pcms::LO;
using pcms::Real;

namespace
{
constexpr Real kGhostSentinel = -1.0e30;
constexpr Real kFillValue = -2.0e30;

// Sentinel used to pre-fill the evaluation buffer. It is deliberately different
// from kFillValue: the FILL policy writes kFillValue into every point that
// localization rejects, so if the buffer were pre-filled with the same value, a
// point the evaluator never wrote and a point it rejected would be
// indistinguishable. With a distinct sentinel the three outcomes separate:
//   kUnresolved -> the evaluator never wrote this point (kernel/copy did not
//   run) kFillValue  -> localization rejected this point and the FILL policy
//   wrote it anything else -> a real interpolated value
constexpr Real kUnresolved = -3.0e30;

KOKKOS_INLINE_FUNCTION Real linear_f(Real x, Real y)
{
  return x + 2.0 * y;
}

template <typename CoordView>
Kokkos::View<Real**, pcms::HostMemorySpace> copy_coordinates_to_host(
  const CoordView& coordinates, LO n, int dim)
{
  using layout = std::conditional_t<
    std::is_same_v<typename CoordView::layout_type, Kokkos::layout_left>,
    Kokkos::LayoutLeft, Kokkos::LayoutRight>;
  Kokkos::View<const Real**, layout, pcms::DeviceMemorySpace,
               Kokkos::MemoryUnmanaged>
    device_coordinates(coordinates.data_handle(), n, dim);
  Kokkos::View<Real**, pcms::HostMemorySpace> host_coordinates("coordinates", n,
                                                               dim);
  pcms::DeepCopyMismatchLayouts(host_coordinates, device_coordinates);
  return host_coordinates;
}

struct DeviceCoordinates
{
  Kokkos::View<Real**, pcms::DeviceMemorySpace> storage;
  pcms::CoordinateView<pcms::DeviceMemorySpace> coordinate_view;
};

DeviceCoordinates make_device_coordinates(const std::vector<Real>& points)
{
  const LO n = static_cast<LO>(points.size() / 2);
  Kokkos::View<Real**, pcms::DeviceMemorySpace> storage("coordinates", n, 2);
  Kokkos::View<Real**, pcms::HostMemorySpace> host("coordinates", n, 2);
  for (LO i = 0; i < n; ++i) {
    host(i, 0) = points[static_cast<size_t>(i) * 2];
    host(i, 1) = points[static_cast<size_t>(i) * 2 + 1];
  }
  pcms::DeepCopyMismatchLayouts(storage, host);
  return {storage,
          {pcms::CoordinateSystem::Cartesian, pcms::MakeRank2View(storage)}};
}

void require(bool condition, const char* message)
{
  if (!condition) {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::fprintf(stderr, "[rank %d] %s\n", rank, message);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

void check_distributed_evaluation(Omega_h::Mesh& mesh,
                                  const pcms::LagrangeFunctionSpace& space)
{
  const auto& layout = *space.GetLayout();
  const LO n_local = layout.GetNumLocalDofHolder();
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  require(n_local > layout.GetNumOwnedDofHolder(),
          "ghosted mesh did not create ghost DOF holders");

  auto coords = copy_coordinates_to_host(
    layout.GetDOFHolderCoordinates().GetValues(), n_local, mesh.dim());
  auto owned = layout.GetOwnedHost();
  std::vector<Real> values(static_cast<size_t>(n_local), kGhostSentinel);
  for (LO i = 0; i < n_local; ++i) {
    if (owned[static_cast<size_t>(i)]) {
      values[static_cast<size_t>(i)] = linear_f(coords(i, 0), coords(i, 1));
    }
  }

  auto field = space.CreateFunction<Real>();
  field.SetDOFHolderDataHost(pcms::Rank2View<const Real, pcms::HostMemorySpace>(
    values.data(), n_local, 1));
  field.SynchronizeGhosts();

  auto data = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());
  for (LO i = 0; i < n_local; ++i) {
    require(std::fabs(data[static_cast<size_t>(i)] -
                      linear_f(coords(i, 0), coords(i, 1))) < 1.0e-12,
            "ghost synchronization did not produce the affine field");
  }

  const LO n_elements = mesh.nents(mesh.dim());
  const int vertices_per_element = mesh.dim() + 1;
  auto element_vertices = Omega_h::HostRead<LO>(mesh.ask_elem_verts());
  auto mesh_coordinates = Omega_h::HostRead<Real>(mesh.coords());
  std::vector<Real> points(static_cast<size_t>(n_elements) * mesh.dim());
  for (LO element = 0; element < n_elements; ++element) {
    for (int dim = 0; dim < mesh.dim(); ++dim) {
      Real centroid = 0.0;
      for (int vertex = 0; vertex < vertices_per_element; ++vertex) {
        const LO vertex_id =
          element_vertices[element * vertices_per_element + vertex];
        centroid += mesh_coordinates[vertex_id * mesh.dim() + dim];
      }
      points[static_cast<size_t>(element) * mesh.dim() + dim] =
        centroid / vertices_per_element;
    }
  }
  auto device_coords = make_device_coordinates(points);
  auto evaluator =
    space.CreatePointEvaluator<Real>(pcms::EvaluationRequest::FromCoordinates(
      device_coords.coordinate_view,
      pcms::OutOfBoundsPolicy{pcms::OutOfBoundsMode::FILL, kFillValue}));

  const LO n_points = n_elements;
  Kokkos::View<Real**, pcms::DeviceMemorySpace> result("result", n_points, 1);
  Kokkos::deep_copy(result, kUnresolved);
  evaluator->Evaluate(field, pcms::MakeRank2View(result));
  auto result_host =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), result);

  // A point still holding kUnresolved was never written by the evaluator (the
  // evaluation kernels and/or the device-to-host copy produced nothing); a
  // point holding kFillValue was written by the FILL policy, i.e. localization
  // rejected it. They are counted separately because the causes differ.
  LO n_unlocalized = 0;
  LO n_untouched = 0;
  LO first_unlocalized = -1;
  LO first_untouched = -1;
  for (LO i = 0; i < n_points; ++i) {
    const Real value = result_host(i, 0);
    if (value == kUnresolved) {
      ++n_untouched;
      if (first_untouched < 0) {
        first_untouched = i;
      }
      continue;
    }
    if (value == kFillValue) {
      ++n_unlocalized;
      if (first_unlocalized < 0) {
        first_unlocalized = i;
      }
      continue;
    }
    const size_t offset = static_cast<size_t>(i) * 2;
    const Real expected = linear_f(points[offset], points[offset + 1]);
    require(std::fabs(value - expected) < 1.0e-12,
            "distributed evaluation returned the wrong value");
  }

  if (n_untouched != 0) {
    std::fprintf(stderr,
                 "[rank %d] %lld/%lld element centroid results were never "
                 "written by the evaluator (still at the pre-fill sentinel "
                 "%.17g, first index %lld): evaluation/copy failure, not a "
                 "localization failure\n",
                 rank, static_cast<long long>(n_untouched),
                 static_cast<long long>(n_points),
                 static_cast<double>(kUnresolved),
                 static_cast<long long>(first_untouched));
  }

  if (n_unlocalized != 0) {
    const size_t offset = static_cast<size_t>(first_unlocalized) * 2;
    std::fprintf(stderr,
                 "[rank %d] %lld/%lld element centroids were not localized "
                 "(first index %lld at point=(%.17g, %.17g))\n",
                 rank, static_cast<long long>(n_unlocalized),
                 static_cast<long long>(n_points),
                 static_cast<long long>(first_unlocalized),
                 static_cast<double>(points[offset]),
                 static_cast<double>(points[offset + 1]));
  }

  // Distinct messages so the failure class is obvious from the abort line
  // alone.
  require(n_unlocalized == 0, "local element centroid was not localized");
  require(n_untouched == 0,
          "element centroid evaluation result was never written");

  if (rank == 0)
    std::printf("distributed field test passed\n");
}
} // namespace

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  {
    Omega_h::Library library(&argc, &argv);
    auto world = library.world();
    require(world->size() >= 2, "test requires at least 2 MPI ranks");

    const int n = 4 * world->size();
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, n, n, 0, false);
    mesh.set_parting(OMEGA_H_GHOSTED, 1, false);
    auto space = pcms::LagrangeFunctionSpace::FromMesh(
      mesh, 1, 1, pcms::CoordinateSystem::Cartesian, "global",
      pcms::LagrangeFunctionSpace::Backend::OmegaH);

    check_distributed_evaluation(mesh, *space);
  }
  MPI_Finalize();
  return 0;
}
