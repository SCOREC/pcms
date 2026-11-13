#ifndef PCMS_TEST_DISTRIBUTED_TEST_UTILS_H
#define PCMS_TEST_DISTRIBUTED_TEST_UTILS_H

// Shared helpers for the MPI tests of distributed fields, search and
// evaluation.
//
// These tests run under mpirun without Catch2, so a failed check aborts every
// rank: a rank that just returned would leave the others blocked in the next
// collective.

#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_mesh.hpp>

#include "pcms/configuration.h"
#include "pcms/field/data/omega_h.h"
#include "pcms/field/data/simple.h"
#ifdef PCMS_ENABLE_MESHFIELDS
#include "pcms/field/data/mesh_fields.h"
#endif
#include "pcms/field/field.h"
#include "pcms/field/value_view.hpp"
#include "pcms/localization/point_search.h"
#include "pcms/utility/mesh_geometry.h"
#include "pcms/utility/mpi_utilities.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace distributed_test
{

using pcms::LO;
using pcms::Real;

/// Abort every rank with `message` unless `condition` holds.
inline void require(bool condition, const char* message)
{
  if (!condition) {
    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::fprintf(stderr, "[rank %d] check failed: %s\n", rank, message);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}

/// Report `what` and abort every rank. For exceptions escaping `main`: a rank
/// that just returned would leave the others blocked in their next collective.
[[noreturn]] inline void AbortAllRanks(const char* what) noexcept
{
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::fprintf(stderr, "[rank %d] uncaught exception: %s\n", rank, what);
  MPI_Abort(MPI_COMM_WORLD, 1);
  std::abort(); // MPI_Abort does not return, but is not declared noreturn
}

// Abort unless every rank holds the same integer value.
inline void require_same_on_all_ranks(Omega_h::CommPtr comm, LO value,
                                      const char* message)
{
  const auto gathered = pcms::AllGatherAllRanks(comm->get_impl(), value);
  for (size_t i = 0; i < gathered.size(); ++i) {
    if (gathered[i] != gathered[0]) {
      std::fprintf(
        stderr, "[rank %d] %s: rank %zu holds %lld, rank 0 holds %lld\n",
        comm->rank(), message, i, static_cast<long long>(gathered[i]),
        static_cast<long long>(gathered[0]));
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}

// Abort unless every rank holds the same real value, to `tolerance`.
inline void require_same_real_on_all_ranks(Omega_h::CommPtr comm, Real value,
                                           Real tolerance, const char* message)
{
  const auto gathered = pcms::AllGatherAllRanks(comm->get_impl(), value);
  for (size_t i = 0; i < gathered.size(); ++i) {
    if (std::fabs(gathered[i] - gathered[0]) > tolerance) {
      std::fprintf(stderr,
                   "[rank %d] %s: rank %zu holds %.17g, rank 0 holds %.17g\n",
                   comm->rank(), message, i, gathered[i], gathered[0]);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}

// The unit square (2D) or cube (3D) as a simplex mesh, partitioned over
// `world` with one ghost layer. The resolution grows with the rank count so
// every rank holds some elements.
//
// `world` must have a power-of-two size: `Omega_h::build_box` partitions by
// recursive bisection. A debug Omega_h asserts on any other size; a release
// one silently produces a broken partition.
template <int Dim>
Omega_h::Mesh BuildGhostedBox(Omega_h::CommPtr world)
{
  const int size = world->size();
  require(size > 0 && (size & (size - 1)) == 0,
          "BuildGhostedBox needs a power-of-two rank count (Omega_h::build_box "
          "partitions by recursive bisection)");
  const LO n = (Dim == 2 ? 4 : 2) * world->size();
  const Real z = Dim == 2 ? 0.0 : 1.0;
  const LO nz = Dim == 2 ? 0 : n;
  auto mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, z, n, n, nz, false);
  mesh.set_parting(OMEGA_H_GHOSTED, 1, false);
  return mesh;
}

// The grid point search for `Dim`, with `n` divisions per axis.
template <int Dim>
struct LocalSearchOf;
template <>
struct LocalSearchOf<2>
{
  using type = pcms::GridPointSearch2D;
};
template <>
struct LocalSearchOf<3>
{
  using type = pcms::GridPointSearch3D;
};

template <int Dim>
typename LocalSearchOf<Dim>::type MakeLocalSearch(Omega_h::Mesh& mesh, LO n)
{
  if constexpr (Dim == 3) {
    return typename LocalSearchOf<Dim>::type(mesh, n, n, n);
  } else {
    return typename LocalSearchOf<Dim>::type(mesh, n, n);
  }
}

// Points gathered from every rank onto every rank, in a deterministic order
// (contributing rank, then local index).
struct GatheredPoints
{
  std::vector<Real> coords;        // [total * dim], interleaved
  std::vector<int> contributor;    // [total], the rank each point came from
  std::vector<pcms::GO> global_id; // [total], when gathered per element
  int total = 0;
};

// Gather `local` (interleaved, `dim` values per point) from every rank.
inline GatheredPoints GatherGlobalPoints(Omega_h::Mesh& mesh,
                                         const std::vector<Real>& local,
                                         int dim)
{
  GatheredPoints out;
  std::vector<int> value_counts;
  out.coords =
    pcms::AllGatherVAllRanks(mesh.comm()->get_impl(), local, &value_counts);
  out.total = static_cast<int>(out.coords.size()) / dim;
  for (size_t r = 0; r < value_counts.size(); ++r) {
    out.contributor.insert(out.contributor.end(),
                           static_cast<size_t>(value_counts[r] / dim),
                           static_cast<int>(r));
  }
  return out;
}

// Every owned element centroid in the decomposition, with its global id.
inline GatheredPoints GatherOwnedCentroids(Omega_h::Mesh& mesh)
{
  const auto dim = mesh.dim();
  const auto owned = Omega_h::HostRead<Omega_h::I8>(mesh.owned(dim));
  const auto centroids =
    Omega_h::HostRead<Real>(pcms::get_entity_centroids(mesh, dim));
  std::vector<Real> local;
  std::vector<pcms::GO> local_globals;
  const bool has_globals = mesh.has_tag(dim, "global");
  const auto globals = has_globals
                         ? Omega_h::HostRead<pcms::GO>(mesh.globals(dim))
                         : Omega_h::HostRead<pcms::GO>();
  for (LO e = 0; e < mesh.nelems(); ++e) {
    if (owned[e] == 0) {
      continue;
    }
    for (int d = 0; d < dim; ++d) {
      local.push_back(centroids[e * dim + d]);
    }
    if (has_globals) {
      local_globals.push_back(globals[e]);
    }
  }
  auto out = GatherGlobalPoints(mesh, local, dim);
  if (has_globals) {
    out.global_id =
      pcms::AllGatherVAllRanks(mesh.comm()->get_impl(), local_globals);
  }
  return out;
}

// Every vertex of every owned element, so the points land on the faces and
// edges elements share with their neighbours.
inline GatheredPoints GatherOwnedElementVertices(Omega_h::Mesh& mesh)
{
  const auto dim = mesh.dim();
  const auto owned = Omega_h::HostRead<Omega_h::I8>(mesh.owned(dim));
  const auto verts = Omega_h::HostRead<LO>(mesh.ask_verts_of(dim));
  const auto coords = Omega_h::HostRead<Real>(mesh.coords());
  std::vector<Real> local;
  for (LO e = 0; e < mesh.nelems(); ++e) {
    if (owned[e] == 0) {
      continue;
    }
    for (int v = 0; v <= dim; ++v) {
      const LO vert = verts[e * (dim + 1) + v];
      for (int d = 0; d < dim; ++d) {
        local.push_back(coords[vert * dim + d]);
      }
    }
  }
  return GatherGlobalPoints(mesh, local, dim);
}

// Interleaved coordinates as a device point view, `[n][Dim]`.
template <int Dim>
Kokkos::View<Real* [Dim]> ToPointView(const std::vector<Real>& interleaved,
                                      const char* label = "points") {
  const auto n = static_cast<LO>(interleaved.size() / Dim);
  Kokkos::View<Real* [Dim]> points(std::string(label), n);
  auto points_h = Kokkos::create_mirror_view(points);
  for (LO i = 0; i < n; ++i) {
    for (int d = 0; d < Dim; ++d) {
      points_h(i, d) = interleaved[static_cast<size_t>(i) * Dim + d];
    }
  }
  Kokkos::deep_copy(points, points_h);
  return points;
}

/// `values` tagged with `field`'s stored value basis: what a field's checked
/// setters and an evaluator's output expect. A Catch2-free counterpart of
/// `pcms::test::TagLike`.
template <typename T, typename MemorySpace, typename Layout>
pcms::ValueView<T, MemorySpace, Layout> TagLike(
  const pcms::Field<std::remove_const_t<T>>& field,
  pcms::Rank2View<T, MemorySpace, Layout> values)
{
  return {field.GetData().GetValueBasis(), values};
}

template <typename T>
pcms::ValueView<T, pcms::DeviceMemorySpace> TagLike(
  const pcms::Field<T>& field,
  Kokkos::View<T**, pcms::DeviceMemorySpace> values)
{
  return {field.GetData().GetValueBasis(), pcms::MakeRank2View(values)};
}

/// The values of every DOF holder of `field`, owned and ghost, in mesh order,
/// copied to the host; ghost values are correct only after Synchronize. Field
/// itself only
/// exposes the owned values, so this reads the backend's storage: the Omega_h
/// tag of an OmegaHFieldData, or the MeshFields backend.
template <typename T>
Kokkos::View<T**, pcms::HostMemorySpace> ValuesWithGhostsHost(
  const pcms::Field<T>& field)
{
  if (const auto* omega_h_data =
        dynamic_cast<const pcms::OmegaHFieldData<T>*>(&field.GetData())) {
    const auto& layout = field.GetLayout();
    const auto nc = static_cast<size_t>(layout.GetNumComponents());
    // The values of every DOF holder are in the field's Omega_h tag.
    using TagValue = typename pcms::OmegaHFieldData<T>::tag_value_type;
    const auto& discretization =
      dynamic_cast<const pcms::OmegaHDiscretization&>(
        *layout.GetDiscretization());
    const auto tag =
      Omega_h::HostRead<TagValue>(discretization.GetMesh().get_array<TagValue>(
        layout.GetDOFHolderEntityDim(), omega_h_data->GetTagName()));
    const size_t n = static_cast<size_t>(tag.size()) / nc;
    Kokkos::View<T**, pcms::HostMemorySpace> host("values_with_ghosts_host", n,
                                                  nc);
    for (size_t i = 0; i < n; ++i) {
      for (size_t comp = 0; comp < nc; ++comp) {
        host(i, comp) = static_cast<T>(tag[i * nc + comp]);
      }
    }
    return host;
  }
#ifdef PCMS_ENABLE_MESHFIELDS
  if (const auto* mesh_fields =
        dynamic_cast<const pcms::MeshFieldsFieldData<T>*>(&field.GetData())) {
    // The MeshFields backend stores every local entity; read it per entity
    // dimension, in the layout's dimension-block order.
    const auto& layout =
      dynamic_cast<const pcms::MeshFieldsAdapterLayout&>(field.GetLayout());
    auto& mesh = layout.GetMesh();
    const auto nodes_per_dim = layout.GetNodesPerDim();
    const auto nc = static_cast<size_t>(layout.GetNumComponents());
    Kokkos::View<T**, pcms::HostMemorySpace> host(
      "values_with_ghosts_host",
      static_cast<size_t>(layout.GetNumDofHolderWithGhosts()), nc);
    size_t row = 0;
    for (int dim = 0; dim <= mesh.dim(); ++dim) {
      if (!nodes_per_dim[dim]) {
        continue;
      }
      const size_t rows = static_cast<size_t>(mesh.nents(dim)) *
                          static_cast<size_t>(nodes_per_dim[dim]);
      Kokkos::View<T*, pcms::DeviceMemorySpace> flat("backend_flat", rows * nc);
      mesh_fields->GetMeshFieldBackend()->GetData(
        pcms::Rank1View<T, pcms::DeviceMemorySpace>(flat.data(), rows * nc),
        nodes_per_dim[dim], nc, dim);
      const auto flat_h =
        Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), flat);
      for (size_t r = 0; r < rows; ++r) {
        for (size_t comp = 0; comp < nc; ++comp) {
          host(row + r, comp) = flat_h(r * nc + comp);
        }
      }
      row += rows;
    }
    return host;
  }
#endif
  throw std::runtime_error("ValuesWithGhostsHost: unsupported field data type");
}

} // namespace distributed_test

#endif // PCMS_TEST_DISTRIBUTED_TEST_UTILS_H
