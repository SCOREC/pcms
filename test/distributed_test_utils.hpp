#ifndef PCMS_TEST_DISTRIBUTED_TEST_UTILS_HPP
#define PCMS_TEST_DISTRIBUTED_TEST_UTILS_HPP

// Shared helpers for the MPI tests of distributed fields.
//
// These tests run under mpirun without Catch2, so a failed check aborts every
// rank: a rank that just returned would leave the others blocked in the next
// collective.

#include <mpi.h>

#include <Kokkos_Core.hpp>

#include "pcms/configuration.hpp"
#include "pcms/backends/omega_h/field_data.hpp"
#include "pcms/field/detail/simple_field_data.hpp"
#ifdef PCMS_ENABLE_OMEGA_H_MESHFIELDS
#include "pcms/backends/omega_h_meshfields/field_data.hpp"
#endif
#include "pcms/field/field.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/value_view.hpp"

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
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
        dynamic_cast<const pcms::omega_h::OmegaHFieldData<T>*>(
          &field.GetData())) {
    const auto& layout = field.GetLayout();
    const auto nc = static_cast<size_t>(layout.GetNumComponents());
    // The values of every DOF holder are in the field's Omega_h tag.
    using TagValue = typename pcms::omega_h::OmegaHFieldData<T>::tag_value_type;
    const auto& geometry =
      dynamic_cast<const pcms::omega_h::OmegaHGeometry&>(*layout.GetGeometry());
    const auto tag =
      Omega_h::HostRead<TagValue>(geometry.GetMesh().get_array<TagValue>(
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
#ifdef PCMS_ENABLE_OMEGA_H_MESHFIELDS
  if (const auto* mesh_fields =
        dynamic_cast<const pcms::omega_h_meshfields::MeshFieldsFieldData<T>*>(
          &field.GetData())) {
    // The MeshFields backend stores every local entity; read it per entity
    // dimension, in the layout's dimension-block order.
    const auto& layout =
      dynamic_cast<const pcms::omega_h_meshfields::MeshFieldsAdapterLayout&>(
        field.GetLayout());
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

/// @return a function linear in the coordinates, a different one per
///         component
inline Real Linear(const Real* x, int dim, int c)
{
  Real value = 1.0 + c;
  for (int d = 0; d < dim; ++d) {
    value += (d + 1.0 + 0.5 * c) * x[d];
  }
  return value;
}

/// @return a function of `space` set at its owned DOF holders to Linear,
///         synchronized
inline pcms::Function<Real> LinearField(const pcms::FunctionSpace& space)
{
  const auto& layout = *space.GetLayout();
  const auto coords = layout.GetDOFHolderCoordinates().GetValues();
  const int dim = static_cast<int>(coords.extent(1));
  Kokkos::View<Real**, pcms::DeviceMemorySpace> device(
    "holder_coordinates", layout.GetNumDofHolder(), dim);
  Kokkos::parallel_for(
    "distributed_test::LinearField",
    Kokkos::MDRangePolicy<pcms::DefaultExecutionSpace, Kokkos::Rank<2>>(
      {0, 0},
      {static_cast<int64_t>(device.extent(0)), static_cast<int64_t>(dim)}),
    KOKKOS_LAMBDA(const int64_t i, const int64_t d) {
      device(i, d) = coords(i, d);
    });
  const auto holders =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, device);
  const int nc = layout.GetNumComponents();
  Kokkos::View<Real**, Kokkos::LayoutRight, pcms::HostMemorySpace> values(
    "linear_values", holders.extent(0), nc);
  std::vector<Real> x(static_cast<std::size_t>(dim));
  for (std::size_t h = 0; h < holders.extent(0); ++h) {
    for (int d = 0; d < dim; ++d) {
      x[static_cast<std::size_t>(d)] = holders(h, d);
    }
    for (int c = 0; c < nc; ++c) {
      values(h, c) = Linear(x.data(), dim, c);
    }
  }
  auto field = space.CreateFunction<Real>();
  field.SetDOFHolderDataHost(pcms::ValueView<const Real, pcms::HostMemorySpace>(
    field.GetData().GetValueBasis(), pcms::MakeConstRank2View(values)));
  field.Synchronize();
  return field;
}

} // namespace distributed_test

#endif // PCMS_TEST_DISTRIBUTED_TEST_UTILS_HPP
