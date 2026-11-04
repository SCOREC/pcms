#ifndef PCMS_TEST_DISTRIBUTED_TEST_UTILS_H
#define PCMS_TEST_DISTRIBUTED_TEST_UTILS_H

// Shared helpers for the MPI tests of distributed fields.
//
// These tests run under mpirun without Catch2, so a failed check aborts every
// rank: a rank that just returned would leave the others blocked in the next
// collective.

#include <mpi.h>

#include <Kokkos_Core.hpp>

#include "pcms/configuration.h"
#include "pcms/field/data/simple.h"
#ifdef PCMS_ENABLE_MESHFIELDS
#include "pcms/field/data/mesh_fields.h"
#endif
#include "pcms/field/field.h"
#include "pcms/field/value_view.hpp"

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <type_traits>

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

/// The values of every DOF holder with ghosts of `field`, owned and ghost,
/// copied to the host. Field itself only exposes the owned values, so this
/// reads the backend's storage.
template <typename T>
Kokkos::View<T**, pcms::HostMemorySpace> ValuesWithGhostsHost(
  const pcms::Field<T>& field)
{
  pcms::Rank2View<const T, pcms::DeviceMemorySpace> local;
  if (const auto* simple =
        dynamic_cast<const pcms::SimpleFieldData<T>*>(&field.GetData())) {
    local = simple->GetValuesWithGhosts();
  }
#ifdef PCMS_ENABLE_MESHFIELDS
  else if (const auto* mesh_fields =
             dynamic_cast<const pcms::MeshFieldsFieldData<T>*>(
               &field.GetData())) {
    local = mesh_fields->GetValuesWithGhosts();
  }
#endif
  else {
    throw std::runtime_error(
      "ValuesWithGhostsHost: unsupported field data type");
  }
  Kokkos::View<T**, pcms::DeviceMemorySpace> staged(
    "local_values_staged", local.extent(0), local.extent(1));
  pcms::CopyDeviceRank2ViewToRank2View(pcms::MakeRank2View(staged), local);
  Kokkos::View<T**, pcms::HostMemorySpace> host(
    "local_values_host", local.extent(0), local.extent(1));
  pcms::DeepCopyMismatchLayouts(host, staged);
  return host;
}

} // namespace distributed_test

#endif // PCMS_TEST_DISTRIBUTED_TEST_UTILS_H
