#ifndef PCMS_MPI_UTILITIES_H
#define PCMS_MPI_UTILITIES_H

// Small wrappers over the MPI collectives used to exchange data between all
// ranks of a communicator.

#include "pcms/utility/mpi_type.h"

#include <mpi.h>

#include <cstddef>
#include <numeric>
#include <type_traits>
#include <vector>

namespace pcms
{

// ---- gathers ----

/**
 * \brief The displacements matching a per-rank count vector, for the `v`
 * collectives (`MPI_Allgatherv`, `MPI_Alltoallv`).
 *
 * \param counts How many entries each rank contributes, indexed by rank.
 * \return The displacement of each rank's block in the concatenated buffer.
 */
[[nodiscard]] inline std::vector<int> DisplacementsOf(
  const std::vector<int>& counts)
{
  std::vector<int> displs(counts.size(), 0);
  for (size_t r = 1; r < counts.size(); ++r) {
    displs[r] = displs[r - 1] + counts[r - 1];
  }
  return displs;
}

/**
 * \brief Gather one scalar from every rank, indexed by rank, on every rank.
 *
 * There is no Omega_h equivalent to reach for here. `Omega_h::Comm::allgather`
 * is a *neighbor* gather: it exchanges with the ranks of the mesh's
 * communication graph (`Comm::sources()`/`destinations()`) rather than with all
 * ranks, and returns one entry per neighbor, ordered by that neighbor list. It
 * is therefore wrong -- and quietly wrong, returning a short buffer -- for
 * anything keyed by rank.
 *
 * `Omega_h::Comm::allreduce` does not have this problem: it is a plain
 * `MPI_Allreduce` on the communicator. So a quantity that is the same on every
 * rank can go through Omega_h, while a per-rank quantity cannot.
 *
 * \param comm The communicator, for example `mesh.comm()->get_impl()`.
 * \param value The value this rank contributes.
 * \return `size`-many values, where element `r` is rank `r`'s contribution.
 */
template <typename T>
[[nodiscard]] std::vector<T> AllGatherAllRanks(MPI_Comm comm, const T& value)
{
  int size = 0;
  MPI_Comm_size(comm, &size);
  std::vector<T> gathered(static_cast<size_t>(size));
  MPI_Allgather(&value, 1, GetMPIType(value), gathered.data(), 1,
                GetMPIType(value), comm);
  return gathered;
}

/**
 * \brief Gather `count` scalars from every rank, each rank's block contiguous.
 *
 * The result is rank-major: rank `r`'s values are at `[r * count, (r+1) *
 * count)`. Use this over the single-value overload when a rank contributes
 * several values, so one exchange replaces several.
 *
 * \param comm The communicator, for example `mesh.comm()->get_impl()`.
 * \param values `count` values this rank contributes.
 * \param count How many values each rank contributes.
 * \return `size * count` values, rank-major.
 */
template <typename T>
[[nodiscard]] std::vector<T> AllGatherAllRanks(MPI_Comm comm, const T* values,
                                               int count)
{
  int size = 0;
  MPI_Comm_size(comm, &size);
  std::vector<T> gathered(static_cast<size_t>(size) *
                          static_cast<size_t>(count));
  MPI_Allgather(values, count, GetMPIType(T{}), gathered.data(), count,
                GetMPIType(T{}), comm);
  return gathered;
}

/**
 * \brief Gather a variable number of scalars from every rank, concatenated in
 * rank order, on every rank.
 *
 * \param comm The communicator, for example `mesh.comm()->get_impl()`.
 * \param values The values this rank contributes; may be empty.
 * \param counts Optional. Filled with the number of values each rank sent.
 * \return All ranks' values, rank 0's first.
 */
template <typename T>
[[nodiscard]] std::vector<T> AllGatherVAllRanks(
  MPI_Comm comm, const std::vector<T>& values,
  std::vector<int>* counts = nullptr)
{
  const auto per_rank =
    AllGatherAllRanks(comm, static_cast<int>(values.size()));
  const auto displs = DisplacementsOf(per_rank);
  std::vector<T> gathered(
    static_cast<size_t>(std::accumulate(per_rank.begin(), per_rank.end(), 0)));
  MPI_Allgatherv(values.data(), static_cast<int>(values.size()),
                 GetMPIType(T{}), gathered.data(), per_rank.data(),
                 displs.data(), GetMPIType(T{}), comm);
  if (counts != nullptr) {
    *counts = per_rank;
  }
  return gathered;
}

// ---- all-to-all exchanges ----

/**
 * \brief Exchange records whose per-rank counts are already known on both
 * sides.
 *
 * `send` holds the records for rank 0 first, then rank 1, and so on, with
 * `send_counts[r]` records for rank `r`. The result is laid out the same way,
 * by source rank, with `recv_counts[r]` records from rank `r`.
 *
 * Records are sent as raw bytes, so `Rec` must be trivially copyable. Every
 * rank runs the same binary, so struct padding only costs a few bytes.
 */
template <typename Rec>
[[nodiscard]] std::vector<Rec> ExchangeRecordsWithKnownCounts(
  MPI_Comm comm, const std::vector<int>& send_counts,
  const std::vector<Rec>& send, const std::vector<int>& recv_counts)
{
  static_assert(std::is_trivially_copyable_v<Rec>,
                "ExchangeRecords sends records as raw bytes");
  const auto to_bytes = [](const std::vector<int>& counts) {
    std::vector<int> bytes(counts.size());
    for (size_t r = 0; r < counts.size(); ++r) {
      bytes[r] = counts[r] * static_cast<int>(sizeof(Rec));
    }
    return bytes;
  };
  const auto send_bytes = to_bytes(send_counts);
  const auto recv_bytes = to_bytes(recv_counts);
  const auto send_displs = DisplacementsOf(send_bytes);
  const auto recv_displs = DisplacementsOf(recv_bytes);

  std::vector<Rec> recv(static_cast<size_t>(
    std::accumulate(recv_counts.begin(), recv_counts.end(), 0)));
  MPI_Alltoallv(send.data(), send_bytes.data(), send_displs.data(), MPI_BYTE,
                recv.data(), recv_bytes.data(), recv_displs.data(), MPI_BYTE,
                comm);
  return recv;
}

/**
 * \brief Exchange records, first telling each rank how many to expect.
 *
 * Same layout as `ExchangeRecordsWithKnownCounts`. `recv_counts` is filled
 * with the number of records received from each rank.
 */
template <typename Rec>
[[nodiscard]] std::vector<Rec> ExchangeRecords(
  MPI_Comm comm, const std::vector<int>& send_counts,
  const std::vector<Rec>& send, std::vector<int>& recv_counts)
{
  recv_counts.assign(send_counts.size(), 0);
  MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT,
               comm);
  return ExchangeRecordsWithKnownCounts(comm, send_counts, send, recv_counts);
}

} // namespace pcms

#endif // PCMS_MPI_UTILITIES_H
