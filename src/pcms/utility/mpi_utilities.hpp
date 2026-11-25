#ifndef PCMS_UTILITY_MPI_UTILITIES_HPP
#define PCMS_UTILITY_MPI_UTILITIES_HPP

#include "pcms/utility/mpi_type.hpp"

#include <mpi.h>

#include <cstddef>
#include <numeric>
#include <type_traits>
#include <vector>

namespace pcms
{

/**
 * The displacements matching per-rank counts, for the `v` collectives
 * (`MPI_Allgatherv`, `MPI_Alltoallv`).
 * \param counts how many entries each rank contributes, indexed by rank
 * \return the offset of each rank's block in the concatenated buffer
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
 * Gathers one scalar from every rank, on every rank. Collective.
 * \param value the value this rank contributes
 * \return one value per rank, indexed by rank
 */
template <typename T>
[[nodiscard]] std::vector<T> AllGatherAllRanks(MPI_Comm comm, const T& value)
{
  int size = 0;
  MPI_Comm_size(comm, &size);
  std::vector<T> gathered(static_cast<size_t>(size));
  MPI_Allgather(&value, 1, utility::GetMPIType(value), gathered.data(), 1,
                utility::GetMPIType(value), comm);
  return gathered;
}

/**
 * Gathers `count` scalars from every rank, on every rank. Collective.
 * \param values the `count` values this rank contributes
 * \param count how many values each rank contributes, the same on every rank
 * \return `count` values per rank, rank-major
 */
template <typename T>
[[nodiscard]] std::vector<T> AllGatherAllRanks(MPI_Comm comm, const T* values,
                                               int count)
{
  int size = 0;
  MPI_Comm_size(comm, &size);
  std::vector<T> gathered(static_cast<size_t>(size) *
                          static_cast<size_t>(count));
  MPI_Allgather(values, count, utility::GetMPIType(T{}), gathered.data(), count,
                utility::GetMPIType(T{}), comm);
  return gathered;
}

/**
 * Gathers a varying number of scalars from every rank, on every rank.
 * Collective.
 * \param values the values this rank contributes; may be empty
 * \param counts if not null, set to the number of values each rank sent
 * \return every rank's values, concatenated in rank order
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
                 utility::GetMPIType(T{}), gathered.data(), per_rank.data(),
                 displs.data(), utility::GetMPIType(T{}), comm);
  if (counts != nullptr) {
    *counts = per_rank;
  }
  return gathered;
}

/**
 * Sends `send_counts[r]` records to each rank `r` and receives
 * `recv_counts[r]` from it, when both sides already know the counts.
 * Collective.
 * \param send the records to send, grouped by destination rank in rank order
 * \return the records received, grouped by source rank in rank order
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
 * As ExchangeRecordsWithKnownCounts, first telling each rank how many records
 * to expect. Collective.
 * \param recv_counts set to the number of records received from each rank
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

#endif // PCMS_UTILITY_MPI_UTILITIES_HPP
