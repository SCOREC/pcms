#ifndef PCMS_COUPLING_XGC_FIELD_ADAPTER_H
#define PCMS_COUPLING_XGC_FIELD_ADAPTER_H
#include "pcms/adapter/omega_h/omega_h_field.h"
#include "pcms/types.h"
#include "pcms/memory_spaces.h"
#include "pcms/field.h"
#include <vector>
#include <redev_variant_tools.h>
#include "xgc_reverse_classification.h"
#include "pcms/assert.h"
#include "pcms/array_mask.h"
#include "pcms/profile.h"

namespace pcms
{
template <typename T, typename CoordinateElementType = Real>
class XGCFieldAdapter
{
public:
  using memory_space = HostMemorySpace;
  using value_type = T;
  using coordinate_element_type = CoordinateElementType;
  /**
   *
   * @param name name of the field
   * @param plane_communicator the communicator of all ranks corresponding to a
   * given XGC plane. This corresponds to sml_plane_comm
   * @param data a view of the data to be used as the field definition
   * @param reverse_classification the reverse classification data for the XGC
   * field
   * @param in_overlap a function describing if an entity defined by the
   * geometric dimension and ID
   */
  XGCFieldAdapter(std::string name, MPI_Comm plane_communicator,
                  Rank1View<T, memory_space> data,
                  const ReverseClassificationVertex& reverse_classification,
                  std::function<int8_t(int, int)> in_overlap)
    : name_(std::move(name)),
      plane_comm_(plane_communicator),
      data_(data),
      gids_(data.size()),
      reverse_classification_(reverse_classification),
      in_overlap_(in_overlap)
  {
    PCMS_FUNCTION_TIMER;
    // PCMS_ALWAYS_ASSERT(reverse_classification.nverts() == data.size());
    MPI_Comm_rank(plane_comm_, &plane_rank_);
    if (RankParticipatesCouplingCommunication()) {
      Kokkos::View<int8_t*, HostMemorySpace> mask("mask", data.size());
      PCMS_ALWAYS_ASSERT((bool)in_overlap);
      for (auto& geom : reverse_classification_) {
        if (in_overlap(geom.first.dim, geom.first.id)) {
          for (auto vert : geom.second) {
            PCMS_ALWAYS_ASSERT(vert < data.size());
            mask(vert) = 1;
          }
        }
      }
      mask_ = ArrayMask<memory_space>{make_const_array_view(mask)};
      PCMS_ALWAYS_ASSERT(!mask_.empty());
      //// XGC meshes are naively ordered in iteration order (full mesh on every
      //// cpu) First ID in XGC is 1!
      std::iota(gids_.begin(), gids_.end(), static_cast<GO>(1));
    }
  }

  int Serialize(Rank1View<T, memory_space> buffer,
                Rank1View<const pcms::LO, memory_space> permutation) const
  {
    PCMS_FUNCTION_TIMER;
    static_assert(std::is_same_v<memory_space, pcms::HostMemorySpace>,
                  "gpu space unhandled\n");
    if (RankParticipatesCouplingCommunication()) {
      auto const_data =
        Rank1View<const T, memory_space>{data_.data_handle(), data_.size()};
      if (buffer.size() > 0) {
        mask_.Apply(const_data, buffer, permutation);
      }
      return mask_.Size();
    }
    return 0;
  }
  void Deserialize(Rank1View<const T, memory_space> buffer,
                   Rank1View<const pcms::LO, memory_space> permutation) const
  {
    PCMS_FUNCTION_TIMER;
    static_assert(std::is_same_v<memory_space, pcms::HostMemorySpace>,
                  "gpu space unhandled\n");
    if (RankParticipatesCouplingCommunication()) {
      mask_.ToFullArray(buffer, data_, permutation);
    }
    // duplicate the data on the root rank of the plane to all other ranks
    MPI_Bcast(data_.data_handle(), data_.size(),
              redev::getMpiType(value_type{}), plane_root_, plane_comm_);
  }

  // REQUIRED
  [[nodiscard]] std::vector<GO> GetGids() const
  {
    PCMS_FUNCTION_TIMER;
    if (RankParticipatesCouplingCommunication()) {
      std::vector<GO> gids(mask_.Size());
      auto v1 = make_array_view(gids_);
      auto v2 = make_array_view(gids);
      mask_.Apply(v1, v2);
      return gids;
    }
    return {};
  }

  // REQUIRED
  [[nodiscard]] ReversePartitionMap GetReversePartitionMap(
    const Partition& partition) const
  {
    PCMS_FUNCTION_TIMER;
    if (RankParticipatesCouplingCommunication()) {

      pcms::ReversePartitionMap reverse_partition;
      // in_overlap_ must contain a function!
      PCMS_ALWAYS_ASSERT(static_cast<bool>(in_overlap_));
      for (const auto& geom : reverse_classification_) {
        // if the geometry is in specified overlap region
        if (in_overlap_(geom.first.dim, geom.first.id)) {

          auto dr = partition.GetDr(geom.first.id, geom.first.dim);
          auto [it, inserted] = reverse_partition.try_emplace(dr);
          // the map gives the local iteration order of the global ids
          auto map = mask_.GetMap();
          std::transform(geom.second.begin(), geom.second.end(),
                         std::back_inserter(it->second), [&map](auto v) {
                           auto idx = map[v];
                           PCMS_ALWAYS_ASSERT(idx > 0);
                           return idx - 1;
                         });
        }
      }

      // Rather than convert an explicit forward classification,
      // we can construct the reverse partitionbased on the geometry
      // and sort the node ids after to get the iteration order correct
      // in XGC the local iteration order maps directly to the global ids
      for (auto& [rank, idxs] : reverse_partition) {
        std::sort(idxs.begin(), idxs.end());
      }
      return reverse_partition;
    }
    return {};
  }
  [[nodiscard]] bool RankParticipatesCouplingCommunication() const noexcept
  {
    PCMS_FUNCTION_TIMER;
    // only do adios communications on 0 rank of the XGC fields
    return (plane_rank_ == plane_root_);
  }

  [[nodiscard]] pcms::mesh_entity_type GetEntityType() const noexcept
  {
    return pcms::mesh_entity_type::VERTEX;
  }

private:
  std::string name_;
  MPI_Comm plane_comm_;
  int plane_rank_;
  Rank1View<T, memory_space> data_;
  std::vector<GO> gids_;
  const ReverseClassificationVertex& reverse_classification_;
  std::function<int8_t(int, int)> in_overlap_;
  ArrayMask<memory_space> mask_;
  static constexpr int plane_root_{0};
};

struct ReadXGCNodeClassificationResult
{
  std::vector<int8_t> dimension;
  std::vector<LO> geometric_id;
};

/**
 *
 *
 * @param in istream input. The input should be in XGC node iteration order
 * (implicit numbering). Each line of the input should have have a dimension and
 * geometric id that the node is classified on.
 *
 * @return the classification for each node in the XGC mesh
 */
[[nodiscard]] ReadXGCNodeClassificationResult ReadXGCNodeClassification(
  std::istream& in);

template <typename T, typename CoordinateElementType>
auto get_nodal_coordinates(
  const XGCFieldAdapter<T, CoordinateElementType>& field)
{
  PCMS_FUNCTION_TIMER;
  Kokkos::View<CoordinateElementType*,
               typename XGCFieldAdapter<T, CoordinateElementType>::memory_space>
    coordinates;
  return coordinates;
}
template <typename T, typename CoordinateElementType, typename MemorySpace>
auto evaluate(const XGCFieldAdapter<T, CoordinateElementType>& field,
              Lagrange<1> /* method */,
              Rank1View<const CoordinateElementType, MemorySpace> coordinates)
  -> Kokkos::View<T*, MemorySpace>
{
  PCMS_FUNCTION_TIMER;
  Kokkos::View<T*, MemorySpace> values("data", coordinates.size() / 2);
  std::cerr << "Evaluation of XGC Field not implemented yet!\n";
  std::abort();
  return values;
}
template <typename T, typename CoordinateElementType, typename MemorySpace>
auto evaluate(const XGCFieldAdapter<T, CoordinateElementType>& field,
              NearestNeighbor /* method */,
              Rank1View<const CoordinateElementType, MemorySpace> coordinates)
  -> Kokkos::View<T*, MemorySpace>
{
  PCMS_FUNCTION_TIMER;
  Kokkos::View<T*, MemorySpace> values("data", coordinates.size() / 2);
  std::cerr << "Evaluation of XGC Field not implemented yet!\n";
  std::abort();
  return values;
}

template <typename T, typename CoordinateElementType, typename U>
auto set_nodal_data(
  const XGCFieldAdapter<T, CoordinateElementType>& field,
  Rank1View<const U,
            typename XGCFieldAdapter<T, CoordinateElementType>::memory_space>
    data) -> void
{
  PCMS_FUNCTION_TIMER;
}

} // namespace pcms

#endif // PCMS_COUPLING_XGC_FIELD_ADAPTER_H
