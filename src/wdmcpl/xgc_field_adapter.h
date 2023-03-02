#ifndef WDM_COUPLING_XGC_FIELD_ADAPTER_H
#define WDM_COUPLING_XGC_FIELD_ADAPTER_H
#include "wdmcpl/types.h"
#include "wdmcpl/memory_spaces.h"
#include "wdmcpl/field.h"
#include <vector>
#include <redev_variant_tools.h>
#include "wdmcpl/xgc_reverse_classification.h"
#include "wdmcpl/assert.h"
#include "wdmcpl/array_mask.h"

namespace wdmcpl
{
namespace detail
{
// Needed since NVHPC doesn't work with overloaded
struct GetRank
{
  using GeomType = DimID;
  GetRank(const GeomType& geom) : geom_(geom) {}
  auto operator()(const redev::ClassPtn& ptn) const
  {
    const auto ent = redev::ClassPtn::ModelEnt({geom_.dim, geom_.id});
    return ptn.GetRank(ent);
  }
  auto operator()(const redev::RCBPtn& /*unused*/) const
  {
    std::cerr << "RCB partition not handled yet\n";
    std::terminate();
    return 0;
  }
  const GeomType& geom_;
};
} // namespace detail

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
                  ScalarArrayView<T, memory_space> data,
                  const ReverseClassificationVertex& reverse_classification,
                  std::function<int8_t(int, int)> in_overlap)
    : name_(std::move(name)),
      plane_comm_(plane_communicator),
      data_(data),
      gids_(data.size()),
      reverse_classification_(reverse_classification),
      in_overlap_(in_overlap)
  {
    // WDMCPL_ALWAYS_ASSERT(reverse_classification.nverts() == data.size());
    Kokkos::View<int8_t*, HostMemorySpace> mask("mask", data.size());
    WDMCPL_ALWAYS_ASSERT((bool)in_overlap);
    for (auto& geom : reverse_classification_) {
      if (in_overlap(geom.first.dim, geom.first.id)) {
        for (auto vert : geom.second) {
          WDMCPL_ALWAYS_ASSERT(vert < data.size());
          mask(vert) = 1;
        }
      }
    }
    mask_ = ArrayMask<memory_space>{make_const_array_view(mask)};
    WDMCPL_ALWAYS_ASSERT(!mask_.empty());
    //// XGC meshes are naively ordered in iteration order (full mesh on every
    //// cpu)
    std::iota(gids_.begin(), gids_.end(), static_cast<GO>(0));
    MPI_Comm_rank(plane_comm_, &plane_rank_);
  }

  int Serialize(
    ScalarArrayView<T, memory_space> buffer,
    ScalarArrayView<const wdmcpl::LO, memory_space> permutation) const
  {
    static_assert(std::is_same_v<memory_space, wdmcpl::HostMemorySpace>,
                  "gpu space unhandled\n");
    if (RankParticipatesCouplingCommunication()) {
      auto const_data = ScalarArrayView<const T, memory_space>{
        data_.data_handle(), data_.size()};
      if (buffer.size() > 0) {
        mask_.Apply(const_data, buffer, permutation);
      }
    }
    return mask_.Size();
  }
  void Deserialize(
    ScalarArrayView<const T, memory_space> buffer,
    ScalarArrayView<const wdmcpl::LO, memory_space> permutation) const
  {
    static_assert(std::is_same_v<memory_space, wdmcpl::HostMemorySpace>,
                  "gpu space unhandled\n");
    if (RankParticipatesCouplingCommunication()) {
      mask_.ToFullArray(buffer, data_, permutation);
    } else {
      // duplicate the data on the root rank of the plane to all other ranks
      MPI_Bcast(data_.data_handle(), data_.size(),
                redev::getMpiType(value_type{}), plane_root_, plane_comm_);
    }
  }

  // REQUIRED
  [[nodiscard]] std::vector<GO> GetGids() const
  {
    std::vector<GO> gids(mask_.Size());
    auto v1 = make_array_view(gids_);
    auto v2 = make_array_view(gids);
    mask_.Apply(v1, v2);
    return gids;
  }

  // REQUIRED
  [[nodiscard]] ReversePartitionMap GetReversePartitionMap(
    const redev::Partition& partition) const
  {
    wdmcpl::ReversePartitionMap reverse_partition;
    // in_overlap_ must contain a function!
    WDMCPL_ALWAYS_ASSERT(static_cast<bool>(in_overlap_));
    for (const auto& geom : reverse_classification_) {
      // if the geometry is in specified overlap region
      if (in_overlap_(geom.first.dim, geom.first.id)) {

        auto dr = std::visit(detail::GetRank{geom.first}, partition);
        auto [it, inserted] = reverse_partition.try_emplace(dr);
        // the map gives the local iteration order of the global ids
        auto map = mask_.GetMap();
        std::transform(geom.second.begin(), geom.second.end(),
                       std::back_inserter(it->second), [&map](auto v) {
                         auto idx = map[v];
                         WDMCPL_ALWAYS_ASSERT(idx > 0);
                         return idx - 1;
                       });
      }
    }
    return reverse_partition;
  }
  [[nodiscard]] bool RankParticipatesCouplingCommunication() const noexcept
  {
    // only do adios communications on 0 rank of the XGC fields
    return (plane_rank_ == plane_root_);
  }

private:
  std::string name_;
  MPI_Comm plane_comm_;
  int plane_rank_;
  ScalarArrayView<T, memory_space> data_;
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
  Kokkos::View<CoordinateElementType*,
               typename XGCFieldAdapter<T, CoordinateElementType>::memory_space>
    coordinates;
  return coordinates;
}
template <typename T, typename CoordinateElementType, typename MemorySpace>
auto evaluate(
  const XGCFieldAdapter<T, CoordinateElementType>& field,
  Lagrange<1> /* method */,
  ScalarArrayView<const CoordinateElementType, MemorySpace> coordinates)
  -> Kokkos::View<T*, MemorySpace>
{
  Kokkos::View<T*, MemorySpace> values("data", coordinates.size() / 2);
  std::cerr << "Evaluation of XGC Field not implemented yet!\n";
  std::abort();
  return values;
}
template <typename T, typename CoordinateElementType, typename MemorySpace>
auto evaluate(
  const XGCFieldAdapter<T, CoordinateElementType>& field,
  NearestNeighbor /* method */,
  ScalarArrayView<const CoordinateElementType, MemorySpace> coordinates)
  -> Kokkos::View<T*, MemorySpace>
{
  Kokkos::View<T*, MemorySpace> values("data", coordinates.size() / 2);
  std::cerr << "Evaluation of XGC Field not implemented yet!\n";
  std::abort();
  return values;
}

template <typename T, typename CoordinateElementType, typename U>
auto set_nodal_data(
  const XGCFieldAdapter<T, CoordinateElementType>& field,
  ScalarArrayView<
    const U, typename XGCFieldAdapter<T, CoordinateElementType>::memory_space>
    data) -> void
{
}

} // namespace wdmcpl

#endif // WDM_COUPLING_XGC_FIELD_ADAPTER_H
