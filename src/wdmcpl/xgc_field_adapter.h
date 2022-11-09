#ifndef WDM_COUPLING_XGC_FIELD_ADAPTER_H
#define WDM_COUPLING_XGC_FIELD_ADAPTER_H
#include "wdmcpl/types.h"
#include "wdmcpl/memory_spaces.h"
#include "wdmcpl/field.h"
#include <vector>
#include <redev_variant_tools.h>

namespace wdmcpl
{

template <typename T, typename CoordinateElementType = Real>
class XGCFieldAdapter
{
public:
  using memory_space = HostMemorySpace;
  using value_type = T;
  using coordinate_element_type = CoordinateElementType;
  XGCFieldAdapter(std::string name,
                  std::vector<int8_t> classification_dimension,
                  std::vector<LO> classification_geometric_id,
                  ScalarArrayView<T, memory_space> data)
    : name_(std::move(name)),
      data_(data),
      gids_(data.size()),
      classification_dimension_(std::move(classification_dimension)),
      classification_geometric_id_(std::move(classification_geometric_id))
  {
    // XGC meshes are naively ordered in iteration order (full mesh on every
    // cpu)
    std::iota(gids_.begin(), gids_.end(), static_cast<GO>(0));
  }

  int Serialize(
    ScalarArrayView<T, memory_space> buffer,
    ScalarArrayView<const wdmcpl::LO, memory_space> permutation) const
  {
    static_assert(std::is_same_v<memory_space, wdmcpl::HostMemorySpace>,
                  "gpu space unhandled\n");
    if (buffer.size() > 0) {
      for (size_t i = 0, j = 0; i < data_.size(); i++) {
        buffer[permutation[j++]] = data_[i];
      }
    }
    return data_.size();
  }
  void Deserialize(
    ScalarArrayView<T, memory_space> buffer,
    ScalarArrayView<const wdmcpl::LO, memory_space> permutation) const
  {
    static_assert(std::is_same_v<memory_space, wdmcpl::HostMemorySpace>,
                  "gpu space unhandled\n");
    REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
    for (int i = 0; i < buffer.size(); ++i) {
      data_[i] = buffer[permutation[i]];
    }
  }

  // REQUIRED
  [[nodiscard]] const std::vector<GO>& GetGids() const { return gids_; }
  // REQUIRED
  [[nodiscard]] ReversePartitionMap GetReversePartitionMap(
    const redev::Partition& partition) const
  {
    wdmcpl::ReversePartitionMap reverse_partition;
    wdmcpl::LO local_index = 0;
    for (size_t i = 0; i < classification_dimension_.size(); ++i) {
      auto& class_dim = classification_dimension_;
      auto& class_id = classification_geometric_id_;
      auto dr =
        std::visit(redev::overloaded{
                     [&class_dim, &class_id, &i](const redev::ClassPtn& ptn) {
                       const auto ent =
                         redev::ClassPtn::ModelEnt({class_dim[i], class_id[i]});
                       return ptn.GetRank(ent);
                     },
                     [](const redev::RCBPtn& ptn) {
                       std::cerr << "RCB partition not handled yet\n";
                       std::terminate();
                       return 0;
                     }},
                   partition);
      reverse_partition[dr].emplace_back(local_index++);
    }
    return reverse_partition;
  }

private:
  std::string name_;
  ScalarArrayView<T, memory_space> data_;
  std::vector<GO> gids_;
  std::vector<int8_t> classification_dimension_;
  std::vector<LO> classification_geometric_id_;
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

} // namespace wdmcpl

#endif // WDM_COUPLING_XGC_FIELD_ADAPTER_H