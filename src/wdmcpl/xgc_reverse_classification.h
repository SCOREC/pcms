#ifndef WDM_COUPLING_XGC_REVERSE_CLASSIFICATION_H
#define WDM_COUPLING_XGC_REVERSE_CLASSIFICATION_H
#include <Kokkos_Core.hpp>
#include <mpi.h>
#include "wdmcpl/types.h"
#include <unordered_map>
#include <set>
#include <wdmcpl/external/mdspan.hpp>
#include <wdmcpl/arrays.h>
#include <wdmcpl/memory_spaces.h>
#include <filesystem>

namespace wdmcpl
{
struct DimID
{
  LO dim;
  LO id;
  bool operator==(const DimID& other) const
  {
    return (dim == other.dim) && (id == other.id);
  }
};
} // namespace wdmcpl
namespace std
{
template <>
struct hash<wdmcpl::DimID>
{
  std::size_t operator()(const wdmcpl::DimID& key) const noexcept
  {
    auto h1 = hash<wdmcpl::LO>{}(key.dim);
    auto h2 = hash<wdmcpl::LO>{}(key.id);
    return h1 ^ (h2 << 1); // hash combine see
                           // https://en.cppreference.com/w/cpp/utility/hash
  }
};
} // namespace std
namespace wdmcpl
{
///
/// This datastructure represents the reverse classification
/// of the mesh verticies on the geometric entities
class ReverseClassificationVertex
{
public:
  // the ordered nature of the set is relied upon to provide the entries for
  // each geometric entity in iteration order (for xgc) where the gids are in
  // ascending order
  using DataMapType = std::unordered_map<DimID, std::set<LO>>;
  void Insert(const DimID& key,
              ScalarArrayView<LO, wdmcpl::HostMemorySpace> data)
  {
    // mdspan doesn't have begin currently. This should be switched
    // to range based for-loop
    for (int i = 0; i < data.extent(0); ++i) {
      Insert(key, data(i));
    }
  }
  void Insert(const DimID& key, LO data) { data_[key].insert(data); }
  std::vector<LO> Serialize();
  void Deserialize(
    ScalarArrayView<LO, wdmcpl::HostMemorySpace> serialized_data);
  bool operator==(const ReverseClassificationVertex& other) const
  {
    return data_ == other.data_;
  }
  const std::set<LO>* Query(const DimID& geometry) const noexcept
  {
    auto it = data_.find(geometry);
    if (it != data_.end()) {
      return &(it->second);
    }
    return nullptr;
  }
  DataMapType::iterator begin() noexcept { return data_.begin(); }
  DataMapType::iterator end() noexcept { return data_.end(); }
  DataMapType::const_iterator begin() const noexcept { return data_.begin(); }
  DataMapType::const_iterator end() const noexcept { return data_.end(); }

private:
  DataMapType data_;
};

ReverseClassificationVertex ReadReverseClassificationVertex(
  std::filesystem::path&);
ReverseClassificationVertex ReadReverseClassificationVertex(std::istream&);
ReverseClassificationVertex ReadReverseClassificationVertex(std::istream&,
                                                            MPI_Comm,
                                                            int root = 0);
ReverseClassificationVertex ReadReverseClassificationVertex(
  std::filesystem::path&, MPI_Comm, int root = 0);

} // namespace wdmcpl

#endif // WDM_COUPLING_XGC_REVERSE_CLASSIFICATION_H
