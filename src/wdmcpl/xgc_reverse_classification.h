#ifndef WDM_COUPLING_XGC_REVERSE_CLASSIFICATION_H
#define WDM_COUPLING_XGC_REVERSE_CLASSIFICATION_H
#include <Kokkos_Core.hpp>
#include <mpi.h>
#include "wdmcpl/types.h"
#include <unordered_map>
#include <set>
#include "wdmcpl/external/mdspan.hpp"
#include "wdmcpl/arrays.h"
#include "wdmcpl/memory_spaces.h"
//#include <filesystem>
#ifdef WDMCPL_HAS_OMEGA_H
#include <Omega_h_mesh.hpp>
#include "wdmcpl/assert.h"
#endif

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
              ScalarArrayView<LO, wdmcpl::HostMemorySpace> data);
  void Insert(const DimID& key, LO data);
  [[nodiscard]] std::vector<LO> Serialize() const;
  void Deserialize(
    ScalarArrayView<LO, wdmcpl::HostMemorySpace> serialized_data);
  [[nodiscard]] bool operator==(const ReverseClassificationVertex& other) const;
  [[nodiscard]] const std::set<LO>* Query(const DimID& geometry) const noexcept;
  [[nodiscard]] DataMapType::iterator begin() noexcept { return data_.begin(); }
  [[nodiscard]] DataMapType::iterator end() noexcept { return data_.end(); }
  [[nodiscard]] DataMapType::const_iterator begin() const noexcept
  {
    return data_.begin();
  }
  [[nodiscard]] DataMapType::const_iterator end() const noexcept
  {
    return data_.end();
  }
  [[nodiscard]] LO GetTotalVerts() const noexcept { return total_verts_; }
  friend std::ostream& operator<<(std::ostream& os,
                                  const ReverseClassificationVertex& v);

private:
  DataMapType data_;
  LO total_verts_{0};
};

ReverseClassificationVertex ReadReverseClassificationVertex(std::string);
ReverseClassificationVertex ReadReverseClassificationVertex(std::istream&);
ReverseClassificationVertex ReadReverseClassificationVertex(std::istream&,
                                                            MPI_Comm,
                                                            int root = 0);
ReverseClassificationVertex ReadReverseClassificationVertex(std::string, MPI_Comm, int root = 0);

#ifdef WDMCPL_HAS_OMEGA_H
enum class IndexBase {
  Zero = 0,
  One = 1
};
template <typename T = Omega_h::LO>
[[nodiscard]] ReverseClassificationVertex ConstructRCFromOmegaHMesh(
  Omega_h::Mesh& mesh, std::string numbering = "simNumbering", IndexBase index_base=IndexBase::One)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // transfer vtx classification to host
  auto classIds_h =
    Omega_h::HostRead<Omega_h::ClassId>(mesh.get_array<Omega_h::ClassId>(0, "class_id"));
  auto classDims_h =
    Omega_h::HostRead<Omega_h::I8>(mesh.get_array<Omega_h::I8>(0, "class_dim"));
  auto vertid = Omega_h::HostRead<T>(mesh.get_array<T>(0, numbering));
  wdmcpl::ReverseClassificationVertex rc;
  WDMCPL_ALWAYS_ASSERT(classDims_h.size() == classIds_h.size());
  for (int i = 0; i < classDims_h.size(); ++i) {
    wdmcpl::DimID geom{classDims_h[i], classIds_h[i]};
    if(index_base == IndexBase::Zero) {
      rc.Insert(geom, vertid[i]);
    } else {
      WDMCPL_ALWAYS_ASSERT(vertid[i] > 0);
      rc.Insert(geom, vertid[i] - 1);
    }
  }
  return rc;
}

#endif
} // namespace wdmcpl

#endif // WDM_COUPLING_XGC_REVERSE_CLASSIFICATION_H
