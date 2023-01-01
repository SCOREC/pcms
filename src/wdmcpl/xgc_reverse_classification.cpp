#include "wdmcpl/xgc_reverse_classification.h"
#include "mpi.h"
#include <fstream>
#include "wdmcpl/assert.h"
#include <string>
namespace wdmcpl
{

std::vector<LO> ReverseClassificationVertex::Serialize() const
{
  std::vector<LO> serialized_data;
  for (auto& geom : data_) {
    serialized_data.push_back(geom.first.dim);
    serialized_data.push_back(geom.first.id);
    serialized_data.push_back(geom.second.size());
    std::copy(geom.second.begin(), geom.second.end(),
              std::back_inserter(serialized_data));
  }
  return serialized_data;
}
void ReverseClassificationVertex::Deserialize(
  ScalarArrayView<LO, wdmcpl::HostMemorySpace> serialized_data)
{
  // expect to deserialize into an empty reverse classification class
  WDMCPL_ALWAYS_ASSERT(data_.empty());
  size_t i = 0;
  while (i < serialized_data.size()) {
    DimID geom{.dim = serialized_data[i], .id = serialized_data[i + 1]};
    auto nverts = serialized_data[i + 2];
    i += 3;
    for (size_t j = i; j < i + nverts; ++j) {
      data_[geom].insert(serialized_data[j]);
    }
    WDMCPL_ALWAYS_ASSERT(data_[geom].size() == static_cast<size_t>(nverts));
    i += nverts;
  }
}
ReverseClassificationVertex ReadReverseClassificationVertex(std::istream& in)
{
  LO total_nverts;
  ReverseClassificationVertex rc;
  in >> total_nverts;
  std::string line;
  while (!in.eof() && in.good()) {
    DimID geometry{};
    in >> geometry.dim;
    // check to make sure that we aren't reading past the end of the file.
    if (in.eof()) {
      break;
    }
    in >> geometry.id;
    getline(in,
            line); // read to the end of the cuurent line with geom_dim/geom_id
    getline(in, line);
    std::stringstream ss{line};
    LO node_id;
    while (ss >> node_id) {
      if (node_id >= 0) {
        // need to subtract 1 since XGC uses 1 based indexing, but we expect 0
        // based indexing in C code
        rc.Insert(geometry, node_id - 1);
      }
    }
  }
  return rc;
}
ReverseClassificationVertex ReadReverseClassificationVertex(std::istream& instr,
                                                            MPI_Comm comm,
                                                            int root)
{
  int rank = -1;
  MPI_Comm_rank(comm, &rank);
  if (rank == root) {
    auto rc = ReadReverseClassificationVertex(instr);
    auto serialized_rc = rc.Serialize();
    static_assert(std::is_same_v<LO, int32_t>,
                  "mpi type hardcoded, must update");
    size_t sz = serialized_rc.size();
    MPI_Bcast(&sz, 1, MPI_INT64_T, root, comm);
    MPI_Bcast(serialized_rc.data(), serialized_rc.size(), MPI_INT32_T, root,
              comm);
    return rc;
  } else {
    size_t sz = 0;
    MPI_Bcast(&sz, 1, MPI_INT64_T, root, comm);
    std::vector<LO> serialized_rc(sz);
    WDMCPL_ALWAYS_ASSERT(serialized_rc.size() == sz);
    MPI_Bcast(serialized_rc.data(), sz, MPI_INT32_T, root, comm);
    wdmcpl::ScalarArrayView<wdmcpl::LO, wdmcpl::HostMemorySpace> av{
      serialized_rc.data(), serialized_rc.size()};
    ReverseClassificationVertex rc;
    rc.Deserialize(av);
    return rc;
  }
}
ReverseClassificationVertex ReadReverseClassificationVertex(
  std::filesystem::path& classification_file)
{
  WDMCPL_ALWAYS_ASSERT(classification_file.has_filename());
  std::ifstream infile(classification_file);
  WDMCPL_ALWAYS_ASSERT(infile.is_open() && infile.good());
  return ReadReverseClassificationVertex(infile);
}

ReverseClassificationVertex ReadReverseClassificationVertex(
  std::filesystem::path& classification_file, MPI_Comm comm, int root)
{
  WDMCPL_ALWAYS_ASSERT(classification_file.has_filename());
  std::ifstream infile(classification_file);
  WDMCPL_ALWAYS_ASSERT(infile.is_open());
  return ReadReverseClassificationVertex(infile, comm, root);
}

void ReverseClassificationVertex::Insert(
  const DimID& key, ScalarArrayView<LO, wdmcpl::HostMemorySpace> data)
{
  // mdspan doesn't have begin currently. This should be switched
  // to range based for-loop
  for (int i = 0; i < data.extent(0); ++i) {
    Insert(key, data(i));
  }
}
void ReverseClassificationVertex::Insert(const DimID& key, LO data)
{
  data_[key].insert(data);
  ++total_verts_;
}

bool ReverseClassificationVertex::operator==(
  const ReverseClassificationVertex& other) const
{
  return data_ == other.data_;
}

const std::set<LO>* ReverseClassificationVertex::Query(
  const DimID& geometry) const noexcept
{
  auto it = data_.find(geometry);
  if (it != data_.end()) {
    return &(it->second);
  }
  return nullptr;
}
std::ostream& operator<<(std::ostream& os, const ReverseClassificationVertex& v)
{
  os << v.GetTotalVerts()<<"\n";
  for (auto& geom : v) {
    os << geom.first.dim << " " << geom.first.id << "\n";
    for (auto& vert : geom.second) {
      os << vert+1 << " ";
    }
    os << "\n";
  }
  return os;
}

} // namespace wdmcpl