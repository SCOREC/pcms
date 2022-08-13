#include <Omega_h_mesh.hpp>
#include <iostream>
#include <wdmcpl.h>
#include <wdmcpl/types.h>
#include <Omega_h_file.hpp>
#include <Omega_h_for.hpp>
#include <redev_variant_tools.h>
#include "test_support.h"

/**
 * return 1 if the specificed model entity is part of the overlap region, 0
 * otherwise
 */
OMEGA_H_DEVICE Omega_h::I8 isModelEntInOverlap(const int dim, const int id)
{
  // the TOMMS generated geometric model has
  // entity IDs that increase with the distance
  // from the magnetic axis
  if (dim == 2 && (id >= 22 && id <= 34)) {
    return 1;
  } else if (dim == 1 && (id >= 21 && id <= 34)) {
    return 1;
  } else if (dim == 0 && (id >= 21 && id <= 34)) {
    return 1;
  }
  return 0;
}

/**
 * Create the tag 'isOverlap' for each mesh vertex whose value is 1 if the
 * vertex is classified on a model entity in the closure of the geometric model
 * faces forming the overlap region; the value is 0 otherwise.
 * OnlyIncludesOverlapping and owned verts
 */
Omega_h::Read<Omega_h::I8> markOverlapMeshEntities(Omega_h::Mesh& mesh)
{
  // transfer vtx classification to host
  auto classIds = mesh.get_array<Omega_h::ClassId>(0, "class_id");
  auto classDims = mesh.get_array<Omega_h::I8>(0, "class_dim");
  auto isOverlap = Omega_h::Write<Omega_h::I8>(classIds.size(), "isOverlap");
  auto markOverlap = OMEGA_H_LAMBDA(int i)
  {
    isOverlap[i] = isModelEntInOverlap(classDims[i], classIds[i]);
  };
  Omega_h::parallel_for(classIds.size(), markOverlap);
  auto isOwned = mesh.owned(0);
  // try masking out to only owned entities
  Omega_h::parallel_for(
    isOverlap.size(),
    OMEGA_H_LAMBDA(int i) { isOverlap[i] = (isOwned[i] && isOverlap[i]); });

  auto isOverlap_r = Omega_h::read(isOverlap);
  mesh.add_tag(0, "isOverlap", 1, isOverlap_r);
  return isOverlap_r;
}
Omega_h::HostRead<Omega_h::I8> markMeshOverlapRegion(Omega_h::Mesh& mesh)
{
  auto isOverlap = markOverlapMeshEntities(mesh);
  return Omega_h::HostRead(isOverlap);
}

struct SerializeOmegaHGids
{
  SerializeOmegaHGids(Omega_h::Mesh& mesh,
                      Omega_h::HostRead<Omega_h::I8> is_overlap_h)
    : mesh_(mesh), is_overlap_h_(is_overlap_h)
  {
  }
  template <typename T>
  int operator()(std::string_view, nonstd::span<T> buffer,
                 nonstd::span<const wdmcpl::LO> permutation) const
  {
    // WDMCPL_ALWAYS_ASSERT(buffer.size() == is_overlap_h_.size());
    auto gids = mesh_.globals(0);
    auto gids_h = Omega_h::HostRead(gids);
    int count = 0;
    for (size_t i = 0, j = 0; i < gids_h.size(); i++) {
      if (is_overlap_h_[i]) {
        if (buffer.size() > 0) {
          buffer[permutation[j++]] = gids_h[i];
        }
        ++count;
      }
    }
    return count;
  }
  Omega_h::Mesh mesh_;
  Omega_h::HostRead<Omega_h::I8> is_overlap_h_;
};

// Serializer is used in a two pass algorithm. Must check that the buffer size
// >0 and return the number of entries.
struct SerializeOmegaH
{
  SerializeOmegaH(Omega_h::Mesh& mesh,
                  Omega_h::HostRead<Omega_h::I8> is_overlap_h)
    : mesh_(mesh), is_overlap_h_(is_overlap_h)
  {
  }
  template <typename T>
  int operator()(std::string_view name, nonstd::span<T> buffer,
                 nonstd::span<const wdmcpl::LO> permutation) const
  {
    // WDMCPL_ALWAYS_ASSERT(buffer.size() == is_overlap_h_.size());
    const auto array = mesh_.get_array<T>(0, std::string(name));
    const auto array_h = Omega_h::HostRead(array);
    int count = 0;
    for (size_t i = 0, j = 0; i < array_h.size(); i++) {
      if (is_overlap_h_[i]) {
        if (buffer.size() > 0) {
          buffer[permutation[j++]] = array_h[i];
        }
        ++count;
      }
    }
    return count;
  }
  Omega_h::Mesh mesh_;
  Omega_h::HostRead<Omega_h::I8> is_overlap_h_;
};
struct DeserializeOmegaH
{
  DeserializeOmegaH(Omega_h::Mesh& mesh,
                    Omega_h::HostRead<Omega_h::I8> is_overlap_h)
    : mesh_(mesh), is_overlap_h_(is_overlap_h)
  {
  }
  template <typename T>
  void operator()(std::string_view, nonstd::span<const T> buffer,
                  nonstd::span<const wdmcpl::LO> permutation) const
  {

    REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
    auto gids = mesh_.globals(0);
    auto gids_h = Omega_h::HostRead(gids);
    std::vector<wdmcpl::GO> global_ids;
    for (size_t i = 0, j = 0; i < gids_h.size(); i++) {
      if (is_overlap_h_[i]) {
        REDEV_ALWAYS_ASSERT(gids_h[i] == buffer[permutation[j++]]);
      }
    }
  }

private:
  Omega_h::Mesh& mesh_;
  Omega_h::HostRead<Omega_h::I8> is_overlap_h_;
};

struct OmegaHGids
{
  OmegaHGids(Omega_h::Mesh& mesh, Omega_h::HostRead<Omega_h::I8> is_overlap_h)
    : mesh_(mesh), is_overlap_h_(is_overlap_h)
  {
  }
  std::vector<wdmcpl::GO> operator()(std::string_view) const
  {
    auto gids = mesh_.globals(0); // GPU
    auto gids_h = Omega_h::HostRead(gids); // CPU
    std::vector<wdmcpl::GO> global_ids;
    for (size_t i = 0; i < gids_h.size(); i++) {
      if (is_overlap_h_[i]) {
        global_ids.push_back(gids_h[i]);
      }
    }
    return global_ids;
  }
  Omega_h::Mesh& mesh_;
  Omega_h::HostRead<Omega_h::I8> is_overlap_h_;
};

struct OmegaHReversePartition
{
  OmegaHReversePartition(Omega_h::Mesh& mesh) : mesh(mesh) {}
  wdmcpl::ReversePartitionMap operator()(std::string_view,
                                         const redev::Partition& partition) const
  {
    auto ohComm = mesh.comm();
    const auto rank = ohComm->rank();
    // transfer vtx classification to host
    auto classIds = mesh.get_array<Omega_h::ClassId>(0, "class_id");
    auto classIds_h = Omega_h::HostRead(classIds);
    auto classDims = mesh.get_array<Omega_h::I8>(0, "class_dim");
    auto classDims_h = Omega_h::HostRead(classDims);
    auto isOverlap =
      mesh.has_tag(0, "isOverlap")
        ? mesh.get_array<Omega_h::I8>(0, "isOverlap")
        : Omega_h::Read<Omega_h::I8>(
            classIds.size(), 1, "isOverlap"); // no mask for overlap vertices
    auto isOverlap_h = Omega_h::HostRead(isOverlap);
    // local_index number of vertices going to each destination process by
    // calling getRank - degree array
    wdmcpl::ReversePartitionMap reverse_partition;
    wdmcpl::LO local_index = 0;
    for (auto i = 0; i < classIds_h.size(); i++) {
      if (isOverlap_h[i]) {
        auto dr = std::visit(
          redev::overloaded{
            [&classDims_h, &classIds_h, &i](const redev::ClassPtn& ptn) {
              const auto ent =
                redev::ClassPtn::ModelEnt({classDims_h[i], classIds_h[i]});
              return ptn.GetRank(ent);
            },
            [](const redev::RCBPtn& ptn) {
              std::cerr << "RCB partition not handled yet\n";
              std::exit(EXIT_FAILURE);
              return 0;
            }},
          partition);
        reverse_partition[dr].emplace_back(local_index++);
      }
    }
    return reverse_partition;
  }
  Omega_h::Mesh& mesh;
};

redev::ClassPtn setupServerPartition(Omega_h::Mesh& mesh,
                                     std::string_view cpnFileName)
{
  namespace ts = test_support;
  auto ohComm = mesh.comm();
  const auto facePartition = !ohComm->rank()
                               ? ts::readClassPartitionFile(cpnFileName)
                               : ts::ClassificationPartition();
  ts::migrateMeshElms(mesh, facePartition);
  auto ptn = ts::CreateClassificationPartition(mesh);
  return redev::ClassPtn(MPI_COMM_WORLD, ptn.ranks, ptn.modelEnts);
}

void xgc_delta_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{


  wdmcpl::Coupler cpl("proxy_couple", wdmcpl::ProcessType::Client, comm,
                      redev::ClassPtn{});
  auto& delta_f = cpl.AddApplication("delta_f");
  auto is_overlap_h = markMeshOverlapRegion(mesh);
  auto& df_gid_field = delta_f.AddField<wdmcpl::GO>(
    "gids", OmegaHGids{mesh, is_overlap_h},
    OmegaHReversePartition{mesh},
    SerializeOmegaHGids{mesh, is_overlap_h},
    DeserializeOmegaH{mesh, is_overlap_h});

  df_gid_field.Send();
  df_gid_field.Receive();
}
void xgc_total_f(MPI_Comm comm, Omega_h::Mesh& mesh)
{
  wdmcpl::Coupler cpl("proxy_couple", wdmcpl::ProcessType::Client, comm,
                      redev::ClassPtn{});
  auto& total_f = cpl.AddApplication("total_f");
  auto is_overlap_h = markMeshOverlapRegion(mesh);
  auto& tf_gid_field = total_f.AddField<wdmcpl::GO>(
    "gids", OmegaHGids{mesh, is_overlap_h}, OmegaHReversePartition{mesh},
    SerializeOmegaHGids{mesh, is_overlap_h},
    DeserializeOmegaH{mesh, is_overlap_h});

  tf_gid_field.Send();
  // get updated field data from coupling server
  tf_gid_field.Receive();
}
struct DeserializeServer
{
  DeserializeServer(std::vector<wdmcpl::GO>& v) : v_(v){};
  template <typename T>
  void operator()(std::string_view name, nonstd::span<const T> buffer,
                  nonstd::span<const wdmcpl::LO> permutation) const
  {
    v_.resize(buffer.size());
    for (int i = 0; i < buffer.size(); ++i) {
      v_[i] = buffer[permutation[i]];
    }
  }

private:
  std::vector<wdmcpl::GO>& v_;
};
struct SerializeServer
{
  SerializeServer(std::vector<wdmcpl::GO>& v) : v_(v){};

  template <typename T>
  int operator()(std::string_view name, nonstd::span<T> buffer,
                 nonstd::span<const wdmcpl::LO> permutation) const
  {
    if (buffer.size() >= 0) {
      for (int i = 0; i < buffer.size(); ++i) {
        buffer[permutation[i]] = v_[i];
      }
    }
    return v_.size();
  }

private:
  std::vector<wdmcpl::GO>& v_;
};
void coupler(MPI_Comm comm, Omega_h::Mesh& mesh, std::string_view cpn_file)
{

  wdmcpl::Coupler cpl("proxy_couple", wdmcpl::ProcessType::Server, comm,
                      setupServerPartition(mesh, cpn_file));
  auto is_overlap_h = markMeshOverlapRegion(mesh);
  std::vector<wdmcpl::GO> delta_f_gids;
  std::vector<wdmcpl::GO> total_f_gids;
  auto& total_f = cpl.AddApplication("total_f");
  auto& delta_f = cpl.AddApplication("delta_f");
  auto& tf_gid_field = total_f.AddField<wdmcpl::GO>(
    "gids", OmegaHGids{mesh, is_overlap_h}, OmegaHReversePartition{mesh},
    SerializeServer{total_f_gids}, DeserializeServer{total_f_gids});

  auto& df_gid_field = delta_f.AddField<wdmcpl::GO>(
    "gids", OmegaHGids{mesh, is_overlap_h}, OmegaHReversePartition{mesh},
    SerializeServer{delta_f_gids}, DeserializeServer{delta_f_gids});

  df_gid_field.Receive();
  tf_gid_field.Receive();
  df_gid_field.Send();
  tf_gid_field.Send();
}

int main(int argc, char** argv)
{
  auto lib = Omega_h::Library(&argc, &argv);
  auto world = lib.world();
  const int rank = world->rank();
  if (argc != 4) {
    if (!rank) {
      std::cerr << "Usage: " << argv[0]
                << " <clientId=-1|0|1> /path/to/omega_h/mesh "
                   "/path/to/partitionFile.cpn\n";
    }
    exit(EXIT_FAILURE);
  }
  OMEGA_H_CHECK(argc == 4);
  const auto clientId = atoi(argv[1]);
  REDEV_ALWAYS_ASSERT(clientId >= -1 && clientId <= 1);
  const auto meshFile = argv[2];
  const auto classPartitionFile = argv[3];
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(meshFile, lib.world(), &mesh);
  MPI_Comm mpi_comm = lib.world()->get_impl();
  const std::string name = "meshVtxIds";
  switch (clientId) {
    case -1: coupler(mpi_comm, mesh, classPartitionFile); break;
    case 0: xgc_delta_f(mpi_comm, mesh); break;
    case 1: xgc_total_f(mpi_comm, mesh); break;
    default:
      std::cerr << "Unhandled client id (should be -1, 0,1)\n";
      exit(EXIT_FAILURE);
  }
  return 0;
}