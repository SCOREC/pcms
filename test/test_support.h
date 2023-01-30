#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H
#include <iostream>
#include <string_view>
#include <chrono>  // steady_clock, duration
#include <numeric> // std::iota
#include <Omega_h_mesh.hpp>
#include <Omega_h_array_ops.hpp>
#include <redev.h>
#include <redev_comm.h>
#include <wdmcpl/external/span.h>
#include <wdmcpl/memory_spaces.h>
#include <wdmcpl/omega_h_field.h>
#include <functional>

namespace test_support
{

struct CSR
{
  redev::GOs off;
  redev::GOs val;
};

struct OutMsg
{
  redev::LOs dest;
  redev::LOs offset;
  redev::LOs permute;
};

struct InMsg
{
  redev::GOs srcRanks;
  redev::GOs offset;
  redev::GOs msgs;
  size_t start;
  size_t count;
};

struct ClassificationPartition
{
  redev::LOs ranks;
  redev::ClassPtn::ModelEntVec modelEnts;
};

void printTime(std::string_view mode, double min, double max, double avg);

void timeMinMaxAvg(double time, double& min, double& max, double& avg);

template <class T>
void getAndPrintTime(T start, std::string_view key, int rank)
{
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  double min, max, avg;
  test_support::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if (!rank)
    printTime(key, min, max, avg);
}

ClassificationPartition readClassPartitionFile(std::string_view cpnFileName);

ClassificationPartition CreateClassificationPartition(Omega_h::Mesh& mesh);

void migrateMeshElms(Omega_h::Mesh& mesh,
                     const ClassificationPartition& partition);

/**
 * Migrate 18 of the mesh elements to rank 1 and return its classification
 * partition on the geometric model. This function is hardcoded for a specific
 * mesh and process count.
 */
ClassificationPartition migrateAndGetPartition(Omega_h::Mesh& mesh);

void writeVtk(Omega_h::Mesh& mesh, std::string_view name, int step);

/**
 * Given the omegah mesh, and the partition object (partition), determine which
 * application vertex should be sending to which rendezvous process, and
 * populate the OutMsg structures dest and offsets array that are required by
 * the Redev::SetOutMessageLayout API. The permutation from the client mesh to
 * the array of message data sent to the server processes is also computed here.
 */
OutMsg prepareAppOutMessage(Omega_h::Mesh& mesh,
                            const redev::ClassPtn& partition);

/**
 * Creates the permutation (rdvPermute) from the input message array of
 * vertex global ids (inGids) to the mesh on this rendezvous process (mesh)
 * such that gids[rdvPermute[i]] == inGids[i].
 */
redev::GOs getRdvPermutation(Omega_h::Mesh& mesh, const redev::GOs& inGids);

/**
 * Creates the rendezvous -> non-rendezvous permutation CSR given inGids and the
 * rdv mesh instance.
 */
CSR getRdvOutPermutation(Omega_h::Mesh& mesh, const redev::GOs& inGids);

/**
 * Construct the meta data for the rendezvous -> non-rendezvous (reverse) send
 * from the meta data associated with the non-rendezvous -> rendezvous (forward)
 * send.
 */
OutMsg prepareRdvOutMessage(Omega_h::Mesh& mesh,
                            const redev::InMessageLayout& in);

/**
 * On the rendezvous processes use the permutation (rdvPermute) from the input
 * message array to the mesh to attach the incoming data (global vertex ids) and
 * check that they match the rendezvous vertex ids.
 */
void checkAndAttachIds(Omega_h::Mesh& mesh, std::string_view name,
                       const redev::GOs& vtxData, const redev::GOs& rdvPermute);

/**
 * Return the index permutation of the input array (v) such that the array is
 * sorted in ascending order.
 * from https://stackoverflow.com/a/12399290
 */
template <typename T>
std::vector<size_t> sortIndexes(const T& v)
{
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  std::stable_sort(idx.begin(), idx.end(),
                   [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
  return idx;
}
/**
 * return 1 if the specificed model entity is part of the overlap region, 0
 * otherwise. Device function must be defined inline
 */
OMEGA_H_DEVICE Omega_h::I8 isModelEntInOverlap(const int dim, const int id)
{
  // the TOMMS generated geometric model has
  // entity IDs that increase with the distance
  // from the magnetic axis
  if ((id >= 22 && id <= 34) && (dim >= 0 && dim <= 2)) {
    return 1;
  }
  return 0;
}
using EntInOverlapFunc = std::function<Omega_h::I8(const int, const int)>;
static_assert(std::is_constructible_v<EntInOverlapFunc, decltype(isModelEntInOverlap)>);
/**
 * On the server we mark the vertices on each process that are in the overlap
 * region and are owned by the process as defined by the Classification
 * Partition GetRank(modelEntity) function.  The ownership of mesh entities on
 * the mesh partition boundary is not guaranteed to match the ownership of the
 * geometric model entities defined by the Classification Partition so we must
 * use GetRank(modelEntity) to ensure consistency with data incoming from
 * the client processes.
 *
 * On the client side we only care that each mesh vertex is sent by exactly one
 * client process (to the server process returned by GetRank(modelEntity)) so
 * using the ownership of mesh entities following the mesh partition ownership
 * is OK. The function markMeshOverlapRegion(...) supports this.
 */
Omega_h::Read<Omega_h::I8> markServerOverlapRegion(
  Omega_h::Mesh& mesh, const redev::ClassPtn& classPtn,
  const EntInOverlapFunc& entInOverlap);
/**
 * Create the tag 'isOverlap' for each mesh vertex whose value is 1 if the
 * vertex is classified on a model entity in the closure of the geometric model
 * faces forming the overlap region; the value is 0 otherwise.
 * OnlyIncludesOverlapping and owned verts
 */
Omega_h::Read<Omega_h::I8> markOverlapMeshEntities(
  Omega_h::Mesh& mesh, const EntInOverlapFunc& entInOverlap);

redev::ClassPtn setupServerPartition(Omega_h::Mesh& mesh,
                                     std::string_view cpnFileName);
template <typename T1, typename T2>
struct Sum
{
  Sum(T1 arr1, T2 arr2) : arr1_(arr1), arr2_(arr2) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const noexcept { arr1_[i] += arr2_[i]; }
  T1 arr1_;
  T2 arr2_;
};
template <typename T1, typename ScalarT>
struct Divide
{
  Divide(T1 arr1, ScalarT scalar) : arr1_(arr1), scalar_(scalar) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const noexcept { arr1_[i] /= scalar_; }
  T1 arr1_;
  ScalarT scalar_;
};

// TODO: move to be internal to wdmcpl
struct MeanCombiner
{
  void operator()(
    const nonstd::span<const std::reference_wrapper<wdmcpl::InternalField>>&
      fields,
    wdmcpl::InternalField& combined_variant) const
  {

    // Internal fields are OmegaHFields
    using execution_space =
      typename wdmcpl::OmegaHMemorySpace::type::execution_space;
    std::visit(
      [&fields](auto&& combined_field) {
        using T = typename std::remove_reference_t<
          decltype(combined_field)>::value_type;
        Omega_h::Write<T> combined_array(combined_field.Size());
        for (auto& field_variant : fields) {
          std::visit(
            [&combined_array, &combined_field](auto&& field) {
              WDMCPL_ALWAYS_ASSERT(field.Size() == combined_array.size());
              auto field_array = get_nodal_data(field);
              auto policy =
                Kokkos::RangePolicy<execution_space>(0, field_array.size());
              Kokkos::parallel_for(policy, Sum{combined_array, field_array});
            },
            field_variant.get());
        }
        auto num_fields = fields.size();
        auto policy =
          Kokkos::RangePolicy<execution_space>(0, combined_array.size());
        Kokkos::parallel_for(policy, Divide{combined_array, num_fields});
        set_nodal_data(combined_field,
                       make_array_view(Omega_h::Read(combined_array)));
      },
      combined_variant);
  }
};

} // namespace test_support
#endif
