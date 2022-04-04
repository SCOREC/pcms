#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H
#include <iostream>
#include <string>
#include <chrono> // steady_clock, duration
#include <numeric> // std::iota
#include <Omega_h_mesh.hpp>
#include <Omega_h_array_ops.hpp>
#include <redev.h>
#include <redev_comm.h>

namespace test_support {

struct OutMsg {
  redev::LOs dest;
  redev::LOs offset;
};

struct InMsg {
  redev::GOs srcRanks;
  redev::GOs offset;
  redev::GOs msgs;
  size_t start;
  size_t count;
};

void printTime(std::string mode, double min, double max, double avg);

void timeMinMaxAvg(double time, double& min, double& max, double& avg);

template <class T>
void getAndPrintTime(T start, std::string key, int rank) {
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  double min, max, avg;
  test_support::timeMinMaxAvg(elapsed_seconds.count(), min, max, avg);
  if(!rank) printTime(key, min, max, avg);
}

/**
 * Compute a partition of the mesh based on its classification
 * on the geometric model. This function is hardcoded for a specific mesh and
 * process count.
 */
void getClassPtn(Omega_h::Mesh& mesh, redev::LOs& ranks, redev::LOs& classIds);


void writeVtk(Omega_h::Mesh& mesh, std::string name, int step);

/**
 * Call the rendezvous Unpack function to receive data and populate the InMsg structure.
 */
void unpack(redev::AdiosComm<redev::GO>& comm, bool knownSizes, InMsg& in);

/**
 * Given the omegah mesh, and the partition object (ptn), determine which application vertex
 * should be sending to which rendezvous process, and populate the OutMsg
 * structures dest and offsets array that are required by the rendezvous Pack API.
 * The permutation from the application mesh to the array of message data sent
 * to the rendezvous processes is also computed here.
 */
void prepareAppOutMessage(Omega_h::Mesh& mesh, const redev::ClassPtn& ptn,
    OutMsg& out, redev::LOs& permute);

/**
 * Creates the permutation (rdvPermute) from the input message array of
 * vertex global ids (inGids) to the mesh on this rendezvous process (mesh)
 * such that gids[rdvPermute[i]] == inGids[i].
 */
void getRdvPermutation(Omega_h::Mesh& mesh, const redev::GOs& inGids, redev::GOs& rdvPermute);

/**
 * On the rendezvous processes use the permutation (rdvPermute) from the input message array
 * to the mesh to attach the incoming data (global vertex ids) and check that they match the
 * rendezvous vertex ids.
 */
void checkAndAttachIds(Omega_h::Mesh& mesh, std::string name, redev::GOs& vtxData, redev::GOs& rdvPermute);

/**
 * Return the index permutation of the input array (v) such that the array is
 * sorted in ascending order.
 * from https://stackoverflow.com/a/12399290
 */
template <typename T>
std::vector<size_t> sort_indexes(const T &v) {
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}

}
#endif
