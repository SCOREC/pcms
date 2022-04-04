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

void getClassPtn(Omega_h::Mesh& mesh, redev::LOs& ranks, redev::LOs& classIds);

void checkAndAttachIds(Omega_h::Mesh& mesh, std::string name, redev::GOs& vtxData, redev::GOs& rdvPermute);

void writeVtk(Omega_h::Mesh& mesh, std::string name, int step);

void unpack(redev::AdiosComm<redev::GO>& comm, bool knownSizes, InMsg& in);

void prepareAppOutMessage(Omega_h::Mesh& mesh, const redev::ClassPtn& ptn,
    OutMsg& out, redev::LOs& permute);

void getRdvPermutation(Omega_h::Mesh& mesh, const redev::GOs& inGids, redev::GOs& rdvPermute);

//from https://stackoverflow.com/a/12399290
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
