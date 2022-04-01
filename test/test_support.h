#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H
#include <iostream>
#include <string>
#include <chrono> //steady_clock, duration
#include <Omega_h_mesh.hpp>
#include <Omega_h_array_ops.hpp>
#include <redev.h>

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

}
#endif
