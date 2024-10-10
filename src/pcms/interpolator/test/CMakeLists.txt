project(mfem-cmake-project C CXX)
cmake_minimum_required(VERSION 3.20)

find_package(Omega_h 10 REQUIRED)

find_package(pcms REQUIRED)
add_executable(test test_polybasis.cpp)
target_link_libraries(test PUBLIC Omega_h::omega_h pcms::pcms)
