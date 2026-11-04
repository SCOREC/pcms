# ########### Use PCMS core without its private dependencies ##################
# A program that uses only the core API needs Kokkos and PCMS: pcms::pcms
# brings neither Omega_h, MPI, MeshFields nor PETSc onto this target's
# include or link line (for a shared PCMS built without the coupler and the
# C/Fortran APIs).
# ##############################################################################

cmake_minimum_required(VERSION 3.23)
project(pcms_core_usage_test CXX)

find_package(pcms REQUIRED)

add_executable(pcms_core_test main.cpp)
target_link_libraries(pcms_core_test PRIVATE pcms::pcms)
