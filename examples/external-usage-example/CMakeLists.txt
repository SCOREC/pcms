# ########### Configure and Link pcms Library ##################################
# This is an example CMakeLists.txt file to configure and link the pcms library
# with an external project.
# ##############################################################################

cmake_minimum_required(VERSION 3.23)
project(pcms_installation_test)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

message(STATUS "MODULE PATH: ${CMAKE_MODULE_PATH}")

# It searches in pcms_ROOT, pcms_DIR, or other standard paths
find_package(pcms REQUIRED)
# The example uses Omega_h (and the Omega_h backend header) and MPI directly;
# PCMS links them privately, so link them here.
find_package(Omega_h REQUIRED)
find_package(MPI REQUIRED)

add_executable(pcms_test main.cpp)

# Link the pcms library to the executable Other target names cound be found in
# your installation's <pcms install dir>/lib/cmake/pcms/pcms-targets.cmake
target_link_libraries(pcms_test PUBLIC pcms::pcms Omega_h::omega_h MPI::MPI_CXX)
