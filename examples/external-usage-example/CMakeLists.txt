# ########### Configure and Link pcms Library ##################################
# This is an example CMakeLists.txt file to configure and link the pcms library
# with an external project.
# ##############################################################################

cmake_minimum_required(VERSION 3.19)
project(pcms_installation_test)

# It searches in pcms_ROOT, pcms_DIR, or other standard paths
find_package(pcms REQUIRED)

add_executable(pcms_test main.cpp)

# Link the pcms library to the executable Other target names cound be found in
# your installation's <pcms install dir>/lib64/cmake/pcms/pcms-targets.cmake
target_link_libraries(pcms_test PUBLIC pcms::pcms)
