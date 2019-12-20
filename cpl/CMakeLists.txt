cmake_minimum_required(VERSION 2.8.12)
	#3.12)
project(wdm_coupling LANGUAGES Fortran C CXX)

# don't look for $ENV{XXX_ROOT} variables
# in find_library/package(xxx) commands
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 OLD)
endif()

cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
message(STATUS "Configuring to build on: ${HOST}")

# RHEL/CENTOS place x86-64 libs in lib64 instead of lib
set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS TRUE)

# prefer static linking over dynamic unless specified
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so")
if(BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

add_subdirectory(src)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()
add_executable(coupling ${CMAKE_SOURCE_DIR}/src/wdmapp.cc)
install(TARGETS coupling DESTINATION bin)
#cpack
#doxy
