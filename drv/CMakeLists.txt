cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
	#2.8)
project(wdm_minapp LANGUAGES Fortran C CXX)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 OLD)
endif()

cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
message(STATUS "Configuring to build on: ${HOST}")

set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS TRUE)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so")
if(BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
endif()

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

add_subdirectory(src)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()

#cpack
#doxy 
