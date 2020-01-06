cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
	#2.8)
project(examples LANGUAGES CXX)

if(POLICY CMP0074)
  cmake_policy(SET CMP0074 OLD)
endif()

cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
message(STATUS "Configuring to build on: ${HOST}")

set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS TRUE)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so")
#if(BUILD_SHARED_LIBS)
#  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
#endif()

add_subdirectory(app1)
#add_subdirectory(app2)
#add_subdirectory(app3)

#if(ENABLE_TESTING)
# add_subdirectory(test)
#endif()

#cpack
#doxy 
