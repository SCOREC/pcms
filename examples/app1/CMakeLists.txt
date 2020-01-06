#cmake_minimum_required(VERSION 2.8)
#project(app CXX)

set(USER_LIB "/users/adesoa/dev")
set(ADIOS2_INCLUDE_DIRS "${USER_LIB}/install/ADIOS2/include")
set(ADIOS2_INC "${USER_LIB}/install/ADIOS2/lib64")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${ADIOS2_INC} )


find_package( MPI REQUIRED)
include_directories( ${MPI_CXX_INCLUDE_PATH} )

find_library(ADIOS2_LIB NAMES adios2 libadios2)
include_directories(${ADIOS2_INCLUDE_DIRS})
if(NOT ADIOS2_LIB)
  message(STATUS "adios2 library not found")
else()
  message(STATUS "adios2 library was found")
endif()



add_executable( app1 app1.cc )
target_link_libraries(  app1 ${MPI_CXX_LIBRARIES} ${ADIOS2_LIB})

INSTALL(TARGETS app1
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION ${PROJECT_BINARY_DIR}
        ARCHIVE DESTINATION ${PROJECT_BINARY_DIR})
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})


