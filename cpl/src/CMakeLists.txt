#cmake_minimum_required(VERSION 2.8)
#project(app CXX)

set(SRCS  
	#coupling.cc
	cpl_scorec.cc
	)
# any headers which need to be public/installed
#  this does not need to be ALL headers
set(HDRS )
set(EXEC cpl_exec)

set(USER_LIB "/users/adesoa/dev")
set(ADIOS2_INCLUDE_DIRS "${USER_LIB}/install/ADIOS2/include")
set(ADIOS2_INC "${USER_LIB}/install/ADIOS2/lib64")

set(CPL_INC "/lore/adesoa/dev/wdmapp_coupling/drv/build")
set(DRV_INC "${CPL_INC}/src")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${ADIOS2_INC} ${CPL_INC} )


find_package( MPI REQUIRED)
include_directories( ${MPI_CXX_INCLUDE_PATH} )

find_library(DRV_LIB NAMES adios2 libadios2)
include_directories(${DRV_INC})
if(NOT DRV_LIB)
  message(STATUS "coupling driver library not found")
else()
  message(STATUS "coupling driver library was found")
endif()

find_library(ADIOS2_LIB NAMES adios2 libadios2)
include_directories(${ADIOS2_INCLUDE_DIRS})
if(NOT ADIOS2_LIB)
  message(STATUS "adios2 library not found")
else()
  message(STATUS "adios2 library was found")
endif()

find_package(ADIOS REQUIRED)
include_directories(${ADIOS_INCLUDE_DIRS})

set(LINK_LIBS ${LINK_LIBS} ${MPI_LIBRARIES} ${ADIOS_LIBRARIES} ${ADIOS2_LIB} ${DRV_LIB})
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${MPI_INCLUDE_PATH} ${ADIOS_INCLUDE_DIRS} ${ADIOS2_INCLUDE_DIRS} ${DRV_INC})

add_executable(${EXEC} ${SRCS} )
target_link_libraries(${EXEC} ${LINK_LIBS})

INSTALL(TARGETS ${EXEC}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION ${PROJECT_BINARY_DIR}
        ARCHIVE DESTINATION ${PROJECT_BINARY_DIR})
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
