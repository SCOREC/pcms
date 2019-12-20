cmake_minimum_required(VERSION 2.8)
# Module load:
#	adios2/2.4.0
#	cray-hdf5-parallel/1.10.2.0
set(CSRC )
set(FSRC adios2_comm_mod.F90
	 adios_comm_mod.F90	
	 new_coupling.F90
	 sim_param.F90)
set(HDRS )
set(MODS new_coupling.mod
	 sim_param.mod
	 adios_comm_mod.mod
	 adios2_comm_mod.mod)

set(LIB_NAME wdm_minapp)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../../cpl/cmake/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/global/common/sw/cray/cnl7/haswell/adios2/2.4.0/intel/19.0.3.199/26fprnf/lib64/cmake/adios2")

#find_package(XGC REQUIRED)
#find_package(GENE REQIRED)

find_package(ADIOS REQUIRED)
include_directories(${ADIOS_INCLUDE_DIRS})

set(USER_LIB_DIR "/project/projectdirs/m499/Software")
set(COUPLING_PLATFORM "cori_haswell")
set(ADIOS2_DIR "${USER_LIB_DIR}/adios2/DEFAULT/${COUPLING_PLATFORM}/DEFAULT")
set(ADIOS2_INCLUDE_DIRS "/global/project/projectdirs/m499/Software/adios2/v2.1.0-3511-g0fa6aeef/cori_haswell/intel/include/adios2/fortran")
#set(ADIOS2_LIB execute_process("${ADIOS2_DIR}/bin/adios2-config --fortran-libs --cxx-libs"))

include_directories(${ADIOS2_INCLUDE_DIRS})

find_package(MPI)
include_directories(${MPI_INCLUDE_DIRS})

find_package(HDF5)
include_directories(${HDF5_INCLUDE_DIRS})

#add_definitions(-DADIOS)
add_definitions(-DADIOS2)

set(LINK_LIBS ${LINK_LIBS} ${ADIOS_LIB} ${ADIOS2_LIB} ${MPI_Fortran_LIBRARIES} ${HDF5_LIBRARIES})
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${ADIOS_INCLUDE_DIRS} ${ADIOS2_INCLUDE_DIRS} ${MPI_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

add_library(${LIB_NAME} ${FSRC} ${CSRC})
target_include_directories(${LIB_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  ${INCLUDE_DIRS}
  )
target_link_libraries(${LIB_NAME} ${LINK_LIBS})

#install(FILES ${MODS} DESTINATION ${PROJECT_BINARY_DIR})
install(TARGETS ${LIB_NAME} LIBRARY DESTINATION ${PROJECT_BINARY_DIR} ARCHIVE DESTINATION ${PROJECT_BINARY_DIR})

