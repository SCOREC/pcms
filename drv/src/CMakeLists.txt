cmake_minimum_required(VERSION 2.8)
	
set(CSRC )
set(FSRC new_coupling_xgc.F90
	 sim_param.F90)
set(HDRS )
set(MODS new_coupling_xgc.mod
	 sim_param.mod
	 adios_comm_mod.mod)

 set(LIB_NAME wdm_minapp)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../../cpl/cmake/")


#find_package(XGC REQUIRED)
#find_package(GENE REQIRED)

find_package(ADIOS REQUIRED)
include_directories(${ADIOS_INCLUDE_DIRS})

find_package(MPI)
include_directories(${MPI_INCLUDE_DIRS})

find_package(HDF5)
include_directories(${HDF5_INCLUDE_DIRS})

set(LINK_LIBS ${LINK_LIBS} ${ADIOS_LIB} ${MPI_Fortran_LIBRARIES} ${HDF5_LIBRARIES})
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${ADIOS_INCLUDE_DIRS} ${MPI_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

add_library(${LIB_NAME} ${FSRC} ${CSRC})
target_include_directories(${LIB_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  ${INCLUDE_DIRS}
  )
target_link_libraries(${LIB_NAME} ${LINK_LIBS})

install(FILES ${MODS} DESTINATION ${PROJECT_BINARY_DIR})
install(TARGETS ${LIB_NAME} LIBRARY DESTINATION ${PROJECT_BINARY_DIR} ARCHIVE DESTINATION ${PROJECT_BINARY_DIR})

