cmake_minimum_required(VERSION 2.8)
	
set(CSRC )
set(FSRC new_coupling_xgc.F90
	 sim_param.F90)
set(HDRS )
set(MODS new_coupling_xgc.mod
	 sim_param.mod
	 adios_read_mod.mod)

 set(LIB_NAME wdm_minapp)
#find_package(XGC REQUIRED)
#find_package(GENE REQIRED)

set(ADIOS_ROOT_DIR /users/adesoa/dev/wdm/wdmapp_coupling/cpl/cmake/)
set(ADIOS_INCLUDE_DIRS ${ADIOS_ROOT_DIR}/adios/)
list(APPEND CMAKE_MODULE_PATH "${ADIOS_ROOT_DIR}")
find_package(ADIOS REQUIRED)
include_directories(${ADIOS_INCLUDE_DIRS})

#--------------------MPI-------------------------
set(MPI_INCLUDE_DIRS /opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/mpich-3.3-diz4f6ieln25ouifyc7ndtqlfksom6nb/include/)
find_package(MPI)
include_directories(${MPI_INCLUDE_DIRS})

set(LINK_LIBS ${LINK_LIBS} ${ADIOS_LIB} ${MPI_Fortran_LIBRARIES})
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${ADIOS_INCLUDE_DIRS} ${MPI_INCLUDE_DIRS})

add_library(${LIB_NAME} ${FSRC} ${CSRC})
target_include_directories(${LIB_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  ${INCLUDE_DIRS}
  )
target_link_libraries(${LIB_NAME} ${LINK_LIBS})

install(FILES ${MODS} DESTINATION ${PROJECT_BINARY_DIR})
install(TARGETS ${LIB_NAME} LIBRARY DESTINATION ${PROJECT_BINARY_DIR} ARCHIVE DESTINATION ${PROJECT_BINARY_DIR})


