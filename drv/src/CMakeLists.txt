set(CSRC )
set(FSRC new_coupling_xgc.F90
	 sim_param.F90)
set(HDRS )

set(LIB_NAME wdm_minapp)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../../cpl/cmake/")
#list(APPEND CMAKE_PROGRAM_PATH "${CMAKE_SOURCE_DIR}/../../WDM_XGC")

#find_package(GENE REQIRED)

find_package(ADIOS REQUIRED)
include_directories(${ADIOS_INCLUDE_DIRS})

find_package(MPI)
include_directories(${MPI_INCLUDE_DIRS})

find_package(HDF5)
include_directories(${HDF5_INCLUDE_DIRS})

set(LINK_LIBS ${LINK_LIBS} ${ADIOS_LIB} ${MPI_Fortran_LIBRARIES} ${HDF5_LIBRARIES})
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${ADIOS_INCLUDE_DIRS} ${MPI_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

# installing XGC
set(XGC_PREFIX "${CMAKE_SOURCE_DIR}/../../WDM_XGC")
find_program(XGC xgc-es-coupling REQUIRED PATHS ${XGC_PREFIX} NO_DEFAULT_PATH)

if(XGC)
	message(STATUS "XGC found")
else()
	message(FATAL_ERROR "XGC not found")
endif()

add_library(${LIB_NAME} ${FSRC} ${CSRC})
target_include_directories(${LIB_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  ${INCLUDE_DIRS}
  )
target_link_libraries(${LIB_NAME} ${LINK_LIBS})

#set(MODS new_coupling_xgc.mod
	#	 sim_param.mod
	 #	 adios_comm_mod.mod)

#install(FILES ${MODS} DESTINATION ${PROJECT_BINARY_DIR})
 install(DIRECTORY "${CMAKE_SOURCE_DIR}/"
	 DESTINATION "include"
	 FILES_MATCHING
	 PATTERN ".mod")

 
install(TARGETS ${LIB_NAME} LIBRARY DESTINATION ${PROJECT_BINARY_DIR} ARCHIVE DESTINATION ${PROJECT_BINARY_DIR})

