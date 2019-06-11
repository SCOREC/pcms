cmake_minimum_required(VERSION 2.8)
set(CSRC )
set(FSRC 
	new_coupling_xgc.F90)

set(HDRS )
set(MODS 
	new_coupling_xgc.mod)

set(LIB_NAME wdm_minapp)

#includes ADIOS dependency
#find_package(wdm_coupler REQUIRED)
#find_package(XGC REQUIRED)
#find_package(GENE REQIRED)

set(ADIOS_ROOT_DIR /users/adesoa/dev/wdm/wdmapp_coupling/cpl/cmake/)
list(APPEND CMAKE_MODULE_PATH "${ADIOS_ROOT_DIR}")

find_package(ADIOS REQUIRED)
include_directories(${ADIOS_INCLUDE_DIRS})
#PATHS /users/adesoa/dev/wdm/wdmapp_coupling/cpl/cmake) # ADIOS_LIBS ADIOS_INCLUDE_DIR

set(LINK_LIBS ${LINK_LIBS} ${ADIOS_LIB})
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${ADIOS_INCLUDE_DIRS} )
set(WDM_ROOT /users/adesoa/dev/install/wdm_coupling/)

add_library(${LIB_NAME} ${FSRC} ${CSRC})
target_include_directories(${LIB_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  ${INCLUDE_DIRS}
  )

target_link_libraries(${LIB_NAME} ${LINK_LIBS})


install(FILES ${MODS} DESTINATION ${WDM_ROOT}/drv/build)
#install(FILES "${HDRS} ${MODS}" DESTINATION include)
#install(TARGETS ${LIB_NAME} LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(TARGETS ${LIB_NAME} LIBRARY DESTINATION ${WDM_ROOT}/drv/build ARCHIVE DESTINATION ${WDM_ROOT}/drv/build)


