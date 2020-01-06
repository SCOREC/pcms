#This is the top XGC CMakeList File for the Build
#Setting Version Number, Project Name

cmake_minimum_required (VERSION 3.3)
project (hello)
#enable_language(CXX)
# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

#unless building shared libs, then select static libs 
# if both static and shared libs are available 
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so") 
if(BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
endif()


#set variables for finding dependent libraries
set(USER_LIB "/users/adesoa/dev")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../../XGC_adios2/cmake/cmake-modules")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/cmake-3.3.1-qqrkgm4rrdmq2cy56k6ngajiu42z2nci/share/cmake-3.3/Modules")
set(ADIOS2_INCLUDE_DIRS "${USER_LIB}/install/ADIOS2/include")
set(PCI_DIR "/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/libpciaccess-0.13.5-gcqhf4ksz5ndlmj5aoydcpexg3sizi4s/lib")
set(XNV_DIR "/usr/lib64")
set(PETSC_DIR "${USER_LIB}/petsc/3.7.7")
set(ADIOS2_INC "${USER_LIB}/install/ADIOS2/lib64")
set(COUPLING_DIR "${USER_LIB}/wdm/wdmapp_coupling/drv/build/src")
set(PERFMOD_INC "${USER_LIB}/XGCM/xgc_scorec/scripts/XGC1/camtimers")
set(PSPLINE_INC "${USER_LIB}/XGCM/xgc_scorec/scripts/XGC1/pspline")
set(NETCDF_DIR "/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/netcdf-4.6.3-n6uagsdaydnxskwuhgjsvnmp6bddod7s")
set(NETCDF_LIB_DIR "${NETCDF_DIR}/lib")
set(NETCDF_INC "${NETCDF_DIR}/include")

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${PERFMOD_INC} ${PSPLINE_INC} ${ADIOS2_INC} ${NETCDF_LIB_DIR} ${COUPLING_DIR} ${PCI_DIR} ${XNV_DIR})

#find dependent packages and libraries
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
if(MPI_CXX_FOUND)
	message(STATUS "${MPI_CXX_INCLUDE_PATH}")
	message(STATUS "${MPI_CXX_LIBRARIES}")
	message(STATUS "${MPI_LINK_FLAGS}")
else()
  message(STATUS "empty variables")
endif()

if(FALSE)
# find libraries that I dont have modules for
find_library(PCI_LIB NAMES pciaccess)
if(NOT PCI_LIB)
  message(STATUS "PCI library not found")
else()
  message(STATUS "PCI library was found")
endif()

find_library(XNV_LIB NAMES XNVCtrl)
if(NOT XNV_LIB)
  message(STATUS "XNV library not found")
else()
  message(STATUS "XNV library was found")
endif()

find_package(X11 REQUIRED)
endif()



find_package(ADIOS)

find_package(PETSc)

find_library(ADIOS2_LIB NAMES adios2 libadios2)
if(NOT ADIOS2_LIB)
  message(STATUS "adios2 library not found")
else()
  message(STATUS "adios2 library was found")
endif()



#comment block starts
if(FALSE)
find_library(COUPLING_LIB NAMES wdm_minapp libwdm_minapp)
if(NOT COUPLING_LIB)
  message(STATUS "coupling library not found")
else()
  message(STATUS "coupling library was found")
endif()

find_library(PERFMOD_LIB NAMES timers libtimers)
if(NOT PERFMOD_LIB)
  message(STATUS "camtimers library not found")
else()
  message(STATUS "camtimers library was found")
endif()

find_library(PSPLINE_LIB NAMES pspline libpspline)
if(NOT PSPLINE_LIB)
  message(STATUS "pspline library not found")
else()
  message(STATUS "pspline library was found")
endif()

find_library(NETCDF_LIB NAMES netcdf libnetcdf)
if(NOT NETCDF_LIB)
  message(STATUS "netcdf library not found")
else()
  message(STATUS "netcdf library was found")
endif()

endif()
#comment block ends


#add_definitions(${MPI_CXX_LINK_FLAGS})
#add_definitions(${-DADIOS2})

#set(DEP_LIBS ${DEP_LIBS} ${ADIOS_LIBRARIES} ${ADIOS2_LIB} ${PSPLINE_LIB} ${PETSC_LIBRARIES} ${PERFMOD_LIB} ${NETCDF_LIB} ${COUPLING_LIB} ${PCI_LIB} ${XNV_LIB} ${MPI_CXX_LIBRARIES})
set(DEP_LIBS ${DEP_LIBS} ${MPI_CXX_LIBRARIES} ${ADIOS_LIBRARIES} ${ADIOS2_LIB})

#set(INCLUDE_DIRS ${INCLUDE_DIRS} ${ADIOS_INCLUDE_DIRS} ${ADIOS2_INCLUDE_DIRS} ${PSPLINE_INC} ${PERMOD_INC} ${PETSC_INCLUDES} ${NETCDF_INC} ${COUPLING_DIR} ${MPI_CXX_INCLUDE_PATH})
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH} ${ADIOS_INCLUDE_DIRS} ${ADIOS2_INCLUDE_DIRS})

add_executable(app1 app1.cc)

target_link_libraries(app1 ${DEP_LIBS})

target_include_directories(app1 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  ${INCLUDE_DIRS}
  )

INSTALL(TARGETS app1
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION ${PROJECT_BINARY_DIR}
	ARCHIVE DESTINATION ${PROJECT_BINARY_DIR})
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
