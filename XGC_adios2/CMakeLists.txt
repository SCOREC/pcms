#This is the top XGC CMakeList File for the Build

#Setting Version Number, Project Name
cmake_minimum_required (VERSION 3.3)
project (xgc)
enable_language(Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message (Fortran_COMPILER_NAME = ${Fortran_COMPILER_NAME})

#unless building shared libs, then select static libs 
# if both static and shared libs are available 
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so") 
if(BUILD_SHARED_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
endif()

#Settings options for testing
enable_testing()
include(CTest)
#This will be set to ON by the CTest driver script (and only by that)
option(ENABLE_TESTING "Build for CTest" OFF)
set(MPIRUN "mpirun"
    CACHE STRING 
    "the mpirun or srun executable")
set(MPIRUN_PROCFLAG "-np"
    CACHE STRING 
    "the command line flag to give process count to MPIRUN")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set(TARGET_LIB_NAME xgc_es)

#set variables for finding dependent libraries
set(USER_LIB "/users/adesoa/dev")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/cmake-modules")
set(ADIOS2_INCLUDE_DIRS "${USER_LIB}/install/ADIOS2/include/adios2/fortran")
set(PETSC_DIR "${USER_LIB}/petsc/3.7.7")
set(ADIOS2_INC "${USER_LIB}/install/ADIOS2/lib64")
set(COUPLING_DIR "${USER_LIB}/wdm/wdmapp_coupling/drv/build/src")
set(PERFMOD_INC "${USER_LIB}/XGCM/xgc_scorec/scripts/XGC1/camtimers")
set(PSPLINE_INC "${USER_LIB}/XGCM/xgc_scorec/scripts/XGC1/pspline")
set(NETCDF_DIR "/opt/scorec/spack/install/linux-rhel7-x86_64/gcc-7.3.0/netcdf-4.6.3-n6uagsdaydnxskwuhgjsvnmp6bddod7s")
set(NETCDF_LIB_DIR "${NETCDF_DIR}/lib")
set(NETCDF_INC "${NETCDF_DIR}/include")

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${PERFMOD_INC} ${PSPLINE_INC} ${ADIOS2_INC} ${NETCDF_LIB_DIR} ${COUPLING_DIR})

#find dependent packages and libraries
find_package(ADIOS)

find_package(PETSc)
#find_package(PETSc COMPONENTS Fortran)

find_library(ADIOS2_LIB NAMES adios2 libadios2)
if(NOT ADIOS2_LIB)
  message(STATUS "adios2 library not found")
else()
  message(STATUS "adios2 library was found")
endif()

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


set(IMSL_OBJ my_imsl.F90)
set(PORT_OBJ bspline90_22.F90 taus88.F90 derf.F90 datanh.F90 pppack.F90 fmin.F90)
set(SER_OBJ mpisingle.F90)
set(PAR_OBJ mpi.F90)
set(EXTRA_OBJ ${PORT_OBJ})
set(MPI_OBJ ${PAR_OBJ})


set(SRC module.F90 search.F90 module_psn.F90 pol_decomp.F90 f0module.F90 adios_comm_mod.F90 adios2_comm_mod.F90 
        elliptics.F90 collisionf.F90 lbal_mod.F90 ${EXTRA_OBJ} psmooth.F90 bicub_mod.F90 one_d_cub_mod.F90 initial_perturbation_GM.F90 
        interfaces.F90 charge.F90 diagnosis.F90 poisson_extra_xgc1.F90 qevaluateandtrapped.F90 setup.F90 read.F90 
	gen_perm.F90 sort_particles.F90 push.F90 pushe.F90 load.F90 efield.F90 interpolation.F90 ${MPI_OBJ} limiter.F90 bounce.F90 
	diagnosis2.F90 collision.F90 collision2.F90 collisionf2.F90 heat.F90 turbulence.F90 neutral.F90 neutral2.F90 
	neutral3.F90 linearsolver.F90 fem2d.F90 fem_ops.F90 poisson_extra_common.F90 petsc_solve.F90) 
#	coupling_core_edge.F90)

set(SRC_ES es_main.F90 es_poisson.F90)	

#link PETSC_KSP_LIB and all other dependencies
set(SOURCES ${SRC} ${SRC_ES})
set( CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Mod )


set(DEP_LIBS ${DEP_LIBS} ${ADIOS_LIBRARIES} ${ADIOS2_LIB} ${PSPLINE_LIB} ${PETSC_LIBRARIES} ${PERFMOD_LIB} ${NETCDF_LIB} ${COUPLING_LIB})
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${ADIOS_INCLUDE_DIRS} ${ADIOS2_INCLUDE_DIRS} ${PSPLINE_INC} ${PERMOD_INC} ${PETSC_INCLUDES} ${NETCDF_INC} ${COUPLING_DIR})
set(FC_FLAG "-ffree-line-length-none -g -O0 ")
add_definitions( ${FC_FLAG}
		  -DPSPLINE 
                  -DUSE_BICUB_MOD 
                  -DUSE_ONE_D_I_CUB_MOD 
                  -UCAM_TIMERS 
                  -DADIOS 
                  -DADIOS2 
                  -DSOLVERLU -llapack -lblas 
                  -DITER_GRID 
                  -DF_USE_MARKER_DEN 
                  -DF_USE_MARKER_DEN2 
                  -DF0_TOR_LINEAR 
                  -DV_PERP 
                  -DDELTAF_MODE2 
                  -DINIT_GENE_PERT 
                  -UPOS2_ORG_HACK 
                  -UFORCE_ZEROCHARGE 
                  -UFORCE_ZEROCHARGE_NEW 
                  -UOPTIM_GYRO_AVG_MAT 
                  -UNEWGYROMATRIX 
                  -DPURE_RK4 
                  -DRK4_ACCURATE_EFIELD 
                  -DSPECIFIC_GENEXGC 
                  -DXGC_COUPLING_CORE_EDGE 
                  -DXGC_COUPLING_CORE_EDGE_FIELD 
                  -DXGC_COUPLING_CORE_EDGE_VARPI2 
                  -UCCE_TEST_1
)

add_library(${TARGET_LIB_NAME} ${SOURCES})
target_link_libraries(${TARGET_LIB_NAME} ${DEP_LIBS})
target_include_directories(${TARGET_LIB_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  ${INCLUDE_DIRS}
  )

INSTALL(TARGETS ${TARGET_LIB_NAME}
	LIBRARY DESTINATION XGC
	ARCHIVE DESTINATION XGC)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})


#binary distribution package
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION "1.0.1")
include(CPack)
