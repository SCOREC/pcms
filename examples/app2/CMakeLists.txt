#cmake_minimum_required(VERSION 2.8)

set(SRCS  
	app2.cc
	)
# any headers which need to be public/installed
#  this does not need to be ALL headers
set(HDRS )
set(EXEC app2)

set(USER_LIB "/users/adesoa/dev")
set(ADIOS2_USER_DIR "/project/projectdirs/m499/Software/adios2/DEFAULT/cori_haswell/DEFAULT")

#For cori build
set(ADIOS2_INCLUDE_DIRS "${ADIOS2_USER_DIR}/include")
set(ADIOS2_LIB_DIR "${ADIOS2_USER_DIR}/lib64")


#For scorec build
#set(ADIOS2_INCLUDE_DIRS "${USER_LIB}/install/ADIOS2/include")
#set(ADIOS2_LIB_DIR "${USER_LIB}/install/ADIOS2/lib64")

set(CPL_INC "/lore/adesoa/dev/wdmapp_coupling/drv/build")
set(DRV_INC "${CPL_INC}/src")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${ADIOS2_LIB_DIR} ${CPL_INC} )


find_package(MPI REQUIRED)
include_directories( ${MPI_CXX_INCLUDE_PATH} )

find_library(ADIOS2_F NAMES adios2_f)
include_directories(${ADIOS2_INCLUDE_DIRS})
if(NOT ADIOS2_F)
  message(STATUS "adios2_f library not found")
else()
  message(STATUS "adios2_f library was found")
endif()

find_library(ADIOS2_LIB NAMES adios2)
if(NOT ADIOS2_LIB)
  message(STATUS "adios2 library not found")
else()
  message(STATUS "adios2 library was found")
endif()

find_library(ADIOS2_SST NAMES adios2_sst)
if(NOT ADIOS2_SST)
	message(STATUS "adios2_sst library not found")
else()
  message(STATUS "adios2_sst library was found")
endif()

find_library(ADIOS2_TAUS NAMES taustubs)
if(NOT ADIOS2_TAUS)
  message(STATUS "taustubs library not found")
else()
  message(STATUS "taustubs library was found")
endif()

<<<<<<< Updated upstream
add_definitions(-DADIOS2_HAVE_MPI)
add_executable(${LIB_NAME} app2_2d.cc)
target_link_libraries(${LIB_NAME} ${MPI_CXX_LIBRARIES} ${ADIOS2_LIB})
=======
find_library(ADIOS2_EV NAMES adios2_evpath)
if(NOT ADIOS2_EV)
	message(STATUS "adios2_evpath library not found")
else()
  message(STATUS "adios2_evpath library was found")
endif()
>>>>>>> Stashed changes

find_library(ADIOS2_FFS NAMES adios2_ffs)
if(NOT ADIOS2_FFS)
  message(STATUS "adios2_ffs library not found")
else()
  message(STATUS "adios2_ffs library was found")
endif()

find_library(ADIOS2_ATL NAMES adios2_atl)
if(NOT ADIOS2_ATL)
  message(STATUS "adios2_atl library not found")
else()
  message(STATUS "adios2_atl library was found")
endif()

find_library(ADIOS2_DILL NAMES adios2_dill)
if(NOT ADIOS2_DILL)
  message(STATUS "adios2_dill library not found")
else()
  message(STATUS "adios2_dill library was found")
endif()

find_library(ADIOS2_ENET NAMES adios2_enet)
if(NOT ADIOS2_ENET)
  message(STATUS "adios2_enet library not found")
else()
  message(STATUS "adios2_enet library was found")
endif()

set(LINK_LIBS ${LINK_LIBS} ${MPI_LIBRARIES} ${ADIOS2_F} ${ADIOS2_LIB} ${ADIOS2_SST} ${ADIOS2_TAUS} ${ADIOS2_EV} ${ADIOS2_FFS} ${ADIOS2_ATL} ${ADIOS2_DILL} ${ADIOS2_ENET}) # ${DRV_LIB})
set(INCLUDE_DIRS ${INCLUDE_DIRS} ${MPI_INCLUDE_PATH} ${ADIOS_INCLUDE_DIRS} ${ADIOS2_INCLUDE_DIRS}) #${DRV_INC})

add_executable(${EXEC} ${SRCS} )
target_link_libraries(${EXEC} ${LINK_LIBS})

INSTALL(TARGETS ${EXEC}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION ${PROJECT_BINARY_DIR}
        ARCHIVE DESTINATION ${PROJECT_BINARY_DIR})
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
