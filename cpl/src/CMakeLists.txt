# source files to compile
set(SRCS cpl_init.cc)
# any headers which need to be public/installed
#  this does not need to be ALL headers
set(HDRS cpl_init.h)

set(LIB_NAME wdm_coupler)

find_package(ADIOS REQUIRED)
set(INCLUDE_DIRS ${ADIOS_INCLUDE_DIRS})
set(LINK_LIBS ${ADIOS_LIBRARIES})

add_library(${LIB_NAME} ${SRCS})
target_include_directories(${LIB_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  ${INCLUDE_DIRS})

target_link_libraries(${LIB_NAME} ${LINK_LIBS})

install(FILES ${HDRS} DESTINATION include)
install(TARGETS ${LIB_NAME} LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
