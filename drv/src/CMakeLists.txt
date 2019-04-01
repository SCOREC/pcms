set(CSRC )
set(FSRC )

set(HDRS )
set(MODS )

set(LIB_NAME wdm_minapp)

find_package(wdm_coupler REQUIRED)

find_package(XGC REQUIRED)
find_package(GENE REQIRED)

add_library(${LIB_NAME} ${FSRC} ${CSRC})
target_include_directories(${LIB_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
  ${INCLUDE_DIRS})
target_link_libraries(${LIB_NAME} ${LINK_LIBS})

install(FILES "${HDRS} ${MODS}" DESTINATION include)
install(TARGETS ${LIB_NAME} LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)


