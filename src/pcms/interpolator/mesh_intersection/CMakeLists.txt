file(GLOB_RECURSE MESH_INTERSECTION_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

file(GLOB_RECURSE MESH_INTERSECTION_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_library(pcms_mesh_intersection
    ${MESH_INTERSECTION_HEADERS}
    ${MESH_INTERSECTION_SOURCES}
)

target_include_directories(pcms_mesh_intersection PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include/pcms/interpolator/mesh_intersection>
)

target_link_libraries(pcms_mesh_intersection PUBLIC
    pcms::core
    Omega_h::omega_h
)

if (PCMS_ENABLE_MESHFIELDS)
    target_link_libraries(pcms_mesh_intersection PUBLIC meshfields::meshfields)
endif()

if (PCMS_ENABLE_PETSC)
	target_include_directories(pcms_mesh_intersection PUBLIC ${PCMS_PETSC_INCLUDES})
	target_link_libraries(pcms_mesh_intersection PUBLIC ${PCMS_PETSC_LIBRARIES})
endif()

target_compile_definitions(pcms_mesh_intersection PUBLIC R3D_USE_CUDA)

install(FILES ${MESH_INTERSECTION_HEADERS}
        DESTINATION include/pcms/interpolator/mesh_intersection)

install(TARGETS pcms_mesh_intersection
        EXPORT meshIntersectionTargets
        INCLUDES DESTINATION include/pcms/interpolator/mesh_intersection
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(EXPORT meshIntersectionTargets
        FILE pcms_mesh_intersection-targets.cmake
        NAMESPACE pcms::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
