find_package(KokkosKernels REQUIRED)

set(PCMS_FIELD_HEADERS
        field.h
        coordinate_system.h
        element_dispatch.h
        evaluation_request.h
        gid_permutation.hpp
        field_layout.h
        field_operator.hpp
        field_metadata.h
        layout/empty.h
        out_of_bounds_policy.h
        field_data.h
        field_factory.h
        function_space.h
        point_evaluator.h
        data/simple.h
        function_space/polynomial_reconstruction.hpp
        layout/point_cloud.h
        data/point_cloud.h
)
set(PCMS_FIELD_SOURCES
        field_layout.cpp
        gid_permutation.cpp
        function_space/polynomial_reconstruction.cpp
        layout/empty.cpp
        layout/point_cloud.cpp
        data/point_cloud.cpp
)

if (PCMS_ENABLE_OMEGA_H)
    list(APPEND PCMS_FIELD_SOURCES
            function_space/lagrange.cpp
            function_space/spline.cpp
            layout/omega_h_entity.cpp
            layout/omega_h_lagrange.cpp
            layout/uniform_grid.cpp
            evaluator/mls_interpolation.cpp
    )
    list(APPEND PCMS_FIELD_HEADERS
            function_space/lagrange.h
            function_space/spline.h
            field_evaluator_factory.h
            layout/omega_h_entity.h
            layout/omega_h_lagrange.h
            evaluator/omega_h_lagrange.h
            layout/uniform_grid.h
            uniform_grid_binary_field.h
            evaluator/uniform_grid.h
            evaluator/uniform_grid_spline.h
            eqdsk_field.h
            evaluator/mls_options.h
            evaluator/mls_interpolation.hpp
            evaluator/mls_interpolation_impl.hpp
            evaluator/pcms_interpolator_aliases.hpp
            evaluator/pcms_interpolator_view_utils.hpp
            evaluator/pcms_interpolator_logger.hpp
            evaluator/mls_point_cloud.h
            evaluator/point_cloud.h
    )
endif ()

if (PCMS_ENABLE_MESHFIELDS)
    list(APPEND PCMS_FIELD_SOURCES
            layout/mesh_fields.cpp
    )
    list(APPEND PCMS_FIELD_HEADERS
            layout/mesh_fields.h
            evaluator/mesh_fields_backend.h
            evaluator/mesh_fields.h
            data/mesh_fields.h
    )
endif ()

if(PCMS_ENABLE_XGC)
    list(APPEND PCMS_FIELD_HEADERS
            data/xgc.h
            function_space/xgc.h
            layout/xgc.h)
    list(APPEND PCMS_FIELD_SOURCES
            layout/xgc.cpp)
endif()

add_library(pcms_field ${PCMS_FIELD_SOURCES})
set_target_properties(pcms_field PROPERTIES
        OUTPUT_NAME pcmsfield
        EXPORT_NAME field)
target_sources(pcms_field PUBLIC
        FILE_SET field
        TYPE HEADERS
        BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/..
        FILES ${PCMS_FIELD_HEADERS})
add_library(pcms::field ALIAS pcms_field)
target_compile_features(pcms_field PUBLIC cxx_std_20)

target_link_libraries(pcms_field PUBLIC
        pcms::utility
        pcms::discretization
        pcms::localization)

if(PCMS_ENABLE_OMEGA_H)
    target_link_libraries(pcms_field PUBLIC Omega_h::omega_h)
    target_link_libraries(pcms_field PRIVATE Kokkos::kokkoskernels)
endif()

if (PCMS_ENABLE_MESHFIELDS)
    target_link_libraries(pcms_field PUBLIC meshfields::meshfields)
endif ()

target_include_directories(pcms_field INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

install(
        TARGETS pcms_field
        EXPORT pcms_field-targets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/field
        FILE_SET field DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

install(
        EXPORT pcms_field-targets
        NAMESPACE pcms::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pcms)
