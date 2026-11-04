find_package(KokkosKernels REQUIRED)

set(PCMS_FIELD_HEADERS
        field.hpp
        coordinate_map.hpp
        value_view.hpp
        element_dispatch.hpp
        evaluation_request.hpp
        gid_permutation.hpp
        field_layout.hpp
        layout/empty.hpp
        out_of_bounds_policy.hpp
        field_data.hpp
        field_factory.hpp
        function_space.hpp
        point_evaluator.hpp
        point_evaluator_factory.hpp
        transformed_source.hpp
        data/simple.hpp
        function_space/polynomial_reconstruction.hpp
        layout/point_cloud.hpp
        data/point_cloud.hpp
)
set(PCMS_FIELD_SOURCES
        field_layout.cpp
        coordinate_map.cpp
        transformed_source.cpp
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
            function_space/lagrange.hpp
            function_space/spline.hpp
            field_evaluator_factory.hpp
            layout/omega_h_entity.hpp
            layout/omega_h_lagrange.hpp
            evaluator/omega_h_lagrange.hpp
            layout/uniform_grid.hpp
            uniform_grid_binary_field.hpp
            evaluator/uniform_grid.hpp
            evaluator/spline_interpolator.hpp
            evaluator/uniform_grid_spline.hpp
            eqdsk_field.hpp
            point_cloud_localization.hpp
            evaluator/mls_interpolation.hpp
            evaluator/mls_interpolation_impl.hpp
            evaluator/pcms_interpolator_aliases.hpp
            evaluator/pcms_interpolator_view_utils.hpp
            evaluator/pcms_interpolator_logger.hpp
            evaluator/mls_point_cloud.hpp
            evaluator/point_cloud.hpp
    )
endif ()

if (PCMS_ENABLE_MESHFIELDS)
    list(APPEND PCMS_FIELD_SOURCES
            ../backends/omega_h_meshfields/layout.cpp
    )
    list(APPEND PCMS_FIELD_HEADERS
            ../backends/omega_h_meshfields/layout.hpp
            ../backends/omega_h_meshfields/evaluator_backend.hpp
            ../backends/omega_h_meshfields/evaluator.hpp
            ../backends/omega_h_meshfields/field_data.hpp
    )
endif ()

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
