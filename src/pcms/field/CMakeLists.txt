find_package(KokkosKernels REQUIRED)

set(PCMS_FIELD_HEADERS
        field.hpp
        coordinate_system.hpp
        coordinate_view.hpp
        coordinate_systems/cartesian.hpp
        coordinate_systems/cylindrical.hpp
        coordinate_map.hpp
        value_view.hpp
        element_dispatch.hpp
        gid_permutation.hpp
        field_layout.hpp
        layout/empty.hpp
        out_of_bounds_policy.hpp
        field_data.hpp
        field_factory.hpp
        function_space.hpp
        point_evaluator.hpp
        point_evaluator_factory.hpp
        query_points.hpp
        differential_operator.hpp
        gradient_source.hpp
        transformed_point_evaluator.hpp
        data/simple.hpp
        function_space/polynomial_reconstruction.hpp
        layout/point_cloud.hpp
        data/point_cloud.hpp
)
set(PCMS_FIELD_SOURCES
        field_layout.cpp
        coordinate_system.cpp
        function_space.cpp
        gradient_source.cpp
        query_points.cpp
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
            evaluator/omega_h_geometry.cpp
    )
    list(APPEND PCMS_FIELD_HEADERS
            function_space/lagrange.hpp
            function_space/spline.hpp
            layout/omega_h_entity.hpp
            layout/omega_h_lagrange.hpp
            evaluator/omega_h_geometry.hpp
            evaluator/omega_h_lagrange.hpp
            evaluator/uniform_grid_geometry.hpp
            evaluator/point_cloud_geometry.hpp
            data/omega_h.hpp
            layout/uniform_grid.hpp
            uniform_grid_binary_field.hpp
            evaluator/uniform_grid.hpp
            evaluator/spline_interpolator.hpp
            evaluator/uniform_grid_spline.hpp
            eqdsk_field.hpp
            evaluator/mls_options.hpp
            evaluator/mls_interpolation.hpp
            evaluator/mls_interpolation_impl.hpp
            evaluator/pcms_interpolator_aliases.hpp
            evaluator/pcms_interpolator_view_utils.hpp
            evaluator/pcms_interpolator_logger.hpp
            evaluator/mls_point_cloud.hpp
            evaluator/point_cloud.hpp
    )
endif ()

if (PCMS_ENABLE_MFEM)
    list(APPEND PCMS_FIELD_SOURCES
            layout/mfem.cpp
            data/mfem.cpp
    )
    list(APPEND PCMS_FIELD_HEADERS
            layout/mfem.hpp
            data/mfem.hpp
            function_space/mfem.hpp
    )
endif ()

if (PCMS_ENABLE_MESHFIELDS)
    list(APPEND PCMS_FIELD_SOURCES
            layout/mesh_fields.cpp
    )
    list(APPEND PCMS_FIELD_HEADERS
            layout/mesh_fields.hpp
            evaluator/mesh_fields_backend.hpp
            evaluator/mesh_fields.hpp
            data/mesh_fields.hpp
    )
endif ()

if(PCMS_ENABLE_XGC)
    list(APPEND PCMS_FIELD_HEADERS
            data/xgc.hpp
            function_space/xgc.hpp
            layout/xgc.hpp
            layout/xgc_reverse_classification.hpp)
    list(APPEND PCMS_FIELD_SOURCES
            layout/xgc.cpp
            layout/xgc_reverse_classification.cpp)
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
        pcms::localization
        pcms::geometry)

if(PCMS_ENABLE_OMEGA_H)
    target_link_libraries(pcms_field PUBLIC Omega_h::omega_h)
    target_link_libraries(pcms_field PRIVATE Kokkos::kokkoskernels)
endif()

if(PCMS_ENABLE_MFEM)
    target_link_libraries(pcms_field PUBLIC mfem)
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
