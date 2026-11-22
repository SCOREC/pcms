set(PCMS_FIELD_HEADERS
        ../../../include/pcms/field/field.hpp
        ../../../include/pcms/field/field_metadata.hpp
        ../../../include/pcms/field/units.hpp
        ../../../include/pcms/field/detail/field_metadata.hpp
        ../../../include/pcms/field/detail/units.hpp
        ../../../include/pcms/field/coordinate_map.hpp
        ../../../include/pcms/field/embedding.hpp
        ../../../include/pcms/field/detail/embedding.hpp
        ../../../include/pcms/field/value_transform.hpp
        ../../../include/pcms/field/value_view.hpp
        element_dispatch.hpp
        gid_permutation.hpp
        point_rows.hpp
        ../../../include/pcms/field/field_layout.hpp
        layout/empty.hpp
        ../../../include/pcms/field/out_of_bounds_policy.hpp
        ../../../include/pcms/field/field_data.hpp
        ../../../include/pcms/field/field_ops.hpp
        ../../../include/pcms/field/detail/field_data_ops.hpp
        ../../../include/pcms/field/detail/field_ops.hpp
        ../../../include/pcms/field/field_factory.hpp
        ../../../include/pcms/field/function_space.hpp
        ../../../include/pcms/field/point_evaluator.hpp
        ../../../include/pcms/field/point_evaluator_factory.hpp
        ../../../include/pcms/field/query_points.hpp
        ../../../include/pcms/field/differential_operator.hpp
        ../../../include/pcms/field/gradient_source.hpp
        ../../../include/pcms/field/detail/transformed_point_evaluator.hpp
        ../../../include/pcms/field/detail/coordinate_map.hpp
        ../../../include/pcms/field/detail/value_transform.hpp
        ../../../include/pcms/field/detail/field_factory.hpp
        ../../../include/pcms/field/detail/field_layout.hpp
        ../../../include/pcms/field/detail/gradient_source.hpp
        ../../../include/pcms/field/detail/point_evaluator.hpp
        ../../../include/pcms/field/detail/point_evaluator_factory.hpp
        ../../../include/pcms/field/detail/value_view.hpp
        ../../../include/pcms/field/detail/simple_field_data.hpp
        ../../../include/pcms/field/function_spaces.hpp
        ../../../include/pcms/field/mls_options.hpp
)
set(PCMS_FIELD_SOURCES
        field_metadata.cpp
        units.cpp
        embedding.cpp
        field_layout.cpp
        function_space.cpp
        function_spaces.cpp
        gradient_source.cpp
        query_points.cpp
        gid_permutation.cpp
        layout/empty.cpp
        function_space/lagrange.cpp
        function_space/spline.cpp
        function_space/periodic.cpp
        function_space/product_geometry.cpp
        function_space/product_space.cpp
)
list(APPEND PCMS_FIELD_HEADERS
)

add_library(pcms_field ${PCMS_FIELD_SOURCES})
set_target_properties(pcms_field PROPERTIES
        OUTPUT_NAME pcmsfield
        EXPORT_NAME field)
# Only the public headers (include/pcms/...) are installed; the private
# headers next to the sources are private sources of the target.
set(PCMS_FIELD_HEADERS_PUBLIC ${PCMS_FIELD_HEADERS})
list(FILTER PCMS_FIELD_HEADERS_PUBLIC INCLUDE REGEX "/include/pcms/")
set(PCMS_FIELD_HEADERS_PRIVATE ${PCMS_FIELD_HEADERS})
list(FILTER PCMS_FIELD_HEADERS_PRIVATE EXCLUDE REGEX "/include/pcms/")
target_sources(pcms_field
  PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include/pcms
         FILES ${PCMS_FIELD_HEADERS_PUBLIC}
  PRIVATE ${PCMS_FIELD_HEADERS_PRIVATE})
add_library(pcms::field ALIAS pcms_field)
target_compile_features(pcms_field PUBLIC cxx_std_20)

target_link_libraries(pcms_field PUBLIC
        pcms::utility
        pcms::localization
        pcms::geometry)
if(PCMS_ENABLE_UNITS)
    target_link_libraries(pcms_field PRIVATE units::units)
endif()

if(PCMS_ENABLE_OMEGA_H)
    # private: core public headers do not include Omega_h. The Omega_h backend
    # header (pcms/backends/omega_h/omega_h.hpp) does; its users link
    # Omega_h::omega_h themselves
    target_link_libraries(pcms_field PRIVATE Omega_h::omega_h)
endif()


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
        DESTINATION ${PCMS_INSTALL_CMAKEDIR})
