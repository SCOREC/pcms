# Omega_h adapter. Public header: include/pcms/backends/omega_h/.

# Adjacency search: part of pcms::localization, which always uses Omega_h.
target_sources(pcms_localization PRIVATE adj_search.cpp adj_search.hpp
                                         queue_visited.hpp)

if(PCMS_ENABLE_OMEGA_H)
  # layouts, field data, Lagrange evaluator, the geometry, the uniform-grid
  # binary field and the public entry points: part of pcms::field
  target_sources(
    pcms_field
    PRIVATE entity_layout.cpp lagrange_layout.cpp geometry.cpp omega_h.cpp
            uniform_grid_binary_field.cpp entity_layout.hpp lagrange_layout.hpp
            geometry.hpp field_data.hpp lagrange_evaluator.hpp
            uniform_grid_binary_field.hpp
    PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
           ${PROJECT_SOURCE_DIR}/include/pcms FILES
           ${PROJECT_SOURCE_DIR}/include/pcms/backends/omega_h/omega_h.hpp)
endif()
