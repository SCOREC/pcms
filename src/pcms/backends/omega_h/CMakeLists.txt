# Omega_h adapter. Public header: include/pcms/backends/omega_h/.

if(PCMS_ENABLE_OMEGA_H)
  # the grid search over the mesh's elements and the mesh-adjacency
  # neighbourhood: part of pcms::localization
  target_sources(
    pcms_localization PRIVATE point_search.cpp adj_search.cpp point_search.hpp
                              adj_search.hpp queue_visited.hpp)

  # layouts, field data, Lagrange evaluator, the geometry, the uniform-grid
  # binary field and the public entry points: part of pcms::field
  target_sources(
    pcms_field
    PRIVATE entity_layout.cpp lagrange_layout.cpp geometry.cpp omega_h.cpp
            uniform_grid_binary_field.cpp entity_layout.hpp lagrange_layout.hpp
            field_data.hpp lagrange_evaluator.hpp
            uniform_grid_binary_field.hpp omega_h_geometry.hpp
    PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
           ${PROJECT_SOURCE_DIR}/include/pcms FILES
           ${PROJECT_SOURCE_DIR}/include/pcms/backends/omega_h/omega_h.hpp
           ${PROJECT_SOURCE_DIR}/include/pcms/backends/omega_h/geometry.hpp)
  # The public header includes Omega_h, which its installed users link
  # themselves; in the build tree pcms::field provides it, so the header-set
  # verification compiles the header.
  target_link_libraries(pcms_field INTERFACE $<BUILD_INTERFACE:Omega_h::omega_h>)
endif()
