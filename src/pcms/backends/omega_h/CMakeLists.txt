# Omega_h adapter.

# Adjacency search: part of pcms::localization, which always uses Omega_h.
# Localization installs its own headers flat, so these are installed here.
target_sources(pcms_localization PRIVATE adj_search.cpp)
install(FILES adj_search.hpp queue_visited.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/backends/omega_h)

if(PCMS_ENABLE_OMEGA_H)

  # layouts, field data, Lagrange evaluator and the uniform-grid binary
  # field: part of pcms::field
  target_sources(pcms_field PRIVATE entity_layout.cpp lagrange_layout.cpp geometry.cpp omega_h.cpp
                                    uniform_grid_binary_field.cpp)
  target_sources(
    pcms_field
    PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
           ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${PROJECT_SOURCE_DIR}/include/pcms FILES entity_layout.hpp
           lagrange_layout.hpp geometry.hpp omega_h.hpp field_data.hpp lagrange_evaluator.hpp
           uniform_grid_binary_field.hpp)
endif()
