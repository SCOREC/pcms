# Omega_h adapter.

# Adjacency search: part of pcms::localization, which always uses Omega_h.
# Localization installs its own headers flat, so these are installed here.
target_sources(pcms_localization PRIVATE adj_search.cpp)
install(FILES adj_search.hpp queue_visited.hpp
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pcms/backends/omega_h)

if(PCMS_ENABLE_OMEGA_H)
  # discretization: part of pcms::discretization
  target_sources(pcms_discretization PRIVATE discretization.cpp omega_h.cpp)
  target_sources(
    pcms_discretization
    PUBLIC FILE_SET field_discretization TYPE HEADERS BASE_DIRS
           ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${PROJECT_SOURCE_DIR}/include/pcms FILES discretization.hpp omega_h.hpp)

  # layouts, Lagrange evaluator and the uniform-grid binary field: part of
  # pcms::field
  target_sources(pcms_field PRIVATE entity_layout.cpp lagrange_layout.cpp)
  target_sources(
    pcms_field
    PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
           ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${PROJECT_SOURCE_DIR}/include/pcms FILES entity_layout.hpp
           lagrange_layout.hpp lagrange_evaluator.hpp
           uniform_grid_binary_field.hpp)
endif()
