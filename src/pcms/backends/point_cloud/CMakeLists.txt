# Point-cloud backend (no external dependency).

# discretization: part of pcms::discretization
target_sources(pcms_discretization PRIVATE discretization.cpp)
target_sources(
  pcms_discretization
  PUBLIC FILE_SET field_discretization TYPE HEADERS BASE_DIRS
         ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${PROJECT_SOURCE_DIR}/include/pcms FILES discretization.hpp)

# layout and field data: part of pcms::field
target_sources(pcms_field PRIVATE layout.cpp field_data.cpp)
target_sources(
  pcms_field
  PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${PROJECT_SOURCE_DIR}/include/pcms
         FILES layout.hpp field_data.hpp)

# MLS evaluation and localization use Omega_h arrays (as before)
if(PCMS_ENABLE_OMEGA_H)
  target_sources(
    pcms_field
    PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
           ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${PROJECT_SOURCE_DIR}/include/pcms FILES localization.hpp
           mls_evaluator.hpp evaluator.hpp)
endif()
