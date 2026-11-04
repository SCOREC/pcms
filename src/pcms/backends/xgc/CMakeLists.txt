# XGC backend.

# The XGC discretization and reverse classification are part of
# pcms::discretization (built regardless of PCMS_ENABLE_XGC, as before).
target_sources(pcms_discretization PRIVATE discretization.cpp
                                           reverse_classification.cpp)
target_sources(
  pcms_discretization
  PUBLIC FILE_SET field_discretization TYPE HEADERS BASE_DIRS
         ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${PROJECT_SOURCE_DIR}/include/pcms FILES discretization.hpp
         ../../../../include/pcms/backends/xgc/reverse_classification.hpp)

# The XGC field layout, field data and field factory are part of pcms::field.
if(PCMS_ENABLE_XGC)
  target_sources(pcms_field PRIVATE layout.cpp xgc.cpp)
  target_sources(
    pcms_field
    PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
           ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${PROJECT_SOURCE_DIR}/include/pcms FILES field_data.hpp
           field_factory.hpp layout.hpp ../../../../include/pcms/backends/xgc/xgc.hpp)
endif()
