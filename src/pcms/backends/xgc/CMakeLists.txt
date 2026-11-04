# XGC backend. Public header: include/pcms/backends/xgc/.

# Discretization and reverse classification: part of pcms::discretization
# (built regardless of PCMS_ENABLE_XGC, as before).
target_sources(
  pcms_discretization
  PRIVATE discretization.cpp reverse_classification.cpp discretization.hpp
  PUBLIC FILE_SET field_discretization TYPE HEADERS BASE_DIRS
         ${PROJECT_SOURCE_DIR}/include/pcms FILES
         ${PROJECT_SOURCE_DIR}/include/pcms/backends/xgc/reverse_classification.hpp)

# Layout, field data and field factory: part of pcms::field.
if(PCMS_ENABLE_XGC)
  target_sources(
    pcms_field
    PRIVATE layout.cpp xgc.cpp field_data.hpp field_factory.hpp layout.hpp
    PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
           ${PROJECT_SOURCE_DIR}/include/pcms FILES
           ${PROJECT_SOURCE_DIR}/include/pcms/backends/xgc/xgc.hpp)
endif()
