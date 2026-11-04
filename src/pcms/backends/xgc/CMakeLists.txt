# XGC backend.

# The XGC layout (with the reverse classification), field data and field
# factory are part of pcms::field.
if(PCMS_ENABLE_XGC)
  target_sources(pcms_field PRIVATE layout.cpp reverse_classification.cpp xgc.cpp)
  target_sources(
    pcms_field
    PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
           ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${PROJECT_SOURCE_DIR}/include/pcms FILES field_data.hpp
           field_factory.hpp layout.hpp ../../../../include/pcms/backends/xgc/reverse_classification.hpp ../../../../include/pcms/backends/xgc/xgc.hpp)
endif()
