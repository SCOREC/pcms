# XGC backend.

# The XGC layout (with the reverse classification), field data and field
# factory are part of pcms::field.
if(PCMS_ENABLE_XGC)
  target_sources(pcms_field PRIVATE layout.cpp reverse_classification.cpp)
  target_sources(
    pcms_field
    PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
           ${CMAKE_CURRENT_SOURCE_DIR}/../.. FILES field_data.hpp
           field_factory.hpp layout.hpp reverse_classification.hpp)
endif()
