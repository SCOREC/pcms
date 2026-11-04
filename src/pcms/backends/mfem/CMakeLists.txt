# MFEM backend: layout, field data and field factory, part of pcms::field.
target_sources(pcms_field PRIVATE layout.cpp field_data.cpp mfem.cpp)
target_sources(
  pcms_field
  PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../.. ${PROJECT_SOURCE_DIR}/include/pcms
         FILES layout.hpp field_data.hpp field_factory.hpp mfem.hpp)
target_link_libraries(pcms_field PUBLIC mfem)
