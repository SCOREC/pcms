# MFEM backend: layout, field data and field factory, part of pcms::field.
# Public header: include/pcms/backends/mfem/.
target_sources(
  pcms_field
  PRIVATE layout.cpp field_data.cpp mfem.cpp layout.hpp field_data.hpp
          field_factory.hpp
  PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS ${PROJECT_SOURCE_DIR}/include/pcms
         FILES ${PROJECT_SOURCE_DIR}/include/pcms/backends/mfem/mfem.hpp)
target_link_libraries(pcms_field PUBLIC mfem)
