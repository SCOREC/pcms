# XGC backend. Public header: include/pcms/backends/xgc/.

# Layout (with the reverse classification), field data and field factory: part
# of pcms::field.
if(PCMS_ENABLE_XGC)
  # the reverse classification is built from an Omega_h mesh and uses MPI
  target_link_libraries(pcms_field PRIVATE MPI::MPI_CXX)
  # The public headers include MPI (and Omega_h, see ../omega_h), which their
  # installed users link themselves; in the build tree pcms::field provides it,
  # so the header-set verification compiles the headers.
  target_link_libraries(pcms_field INTERFACE $<BUILD_INTERFACE:MPI::MPI_CXX>)
  target_sources(
    pcms_field
    PRIVATE layout.cpp reverse_classification.cpp xgc.cpp field_data.hpp
            field_factory.hpp layout.hpp
    PUBLIC FILE_SET field TYPE HEADERS BASE_DIRS
           ${PROJECT_SOURCE_DIR}/include/pcms FILES
           ${PROJECT_SOURCE_DIR}/include/pcms/backends/xgc/reverse_classification.hpp
           ${PROJECT_SOURCE_DIR}/include/pcms/backends/xgc/detail/reverse_classification.hpp
           ${PROJECT_SOURCE_DIR}/include/pcms/backends/xgc/xgc.hpp)
endif()
