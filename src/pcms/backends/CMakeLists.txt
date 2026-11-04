# Backend adapters. A backend is a directory and a namespace
# (pcms::backend::<name>), not a separate library: each backend adds its
# sources and headers to the core target that uses it.
add_subdirectory(xgc)
if(PCMS_ENABLE_MFEM)
  add_subdirectory(mfem)
endif()
