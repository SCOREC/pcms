# Backend adapters. A backend is a directory and a namespace
# (pcms::backend::<name>), not a separate library: each backend adds its
# sources and headers to the core target that uses it.
add_subdirectory(xgc)
add_subdirectory(point_cloud)
add_subdirectory(uniform_grid)
if(PCMS_ENABLE_MFEM)
  add_subdirectory(mfem)
endif()
# bridges: built only when both backends are enabled
if(PCMS_ENABLE_OMEGA_H AND PCMS_ENABLE_MESHFIELDS)
  add_subdirectory(omega_h_meshfields)
endif()
