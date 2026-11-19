# Backend adapters. A backend is a directory and a namespace
# (pcms::<name>), not a separate library: each backend adds its
# sources and headers to the core target that uses it.
add_subdirectory(omega_h)
add_subdirectory(xgc)
add_subdirectory(point_cloud)
add_subdirectory(uniform_grid)
if(PCMS_ENABLE_MFEM)
  add_subdirectory(mfem)
endif()
# the MeshFields spaces on Omega_h meshes (PCMS_ENABLE_OMEGA_H_MESHFIELDS,
# which needs PCMS_ENABLE_OMEGA_H)
if(PCMS_ENABLE_OMEGA_H_MESHFIELDS)
  add_subdirectory(omega_h_meshfields)
endif()
