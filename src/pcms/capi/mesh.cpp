#include "mesh.h"

#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

PcmsOmegaHLibraryHandle pcms_create_omega_h_library()
{
#ifdef Omega_H_USE_MPI
  Omega_h::Library* oh_lib =
    new Omega_h::Library(nullptr, nullptr, MPI_COMM_SELF);
#else
  Omega_h::Library* oh_lib = new Omega_h::Library(nullptr, nullptr);
#endif

  return {reinterpret_cast<void*>(oh_lib)};
}

PcmsOmegaHMeshHandle pcms_create_omega_h_mesh(
  const char* filename, const PcmsOmegaHLibraryHandle oh_lib_handle)
{
  auto fname = std::string(filename);
  // trim the filename since it is coming from c or fortran api which may have
  // extra spaces at the end
  fname.erase(fname.find_last_not_of(" \n\r\t") + 1);
  Omega_h::Library* oh_lib =
    reinterpret_cast<Omega_h::Library*>(oh_lib_handle.lib_handle);
  auto* mesh = new Omega_h::Mesh(Omega_h::binary::read(fname, oh_lib->world()));

  return {reinterpret_cast<void*>(mesh)};
}

void pcms_destroy_omega_h_mesh(PcmsOmegaHMeshHandle oh_mesh)
{
  if (oh_mesh.mesh_handle != nullptr) {
    delete reinterpret_cast<Omega_h::Mesh*>(oh_mesh.mesh_handle);
  }
}

void pcms_destroy_omega_h_library(PcmsOmegaHLibraryHandle oh_lib_handle)
{
  if (oh_lib_handle.lib_handle != nullptr) {
    delete reinterpret_cast<Omega_h::Library*>(oh_lib_handle.lib_handle);
  }
}
