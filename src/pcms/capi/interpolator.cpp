//
// Created by hasanm4 on 2/17/25.
//
#include <pcms/capi/kokkos.h>
#include "interpolator.h"
#include <interpolation_base.h>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>


//[[nodiscard]]
PcmsInterpolatorHandle pcms_create_interpolator(PcmsInterpolatorOHMeshHandle oh_mesh, double radius)
{
  auto* source_mesh = reinterpret_cast<Omega_h::Mesh*>(oh_mesh.mesh_handle);
  auto* interpolator = new MLSInterpolationHandler(*source_mesh, radius);
  return {reinterpret_cast<void*>(interpolator)};
}

void pcms_destroy_interpolator(PcmsInterpolatorHandle interpolator)
{
  if (interpolator.pointer != nullptr) {
    delete reinterpret_cast<MLSInterpolationHandler*>(interpolator.pointer);
  }
}

PcmsInterpolatorOHMeshHandle read_oh_mesh(const char * filename){
  auto fname = std::string(filename);
  // trim the filename since it is coming from c or fortran api which may have extra spaces at the end
  fname.erase(fname.find_last_not_of(" \n\r\t")+1);
  auto* oh_lib = new Omega_h::Library();
  auto* mesh = new Omega_h::Mesh(Omega_h::binary::read(fname, oh_lib->world()));

  return {reinterpret_cast<void*>(mesh), reinterpret_cast<void*>(oh_lib)};
}

void release_oh_mesh(PcmsInterpolatorOHMeshHandle oh_mesh){
  if (oh_mesh.mesh_handle != nullptr) {
    assert(oh_mesh.lib_handle != nullptr);
    delete reinterpret_cast<Omega_h::Mesh*>(oh_mesh.mesh_handle);
    delete reinterpret_cast<Omega_h::Library*>(oh_mesh.lib_handle);
  }
}

void pcms_interpolate(PcmsInterpolatorHandle interpolator, double* input, int input_size, double* output, int output_size) {
  auto* mls_interpolator = reinterpret_cast<MLSInterpolationHandler*>(interpolator.pointer);

  OMEGA_H_CHECK_PRINTF(input_size == mls_interpolator->getSourceSize(), "Input array size does not match the source size %d != %d\n", input_size, mls_interpolator->getSourceSize());
  OMEGA_H_CHECK_PRINTF(output_size == mls_interpolator->getTargetSize(), "Output array size does not match the target size %d != %d\n", output_size, mls_interpolator->getTargetSize());

  pcms::ScalarArrayView<double, pcms::HostMemorySpace> input_array(input, input_size);
  pcms::ScalarArrayView<double, pcms::HostMemorySpace> output_array(output, output_size);

  mls_interpolator->eval(input_array, output_array);
}
