//
// Created by hasanm4 on 2/17/25.
//
#include <pcms/capi/kokkos.h>
#include "interpolator.h"
#include <interpolation_base.h>


//[[nodiscard]]
PcmsInterpolatorHandle pcms_create_interpolator(PcmsInterpolatorOHMeshHandle oh_mesh, double radius)
{
  auto* source_mesh = reinterpret_cast<Omega_h::Mesh*>(oh_mesh.pointer);
  auto* interpolator = new MLSInterpolationHandler(*source_mesh, radius);
  return {reinterpret_cast<void*>(interpolator)};
}

void pcms_destroy_interpolator(PcmsInterpolatorHandle interpolator)
{
  if (interpolator.pointer != nullptr) {
    delete reinterpret_cast<MLSInterpolationHandler*>(interpolator.pointer);
  }
}
