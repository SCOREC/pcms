#include "pcms/adapter/omega_h/omega_h_field_layout.h"
#include "omega_h_field_layout.h"

namespace pcms
{
/*
 * Field Layout
 */

OmegaHFieldLayout::OmegaHFieldLayout(Omega_h::Mesh& mesh,
                                     OmegaHFieldLayoutLocation location,
                                     int num_components)
  : mesh_(mesh), location_(location), num_components_(num_components)
{
}

int OmegaHFieldLayout::GetNumComponents() const
{
  return num_components_;
}
// nodes for standard lagrange FEM
LO OmegaHFieldLayout::GetNumOwnedDofHolder() const
{
  switch (location_) {
    case OmegaHFieldLayoutLocation::PieceWise:
      // global ids are faces in 2D and regions in 3D
      return mesh_.nents(mesh_.dim());
      break;
    case OmegaHFieldLayoutLocation::Linear:
      // global ids are vertices
      return mesh_.nents(0);
      break;
  }
  // should never reach here, but return 0 to avoid compiler warning
  return 0;
}

GO OmegaHFieldLayout::GetNumGlobalDofHolder() const
{
  switch (location_) {
    case OmegaHFieldLayoutLocation::PieceWise:
      // global ids are faces in 2D and regions in 3D
      return mesh_.nglobal_ents(mesh_.dim());
      break;
    case OmegaHFieldLayoutLocation::Linear:
      // global ids are vertices
      return mesh_.nglobal_ents(0);
      break;
  }
  // should never reach here, but return 0 to avoid compiler warning
  return 0;
}

GlobalIDView<HostMemorySpace> OmegaHFieldLayout::GetOwnedGids()
{
  int dim = 0;
  switch (location_) {
    case OmegaHFieldLayoutLocation::PieceWise:
      // global ids are faces in 2D and regions in 3D
      dim = mesh_.dim();
      break;
    case OmegaHFieldLayoutLocation::Linear:
      // ids are just global ids of the vertex
      dim = 0;
      break;
  }
  auto gids = mesh_.globals(dim);
  static_assert(std::is_same_v<HostMemorySpace, DefaultExecutionSpace::memory_space>, "types must match");
  return GlobalIDView<HostMemorySpace>(gids.data(), gids.size());
}

OmegaHFieldLayoutLocation OmegaHFieldLayout::GetLocation() const
{
  return location_;
}

bool OmegaHFieldLayout::IsDistributed()
{
  return true;
}

} // namespace pcms
