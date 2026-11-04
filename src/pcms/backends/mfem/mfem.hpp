#ifndef PCMS_BACKENDS_MFEM_MFEM_HPP
#define PCMS_BACKENDS_MFEM_MFEM_HPP

// Public entry points of the MFEM backend: wrap MFEM grid functions in place.

#include "pcms/field/field_factory.hpp"
#include "pcms/utility/coordinate_system.hpp"
#include "pcms/utility/memory_spaces.hpp"

#include <Kokkos_Core.hpp>
#include <mfem.hpp>

#include <memory>
#include <string>

namespace pcms::backend::mfem
{

/// Field factory for order-1 H1 (vertex) fields of pfes on pmesh whose values
/// live in gf (not copied). The MFEM objects are borrowed and must outlive
/// the factory and its fields. layout_name identifies the layout in the
/// coupler.
std::shared_ptr<FieldFactory> CreateFieldFactory(
  ::mfem::ParMesh& pmesh, ::mfem::ParFiniteElementSpace& pfes,
  ::mfem::ParGridFunction& gf,
  std::shared_ptr<const CoordinateSystem> coordinate_system,
  std::string layout_name = "");

/// Per-vertex overlap mask (true where a vertex belongs to an element with
/// the given attribute), for ApplicationComm::SetLayoutOverlapMask.
Kokkos::View<bool*, HostMemorySpace> OverlapMaskFromAttribute(
  ::mfem::ParMesh& pmesh, int attribute);

} // namespace pcms::backend::mfem

#endif // PCMS_BACKENDS_MFEM_MFEM_HPP
