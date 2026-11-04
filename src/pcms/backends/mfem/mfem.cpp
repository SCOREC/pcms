#include "pcms/backends/mfem/mfem.hpp"
#include "pcms/backends/mfem/field_factory.hpp"
#include "pcms/backends/mfem/layout.hpp"

namespace pcms::backend::mfem
{

std::shared_ptr<FieldFactory> CreateFieldFactory(
  ::mfem::ParMesh& pmesh, ::mfem::ParFiniteElementSpace& pfes,
  ::mfem::ParGridFunction& gf,
  std::shared_ptr<const CoordinateSystem> coordinate_system,
  std::string layout_name)
{
  return std::make_shared<MFEMFieldFactory>(
    pmesh, pfes, gf, std::move(coordinate_system), std::move(layout_name));
}

Kokkos::View<bool*, HostMemorySpace> OverlapMaskFromAttribute(
  ::mfem::ParMesh& pmesh, int attribute)
{
  return MFEMLayout::OverlapMaskFromAttribute(pmesh, attribute);
}

} // namespace pcms::backend::mfem
