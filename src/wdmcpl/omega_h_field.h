#ifndef WDM_COUPLING_OMEGA_H_FIELD_H
#define WDM_COUPLING_OMEGA_H_FIELD_H
#include "wdmcpl/types.h"
#include "wdmcpl/external/span.h"
#include <Omega_h_mesh.hpp>
#include "wdmcpl/field.h"
#include "wdmcpl/coordinate_systems.h"
#include <Kokkos_core.hpp>
#include <wdmcpl/assert.h>
#include <Omega_h_for.hpp>
#include "wdmcpl/arrays.h"
#include "wdmcpl/field_evaluation_methods.h"

// FIXME add executtion spaces (don't use kokkos exe spaces directly)

namespace wdmcpl
{

// TODO different types dependent on active OmegaHBackend
struct OmegaHMemorySpace
{
  using type = typename Kokkos::DefaultExecutionSpace::memory_space;
};

template <typename T,
          typename CoordinateElementType =
            Real> // CoordinateElement<Cartesian, Real>>
class OmegaHField
{
public:
  // using memory_space = OmegaHMemorySpace::type;
  // using value_type = T;
  //  using coordinate_element_type = CoordinateElement<coordinate_system,
  //  Real>;
  // using coordinate_element_type = CoordinateElementType;

  OmegaHField(std::string name, Omega_h::Mesh& mesh)
    : name_(std::move(name)), mesh_(mesh)
  {
  }
  [[nodiscard]] const std::string& GetName() const noexcept { return name_; }
  [[nodiscard]] Omega_h::Mesh& GetMesh() const noexcept { return mesh_; }

private:
  std::string name_;
  Omega_h::Mesh& mesh_;
};

template <typename T, typename CoordinateElementType>
auto get_nodal_data(const OmegaHField<T, CoordinateElementType>& field)
  -> Omega_h::Read<T>
{
  return field.GetMesh().template get_array<T>(0, field.GetName());
}

// TODO since Omega_h owns coordinate data, we could potentially
// return a view of the data without lifetime issues.
template <typename T, typename CoordinateElementType>
auto get_nodal_coordinates(const OmegaHField<T, CoordinateElementType>& field)
{
  if constexpr (detail::HasCoordinateSystem<CoordinateElementType>::value) {
    const auto coords = field.GetMesh().coords();
    return MDArray<CoordinateElementType>{};
    // FIXME implement copy to
    throw;
  } else {
    return Omega_h::Reals{field.GetMesh().coords()};
  }
}

/**
 * Sets the data on the entire mesh
 */
template <typename T, typename CoordinateElementType, typename U>
auto set(const OmegaHField<T, CoordinateElementType>& field,
         ScalarArrayView<const U, OmegaHMemorySpace::type> data) -> void
{

  WDMCPL_ALWAYS_ASSERT(data.extent(0) == field.GetMesh().nents(0));
  auto& mesh = field.GetMesh();
  Omega_h::Write<U> array(data.extent(0), 0);
  Omega_h::parallel_for(
    data.extent(0), OMEGA_H_LAMBDA(size_t i) { array[i] = data(i); });
  if (mesh.has_tag(0, field.GetName())) {
    mesh.set_tag(0, field.GetName(), Omega_h::Read(array));
  } else {
    mesh.add_tag(0, field.GetName(), 1, Omega_h::Read(array));
  }
}

template <typename T, typename CoordinateElementType, int o>
auto evaluate(
  const OmegaHField<T, CoordinateElementType>& field, Lagrange<o> method,
  ScalarArrayView<const CoordinateElementType, OmegaHMemorySpace::type>
    coordinates) -> Omega_h::Reals
{
  // FIXME implement
}

} // namespace wdmcpl
namespace Omega_h
{
template <typename T>
auto make_array_view(const Omega_h::Read<T>& array)
  -> wdmcpl::ScalarArrayView<const T, typename wdmcpl::OmegaHMemorySpace::type>
{
  wdmcpl::ScalarArrayView<const T, typename wdmcpl::OmegaHMemorySpace::type>
    view(array.data(), array.size());
  return view;
}
} // namespace Omega_h

#endif // WDM_COUPLING_OMEGA_H_FIELD_H
