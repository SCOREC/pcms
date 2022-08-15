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

// FIXME add executtion spaces (don't use kokkos exe spaces directly)

namespace wdmcpl
{

// TODO different types dependent on active OmegaHBackend
struct OmegaHMemorySpace
{
  using type = Kokkos::DefaultExecutionSpace;
};

template <typename T>
class OmegaHField
{
public:
  using node_handle_type = Omega_h::GO;
  using node_handle_array = Omega_h::GOs;
  using coordinate_system = Cartesian;
  using memory_space = OmegaHMemorySpace::type;
  using value_type = T;
  using coordinate_element_type = CoordinateElement<coordinate_system, Real>;

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

/*
namespace detail {
  constexpr void check_omega_h_field() {
    check_field<Omega_h_field>();
  }
}
*/

// template <>
// nonstd::span<const Omega_h::Real> GetSpan(const Omega_h::Reals& array)
//{
//   return {array.data(), static_cast<size_t>(array.size())};
// }

/*
template <>
typename Omega_h_field::NodeHandleArray get_node_handles(
  const Omega_h_field& field)
{
  return field.mesh.globals(0);
}

template <>
Omega_h_field::CoordinateArrayType get_nodal_coordinates(
  const Omega_h_field& field)
{
  return {field.mesh.coords()};
}

// set the data on the entire mesh
template <>
void set(Omega_h_field & field,
         ScalarArrayView<const typename Omega_h_field::DataType,
           typename Omega_h_field::ExecutionSpace> data)
 {
  WDMCPL_ALWAYS_ASSERT(data.extent(0) == field.mesh.nents(0));
  auto& mesh = field.mesh;
  Omega_h::Write array(data.extent(0),0);
  Omega_h::parallel_for(data.extent(0),OMEGA_H_LAMBDA(size_t i){
                                            array[i] = data[i];
                                          });
  mesh.set_tag(0,field.name,Omega_h::Read(array));
}
*/

template <typename T>
auto get_nodal_data(const OmegaHField<T>& field)
  -> Omega_h::Read<T>
{
  return field.GetMesh().template get_array<T>(0, field.GetName());
}

/**
 * Sets the data on the entire mesh
 */
template <typename ElementType>
auto set(const OmegaHField<ElementType>& field,
         ScalarArrayView<const ElementType, OmegaHMemorySpace::type> data) -> void
{

  WDMCPL_ALWAYS_ASSERT(data.extent(0) == field.GetMesh().nents(0));
  auto& mesh = field.GetMesh();
  Omega_h::Write<ElementType> array(data.extent(0), 0);
  Omega_h::parallel_for(
    data.extent(0), OMEGA_H_LAMBDA(size_t i) { array[i] = data(i); });
  if (mesh.has_tag(0, field.GetName())) {
    mesh.set_tag(0, field.GetName(), Omega_h::Read(array));
  } else {
    mesh.add_tag(0, field.GetName(), 1, Omega_h::Read(array));
  }
}

} // namespace wdmcpl
namespace Omega_h
{
template <typename T>
wdmcpl::ScalarArrayView<const T, typename wdmcpl::OmegaHMemorySpace::type>
make_array_view(const Omega_h::Read<T>& array)
{
  wdmcpl::ScalarArrayView<const T, typename wdmcpl::OmegaHMemorySpace::type>
    view(array.data(), array.size());
  return view;
}
} // namespace Omega_h

#endif // WDM_COUPLING_OMEGA_H_FIELD_H
