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

class OmegaHField
{
public:
  using node_handle_type = Omega_h::GO;
  using node_handle_array = Omega_h::GOs;
  using coordinate_system = Cartesian;
  using memory_space = OmegaHMemorySpace::type;
  using value_type = LO;
  using coordinate_element_type = CoordinateElement<coordinate_system, Real>;

  OmegaHField(std::string name, Omega_h::Mesh& mesh)
    : name_(std::move(name)), mesh_(mesh)
  {
  }
  // using DataType = wdmcpl::Real;
  // using ArrayType = Omega_h::Reals;

  // using CoordinateArrayType =
  //   CoordinateArray<CoordinateSystem, ArrayType, DataType, ExecutionSpace>;
  // using ScalarArrayType =
  //   ScalarArray<ArrayType, ExecutionSpace>;
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
struct ArrayToViewTypesMap<Omega_h::Read<T>> {
  using memory_space = typename OmegaHMemorySpace::type;
  using value_type = const T;
};

// template <typename T>
template <>
ScalarArrayView<const Real, typename OmegaHMemorySpace::type>
make_array_view(const Omega_h::Read<Real>& array)
{
  ScalarArrayView<const Real, typename OmegaHMemorySpace::type> view(array.data(),array.size());
  return view;
}
template <>
ScalarArrayView<const LO, typename OmegaHMemorySpace::type>
make_array_view(const Omega_h::Read<LO>& array)
{
  ScalarArrayView<const LO, typename OmegaHMemorySpace::type> view(array.data(),array.size());
  return view;
}

auto get_nodal_data(const OmegaHField& field)
  -> Omega_h::Read<OmegaHField::value_type>
{
  return field.GetMesh().get_array<OmegaHField::value_type>(0, field.GetName());
}

/**
 * Sets the data on the entire mesh
 */
template <typename ElementType>
auto set(const OmegaHField& field,
         ScalarArrayView<ElementType, OmegaHMemorySpace::type> data) -> void
{

  WDMCPL_ALWAYS_ASSERT(data.extent(0) == field.GetMesh().nents(0));
  auto& mesh = field.GetMesh();
  Omega_h::Write<std::remove_cv_t<ElementType>> array(data.extent(0), 0);
  Omega_h::parallel_for(
    data.extent(0), OMEGA_H_LAMBDA(size_t i) { array[i] = data(i); });
  if(mesh.has_tag(0, field.GetName())) {
    mesh.set_tag(0, field.GetName(), Omega_h::Read(array));
  }
  else {
    mesh.add_tag(0, field.GetName(), 1, Omega_h::Read(array));
  }
}

} // namespace wdmcpl

#endif // WDM_COUPLING_OMEGA_H_FIELD_H
