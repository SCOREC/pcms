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

// FIXME add executtion spaces (don't use kokkos exe spaces directly)

namespace wdmcpl
{

struct Omega_h_field
{
  using NodeHandleType = Omega_h::GO;
  using NodeHandleArray = Omega_h::GOs;
  using CoordinateSystem = Cartesian;
  using ExecutionSpace = Kokkos::DefaultExecutionSpace;
  using DataType = wdmcpl::Real;
  using ArrayType = Omega_h::Reals;

  using CoordinateArrayType =
    CoordinateArray<CoordinateSystem, ArrayType, DataType, ExecutionSpace>;
  using ScalarArrayType =
    ScalarArray<ArrayType, DataType, ExecutionSpace>;

  std::string name;
  Omega_h::Mesh& mesh;
};

namespace detail {
  constexpr void check_omega_h_field() {
    check_field<Omega_h_field>();
  }
}

template <>
nonstd::span<const Omega_h::Real> GetSpan(const Omega_h::Reals& array)
{
  return {array.data(), static_cast<size_t>(array.size())};
}


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
void set(Omega_h_field & field, const typename Omega_h_field::ScalarArrayType & data) {
  const auto data_view = data.GetSpan();
  WDMCPL_ALWAYS_ASSERT(data_view.size() == field.mesh.nents(0));
  auto& mesh = field.mesh;
  Omega_h::Write array(data_view.size(),0);
  Omega_h::parallel_for(data_view.size(),OMEGA_H_LAMBDA(size_t i){
                                            array[i] = data_view[i];
                                          });
  mesh.set_tag(0,field.name,Omega_h::Read(array));
}


} // namespace wdmcpl

#endif // WDM_COUPLING_OMEGA_H_FIELD_H
