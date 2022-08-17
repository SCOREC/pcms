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
#include "wdmcpl/point_search.h"

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
  OmegaHField(std::string name, Omega_h::Mesh& mesh, int search_nx = 10,
              int search_ny = 10)
    : name_(std::move(name)), mesh_(mesh), search_{mesh, search_nx, search_ny}
  {
  }
  [[nodiscard]] const std::string& GetName() const noexcept { return name_; }
  [[nodiscard]] Omega_h::Mesh& GetMesh() const noexcept { return mesh_; }
  // pass through to search function
  template <typename... Ts>
  auto Search(Ts... args) const
  {
    return search_(std::forward<Ts>(args)...);
  }

private:
  std::string name_;
  Omega_h::Mesh& mesh_;
  // FIXME GridPointSearch take mdspan
  GridPointSearch search_;
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

// TODO abstract out repeat parts of lagrange/nearest neighbor evaluation
template <typename T, typename CoordinateElementType>
auto evaluate(
  const OmegaHField<T, CoordinateElementType>& field, Lagrange<1> method,
  ScalarArrayView<const CoordinateElementType, OmegaHMemorySpace::type>
    coordinates) -> Omega_h::Read<T>
{
  Omega_h::Write<T> values(coordinates.size() / 2);
  auto tris2verts = field.GetMesh().ask_elem_verts();
  auto field_values = field.GetMesh().template get_array<T>(0, field.GetName());
  // FIXME field search operation needs to be on GPU
  // FIXME field search takes mdspan
  // FIXME coordinate array view as nx2 span
  Omega_h::parallel_for(
    coordinates.size() / 2, OMEGA_H_LAMBDA(LO i) {
      auto [elem_idx, coord] = field.Search(
        Omega_h::Vector<2>{coordinates(2 * i), coordinates(2 * i + 1)});
      // TODO deal with case for elem_idx < 0 (point outside of mesh)
      KOKKOS_ASSERT(elem_idx >= 0);
      const auto elem_tri2verts =
        Omega_h::gather_verts<3>(tris2verts, elem_idx);
      Real val = 0;
      for (int j = 0; j < 3; ++j) {
        val += field_values[elem_tri2verts[j]] * coord[j];
      }
      if constexpr (std::is_integral_v<T>) {
        val = std::round(val);
      }
      values[i] = val;
    });
  return values;
}

template <typename T, typename CoordinateElementType>
auto evaluate(
  const OmegaHField<T, CoordinateElementType>& field, NearestNeighbor method,
  ScalarArrayView<const CoordinateElementType, OmegaHMemorySpace::type>
    coordinates) -> Omega_h::Read<T>
{
  Omega_h::Write<T> values(coordinates.size() / 2);
  auto tris2verts = field.GetMesh().ask_elem_verts();
  auto field_values = field.GetMesh().template get_array<T>(0, field.GetName());
  // FIXME field search operation needs to be on GPU
  // FIXME field search takes mdspan
  // FIXME coordinate array view as nx2 span
  Omega_h::parallel_for(
    coordinates.size() / 2, OMEGA_H_LAMBDA(LO i) {
      auto [elem_idx, coord] = field.Search(
        Omega_h::Vector<2>{coordinates(2 * i), coordinates(2 * i + 1)});
      // TODO deal with case for elem_idx < 0 (point outside of mesh)
      KOKKOS_ASSERT(elem_idx >= 0);
      const auto elem_tri2verts =
        Omega_h::gather_verts<3>(tris2verts, elem_idx);
      // value is closest to point has the largest coordinate
      int vert = 0;
      auto max_val = coord[0];
      for (int j = 1; j <= 2; ++j) {
        auto next_val = coord[j];
        if (next_val > max_val) {
          max_val = next_val;
          vert = j;
        }
      }
      values[i] = field_values[elem_tri2verts[vert]];
    });
  return values;
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
