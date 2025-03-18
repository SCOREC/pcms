#include "omega_h_field2.h"
#include "omega_h_field.h"
#include "pcms/adapter/omega_h/omega_h_field_layout.h"
#include "pcms/adapter/omega_h/omega_h_field2.h"

namespace pcms
{
/*
 * Field
 */
OmegaHField2::OmegaHField2(std::string name, CoordinateSystem coordinate_system,
                         const OmegaHFieldLayout& layout, Omega_h::Mesh& mesh)
  : name_(name), coordinate_system_(coordinate_system), layout_(layout),
    mesh_(mesh)
{
}

std::string& OmegaHField2::GetName() const
{
  return name_;
}

mesh_entity_type OmegaHField2::GetEntityType() const
{
  return mesh_entity_type::VERTEX;
}

CoordinateSystem OmegaHField2::GetCoordinateSystem() const
{
  return coordinate_system_;
}

Kokkos::View<const Real*> OmegaHField2::GetNodalData() const
{
  return mesh_.template get_array<Real>(0, name_).view();
}

void OmegaHField2::SetNodalData(Kokkos::View<const Real*> data)
{
  const auto has_tag = mesh_.has_tag(0, name_);
  PCMS_ALWAYS_ASSERT(static_cast<LO>(data.size()) == mesh_.nents(0));
  Omega_h::Write<T> array(data.size());
  Omega_h::parallel_for(
    data.size(), OMEGA_H_LAMBDA(size_t i) { array[i] = data(i); });
  if (has_tag) {
    mesh_.set_tag(0, name_, Omega_h::Read<T>(array));
  } else {
    mesh_.add_tag(0, name_, 1, Omega_h::Read<T>(array));
  }
}

void OmegaHField2::SetEvaluationCoordinates(LocalizationHint hint)
{
  hint_ = (Kokkos::View<GridPointSearch::Result*> *) hint.data;
}

LocalizationHint OmegaHField2::GetLocalizationHint(
  CoordinateView<HostMemorySpace> coordinate_view)
{
  // TODO decide if we want to implicitly perform the coordinate transformations
  // when possible
  if (coordinates.GetCoordinateSystem() != coordinate_system_) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Coordinate system mismatch");
  }

  auto coordinates = coordinate_view.GetCoordinates();
  Kokkos::View<Real* [2]> coords("coords", coordinates.size() / 2);
  Kokkos::parallel_for(
    coordinates.size() / 2, KOKKOS_LAMBDA(LO i) {
      coords(i, 0) = coordinates(2 * i);
      coords(i, 1) = coordinates(2 * i + 1);
    });
  auto results = field.Search(coords);

  return LocalizationHint { results };
}

void OmegaHField2::Evaluate(FieldDataView<double, HostMemorySpace> results)
{
  // TODO decide if we want to implicitly perform the coordinate transformations
  // when possible
  if (results.GetCoordinateSystem() != coordinate_system_) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Coordinate system mismatch");
  }

  if (!hint_) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Evaluation coordinates not set");
  }

  auto tris2verts = mesh_.ask_elem_verts();
  auto field_values = mesh_.template get_array<Real>(0, name_);

  auto search_results = *hint_;
  Omega_h::Write<T> values(search_results.size());

  switch (layout_.location) {
    case OmegaHFieldLayoutLocation::Linear:
      Kokkos::parallel_for(
        search_results.size(), KOKKOS_LAMBDA(LO i) {
          auto [dim, elem_idx, coord] = search_results(i);
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
      break;

    case OmegaHFieldLayoutLocation::Piecewise:
      Kokkos::parallel_for(
        results.size(), KOKKOS_LAMBDA(LO i) {
          auto [dim, elem_idx, coord] = results(i);
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
      break;
  }

  // if field is piecewise
  // 1. use localization to get the appropriate element
  // 2. get the face/region value in the element
  // OPEN Question: what to do if point is on the interface between two regions?
  // a) return mean?, b) return 2 values?, c) pick one at random

  // if the field is linear:
  // 1. use localization to get the appropriate element
  // 2. get vertexes of the element (e2v map)
  // 3. get the values of the vertexes
  // 4. get barycentric coordinates of the point in the element
  // 5. interpolate the values of the vertexes to the point
}

void OmegaHField2::EvaluateGradient(FieldDataView<double, HostMemorySpace> results)
{
  // TODO when moved to PCMS throw PCMS exception
  throw std::runtime_error("Not implemented");
}

const FieldLayout& OmegaHField2::GetLayout()
{
  return layout_;
}

bool OmegaHField2::CanEvaluateGradient()
{
  // TODO compute the gradient field using element shape functions
  return false;
}

int OmegaHField2::Serialize(
  Rank1View<double, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const
{
  return 1;
}

void OmegaHField2::Deserialize(
  Rank1View<const double, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const
{
}

void OmegaHField2::Copy(OmegaHField2& target) const
{
  auto own_field = mesh_.template get_array<T>(0, name_);
  auto mesh = target.mesh_;
  std::string& name = target.name_;

  const auto has_tag = mesh.has_tag(0, name);
  PCMS_ALWAYS_ASSERT(static_cast<LO>(own_field.size()) == mesh.nents(0));
  Omega_h::Write<T> array(own_field.size());
  Omega_h::parallel_for(
    own_field.size(), OMEGA_H_LAMBDA(size_t i) { array[i] = own_field(i); });
  if (has_tag) {
    mesh.set_tag(0, name, Omega_h::Read<T>(array));
  } else {
    mesh.add_tag(0, name, 1, Omega_h::Read<T>(array));
  }
}

void OmegaHField2::Interpolate(OmegaHField2& target) const {

}

} // namespace pcms
