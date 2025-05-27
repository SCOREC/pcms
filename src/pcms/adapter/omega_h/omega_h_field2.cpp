#include "omega_h_field2.h"
#include "omega_h_field.h"
#include "pcms/adapter/omega_h/omega_h_field2.h"
#include <memory>

namespace pcms
{
OmegaHFIeld2LocalizationHint::OmegaHFIeld2LocalizationHint() {}

OmegaHFIeld2LocalizationHint::OmegaHFIeld2LocalizationHint(
  Kokkos::View<GridPointSearch::Result*> search_results)
  : search_results_(search_results)
{
}

Kokkos::View<GridPointSearch::Result*> OmegaHFIeld2LocalizationHint::GetResults()
{
  return search_results_;
}

/*
 * Field
 */
OmegaHField2::OmegaHField2(std::string name, CoordinateSystem coordinate_system,
                         const OmegaHFieldLayout& layout, Omega_h::Mesh& mesh)
  : name_(name), coordinate_system_(coordinate_system), layout_(layout),
    mesh_(mesh), search_(mesh, 10, 10)
{
}

const std::string& OmegaHField2::GetName() const
{
  return name_;
}

CoordinateSystem OmegaHField2::GetCoordinateSystem() const
{
  return coordinate_system_;
}

FieldDataView<const Real, HostMemorySpace> OmegaHField2::GetDOFHolderData() const
{
  PCMS_FUNCTION_TIMER;
  auto array = mesh_.template get_array<Real>(0, name_);
  Rank1View<const Real, pcms::HostMemorySpace> array_view{std::data(array), std::size(array)};
  FieldDataView<const Real, HostMemorySpace> data_view{array_view, GetCoordinateSystem()};
  return data_view;
};

CoordinateView<HostMemorySpace> OmegaHField2::GetDOFHolderCoordinates() const
{
  PCMS_FUNCTION_TIMER;
  auto coords = Omega_h::get_ent_centroids(mesh_, 0);
  Rank2View<const double, HostMemorySpace> coords_view(coords.data(), coords.size() / 2, 2);
  return CoordinateView<HostMemorySpace>{GetCoordinateSystem(), coords_view};
}

void OmegaHField2::SetDOFHolderData(FieldDataView<const Real, HostMemorySpace> data) {
  PCMS_FUNCTION_TIMER;
  if (data.GetCoordinateSystem() != coordinate_system_) {
    throw std::runtime_error("Coordinate system mismatch");
  }
  const auto has_tag = mesh_.has_tag(0, name_);
  Rank1View<const double, HostMemorySpace> values = data.GetValues();
  PCMS_ALWAYS_ASSERT(static_cast<LO>(values.size()) == mesh_.nents(0));
  Omega_h::Write<Real> array(values.size());
  Omega_h::parallel_for(
    values.size(), OMEGA_H_LAMBDA(size_t i) { array[i] = values[i]; });
  if (has_tag) {
    mesh_.set_tag(0, name_, Omega_h::Read<Real>(array));
  } else {
    mesh_.add_tag(0, name_, 1, Omega_h::Read<Real>(array));
  }
}

LocalizationHint OmegaHField2::GetLocalizationHint(
  CoordinateView<HostMemorySpace> coordinate_view) const
{
  PCMS_FUNCTION_TIMER;
  // TODO decide if we want to implicitly perform the coordinate transformations
  // when possible
  if (coordinate_view.GetCoordinateSystem() != coordinate_system_) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Coordinate system mismatch");
  }

  auto coordinates = coordinate_view.GetCoordinates();
  Kokkos::View<Real* [2]> coords("coords", coordinates.size() / 2);
  Kokkos::parallel_for(
    coordinates.size() / 2, KOKKOS_LAMBDA(LO i) {
      coords(i, 0) = coordinates(i, 0);
      coords(i, 1) = coordinates(i, 1);
    });
  auto results = search_(coords);
  auto hint = std::make_shared<OmegaHFIeld2LocalizationHint>(results);

  return LocalizationHint{hint};
}

void OmegaHField2::Evaluate(LocalizationHint location,
                            FieldDataView<double, HostMemorySpace> results) const
{
  PCMS_FUNCTION_TIMER;
  // TODO decide if we want to implicitly perform the coordinate transformations
  // when possible
  if (results.GetCoordinateSystem() != coordinate_system_) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Coordinate system mismatch");
  }

  OmegaHFIeld2LocalizationHint hint =
    *((OmegaHFIeld2LocalizationHint*)location.data.get());

  auto search_results = hint.GetResults();

  if (!search_results.is_allocated()) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Evaluation coordinates not set");
  }

  auto tris2verts = mesh_.ask_elem_verts();
  auto field_values = mesh_.template get_array<Real>(0, name_);

  Rank1View<double, HostMemorySpace> values = results.GetValues();

  switch (layout_.GetLocation()) {
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
          values[i] = val;
        });
      break;

    case OmegaHFieldLayoutLocation::PieceWise:
      Kokkos::parallel_for(
        search_results.size(), KOKKOS_LAMBDA(LO i) {
          auto [dim, elem_idx, coord] = search_results(i);
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
}

void OmegaHField2::EvaluateGradient(FieldDataView<double, HostMemorySpace> results)
{
  // TODO when moved to PCMS throw PCMS exception
  throw std::runtime_error("Not implemented");
}

const FieldLayout& OmegaHField2::GetLayout() const
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
  PCMS_FUNCTION_TIMER;
  // host copy of filtered field data array
  const auto array_h = GetDOFHolderData().GetValues();
  if (buffer.size() > 0) {
    for (LO i = 0; i < array_h.size(); i++) {
      buffer[i] = array_h[permutation[i]];
    }
  }
  return array_h.size();
}

void OmegaHField2::Deserialize(
  Rank1View<const double, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
{
  PCMS_FUNCTION_TIMER;
  REDEV_ALWAYS_ASSERT(buffer.size() == permutation.size());
  Omega_h::HostWrite<Real> sorted_buffer(buffer.size());
  for (size_t i = 0; i < buffer.size(); ++i) {
    sorted_buffer[permutation[i]] = buffer[i];
  }
  const auto sorted_buffer_d = Omega_h::Read<Real>(sorted_buffer);
  Rank1View<const Real, HostMemorySpace> sorted_buffer_view{std::data(sorted_buffer_d), std::size(sorted_buffer_d)};
  FieldDataView<const Real, HostMemorySpace> field_view{sorted_buffer_view, GetCoordinateSystem()};
  SetDOFHolderData(field_view);
}
} // namespace pcms
