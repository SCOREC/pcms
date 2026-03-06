#include "uniform_grid_field.h"
#include "pcms/utility/profile.h"
#include "pcms/interpolator/linear_interpolant.hpp"
#include "pcms/interpolator/multidimarray.hpp"
#include "pcms/utility/assert.h"
#include <memory>

namespace pcms
{

// Localization hint structure for UniformGrid
template <unsigned Dim>
struct UniformGridFieldLocalizationHint
{
  UniformGridFieldLocalizationHint(
    Kokkos::View<LO*, HostMemorySpace> cell_indices,
    Kokkos::View<Real**, HostMemorySpace> coordinates, OutOfBoundsMode mode,
    Kokkos::View<bool*, HostMemorySpace> is_out_of_bounds,
    size_t num_out_of_bounds)
    : cell_indices_(cell_indices),
      coordinates_(coordinates),
      mode_(mode),
      is_out_of_bounds_(is_out_of_bounds),
      num_out_of_bounds_(num_out_of_bounds)
  {
  }

  // Cell indices for each point as 1d array
  Kokkos::View<LO*, HostMemorySpace> cell_indices_;
  // Coordinates of each point
  Kokkos::View<Real**, HostMemorySpace> coordinates_;
  // Out of bounds handling
  OutOfBoundsMode mode_;
  Kokkos::View<bool*, HostMemorySpace> is_out_of_bounds_;
  size_t num_out_of_bounds_;
};

template <unsigned Dim>
UniformGridField<Dim>::UniformGridField(
  const UniformGridFieldLayout<Dim>& layout)
  : layout_(layout),
    grid_(layout.GetGrid()),
    dof_holder_data_("dof_holder_data", static_cast<size_t>(layout.OwnedSize()))
{
  PCMS_FUNCTION_TIMER;
  // Default to NEAREST_BOUNDARY for uniform grid fields
  out_of_bounds_mode_ = OutOfBoundsMode::NEAREST_BOUNDARY;
}

template <unsigned Dim>
Rank1View<const Real, HostMemorySpace> UniformGridField<Dim>::GetDOFHolderData()
  const
{
  PCMS_FUNCTION_TIMER;
  return make_const_array_view(dof_holder_data_);
}

template <unsigned Dim>
void UniformGridField<Dim>::SetDOFHolderData(
  Rank1View<const Real, HostMemorySpace> data)
{
  PCMS_FUNCTION_TIMER;
  PCMS_ALWAYS_ASSERT(data.size() == dof_holder_data_.size());

  for (size_t i = 0; i < data.size(); ++i) {
    dof_holder_data_[i] = data[i];
  }
}

template <unsigned Dim>
View<Dim, Real, HostMemorySpace> UniformGridField<Dim>::to_mdspan()
{
  PCMS_FUNCTION_TIMER;

  if constexpr (Dim == 2) {
    return View<Dim, Real, HostMemorySpace>(
      dof_holder_data_.data(), grid_.divisions[0] + 1, grid_.divisions[1] + 1);
  } else if constexpr (Dim == 3) {
    return View<Dim, Real, HostMemorySpace>(
      dof_holder_data_.data(), grid_.divisions[0] + 1, grid_.divisions[1] + 1,
      grid_.divisions[2] + 1);
  } else {
    static_assert(Dim == 2 || Dim == 3,
                  "to_mdspan only supports 2D or 3D uniform grids");
  }
}

template <unsigned Dim>
View<Dim, const Real, HostMemorySpace> UniformGridField<Dim>::to_mdspan() const
{
  PCMS_FUNCTION_TIMER;

  if constexpr (Dim == 2) {
    return View<Dim, const Real, HostMemorySpace>(
      dof_holder_data_.data(), grid_.divisions[0] + 1, grid_.divisions[1] + 1);
  } else if constexpr (Dim == 3) {
    return View<Dim, const Real, HostMemorySpace>(
      dof_holder_data_.data(), grid_.divisions[0] + 1, grid_.divisions[1] + 1,
      grid_.divisions[2] + 1);
  } else {
    static_assert(Dim == 2 || Dim == 3,
                  "to_mdspan only supports 2D or 3D uniform grids");
  }
}

template <unsigned Dim>
LocalizationHint UniformGridField<Dim>::GetLocalizationHint(
  CoordinateView<HostMemorySpace> coordinate_view) const
{
  PCMS_FUNCTION_TIMER;

  if (coordinate_view.GetCoordinateSystem() !=
      layout_.GetDOFHolderCoordinates().GetCoordinateSystem()) {
    throw std::runtime_error(
      "Coordinate system mismatch in GetLocalizationHint");
  }

  auto coordinates = coordinate_view.GetCoordinates();
  LO num_points = coordinates.extent(0);

  Kokkos::View<LO*, HostMemorySpace> cell_indices("cell_indices", num_points);
  Kokkos::View<Real**, HostMemorySpace> coords_copy("coords_copy", num_points,
                                                    Dim);
  Kokkos::View<bool*, HostMemorySpace> is_out_of_bounds("is_out_of_bounds",
                                                        num_points);
  size_t num_out_of_bounds = 0;

  // Find which cell each point belongs to and detect out-of-bounds
  for (LO i = 0; i < num_points; ++i) {
    Omega_h::Vector<Dim> point;
    for (unsigned d = 0; d < Dim; ++d) {
      point[d] = coordinates(i, d);
      coords_copy(i, d) = coordinates(i, d);
    }

    // Check if point is within grid bounds
    bool out_of_bounds = !grid_.IsPointInBounds(point);

    is_out_of_bounds[i] = out_of_bounds;
    if (out_of_bounds) {
      num_out_of_bounds++;
      if (out_of_bounds_mode_ == OutOfBoundsMode::ERROR) {
        PCMS_ALWAYS_ASSERT(false && "Point found outside uniform grid domain");
      }
    }

    cell_indices[i] = grid_.ClosestCellID(point);
  }

  auto hint = std::make_shared<UniformGridFieldLocalizationHint<Dim>>(
    cell_indices, coords_copy, out_of_bounds_mode_, is_out_of_bounds,
    num_out_of_bounds);
  return LocalizationHint{hint};
}

template <unsigned Dim>
void UniformGridField<Dim>::Evaluate(
  LocalizationHint location, FieldDataView<Real, HostMemorySpace> results) const
{
  PCMS_FUNCTION_TIMER;

  if (results.GetCoordinateSystem() !=
      layout_.GetDOFHolderCoordinates().GetCoordinateSystem()) {
    throw std::runtime_error("Coordinate system mismatch in Evaluate");
  }

  UniformGridFieldLocalizationHint<Dim> hint =
    *reinterpret_cast<UniformGridFieldLocalizationHint<Dim>*>(
      location.data.get());

  auto values = dof_holder_data_;
  auto coordinates = hint.coordinates_;
  auto cell_indices = hint.cell_indices_;
  LO num_points = coordinates.extent(0);

  Kokkos::View<LO**, HostMemorySpace> cell_dimensioned_indices(
    "cell_dimensioned_indices", num_points, Dim);
  // Convert cell indices to multi-dimensional indices
  for (LO i = 0; i < num_points; ++i) {
    auto dim_indices = grid_.GetDimensionedIndex(cell_indices[i]);
    for (unsigned d = 0; d < Dim; ++d) {
      cell_dimensioned_indices(i, d) = dim_indices[d];
    }
  }

  // Dimensions for vertex grid (m+1 vertices per dimension for m cells)
  auto cell_divisions = layout_.GetGrid().divisions;
  IntVecView dimensions_view("dimensions", Dim);
  auto dimensions_view_host = Kokkos::create_mirror_view(dimensions_view);
  for (unsigned d = 0; d < Dim; ++d) {
    dimensions_view_host(d) = cell_divisions[d] + 1;
  }
  Kokkos::deep_copy(dimensions_view, dimensions_view_host);

  // Compute parametric coordinates for each point
  RealMatView parametric_coords("parametric_coords", num_points, Dim);
  auto parametric_coords_host = Kokkos::create_mirror_view(parametric_coords);

  for (LO i = 0; i < num_points; ++i) {
    auto cell_bbox = grid_.GetCellBBOX(cell_indices[i]);
    for (unsigned d = 0; d < Dim; ++d) {
      Real coord = coordinates(i, d);
      Real cell_min = cell_bbox.center[d] - cell_bbox.half_width[d];
      Real cell_max = cell_bbox.center[d] + cell_bbox.half_width[d];
      parametric_coords_host(i, d) = (coord - cell_min) / (cell_max - cell_min);
    }
  }

  Kokkos::deep_copy(parametric_coords, parametric_coords_host);

  // Convert cell indices and values to the required view types
  IntMatView cell_indices_interp("cell_indices_interp", num_points, Dim);
  auto cell_indices_interp_host =
    Kokkos::create_mirror_view(cell_indices_interp);
  for (LO i = 0; i < num_points; ++i) {
    for (unsigned d = 0; d < Dim; ++d) {
      cell_indices_interp_host(i, d) = cell_dimensioned_indices(i, d);
    }
  }
  Kokkos::deep_copy(cell_indices_interp, cell_indices_interp_host);

  RealVecView values_interp("values_interp", values.extent(0));
  Kokkos::deep_copy(values_interp, values);

  auto interpolator = RegularGridInterpolator(
    parametric_coords, values_interp, cell_indices_interp, dimensions_view);
  auto evaluated_values = results.GetValues();
  auto interpolated_values = interpolator.linear_interpolation();
  auto interpolated_values_host =
    Kokkos::create_mirror_view(interpolated_values);
  Kokkos::deep_copy(interpolated_values_host, interpolated_values);

  // Copy interpolated values and handle out-of-bounds points
  for (LO i = 0; i < evaluated_values.size(); ++i) {
    if (hint.is_out_of_bounds_[i] && hint.mode_ == OutOfBoundsMode::FILL) {
      // Fill out-of-bounds points with fill value
      evaluated_values[i] = fill_value_;
    } else {
      // Use interpolated value (for in-bounds or NEAREST_BOUNDARY mode)
      evaluated_values[i] = interpolated_values_host[i];
    }
  }
}

template <unsigned Dim>
void UniformGridField<Dim>::EvaluateGradient(
  FieldDataView<Real, HostMemorySpace>)
{
  throw std::runtime_error("Not implemented");
}

template <unsigned Dim>
const FieldLayout& UniformGridField<Dim>::GetLayout() const
{
  return layout_;
}

template <unsigned Dim>
bool UniformGridField<Dim>::CanEvaluateGradient()
{
  return false;
}

template <unsigned Dim>
int UniformGridField<Dim>::Serialize(
  Rank1View<Real, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const
{
  PCMS_FUNCTION_TIMER;

  const auto array_h = GetDOFHolderData();
  if (buffer.size() > 0) {
    PCMS_ALWAYS_ASSERT(buffer.size() == array_h.size());
    for (LO i = 0; i < array_h.size(); ++i) {
      buffer[permutation[i]] = array_h[i];
    }
  }
  return array_h.size();
}

template <unsigned Dim>
void UniformGridField<Dim>::Deserialize(
  Rank1View<const Real, pcms::HostMemorySpace> buffer,
  Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation)
{
  PCMS_FUNCTION_TIMER;

  Kokkos::View<Real*, HostMemorySpace> sorted_buffer("sorted_buffer",
                                                     permutation.size());
  auto owned = layout_.GetOwned();

  for (LO i = 0; i < sorted_buffer.size(); ++i) {
    PCMS_ALWAYS_ASSERT(owned[i]);
    sorted_buffer[i] = buffer[permutation[i]];
  }

  SetDOFHolderData(pcms::make_const_array_view(sorted_buffer));
}

// Explicit template instantiations
template class UniformGridField<2>;
template class UniformGridField<3>;

} // namespace pcms
