//
// Created by hasanm4 on 1/17/25.
//

#include "interpolation_base.h"
#include "interpolation_helpers.h"

#include <execution>

Omega_h::Reals getCentroids(Omega_h::Mesh& mesh)
{
  OMEGA_H_CHECK_PRINTF(
    mesh.dim() == 2, "Only 2D meshes are supported but found %d\n", mesh.dim());

  const auto& coords = mesh.coords();
  Omega_h::Write<Omega_h::Real> centroids(mesh.nfaces() * mesh.dim(), 0.0);

  auto face2node = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  Omega_h::parallel_for(
    mesh.nfaces(), OMEGA_H_LAMBDA(Omega_h::LO face) {
      auto nodes = Omega_h::gather_verts<3>(face2node, face);
      Omega_h::Few<Omega_h::Vector<2>, 3> face_coords =
        Omega_h::gather_vectors<3, 2>(coords, nodes);
      Omega_h::Vector<2> centroid = Omega_h::average(face_coords);
      centroids[2 * face + 0] = centroid[0];
      centroids[2 * face + 1] = centroid[1];
    });

  return {centroids};
}

MLSInterpolationHandler::MLSInterpolationHandler(Omega_h::Mesh& source_mesh,
                                                 double radius,
                                                 uint min_req_support,
                                                 uint degree, bool adapt_radius)
  : source_mesh_(source_mesh),
    target_mesh_(source_mesh),
    radius_(radius),
    min_req_supports_(min_req_support),
    degree_(degree),
    adapt_radius_(adapt_radius)
{
  single_mesh_ = true;
  target_coords_ = source_mesh_.coords();
  source_coords_ = getCentroids(source_mesh_);

  OMEGA_H_CHECK_PRINTF(source_mesh_.dim() == 2,
                       "Only 2D meshes are supported but found %d\n",
                       source_mesh_.dim());

  source_field_ =
    Omega_h::HostWrite<Omega_h::Real>(source_mesh_.nfaces(), "source field");

  target_field_ =
    Omega_h::HostWrite<Omega_h::Real>(source_mesh_.nverts(), "target field");

  find_supports(min_req_supports_);
}

MLSInterpolationHandler::MLSInterpolationHandler(
  Omega_h::Mesh& source_mesh, Omega_h::Mesh& target_mesh, const double radius,
  uint min_req_support, uint degree, const bool adapt_radius)
  : source_mesh_(source_mesh),
    target_mesh_(target_mesh),
    radius_(radius),
    min_req_supports_(min_req_support),
    degree_(degree),
    adapt_radius_(adapt_radius)
{
  OMEGA_H_CHECK_PRINTF(source_mesh_.dim() == 2 && target_mesh_.dim() == 2,
                       "Only 2D meshes are supported but found %d, %d\n",
                       source_mesh_.dim(), target_mesh_.dim());

  source_coords_ = source_mesh_.coords();
  target_coords_ = target_mesh_.coords();

  source_field_ =
    Omega_h::HostWrite<Omega_h::Real>(source_mesh_.nverts(), "source field");
  target_field_ =
    Omega_h::HostWrite<Omega_h::Real>(target_mesh_.nverts(), "target field");

  find_supports(min_req_supports_);
}

KOKKOS_INLINE_FUNCTION
double pointDistanceSquared(const double x1, const double y1, const double z1,
                            const double x2, const double y2, const double z2)
{
  return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
}

// replace with Kokkos::minmax_element when out of experimental
// https://kokkos.org/kokkos-core-wiki/API/algorithms/std-algorithms/all/StdMinMaxElement.html
void minmax(Omega_h::Read<Omega_h::LO> num_supports, uint& min_supports_found,
            uint& max_supports_found)
{
  using minMaxReducerType = Kokkos::MinMax<uint, Kokkos::HostSpace>;
  using minMaxValueType = minMaxReducerType::value_type;
  minMaxValueType minmax;
  Kokkos::parallel_reduce(
    num_supports.size(),
    KOKKOS_LAMBDA(int i, minMaxValueType& update) {
      if (num_supports[i] < update.min_val)
        update.min_val = num_supports[i];
      if (num_supports[i] > update.max_val)
        update.max_val = num_supports[i];
    },
    minMaxReducerType(minmax));
  Kokkos::fence();
  min_supports_found = minmax.min_val;
  max_supports_found = minmax.max_val;
}

void adapt_radii(uint min_req_supports, uint max_allowed_supports,
                 Omega_h::LO n_targets, Omega_h::Write<Omega_h::Real> radii2_l,
                 Omega_h::Write<Omega_h::LO> num_supports)
{
  Omega_h::parallel_for(
    "increase radius", n_targets, OMEGA_H_LAMBDA(const int& i) {
      Omega_h::LO nsupports = num_supports[i];
      if (nsupports < min_req_supports) {
        double factor =
          Omega_h::Real(min_req_supports) / Omega_h::Real(nsupports);
        OMEGA_H_CHECK_PRINTF(factor > 1.0,
                             "Factor should be more than 1.0: %f\n", factor);
        factor = (nsupports == 0 || factor > 1.5) ? 1.5 : factor;
        radii2_l[i] *= factor;
      } else if (nsupports > max_allowed_supports) { // if too many supports
        double factor =
          Omega_h::Real(min_req_supports) / Omega_h::Real(nsupports);
        OMEGA_H_CHECK_PRINTF(factor < 1.0,
                             "Factor should be less than 1.0: %f\n", factor);
        factor = (factor < 0.1) ? 0.33 : factor;
        radii2_l[i] *= factor;
      }
      num_supports[i] = 0; // reset the support pointer
    });
  Kokkos::fence();
}

// TODO Merge this with distance based search like NeighborSearch for
// consistency
void MLSPointCloudInterpolation::fill_support_structure(
  Omega_h::Write<Omega_h::Real> radii2_l,
  Omega_h::Write<Omega_h::LO> num_supports)
{
  // parallel scan for fill the support index with cumulative sum
  auto support_ptr_l = Omega_h::Write<Omega_h::LO>(n_targets_ + 1, 0);
  uint total_supports = 0;
  Kokkos::fence();
  Kokkos::parallel_scan(
    "scan", n_targets_,
    KOKKOS_LAMBDA(const int& i, uint& update, const bool final) {
      update += num_supports[i];
      if (final) {
        support_ptr_l[i + 1] = update;
      }
    },
    total_supports);

  pcms::printInfo("Total supports found: %d\n", total_supports);
  // resize the support index
  // supports_.supports_idx = Omega_h::Write<Omega_h::LO>(total_supports, 0);
  auto support_idx_l = Omega_h::Write<Omega_h::LO>(total_supports, 0);

  // fill the support index
  const auto dim = dim_;
  const auto target_coords_l = target_coords_;
  const auto source_coords_l = source_coords_;
  const auto n_sources = n_sources_;
  Kokkos::parallel_for(
    "fill support index", n_targets_, KOKKOS_LAMBDA(const int& target_id) {
      auto target_radius2 = radii2_l[target_id];
      auto target_coord = Omega_h::Vector<3>{0, 0, 0};
      for (int d = 0; d < dim; ++d) {
        target_coord[d] = target_coords_l[target_id * dim + d];
      }

      auto start_ptr = support_ptr_l[target_id];
      auto end_ptr = support_ptr_l[target_id + 1];

      for (int source_id = 0; source_id < n_sources; source_id++) {
        auto source_coord = Omega_h::Vector<3>{0, 0, 0};
        for (int d = 0; d < dim; ++d) {
          source_coord[d] = source_coords_l[source_id * dim + d];
        }
        auto dist2 = pointDistanceSquared(source_coord[0], source_coord[1],
                                          source_coord[2], target_coord[0],
                                          target_coord[1], target_coord[2]);
        if (dist2 <= target_radius2) {
          support_idx_l[start_ptr] = source_id;
          start_ptr++;
          OMEGA_H_CHECK_PRINTF(
            start_ptr <= end_ptr,
            "Support index out of bounds:start %d end %d target_id %d\n",
            start_ptr, end_ptr, target_id);
        }
      }
    });
  Kokkos::fence();

  // copy the support index to the supports
  supports_.radii2 = radii2_l;
  supports_.supports_ptr = Omega_h::LOs(support_ptr_l);
  supports_.supports_idx = Omega_h::LOs(support_idx_l);
}

// use uniform grid based point search when available
void MLSPointCloudInterpolation::distance_based_pointcloud_search(
  Omega_h::Write<Omega_h::Real> radii2_l,
  Omega_h::Write<Omega_h::LO> num_supports) const
{
  // n^2 search, compare each point with all other points
  const auto dim = dim_;
  const auto target_coords_l = target_coords_;
  const auto source_coords_l = source_coords_;
  const auto n_sources = n_sources_;
  const auto n_targets = n_targets_;
  Kokkos::parallel_for(
    "n^2 search", n_targets, KOKKOS_LAMBDA(const int& target_id) {
      auto target_coord = Omega_h::Vector<3>{0, 0, 0};
      for (int d = 0; d < dim; ++d) {
        target_coord[d] = target_coords_l[target_id * dim + d];
      }
      auto target_radius2 = radii2_l[target_id];

      // TODO: parallel with kokkos parallel_for
      for (int i = 0; i < n_sources; ++i) {
        auto source_coord = Omega_h::Vector<3>{0, 0, 0};
        for (int d = 0; d < dim; ++d) {
          source_coord[d] = source_coords_l[i * dim + d];
        }
        auto dist2 = pointDistanceSquared(source_coord[0], source_coord[1],
                                          source_coord[2], target_coord[0],
                                          target_coord[1], target_coord[2]);
        if (dist2 <= target_radius2) {
          num_supports[target_id]++; // only one thread is updating
        }
      }
    });
  Kokkos::fence();
}

void MLSPointCloudInterpolation::find_supports(uint min_req_supports,
                                               uint max_allowed_supports,
                                               uint max_count)
{
  pcms::printDebugInfo("First 10 Target Points with %d points:\n", n_targets_);
  Omega_h::parallel_for(
    "print target points", 10, OMEGA_H_LAMBDA(const int& i) {
      pcms::printDebugInfo("Target Point %d: (%f, %f)\n", i,
                           target_coords_[i * 2 + 0],
                           target_coords_[i * 2 + 1]);
    });
  pcms::printDebugInfo("First 10 Source Points with %d points:\n", n_sources_);
  Omega_h::parallel_for(
    "print source points", 10, OMEGA_H_LAMBDA(const int& i) {
      pcms::printDebugInfo("Source Point %d: (%f, %f)\n", i,
                           source_coords_[i * 2 + 0],
                           source_coords_[i * 2 + 1]);
    });

  auto radii2_l = Omega_h::Write<Omega_h::Real>(n_targets_, radius_);
  auto num_supports = Omega_h::Write<Omega_h::LO>(n_targets_, 0);

  uint min_supports_found = 0;
  uint max_supports_found = 0;
  // radius adjustment loop
  int loop_count = 0;
  while (!within_number_of_support_range(min_supports_found, max_supports_found,
                                         min_req_supports,
                                         max_allowed_supports)) {
    distance_based_pointcloud_search(radii2_l, num_supports);

    loop_count++;
    if (!adapt_radius_) {
      break;
    }

    Kokkos::fence();
    if (loop_count > max_count) {
      pcms::printError(
        "Loop count exceeded 100 and still not converged.\n"
        "Manually check if the number of minimum and maximum supports are "
        "reasonable.\n"
        "There are situations when it may not converge.\n");
      break;
    }

    // find minimum and maximum number of supports found
    minmax(num_supports, min_supports_found, max_supports_found);

    // increase radius if not enough supports or too many supports
    if (!within_number_of_support_range(min_supports_found, max_supports_found,
                                        min_req_supports,
                                        max_allowed_supports)) {
      pcms::printInfo(
        "Adjusting radius: Iter %d:(min: %d max: %d) Min Req: %d Max "
        "Allowed %d\n",
        loop_count, min_supports_found, max_supports_found, min_req_supports,
        max_allowed_supports);
      adapt_radii(min_req_supports, max_allowed_supports, n_targets_, radii2_l,
                  num_supports);
    }
  }
  pcms::printInfo("Searched %d times and supports found: min: %d max: %d\n",
                  loop_count, min_supports_found, max_supports_found);

  fill_support_structure(radii2_l, num_supports);
}

// TODO : find way to avoid this copy
void copyHostScalarArrayView2HostWrite(
  pcms::Rank1View<double, pcms::HostMemorySpace> source,
  Omega_h::HostWrite<Omega_h::Real>& target)
{
  OMEGA_H_CHECK_PRINTF(
    source.size() == target.size(),
    "Size mismatch in copy_data_from_ScalarArray_to_HostWrite: %zu %d\n",
    source.size(), target.size());
  OMEGA_H_CHECK_PRINTF(source.data_handle() != target.data(),
                       "Source and Target contain the same pointer %p\n",
                       source.data_handle());

  for (int i = 0; i < source.size(); ++i) {
    target[i] = source[i];
  }
}
void copyHostWrite2ScalarArrayView(
  const Omega_h::HostWrite<Omega_h::Real>& source,
  pcms::Rank1View<double, pcms::HostMemorySpace> target)
{
  OMEGA_H_CHECK_PRINTF(
    source.size() == target.size(),
    "Size mismatch in copy_data_from_HostWrite_to_ScalarArray: %d %zu\n",
    source.size(), target.size());
  OMEGA_H_CHECK_PRINTF(source.data() != target.data_handle(),
                       "Source and Target contain the same pointer %p\n",
                       source.data());

  for (int i = 0; i < source.size(); ++i) {
    target[i] = source[i];
  }
}

void MLSPointCloudInterpolation::eval(
  pcms::Rank1View<double, pcms::HostMemorySpace> source_field,
  pcms::Rank1View<double, pcms::HostMemorySpace> target_field)
{
  OMEGA_H_CHECK_PRINTF(target_field.size() == target_coords_.size() / dim_,
                       "Source Data and Source Points size mismatch: %zu %d\n",
                       target_field.size(), target_coords_.size() / dim_);

  OMEGA_H_CHECK_PRINTF(source_field.size() == source_coords_.size() / dim_,
                       "Target Data and Target Points size mismatch: %zu %d\n",
                       source_field.size(), source_coords_.size() / dim_);
  for (int i = 0; i < 10; i++) {
    pcms::printDebugInfo("i = %d	field = %f\n", i, source_field[i]);
  }

  copyHostScalarArrayView2HostWrite(source_field, source_field_);

  // TODO: make the basis function a template or pass it as a parameter
  auto target_field_write = mls_interpolation(
    Omega_h::Reals(source_field_), source_coords_, target_coords_, supports_, 2,
    degree_, pcms::RadialBasisFunction::RBF_GAUSSIAN, lambda_, 1e-6,
    decay_factor_);

  target_field_ = Omega_h::HostWrite<Omega_h::Real>(target_field_write);
  copyHostWrite2ScalarArrayView(target_field_, target_field);
}

void MLSInterpolationHandler::eval(
  pcms::Rank1View<double, pcms::HostMemorySpace> source_field,
  pcms::Rank1View<double, pcms::HostMemorySpace> target_field)
{
  OMEGA_H_CHECK_PRINTF(
    target_field.size() == target_coords_.size() / target_mesh_.dim(),
    "Source Data and Source Points size mismatch: %zu %d\n",
    target_field.size(), target_coords_.size() / target_mesh_.dim());

  OMEGA_H_CHECK_PRINTF(
    source_field.size() == source_coords_.size() / source_mesh_.dim(),
    "Target Data and Target Points size mismatch: %zu %d\n",
    source_field.size(), source_coords_.size() / source_mesh_.dim());

  copyHostScalarArrayView2HostWrite(source_field, source_field_);

  // TODO: make the basis function a template or pass it as a parameter
  auto target_field_write = mls_interpolation(
    Omega_h::Reals(source_field_), source_coords_, target_coords_, supports_, 2,
    degree_, pcms::RadialBasisFunction::RBF_GAUSSIAN, 0, 1e-6, 5);

  target_field_ = Omega_h::HostWrite<Omega_h::Real>(target_field_write);
  copyHostWrite2ScalarArrayView(target_field_, target_field);
}

void MLSInterpolationHandler::find_supports(const uint min_req_support)
{
  if (single_mesh_) {
    supports_ =
      searchNeighbors(source_mesh_, radius_, min_req_support, adapt_radius_);
  } else { // two mesh : vert to vert
    supports_ = searchNeighbors(source_mesh_, target_mesh_, radius_,
                                min_req_support, adapt_radius_);
  }

#ifndef NDEBUG
  Omega_h::HostRead<Omega_h::Real> hostRadii2(supports_.radii2);
  for (size_t i = 0; i < hostRadii2.size(); ++i) {
    OMEGA_H_CHECK_PRINTF(
      hostRadii2[i] > 1e-10,
      "Radius squared has to be more than zero. Found [%zu] = %f\n", i,
      hostRadii2[i]);
  }
#endif
}

size_t MLSInterpolationHandler::getSourceSize()
{
  if (single_mesh_) {
    return source_mesh_.nfaces();
  } else {
    return source_mesh_.nverts();
  }
}

size_t MLSInterpolationHandler::getTargetSize()
{
  if (single_mesh_) {
    return source_mesh_.nverts();
  } else {
    return target_mesh_.nverts();
  }
}
