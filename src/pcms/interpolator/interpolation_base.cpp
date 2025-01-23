//
// Created by hasanm4 on 1/17/25.
//

#include "interpolation_base.h"
#include "MLS_rbf_options.hpp"

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

// TODO : find way to avoid this copy
void copyHostScalarArrayView2HostWrite(
  pcms::ScalarArrayView<double, pcms::HostMemorySpace> source,
  Omega_h::HostWrite<Omega_h::Real>& target)
{
  OMEGA_H_CHECK_PRINTF(
    source.size() == target.size(),
    "Size mismatch in copy_data_from_ScalarArray_to_HostWrite: %zu %d\n",
    source.size(), target.size());

  for (int i = 0; i < source.size(); ++i) {
    target[i] = source[i];
  }
}
void copyHostWrite2ScalarArrayView(
  const Omega_h::HostWrite<Omega_h::Real>& source,
  pcms::ScalarArrayView<double, pcms::HostMemorySpace> target)
{
  OMEGA_H_CHECK_PRINTF(
    source.size() == target.size(),
    "Size mismatch in copy_data_from_HostWrite_to_ScalarArray: %d %zu\n",
    source.size(), target.size());

  for (int i = 0; i < source.size(); ++i) {
    target[i] = source[i];
  }
}

void MLSInterpolationHandler::eval(
  pcms::ScalarArrayView<double, pcms::HostMemorySpace> source_field,
  pcms::ScalarArrayView<double, pcms::HostMemorySpace> target_field)
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
    degree_, supports_.radii2, RadialBasisFunction::RBF_GAUSSIAN);

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
      "Radius squared has to be more than zero found found [%zu] = %f\n", i,
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