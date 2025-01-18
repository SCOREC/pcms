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
  Omega_h::Write<Omega_h::Real> centroids(mesh.nfaces() * 2, 0.0);

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

MLSInterpolationHandler::MLSInterpolationHandler(
  const std::string& source_mesh_fname, double radius, bool adapt_radius)
  : radius_(radius), adapt_radius_(adapt_radius)
{
  single_mesh_ = true;
  library_ = Omega_h::Library(nullptr, nullptr);

  source_mesh_ = Omega_h::binary::read(source_mesh_fname, &library_);
  source_coords_ = source_mesh_.coords();
  target_coords_ = getCentroids(source_mesh_);

  OMEGA_H_CHECK_PRINTF(source_mesh_.dim() == 2,
                       "Only 2D meshes are supported but found %d\n",
                       source_mesh_.dim());

  source_field_ =
    Omega_h::HostWrite<Omega_h::Real>(source_mesh_.nverts(), "source field");

  // used nfaces to interpolate to the centroid of the faces
  target_field_ =
    Omega_h::HostWrite<Omega_h::Real>(source_mesh_.nfaces(), "target field");

  find_supports();
}

MLSInterpolationHandler::MLSInterpolationHandler(Omega_h::Mesh source_mesh,
                                                 double radius,
                                                 bool adapt_radius)
  : radius_(radius), adapt_radius_(adapt_radius)
{
  single_mesh_ = true;
  source_mesh_ = std::move(source_mesh);
  source_coords_ = source_mesh_.coords();
  target_coords_ = getCentroids(source_mesh_);

  OMEGA_H_CHECK_PRINTF(source_mesh_.dim() == 2,
                       "Only 2D meshes are supported but found %d\n",
                       source_mesh_.dim());

  source_field_ =
    Omega_h::HostWrite<Omega_h::Real>(source_mesh_.nverts(), "source field");

  target_field_ =
    Omega_h::HostWrite<Omega_h::Real>(source_mesh_.nfaces(), "target field");
}

// TODO: delegate this constructor to a single one
MLSInterpolationHandler::MLSInterpolationHandler(
  const std::string& source_mesh_fname, const std::string& target_mesh_fname,
  const double radius, const bool adapt_radius)
  : radius_(radius), adapt_radius_(adapt_radius)
{
  // TODO: take argc and argv
  library_ = Omega_h::Library(nullptr, nullptr);

  source_mesh_ = Omega_h::binary::read(source_mesh_fname, &library_);
  target_mesh_ = Omega_h::binary::read(target_mesh_fname, &library_);

  OMEGA_H_CHECK_PRINTF(source_mesh_.dim() == 2 && target_mesh_.dim() == 2,
                       "Only 2D meshes are supported but found %d, %d\n",
                       source_mesh_.dim(), target_mesh_.dim());

  source_coords_ = source_mesh_.coords();
  target_coords_ = target_mesh_.coords();

  source_field_ =
    Omega_h::HostWrite<Omega_h::Real>(source_mesh_.nverts(), "source field");
  target_field_ =
    Omega_h::HostWrite<Omega_h::Real>(target_mesh_.nverts(), "target field");

  find_supports();
}

MLSInterpolationHandler::MLSInterpolationHandler(Omega_h::Mesh source_mesh,
                                                 Omega_h::Mesh target_mesh,
                                                 const double radius,
                                                 const bool adapt_radius)
  : radius_(radius), adapt_radius_(adapt_radius)
{
  // TODO: check if move works for mesh
  source_mesh_ = std::move(source_mesh);
  target_mesh_ = std::move(target_mesh);

  OMEGA_H_CHECK_PRINTF(source_mesh_.dim() == 2 && target_mesh_.dim() == 2,
                       "Only 2D meshes are supported but found %d, %d\n",
                       source_mesh_.dim(), target_mesh_.dim());

  source_coords_ = source_mesh_.coords();
  target_coords_ = target_mesh_.coords();

  source_field_ =
    Omega_h::HostWrite<Omega_h::Real>(source_mesh_.nverts(), "source field");
  target_field_ =
    Omega_h::HostWrite<Omega_h::Real>(target_mesh_.nverts(), "target field");

  find_supports();
}

// TODO : find way to avoid this copy
void copyHostScalarArrayView2HostWrite(
  const pcms::ScalarArrayView<double, pcms::HostMemorySpace>& source,
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
  pcms::ScalarArrayView<double, pcms::HostMemorySpace>& target)
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
  const pcms::ScalarArrayView<double, pcms::HostMemorySpace>& source_field,
  pcms::ScalarArrayView<double, pcms::HostMemorySpace>& target_field)
{
  OMEGA_H_CHECK_PRINTF(source_field.size() ==
                         source_coords_.size() / source_mesh_.dim(),
                       "Source Data and Source Points size mismatch: %zu %d\n",
                       source_field.size(), source_coords_.size());
  if (single_mesh_) {
    OMEGA_H_CHECK_PRINTF(
      target_field.size() == source_mesh_.nfaces(),
      "Target Data and Target Points size mismatch: %zu %d\n",
      target_field.size(), target_coords_.size());
  } else {
    OMEGA_H_CHECK_PRINTF(
      target_field.size() == target_coords_.size() / target_mesh_.dim(),
      "Target Data and Target Points size mismatch: %zu %d\n",
      target_field.size(), target_coords_.size());
  }

  copyHostScalarArrayView2HostWrite(source_field, source_field_);

  auto target_field_write = mls_interpolation(
    Omega_h::Reals(source_field_), source_coords_, target_coords_, supports_, 2,
    4, supports_.radii2, RadialBasisFunction::RBF_GAUSSIAN);

  target_field_ = Omega_h::HostWrite<Omega_h::Real>(target_field_write);
  copyHostWrite2ScalarArrayView(target_field_, target_field);
}

void MLSInterpolationHandler::find_supports(
  const int min_req_support) // ? Needed?
{
  if (single_mesh_) {
    supports_ =
      searchNeighbors(source_mesh_, radius_, min_req_support, adapt_radius_);
  } else { // not single mesh
    supports_ = searchNeighbors(source_mesh_, target_mesh_, radius_,
                                min_req_support, adapt_radius_);
  }

  OMEGA_H_CHECK_PRINTF(
    supports_.radii2[0] < 1e-10,
    "Radius squared has to be more than zero found found [0] = %f\n",
    supports_.radii2[0]);

  for (int i = 0; i < supports_.radii2.size(); ++i) {
    printf("%f %d\n", supports_.radii2[i], i);
  }
}