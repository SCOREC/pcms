#ifndef PCMS_INTERPOLATION_BASE_H
#define PCMS_INTERPOLATION_BASE_H
//
// Created by Fuad Hasan on 1/12/25.
//
#include "mls_interpolation.hpp"
#include "adj_search.hpp"
#include "interpolation_helpers.h"
#include <Omega_h_file.hpp>
// #include <Omega_h_library.hpp>
#include <pcms/arrays.h>
#include <string>

class InterpolationBase
{
  /**
   * @brief Evaluate the interpolation
   * @param source_field The field to interpolate from
   * @param target_field The field to interpolate to
   */
  virtual void eval(
    pcms::Rank1View<double, pcms::HostMemorySpace> source_field,
    pcms::Rank1View<double, pcms::HostMemorySpace> target_field) = 0;
};

/**
 *@brief Meshless interpolation using MLS
 */
class MLSPointCloudInterpolation : public InterpolationBase
{
public:
  MLSPointCloudInterpolation(
    pcms::Rank1View<double, pcms::HostMemorySpace> source_points,
    pcms::Rank1View<double, pcms::HostMemorySpace> target_points, int dim,
    double radius, uint min_req_supports = 10, uint degree = 3,
    bool adapt_radius = true);

  void eval(
    pcms::Rank1View<double, pcms::HostMemorySpace> source_field,
    pcms::Rank1View<double, pcms::HostMemorySpace> target_field) override;

  SupportResults getSupports() { return supports_; }
  size_t getSourceSize() const { return source_coords_.size() / dim_; }
  size_t getTargetSize() const { return target_coords_.size() / dim_; }

private:
  int dim_;
  double radius_;
  bool adapt_radius_;
  // bool single_mesh_ = false;
  uint degree_;
  uint min_req_supports_;

  // InterpolationType interpolation_type_;
  Omega_h::LO n_sources_ = 0;
  Omega_h::LO n_targets_ = 0;
  Omega_h::Reals source_coords_;
  Omega_h::Reals target_coords_;

  SupportResults supports_;

  Omega_h::HostWrite<Omega_h::Real> target_field_;
  Omega_h::HostWrite<Omega_h::Real> source_field_;

  void fill_support_structure(Omega_h::Write<Omega_h::Real> radii2_l,
                              Omega_h::Write<Omega_h::LO> num_supports);
  void distance_based_pointcloud_search(
    Omega_h::Write<Omega_h::Real> radii2_l,
    Omega_h::Write<Omega_h::LO> num_supports) const;
  void find_supports(uint min_req_supports = 10, uint max_allowed_supports = 30,
                     uint max_count = 100);
};

/**
 * @brief Moving Least Square Radial Basis Function Interpolation
 */
class MLSInterpolationHandler : public InterpolationBase
{

public:
  void eval(
    pcms::Rank1View<double, pcms::HostMemorySpace> source_field,
    pcms::Rank1View<double, pcms::HostMemorySpace> target_field) override;

  /**
   * @brief Vertex to Vertex interpolation for two given meshes
   * @param source_mesh The source mesh
   * @param target_mesh The target mesh
   * @param radius The cutoff radius for the MLS interpolation
   * @param min_req_supports The minimum number of source locations required for
   * interpolation
   * @param degree The degree of the polynomial used in the MLS interpolation
   * @param adapt_radius Whether to adapt the radius based on the local density
   */
  MLSInterpolationHandler(Omega_h::Mesh& source_mesh,
                          Omega_h::Mesh& target_mesh, double radius,
                          uint min_req_supports = 10, uint degree = 3,
                          bool adapt_radius = true);

  /**
   * @brief Centroids to Vertices interpolation for a single mesh
   * @param source_mesh The source mesh
   * @param radius The cutoff radius for the MLS interpolation
   * @param adapt_radius Whether to adapt the radius based on the local density
   * @param min_req_supports Min number of source locations required
   * @param degree The degree of the polynomial used in the MLS interpolation
   */
  MLSInterpolationHandler(Omega_h::Mesh& source_mesh, double radius,
                          uint min_req_supports = 10, uint degree = 3,
                          bool adapt_radius = true);

  size_t getSourceSize();
  size_t getTargetSize();

  SupportResults getSupports() { return supports_; }

private:
  double radius_;
  bool adapt_radius_;
  bool single_mesh_ = false;
  uint degree_;
  uint min_req_supports_;

  // InterpolationType interpolation_type_;

  Omega_h::Mesh& source_mesh_;
  // TODO: handle what to do with this when only 1 mesh is provided
  Omega_h::Mesh& target_mesh_;
  Omega_h::Reals source_coords_;
  Omega_h::Reals target_coords_;

  SupportResults supports_;

  Omega_h::HostWrite<Omega_h::Real> target_field_;
  Omega_h::HostWrite<Omega_h::Real> source_field_;

  void find_supports(uint min_req_supports = 10);
};

#endif // PCMS_INTERPOLATION_BASE_H
