#ifndef PCMS_INTERPOLATION_BASE_H
#define PCMS_INTERPOLATION_BASE_H
//
// Created by Fuad Hasan on 1/12/25.
//
#include "mls_interpolation.hpp"
#include "adj_search.hpp"
#include "interpolation_helpers.h"
#include <Omega_h_file.hpp>
#include <pcms/arrays.h>
#include <string>

class InterpolationBase
{
public:
  virtual ~InterpolationBase() = default;
  /**
   * @brief Evaluate the interpolation
   * @param source_field The field to interpolate from
   * @param target_field The field to interpolate to
   */
  virtual void eval(
    // TODO: Should these be templated to support different types?
    pcms::Rank1View<double, pcms::HostMemorySpace> source_field,
    pcms::Rank1View<double, pcms::HostMemorySpace> target_field) = 0;

  /**
   * @brief Get the size of the source field
   * @return size_t The size of the source field
   */
  virtual size_t getSourceSize() const = 0;

  /**
   * @brief Get the size of the target field
   * @return size_t The size of the target field
   */
  virtual size_t getTargetSize() const = 0;
};

/**
 *@brief Meshless Point-Cloud Based Interpolation Using MLS
 */
class MLSPointCloudInterpolation final : public InterpolationBase
{
public:
  template <typename SourceType, typename TargetType>
  MLSPointCloudInterpolation(
    pcms::Rank1View<SourceType, pcms::HostMemorySpace> source_points,
    pcms::Rank1View<TargetType, pcms::HostMemorySpace> target_points, int dim,
    double radius, uint min_req_supports = 10, uint degree = 3,
    bool adapt_radius = true, double lambda = 0.0, double decay_factor = 5.0)
    : dim_(dim),
      radius_(radius),
      adapt_radius_(adapt_radius),
      degree_(degree),
      min_req_supports_(min_req_supports),
      n_targets_(target_points.size() / dim),
      n_sources_(source_points.size() / dim),
      lambda_(lambda),
      decay_factor_(decay_factor)
  {
    // check if source and target types can fall back to const double
    static_assert(std::is_convertible_v<SourceType, const Omega_h::Real>,
                  "SourceType must be convertible to const double");
    static_assert(std::is_convertible_v<TargetType, const Omega_h::Real>,
                  "TargetType must be convertible to const double");

    source_field_ = Omega_h::HostWrite<Omega_h::Real>(
      source_points.size() / dim_, "source field");
    target_field_ = Omega_h::HostWrite<Omega_h::Real>(
      target_points.size() / dim_, "target field");

    Omega_h::HostWrite<Omega_h::Real> source_coords_host(source_points.size(),
                                                         "source points");
    Omega_h::HostWrite<Omega_h::Real> target_coords_host(target_points.size(),
                                                         "target points");
    // TODO Remove these copies
    for (int i = 0; i < source_points.size(); ++i) {
      source_coords_host[i] = source_points[i];
    }
    for (int i = 0; i < target_points.size(); ++i) {
      target_coords_host[i] = target_points[i];
    }
    source_coords_ = Omega_h::Reals(source_coords_host);
    target_coords_ = Omega_h::Reals(target_coords_host);

    find_supports(min_req_supports_, 3 * min_req_supports_, 100);
  }

  void eval(
    pcms::Rank1View<double, pcms::HostMemorySpace> source_field,
    pcms::Rank1View<double, pcms::HostMemorySpace> target_field) override;

  SupportResults getSupports() { return supports_; }
  size_t getSourceSize() const override { return source_coords_.size() / dim_; }
  size_t getTargetSize() const override { return target_coords_.size() / dim_; }

private:
  int dim_;
  double radius_;
  bool adapt_radius_;
  uint degree_;
  uint min_req_supports_;
  double lambda_;
  double decay_factor_;

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
class MLSMeshInterpolation final : public InterpolationBase
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
  MLSMeshInterpolation(Omega_h::Mesh& source_mesh, Omega_h::Mesh& target_mesh,
                       double radius, uint min_req_supports = 10,
                       uint degree = 3, bool adapt_radius = true);

  /**
   * @brief Centroids to Vertices interpolation for a single mesh
   * @param source_mesh The source mesh
   * @param radius The cutoff radius for the MLS interpolation
   * @param adapt_radius Whether to adapt the radius based on the local density
   * @param min_req_supports Min number of source locations required
   * @param degree The degree of the polynomial used in the MLS interpolation
   */
  MLSMeshInterpolation(Omega_h::Mesh& source_mesh, double radius,
                       uint min_req_supports = 10, uint degree = 3,
                       bool adapt_radius = true);

  size_t getSourceSize() const override;
  size_t getTargetSize() const override;

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
