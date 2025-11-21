/**
 * @file interpolation_base.h
 * @brief Interpolation classes to wrap and streamline various interpolation
 * methods
 *
 * These classes provide a unified interface for different interpolation
 * techniques, hold necessary data structures, and manage resources effectively.
 * This also provides the foundation of the C and Fortran bindings.
 *
 */

#ifndef PCMS_INTERPOLATION_BASE_H
#define PCMS_INTERPOLATION_BASE_H

#include "mls_interpolation.hpp"
#include "adj_search.hpp"
#include "interpolation_helpers.h"
#include <Omega_h_file.hpp>
#include <pcms/arrays.h>
#include <string>

/**
 * @brief Pure virtual base class for interpolation methods
 * @details Provides external interface for interpolation methods.
 */
class InterpolationBase
{
public:
  virtual ~InterpolationBase() = default;
  /**
   * @brief Evaluate the interpolation
   * @param source_field The field to interpolate from
   * @param target_field The field to interpolate to
   * @note User is responsible for ensuring that the source and target fields
   * are sustained during the evaluation. The size of the source and target
   * fields must match the sizes returned by getSourceSize() and getTargetSize()
   * respectively.
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
 *@brief Meshless Point-Cloud Based Moving Least Square (MLS) Interpolation
 *@details If two point clouds are provided, mesh based accelerated adjacency
 *search cannot be used. Instead, each target points search the whole source
 *point cloud for neighbors within a given radius.
 */
class MLSPointCloudInterpolation final : public InterpolationBase
{
public:
  /**
   * @brief Constructor for point-cloud based MLS interpolation
   * @tparam SourceType Source-point data type: any type convertible to const
   * double
   * @tparam TargetType Target-point data type: any type convertible to const
   * double
   * @param source_points The source point coordinates
   * @param target_points The target point coordinates
   * @param dim The spatial dimension of the points
   * @param radius The cutoff radius for the MLS interpolation
   * @param min_req_supports The minimum number of source locations required for
   * interpolation
   * @param degree The degree of the polynomial used in the MLS interpolation
   * @param adapt_radius Whether to adapt the radius to satisfy the minimum and
   * maximum number of supports (maximum is set to 3 times the minimum)
   * @param lambda Regularization parameter for the MLS interpolation
   * @param decay_factor Decay factor for the weight function in the MLS
   * interpolation
   *
   * For more details about MLS interpolation parameters, refer to
   * the documentation of mls_interpolation and RadialBasisFunction.
   *
   * - Source and target point coordinates are expected to be in a flattened
   * array format. For example, for 3D points, the coordinates should be
   * provided as [x1, y1, z1, x2, y2, z2, ..., xn, yn, zn] with n being the
   * number of points.
   * - It constructors the support structure by performing \f$(N_{targets}
   * \times N_{sources})\f$ searches.
   *
   * @see mls_interpolation, RadialBasisFunction
   */
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
  int dim_;           /*!< Spatial dimension of the point clouds */
  double radius_;     /*!< Cutoff radius for the MLS interpolation */
  bool adapt_radius_; /*!< Whether to adapt the radius based on local density */
  uint degree_; /*!< Degree of the polynomial used in the MLS interpolation */
  uint min_req_supports_; /*!< Minimum number of source locations required */
  double lambda_; /*!< Regularization parameter for the MLS interpolation */
  double decay_factor_; /*!< Decay factor for the weight function in the MLS
                         interpolation */

  // InterpolationType interpolation_type_;
  Omega_h::LO n_sources_ = 0;    /*!< Number of source points */
  Omega_h::LO n_targets_ = 0;    /*!< Number of target points */
  Omega_h::Reals source_coords_; /*!< Source point coordinates */
  Omega_h::Reals target_coords_; /*!< Target point coordinates */

  SupportResults supports_; /*!< Support structure for MLS interpolation */

  Omega_h::HostWrite<Omega_h::Real> target_field_; /*!< Target field storage */
  Omega_h::HostWrite<Omega_h::Real> source_field_; /*!< Provided source field */

  void fill_support_structure(Omega_h::Write<Omega_h::Real> radii2_l,
                              Omega_h::Write<Omega_h::LO> num_supports);

  /**
   * @brief Perform distance-based search with the given radii
   * @param radii2_l Squared radii for each target point
   * @param num_supports Number of supports found for each target point
   */
  void distance_based_pointcloud_search(
    Omega_h::Write<Omega_h::Real> radii2_l,
    Omega_h::Write<Omega_h::LO> num_supports) const;

  /**
   * @brief Find supports for each target point
   * @param min_req_supports Minimum required supports
   * @param max_allowed_supports Maximum allowed supports
   * @param max_count Maximum number of iterations to adjust radius
   */
  void find_supports(uint min_req_supports = 10, uint max_allowed_supports = 30,
                     uint max_count = 100);
};

/**
 * @brief Moving Least Square (MLS) Interpolation with adjacency search
 * @details Supports two modes:
 * - Vertex to Vertex interpolation between two meshes
 * - Centroid to Vertex interpolation within a single mesh
 */
class MLSMeshInterpolation final : public InterpolationBase
{

public:
  void eval(
    pcms::Rank1View<double, pcms::HostMemorySpace> source_field,
    pcms::Rank1View<double, pcms::HostMemorySpace> target_field) override;

  /**
   * @brief Vertex to Vertex interpolation between two meshes
   * @param source_mesh Source mesh
   * @param target_mesh Target mesh
   * @param radius The cutoff radius for the MLS interpolation
   * @param min_req_supports Min number of source locations required (max is set
   * to 3x min)
   * @param degree The degree of the polynomial used in the MLS interpolation
   * @param adapt_radius Whether to adapt the radius based on the local density
   * @param lambda Regularization parameter for the MLS interpolation
   * @param decay_factor Decay factor for the weight function in the MLS
   * interpolation
   *
   * @details For more details about MLS interpolation parameters, refer to
   * the documentation of mls_interpolation and RadialBasisFunction.
   *
   * @Note Both source and target meshes must have the same spatial dimension.
   *
   * @see mls_interpolation, RadialBasisFunction
   */
  MLSMeshInterpolation(Omega_h::Mesh& source_mesh, Omega_h::Mesh& target_mesh,
                       double radius, uint min_req_supports = 10,
                       uint degree = 3, bool adapt_radius = true,
                       double lambda = 0.0, double decay_factor = 5.0);

  /**
   * @brief Centroid to Vertex interpolation within a single mesh
   * @param source_mesh The mesh to interpolate within
   * @param radius The cutoff radius for the MLS interpolation
   * @param min_req_supports Min number of source locations required (max is set
   * to 3x min)
   * @param degree The degree of the polynomial used in the MLS interpolation
   * @param adapt_radius Whether to adapt the radius based on the local density
   * @param lambda Regularization parameter for the MLS interpolation
   * @param decay_factor Decay factor for the weight function in the MLS
   * interpolation
   *
   * @details For more details about MLS interpolation parameters, refer to
   * the documentation of mls_interpolation and RadialBasisFunction.
   *
   * @see mls_interpolation, RadialBasisFunction
   */
  MLSMeshInterpolation(Omega_h::Mesh& source_mesh, double radius,
                       uint min_req_supports = 10, uint degree = 3,
                       bool adapt_radius = true, double lambda = 0.0,
                       double decay_factor = 5.0);

  size_t getSourceSize() const override;
  size_t getTargetSize() const override;

  SupportResults getSupports() { return supports_; }

private:
  double radius_; /*!< Cutoff radius for the MLS interpolation */
  double lambda_; /*!< Regularization parameter for the MLS interpolation */
  double decay_factor_; /*!< Decay factor for the weight function in the MLS
                     interpolation */
  bool adapt_radius_; /*!< Whether to adapt the radius based on local density */
  bool single_mesh_ = false; /*!< Whether single mesh mode is used */
  uint degree_; /*!< Degree of the polynomial used in the MLS interpolation */
  uint min_req_supports_; /*!< Minimum number of source locations required */

  // InterpolationType interpolation_type_;

  Omega_h::Mesh& source_mesh_; /*!< Reference to the source mesh */
  // TODO: handle what to do with this when only 1 mesh is provided
  Omega_h::Mesh& target_mesh_;   /*!< Reference to the target mesh */
  Omega_h::Reals source_coords_; /*!< Source point coordinates */
  Omega_h::Reals target_coords_; /*!< Target point coordinates */

  SupportResults supports_; /*!< Support structure for MLS interpolation */

  Omega_h::HostWrite<Omega_h::Real> target_field_; /*!< Target field storage */
  Omega_h::HostWrite<Omega_h::Real> source_field_; /*!< Provided source field */

  /**
   * @brief Adjacency-based search to find supports
   * @param min_req_supports Minimum required supports
   * @param max_allowed_supports Maximum allowed supports
   *
   * @see searchNeighbors
   */
  void find_supports(uint min_req_supports = 10,
                     uint max_allowed_supports = 30);
};

#endif // PCMS_INTERPOLATION_BASE_H
