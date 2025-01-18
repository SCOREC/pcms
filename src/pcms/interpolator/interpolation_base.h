//
// Created by Fuad Hasan on 1/12/25.
//
#include "MLSInterpolation.hpp"
#include "adj_search.hpp"
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <pcms/arrays.h>
#include <string>

#ifndef PCMS_INTERPOLATION_BASE_H
#define PCMS_INTERPOLATION_BASE_H

class InterpolationBase
{
  /**
   * @brief Evaluate the interpolation
   * @param source_field The field to interpolate from
   * @param target_field The field to interpolate to
   */
  virtual void eval(
    const pcms::ScalarArrayView<double, pcms::HostMemorySpace>& source_field,
    pcms::ScalarArrayView<double, pcms::HostMemorySpace>& target_field) = 0;
};

/**
 * @brief Moving Least Square Radial Basis Function Interpolation
 */
class MLSInterpolationHandler : public InterpolationBase
{

public:
  void eval(
    const pcms::ScalarArrayView<double, pcms::HostMemorySpace>& source_field,
    pcms::ScalarArrayView<double, pcms::HostMemorySpace>& target_field)
    override;

  MLSInterpolationHandler(const std::string& source_mesh_fname,
                   const std::string& target_mesh_fname, double radius,
                   bool adapt_radius = true);

  MLSInterpolationHandler(Omega_h::Mesh source_mesh, Omega_h::Mesh target_mesh,
                   double radius, bool adapt_radius = true);

  MLSInterpolationHandler(const std::string& source_mesh_fname, double radius,
                   bool adapt_radius = true);

  MLSInterpolationHandler(Omega_h::Mesh source_mesh, double radius,
                   bool adapt_radius = true);

private:
  double radius_;
  bool adapt_radius_;
  bool single_mesh_ = false;

  std::string interpolation_type_;

  Omega_h::Library library_;
  Omega_h::Mesh source_mesh_;
  // TODO: handle what to do with this when only 1 mesh is provided
  Omega_h::Mesh target_mesh_;
  Omega_h::Reals source_coords_;
  Omega_h::Reals target_coords_;

  SupportResults supports_;

  Omega_h::HostWrite<Omega_h::Real> target_field_;
  Omega_h::HostWrite<Omega_h::Real> source_field_;

  void find_supports(int min_req_supports = 12);
};

#endif // PCMS_INTERPOLATION_BASE_H
