#include <Omega_h_bbox.cpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_MathematicalFunctions.hpp>
#include <Kokkos_Random.hpp>
#include <MeshField_Shape.hpp>
#include <Omega_h_array.hpp>
#include <Omega_h_for.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_shape.hpp>
#include <fstream>
#include <iostream>
#include <pcms/localization/point_search.h>
#include <sstream>
#include <vector>
#include <cstdint>
#include <stdexcept>

constexpr int MAX_POINTS = 1600;
/**
 * @brief Read precomputed Sobol barycentric samples from file.
 *
 * For now this routine reads barycentric Sobol samples from a text file with a
 * header row. TODO: replace this with a function that directly generates Sobol
 * sequence samples and maps them to barycentric coordinates.
 *
 * @param[in] file_path Path to the sample file.
 * @return Device view of shape `(nsamples, 3)` containing barycentric samples.
 */
Kokkos::View<MeshField::Real* [3]> read_sobol_barycentric_samples_from_file(
  std::string file_path);

/**
 * @brief Generate uniform random barycentric coordinates for triangle sampling.
 *
 * Uses the standard square-to-triangle transform to produce barycentric
 * coordinates uniformly distributed over triangle area.
 *
 * @param[in] npoints_each_tri Number of samples to generate.
 * @return Device view of shape `(npoints_each_tri, 3)`.
 */
Kokkos::View<MeshField::Real* [3]> generate_uniform_random_barycentric_coords(
  const int npoints_each_tri);

/**
 * @brief Compute global sample coordinates in each element from reference
 *        barycentric sample coordinates.
 *
 * Reuses the same reference-triangle barycentric sample pattern in every target
 * element and maps each barycentric sample to its corresponding physical
 * coordinate using the element vertex coordinates.
 *
 * @param[in] target_mesh Target triangular mesh.
 * @param[in] ref_barycentric_coords Reference barycentric coordinates of shape
 *                                   `(npoints_each_tri, 3)`.
 * @return Device view of shape `(nelems * npoints_each_tri, 2)` containing the
 *         global sample coordinates in all elements.
 */
Kokkos::View<pcms::Real* [2]> global_coords_from_ref_barycentric_coords(
  Omega_h::Mesh& target_mesh,
  const Kokkos::View<MeshField::Real* [3]>& ref_barycentric_coords);

/**
 * @brief Locate query points in a 2D triangular mesh.
 *
 * This routine performs point localization for a set of query points in the
 * given 2D mesh using a structured grid search. For each query
 * point, it returns the containing element id together with the associated
 * parametric coordinates stored in `pcms::GridPointSearch::Result`.
 *
 * @param[in] mesh Input 2D mesh in which the query points are to be localized.
 * @param[in] points Query point coordinates of shape `(npoints, 2)`.
 *
 * @return Search results for all query points.
 *
 * @note This routine assumes that `mesh` is two-dimensional.
 */
Kokkos::View<pcms::GridPointSearch::Result*> localize_points_in_mesh(
  Omega_h::Mesh& mesh, const Kokkos::View<pcms::Real* [2]>& points);

/**
 * @brief Evaluate a nodal field at localized query points.
 *
 * Uses the containing element id and parametric or barycentric coordinates from
 * a point-localization step to interpolate the nodal field values at the query
 * points.
 *
 * @param[in] mesh Input triangular mesh.
 * @param[in] nodal_field_values Field values stored at mesh vertices.
 * @param[in] results Point-localization results for the query points.
 * @return Field values evaluated at the query points.
 *
 * @note Points with invalid element ids retain the default value `0.0`.
 */
Omega_h::Reals evaluate_field_from_point_localization(
  Omega_h::Mesh& mesh, const Omega_h::Reals& nodal_field_values,
  const Kokkos::View<pcms::Real* [2]>& points,
  const Kokkos::View<pcms::GridPointSearch::Result*>& results);

/**
 * @brief Compute a Monte Carlo estimate of an element integral.
 *
 * Approximates the integral of the product of shape-function values and source
 * values over an element using uniformly distributed sample points.
 *
 * @param[in] shape_func_values_at_points Shape-function values at sample
 * points.
 * @param[in] src_values_at_points Source-field values at the same sample
 * points.
 * @param[in] npoints_each_tri Number of sample points.
 * @param[in] volume Element volume.
 * @return Estimated integral value.
 */
KOKKOS_INLINE_FUNCTION
double monte_carlo_integral(const Omega_h::Real* shape_func_values_at_points,
                            const Omega_h::Real* src_values_at_points,
                            const int npoints_each_tri,
                            const Omega_h::Real volume)
{

  if (npoints_each_tri <= 0)
    return 0.0;
  Omega_h::Real sum = 0;

  for (int i = 0; i < npoints_each_tri; ++i) {
    sum += shape_func_values_at_points[i] * src_values_at_points[i];
  }

  sum *= volume;
  return sum / npoints_each_tri;
}

/**
 * @brief Compute element-local load vectors using Monte Carlo integration.
 *
 * For each target triangle \f$\Omega_t\f$, this routine approximates the local
 * load-vector entries by uniform Monte Carlo sampling:
 * \f[
 * \widehat b_k^{\,t}
 * = \frac{|\Omega_t|}{N}\sum_{i=1}^N
 * f^s(\mathbf X_{t,i})\,\psi_k(\mathbf X_{t,i}),
 * \f]
 * where \f$N\f$ is the number of sample points in the element,
 * \f$f^s(\mathbf X_{t,i})\f$ is the source-field value at the sampled point,
 * and \f$\psi_k\f$ is the local target basis function.
 *
 * The sampling points are generated on the reference triangle in barycentric
 * coordinates and the same reference sample set is reused for all target
 * elements. For linear triangular elements, the barycentric coordinates are
 * equal to the local shape-function values, so they are used directly in the
 * Monte Carlo estimator.
 *
 * @param[in] target_mesh Target 2D triangular mesh.
 * @param[in] field_values_at_points Source-field values evaluated at the
 * sampled physical points, stored element-by-element.
 * @param[in] npoints_each_tri Number of Monte Carlo sample points per target
 *                             triangle.
 * @param[in] method Sampling method used to generate the reference barycentric
 *                   coordinates.
 * @param[in] sobol_filename File containing precomputed Sobol barycentric
 *                           samples when Sobol sampling is selected.
 *
 * @return Flattened element-local load vectors of size
 *         `3 * target_mesh.nelems()`.
 *
 * @note This routine computes element-local contributions only; it does not
 *       assemble a global load vector.
 * @note This routine assumes that `field_values_at_points` has already been
 *       evaluated at the sampled physical points.
 */

Kokkos::View<MeshField::Real*> loadVectorMCIntegrator(
  Omega_h::Mesh& target_mesh, const Omega_h::Reals& field_values_at_points,
  const int npoints_each_tri, SamplingMethod method,
  const std::string& sobol_filename)
{

  Kokkos::View<MeshField::Real* [3]> ref_barycentric_coords;
  if (method == SamplingMethod::SOBOL) {
    if (sobol_filename.empty()) {
                throw std::runtime_error("Could not open sobol sample file : ";
    }
    ref_barycentric_coords =
      read_sobol_barycentric_samples_from_file(sobol_filename);
  } else {
    ref_barycentric_coords =
      generate_uniform_random_barycentric_coords(npoints_each_tri);
  }

  int subVectorSize = 3;

  Omega_h::Reals elementsArea;
  elementsArea = Omega_h::measure_elements_real(&target_mesh);

  Kokkos::View<MeshField::Real*> elmLoadVector("elmLoadVector",
                                               target_mesh.nelems() * 3);

  Kokkos::parallel_for(
    "eval load vector using MC", target_mesh.nelems(),
    KOKKOS_LAMBDA(const int elm) {
      Omega_h::Real N0[MAX_POINTS] = {};
      Omega_h::Real N1[MAX_POINTS] = {};
      Omega_h::Real N2[MAX_POINTS] = {};
      Omega_h::Real src_field_values[MAX_POINTS] = {};

      int base_idx_src_values = elm * npoints_each_tri;

      for (int i = 0; i < npoints_each_tri; ++i) {
        N0[i] = ref_barycentric_coords(i, 0);
        N1[i] = ref_barycentric_coords(i, 1);
        N2[i] = sampled_barycentric_coords(i, 2);
        src_field_values[i] = field_values_at_points[base_idx_src_values + i];
      }

      Omega_h::Vector<3> result;
      result[0] = montecarlo_integration(N0, src_field_values, npoints_each_tri,
                                         elementsArea[elm]);
      result[1] = montecarlo_integration(N1, src_field_values, npoints_each_tri,
                                         elementsArea[elm]);
      result[2] = montecarlo_integration(N2, src_field_values, npoints_each_tri,
                                         elementsArea[elm]);

      for (int i = 0; i < 3; ++i) {
        elmLoadVector(elm * subVectorSize + i) = result[i];
      }
    });

  return elmLoadVector;
}
