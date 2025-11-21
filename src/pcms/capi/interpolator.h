/**
 * @file interpolator.h
 * @brief PCMS Interpolator C API Header
 *
 * This header defines the C API for the PCMS Interpolator library, which
 * provides functionality for creating and using interpolators based on Moving
 * Least Squares (MLS) methods. For now, it only supports 2D interpolation and
 * hard-coded for RBF, Gaussian weight function.
 */

#ifndef PCMS_INTERPOLATOR_CAPI_H
#define PCMS_INTERPOLATOR_CAPI_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Holds a void pointer of InterpolationBase object
 * @see InterpolationBase, MLSMeshInterpolation, MLSPointCloudInterpolation
 */
struct PcmsInterpolatorHandle
{
  void* pointer;
};

/**
 * @brief Typedef for PcmsInterpolatorHandle struct
 * @copydetails PcmsInterpolatorHandle
 * @see PcmsInterpolatorHandle, InterpolationBase, MLSMeshInterpolation,
 * MLSPointCloudInterpolation
 */
typedef struct PcmsInterpolatorHandle PcmsInterpolatorHandle;

/**
 * @brief Holds void pointers of Omega_h Mesh object and the library
 * (communicator) handle
 */
struct PcmsInterpolatorOHMeshHandle
{
  void* mesh_handle;
  void* lib_handle;
};

/**
 * @brief Typedef for PcmsInterpolatorOHMeshHandle struct
 * @copydetails PcmsInterpolatorOHMeshHandle
 * @see PcmsInterpolatorOHMeshHandle
 */
typedef struct PcmsInterpolatorOHMeshHandle PcmsInterpolatorOHMeshHandle;

/**
 * @brief Create centroid to node interpolator
 * @param oh_mesh Omega_h mesh handle
 * @param radius Starting radius for support search
 * @return PcmsInterpolatorHandle Handle to the created interpolator
 *
 * @details Creates an interpolator that maps values from element centroids to
 * mesh nodes It uses MLS interpolation and uses the radius to start adaptive
 * support search.
 *
 * @note Uses default parameters of MLSMeshInterpolation.
 * @see MLSMeshInterpolation
 */
PcmsInterpolatorHandle pcms_create_interpolator(
  PcmsInterpolatorOHMeshHandle oh_mesh, double radius);

/**
 * @brief Create 2D point-based MLS interpolator for RBF interpolation
 * @param source_points Void pointer to the first source point of double of size
 * source_points_size
 * @param source_points_size Size of source points array \f$
 * Number\ of\ source\ points \times 2 \f$
 * @param target_points Void pointer to the first target point of double of size
 * target_points_size
 * @param target_points_size Size of target points array \f$
 * Number\ of\ target\ points \times 2 \f$
 * @param radius Starting radius for support search (always adaptive)
 * @param degree Degree of the MLS basis functions
 * @param min_req_supports Minimum required supports for each target point
 * (maximum allowed supports is three times this value)
 * @param lambda Regularization parameter
 * @param decay_factor Decay factor for weight function
 * @return PcmsInterpolatorHandle Handle to the created interpolator
 *
 * @details To better select the interpolation parameters, look at their
 * explanation in mls_interpolation function documentation. Call
 * pcms_interpolate function with this interpolator to perform interpolation
 * from source points to target points. Remember to delete the created
 * interpolator using pcms_destroy_interpolator function after use to avoid
 * memory leaks. The user is responsible for ensuring that the source and target
 * points are sustained in memory during the call of this constructor. They are
 * no longer needed after the interpolator is created as they are copied
 * internally.
 *
 * @see mls_interpolation, InterpolationBase, MLSPointCloudInterpolation
 */
PcmsInterpolatorHandle pcms_create_point_based_interpolator(
  void* source_points, int source_points_size, void* target_points,
  int target_points_size, double radius, int degree, int min_req_supports,
  double lambda, double decay_factor);

/**
 * @brief Create 2D point-based interpolator to map from Degas2 mesh element
 * centroids (source) to XGC nodes (target)
 * @param target_points Void pointer to the first target point of double of size
 * target_points_size
 * @param target_points_size Size of target points array \f$
 * Number\ of\ target\ points \times 2 \f$
 * @param dg2_mesh_filename C-string of Degas2 mesh filename (Omega_h (.osh)
 * format)
 * @param radius Starting radius for support search (always adaptive)
 * @param dg2_elem_count Void pointer to integer to write the number of Degas2
 * elements
 * @param degree Degree of the MLS basis functions
 * @param min_req_supports Minimum required supports for each target point
 * (maximum allowed supports is three times this value)
 * @param lambda Regularization parameter of MLS
 * @param decay_factor Decay factor for weight function of MLS
 * @return PcmsInterpolatorHandle Handle to the created interpolator
 *
 * @note This is used to interpolate after the Degas2 step in a coupled
 * Degas2-XGC simulation.
 * @copydetails pcms_create_point_based_interpolator
 * @see mls_interpolation
 */
PcmsInterpolatorHandle pcms_create_degas2xgcnode_interpolator(
  void* target_points, int target_points_size, const char* dg2_mesh_filename,
  double radius, void* dg2_elem_count, int degree, int min_req_supports,
  double lambda, double decay_factor);

/**
 * @brief Create 2D point-based interpolator to map from XGC nodes (source) to
 * Degas2 mesh element centroids (target)
 * @param dg2_mesh_filename C-string of Degas2 mesh filename (Omega_h (.osh)
 * format)
 * @param source_points Void pointer to the first source point of double of size
 * source_points_size
 * @param source_points_size Size of source points array \f$
 * Number\ of\ source\ points \times 2 \f$
 * @param radius Starting radius for support search (always adaptive)
 * @param dg2_elem_count Void pointer to integer to write the number of Degas2
 * elements
 * @param degree Degree of the MLS basis functions
 * @param min_req_supports Minimum required supports for each target point
 * (maximum allowed supports is three times this value)
 * @param lambda Regularization parameter of MLS
 * @param decay_factor Decay factor for weight function of MLS
 * @return PcmsInterpolatorHandle Handle to the created interpolator
 *
 * @note This is used to interpolate before the Degas2 step in a coupled
 * Degas2-XGC simulation.
 * @copydetails pcms_create_point_based_interpolator
 * @see mls_interpolation
 */
PcmsInterpolatorHandle pcms_create_xgcnodedegas2_interpolator(
  const char* dg2_mesh_filename, void* source_points, int source_points_size,
  double radius, void* dg2_elem_count, int degree, int min_req_supports,
  double lambda, double decay_factor);

/**
 * @brief Destroy interpolator
 * @param interpolator Handle to the created interpolator
 *
 * @details Call this function to delete the created interpolator using
 * pcms_create_interpolator or pcms_create_point_based_interpolator functions
 * after use to avoid memory leaks.
 */
void pcms_destroy_interpolator(PcmsInterpolatorHandle interpolator);

/**
 * @brief Read Omega_h mesh from file
 * @param filename C-string of Omega_h mesh filename
 * @return PcmsInterpolatorOHMeshHandle Handle to the read Omega_h mesh
 *
 * @details Reads an Omega_h mesh from the given filename and returns a handle
 * to it. This handle can be used to create interpolators that operate on the
 * mesh.
 */
PcmsInterpolatorOHMeshHandle read_oh_mesh(const char* filename);

/**
 * @brief Release Omega_h mesh
 * @param oh_mesh_handle Handle to the Omega_h mesh to be released
 *
 * @details Releases the Omega_h mesh associated with the given handle to
 * free up resources.
 */
void release_oh_mesh(PcmsInterpolatorOHMeshHandle oh_mesh_handle);

/**
 * @brief Perform interpolation
 * @param interpolator Handle to the created interpolator
 * @param input Void pointer to the first element of the input data array of
 * type double and size of input_size
 * @param input_size Size of the input data array which is the same as the
 * number of source points/elements
 * @param output Void pointer to the first element of the output data array of
 * type double and size of output_size
 * @param output_size Size of the output data array which is the same as the
 * number of target points/elements
 *
 * @details Uses the given interpolator to perform interpolation using the input
 * data and writes the results to the output data array.
 *
 * @see pcms_create_interpolator, pcms_create_point_based_interpolator,
 * mls_interpolation, InterpolationBase
 */
void pcms_interpolate(PcmsInterpolatorHandle interpolator, void* input,
                      int input_size, void* output, int output_size);

#ifdef __cplusplus
}
#endif

#endif // PCMS_INTERPOLATOR_CAPI_H
