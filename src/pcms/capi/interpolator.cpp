//
// Created by hasanm4 on 2/17/25.
//
#include <pcms/capi/kokkos.h>
#include "interpolator.h"
#include <interpolation_base.h>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <pcms/print.h>

//[[nodiscard]]
PcmsInterpolatorHandle pcms_create_interpolator(
  PcmsInterpolatorOHMeshHandle oh_mesh, double radius)
{
  auto* source_mesh = reinterpret_cast<Omega_h::Mesh*>(oh_mesh.mesh_handle);
  auto* interpolator = new MLSInterpolationHandler(*source_mesh, radius);
  return {reinterpret_cast<void*>(interpolator)};
}

PcmsPointBasedInterpolatorHandle pcms_create_point_based_interpolator(
  void* source_points, int source_points_size, void* target_points,
  int target_points_size, double radius)
{

  auto source_points_view = pcms::Rank1View<double, pcms::HostMemorySpace>(
    reinterpret_cast<double*>(source_points), source_points_size);
  auto target_points_view = pcms::Rank1View<double, pcms::HostMemorySpace>(
    reinterpret_cast<double*>(target_points), target_points_size);
  auto* interpolator = new MLSPointCloudInterpolation(
    source_points_view, target_points_view, 2, radius, 15, 3, true);
  return {reinterpret_cast<void*>(interpolator)};
}

PcmsPointBasedInterpolatorHandle pcms_create_degas2xgc_interpolator(
  const char* xgc_mesh_filename, const char* dg2_mesh_filename, double radius)
{
  // use the xgc_mesh_nodes as source points and
  // dg2 element centroids as target points
  // then, create a point based interpolator like above
  auto xgc_fname = std::string(xgc_mesh_filename);
  auto dg2_fname = std::string(dg2_mesh_filename);

  // trim the filenames since they are coming from c or fortran api which may
  // have extra spaces at the end
  xgc_fname = xgc_fname.erase(xgc_fname.find_last_not_of(" \n\r\t") + 1);
  dg2_fname = dg2_fname.erase(dg2_fname.find_last_not_of(" \n\r\t") + 1);
  pcms::printInfo(
    "The interpolator got xgc mesh file: %s and dg2 mesh file: %s\n",
    xgc_fname.c_str(), dg2_fname.c_str());

  // read the meshes
  auto xgc_mesh_lib = Omega_h::Library();
  auto xgc_mesh = Omega_h::binary::read(xgc_fname, xgc_mesh_lib.world());
  auto dg2_mesh_lib = Omega_h::Library();
  auto dg2_mesh = Omega_h::binary::read(dg2_fname, dg2_mesh_lib.world());

  auto xgc_nodes = xgc_mesh.coords();
  auto xgc_num_nodes = xgc_mesh.nverts();
  OMEGA_H_CHECK_PRINTF(xgc_mesh.dim() == 2, "XGC mesh dimension is not 2D %d\n",
                       xgc_mesh.dim());
  OMEGA_H_CHECK_PRINTF(
    xgc_num_nodes * 2 == xgc_nodes.size(),
    "XGC mesh nodes size does not match the number of vertices %d != %d\n",
    xgc_num_nodes * 2, xgc_nodes.size());
  OMEGA_H_CHECK_PRINTF(dg2_mesh.dim() == 2, "DG2 mesh dimension is not 2D %d\n",
                       dg2_mesh.dim());

  auto dg2_num_elems = dg2_mesh.nelems();
  auto dg2_elem_centroids = getCentroids(dg2_mesh);
  OMEGA_H_CHECK_PRINTF(dg2_num_elems * 2 == dg2_elem_centroids.size(),
                       "DG2 mesh element centroids size does not match the "
                       "number of elements %d != %d\n",
                       dg2_num_elems * 2, dg2_elem_centroids.size());

  Omega_h::HostRead<Omega_h::Real> xgc_nodes_host(xgc_nodes);
  Omega_h::HostRead<Omega_h::Real> dg2_elem_centroids_host(dg2_elem_centroids);

  return pcms_create_point_based_interpolator(
    (void*)dg2_elem_centroids_host.data(), dg2_elem_centroids.size(),
    (void*)xgc_nodes_host.data(), xgc_nodes.size(), radius);
}

Omega_h::HostRead<Omega_h::Real> read_mesh_centroids(const char* mesh_filename,
                                                     int& num_elements)
{
  auto fname = std::string(mesh_filename);
  fname = fname.erase(fname.find_last_not_of(" \n\r\t") + 1);
  pcms::printInfo("The interpolator got dg2 mesh file: %s\n", fname.c_str());
  auto mesh_lib = Omega_h::Library(nullptr, nullptr, MPI_COMM_SELF);
  auto mesh = Omega_h::binary::read(fname, mesh_lib.world());
  auto elem_centroids = getCentroids(mesh);
  num_elements = mesh.nelems();
  OMEGA_H_CHECK_PRINTF(num_elements * 2 == elem_centroids.size(),
                       "Mesh element centroids size does not match the number "
                       "of elements %d != %d\n",
                       num_elements * 2, elem_centroids.size());

  pcms::printInfo("Number of element centroids: %d\n",
                  elem_centroids.size() / 2);
  OMEGA_H_CHECK_PRINTF(mesh.dim() == 2, "Mesh dimension is not 2D %d\n",
                       mesh.dim());

  return {elem_centroids};
}

void write_void_int_pointer(void* pointer, int value)
{
  if (pointer) {
    int* dg2_elem_count_int = reinterpret_cast<int*>(pointer);
    *dg2_elem_count_int = value;
  } else {
    pcms::printError("Error: NULL pointer provided to write integer value\n");
  }
}

PcmsPointBasedInterpolatorHandle pcms_create_degas2xgcnode_interpolator(
  void* target_points, int target_points_size, const char* dg2_mesh_filename,
  double radius, void* dg2_elem_count)
{
  // same as above pcms_create_degas2xgc_interpolator but the target points are
  // provided by the user this is useful when the corresponding xgc mesh is not
  // available

  int dg2_num_elems = 0;
  auto dg2_elem_centroids_host =
    read_mesh_centroids(dg2_mesh_filename, dg2_num_elems);
  write_void_int_pointer(dg2_elem_count, dg2_num_elems);

  return pcms_create_point_based_interpolator(
    (void*)dg2_elem_centroids_host.data(), dg2_elem_centroids_host.size(),
    target_points, target_points_size, radius);
}

PcmsPointBasedInterpolatorHandle pcms_create_xgcnodedegas2_interpolator(
  const char* dg2_mesh_filename, void* source_points, int source_points_size,
  double radius, void* dg2_elem_count)
{
  int dg2_num_elems = 0;
  auto dg2_elem_centroids_host =
    read_mesh_centroids(dg2_mesh_filename, dg2_num_elems);
  write_void_int_pointer(dg2_elem_count, dg2_num_elems);

  return pcms_create_point_based_interpolator(
    source_points, source_points_size, (void*)dg2_elem_centroids_host.data(),
    dg2_elem_centroids_host.size(), radius);
}

void pcms_destroy_point_based_interpolator(
  PcmsPointBasedInterpolatorHandle interpolator)
{
  if (interpolator.pointer != nullptr) {
    delete reinterpret_cast<MLSPointCloudInterpolation*>(interpolator.pointer);
  }
}

void pcms_destroy_interpolator(PcmsInterpolatorHandle interpolator)
{
  if (interpolator.pointer != nullptr) {
    delete reinterpret_cast<MLSInterpolationHandler*>(interpolator.pointer);
  }
}

PcmsInterpolatorOHMeshHandle read_oh_mesh(const char* filename)
{
  auto fname = std::string(filename);
  // trim the filename since it is coming from c or fortran api which may have
  // extra spaces at the end
  fname.erase(fname.find_last_not_of(" \n\r\t") + 1);
  auto* oh_lib = new Omega_h::Library();
  auto* mesh = new Omega_h::Mesh(Omega_h::binary::read(fname, oh_lib->world()));

  return {reinterpret_cast<void*>(mesh), reinterpret_cast<void*>(oh_lib)};
}

void release_oh_mesh(PcmsInterpolatorOHMeshHandle oh_mesh)
{
  if (oh_mesh.mesh_handle != nullptr) {
    assert(oh_mesh.lib_handle != nullptr);
    delete reinterpret_cast<Omega_h::Mesh*>(oh_mesh.mesh_handle);
    delete reinterpret_cast<Omega_h::Library*>(oh_mesh.lib_handle);
  }
}

void pcms_interpolate(PcmsInterpolatorHandle interpolator, void* input,
                      int input_size, void* output, int output_size)
{
  auto* mls_interpolator =
    reinterpret_cast<MLSInterpolationHandler*>(interpolator.pointer);

  OMEGA_H_CHECK_PRINTF(
    input_size == mls_interpolator->getSourceSize(),
    "Input array size does not match the source size %d != %d\n", input_size,
    mls_interpolator->getSourceSize());
  OMEGA_H_CHECK_PRINTF(
    output_size == mls_interpolator->getTargetSize(),
    "Output array size does not match the target size %d != %d\n", output_size,
    mls_interpolator->getTargetSize());

  pcms::Rank1View<double, pcms::HostMemorySpace> input_array(
    reinterpret_cast<double*>(input), input_size);
  pcms::Rank1View<double, pcms::HostMemorySpace> output_array(
    reinterpret_cast<double*>(output), output_size);

  mls_interpolator->eval(input_array, output_array);
}

void pcms_interpolate_point_based(PcmsPointBasedInterpolatorHandle interpolator,
                                  void* input, int input_size, void* output,
                                  int output_size)
{
  auto* mls_interpolator =
    reinterpret_cast<MLSPointCloudInterpolation*>(interpolator.pointer);

  OMEGA_H_CHECK_PRINTF(
    input_size == mls_interpolator->getSourceSize(),
    "Input array size does not match the source size %d != %d\n", input_size,
    mls_interpolator->getSourceSize());
  OMEGA_H_CHECK_PRINTF(
    output_size == mls_interpolator->getTargetSize(),
    "Output array size does not match the target size %d != %d\n", output_size,
    mls_interpolator->getTargetSize());

  pcms::Rank1View<double, pcms::HostMemorySpace> input_array(
    reinterpret_cast<double*>(input), input_size);
  pcms::Rank1View<double, pcms::HostMemorySpace> output_array(
    reinterpret_cast<double*>(output), output_size);

  mls_interpolator->eval(input_array, output_array);
}
