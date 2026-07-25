//
// Created by hasanm4 on 2/17/25.
//
#include <pcms/capi/kokkos.h>
#include <pcms/capi/interpolator.h>
#include <pcms/transfer/interpolation_base.h>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <pcms/utility/mesh_geometry.h>
#include <pcms/utility/print.h>

//[[nodiscard]]
PcmsInterpolatorHandle pcms_create_interpolator(PcmsOmegaHMeshHandle oh_mesh,
                                                double radius)
{
  auto* source_mesh = reinterpret_cast<Omega_h::Mesh*>(oh_mesh.mesh_handle);
  auto* interpolator = new pcms::MLSMeshInterpolation(*source_mesh, radius);
  return {reinterpret_cast<void*>(interpolator)};
}

PcmsInterpolatorHandle pcms_create_point_based_interpolator(
  void* source_points, int source_points_size, void* target_points,
  int target_points_size, double radius, int degree, int min_req_supports,
  double lambda, double decay_factor)
{

  auto source_points_view = pcms::Rank1View<double, pcms::HostMemorySpace>(
    reinterpret_cast<double*>(source_points), source_points_size);
  auto target_points_view = pcms::Rank1View<double, pcms::HostMemorySpace>(
    reinterpret_cast<double*>(target_points), target_points_size);
  auto* interpolator = new pcms::MLSPointCloudInterpolation(
    source_points_view, target_points_view, 2, radius, min_req_supports, degree,
    true, lambda, decay_factor);
  return {reinterpret_cast<void*>(interpolator)};
}

Omega_h::HostRead<Omega_h::Real> read_mesh_centroids(const char* mesh_filename,
                                                     int& num_elements)
{
  auto fname = std::string(mesh_filename);
  fname = fname.erase(fname.find_last_not_of(" \n\r\t") + 1);
  pcms::printInfo("The interpolator got dg2 mesh file: %s\n", fname.c_str());
  auto mesh_lib = Omega_h::Library(nullptr, nullptr, MPI_COMM_SELF);
  auto mesh = Omega_h::binary::read(fname, mesh_lib.world());
  OMEGA_H_CHECK_PRINTF(mesh.dim() == 2, "Mesh dimension is not 2D %d\n",
                       mesh.dim());
  auto elem_centroids = pcms::get_entity_centroids(mesh, Omega_h::FACE);
  num_elements = mesh.nelems();
  OMEGA_H_CHECK_PRINTF(num_elements * 2 == elem_centroids.size(),
                       "Mesh element centroids size does not match the number "
                       "of elements %d != %d\n",
                       num_elements * 2, elem_centroids.size());

  pcms::printInfo("Number of element centroids: %d\n",
                  elem_centroids.size() / 2);

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

PcmsInterpolatorHandle pcms_create_degas2xgcnode_interpolator(
  void* target_points, int target_points_size, const char* dg2_mesh_filename,
  double radius, void* dg2_elem_count, int degree, int min_req_supports,
  double lambda, double decay_factor)
{
  // same as above pcms_create_degas2xgc_interpolator but the target points are
  // provided by the user this is useful when the corresponding xgc mesh is not
  // available

  int dg2_num_elems = 0;
  Omega_h::HostRead<Omega_h::Real> dg2_elem_centroids_host =
    read_mesh_centroids(dg2_mesh_filename, dg2_num_elems);
  write_void_int_pointer(dg2_elem_count, dg2_num_elems);

  return pcms_create_point_based_interpolator(
    (void*)dg2_elem_centroids_host.data(), dg2_elem_centroids_host.size(),
    target_points, target_points_size, radius, degree, min_req_supports, lambda,
    decay_factor);
}

PcmsInterpolatorHandle pcms_create_xgcnodedegas2_interpolator(
  const char* dg2_mesh_filename, void* source_points, int source_points_size,
  double radius, void* dg2_elem_count, int degree, int min_req_supports,
  double lambda, double decay_factor)
{
  int dg2_num_elems = 0;
  Omega_h::HostRead<Omega_h::Real> dg2_elem_centroids_host =
    read_mesh_centroids(dg2_mesh_filename, dg2_num_elems);
  write_void_int_pointer(dg2_elem_count, dg2_num_elems);

  return pcms_create_point_based_interpolator(
    source_points, source_points_size, (void*)dg2_elem_centroids_host.data(),
    dg2_elem_centroids_host.size(), radius, degree, min_req_supports, lambda,
    decay_factor);
}

void pcms_destroy_interpolator(PcmsInterpolatorHandle interpolator)
{
  if (interpolator.pointer != nullptr) {
    delete reinterpret_cast<pcms::InterpolationBase*>(interpolator.pointer);
  }
}

void pcms_interpolate(PcmsInterpolatorHandle interpolator, void* input,
                      int input_size, void* output, int output_size)
{
  auto* mls_interpolator =
    reinterpret_cast<pcms::InterpolationBase*>(interpolator.pointer);

  OMEGA_H_CHECK_PRINTF(
    input_size == mls_interpolator->getSourceSize(),
    "Input array size does not match the source size %d != %zu\n", input_size,
    mls_interpolator->getSourceSize());
  OMEGA_H_CHECK_PRINTF(
    output_size == mls_interpolator->getTargetSize(),
    "Output array size does not match the target size %d != %zu\n", output_size,
    mls_interpolator->getTargetSize());

  pcms::Rank1View<double, pcms::HostMemorySpace> input_array(
    reinterpret_cast<double*>(input), input_size);
  pcms::Rank1View<double, pcms::HostMemorySpace> output_array(
    reinterpret_cast<double*>(output), output_size);

  mls_interpolator->eval(input_array, output_array);
}

// ---------------------------------------------------------------------------
// Conservative Projection (mesh-intersection-based)
// ---------------------------------------------------------------------------

#if defined(PCMS_ENABLE_PETSC) && defined(PCMS_ENABLE_MESHFIELDS)
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <pcms/field/field.h>
#include <pcms/field/function_space/lagrange.h>
#include <pcms/transfer/omega_h_conservative_projection.hpp>
#include <pcms/utility/arrays.h>
#include <cstring>
#include <memory>
#include <string>

namespace
{

// Trims trailing whitespace from a C string (Fortran strings may be padded).
std::string trim_filename(const char* filename)
{
  auto fname = std::string(filename);
  fname.erase(fname.find_last_not_of(" \n\r\t") + 1);
  return fname;
}

struct ConservativeProjectionContext
{
  // Library must be declared before meshes (Omega_h::Mesh holds a pointer to it).
  Omega_h::Library library;
  Omega_h::Mesh source_mesh;
  Omega_h::Mesh target_mesh;
  pcms::LagrangeFunctionSpace source_space;
  pcms::LagrangeFunctionSpace target_space;
  std::unique_ptr<pcms::OmegaHConservativeProjection> projection;
  pcms::Field<pcms::Real> source_field;
  pcms::Field<pcms::Real> target_field;

  ConservativeProjectionContext(const char* source_mesh_name,
                                const char* target_mesh_name, int src_order,
                                int tgt_order)
    : library(nullptr, nullptr, MPI_COMM_SELF),
      source_mesh(Omega_h::binary::read(trim_filename(source_mesh_name),
                                        library.world())),
      target_mesh(Omega_h::binary::read(trim_filename(target_mesh_name),
                                        library.world())),
      source_space(pcms::LagrangeFunctionSpace::FromMesh(
        source_mesh, src_order, 1, pcms::CoordinateSystem::Cartesian, "global",
        pcms::LagrangeFunctionSpace::Backend::OmegaH)),
      target_space(pcms::LagrangeFunctionSpace::FromMesh(
        target_mesh, tgt_order, 1, pcms::CoordinateSystem::Cartesian, "global",
        pcms::LagrangeFunctionSpace::Backend::OmegaH)),
      projection(std::make_unique<pcms::OmegaHConservativeProjection>(
        source_space, target_space)),
      source_field(source_space.CreateField<pcms::Real>()),
      target_field(target_space.CreateField<pcms::Real>())
  {
  }
};

} // namespace

PcmsConservativeProjectionHandle pcms_create_conservative_projection(
  const char* source_mesh_name, int source_order,
  const char* target_mesh_name, int target_order)
{
  auto* ctx = new ConservativeProjectionContext(source_mesh_name,
                                                target_mesh_name, source_order,
                                                target_order);
  return {reinterpret_cast<void*>(ctx)};
}

int pcms_conservative_projection_get_source_size(
  PcmsConservativeProjectionHandle projection)
{
  auto* ctx =
    reinterpret_cast<ConservativeProjectionContext*>(projection.pointer);
  return ctx->source_space.GetLayout()->GetNumOwnedDofHolder();
}

int pcms_conservative_projection_get_target_size(
  PcmsConservativeProjectionHandle projection)
{
  auto* ctx =
    reinterpret_cast<ConservativeProjectionContext*>(projection.pointer);
  return ctx->target_space.GetLayout()->GetNumOwnedDofHolder();
}

void pcms_conservative_projection_apply(
  PcmsConservativeProjectionHandle projection, void* source_data,
  int source_size, void* target_data, int target_size)
{
  auto* ctx =
    reinterpret_cast<ConservativeProjectionContext*>(projection.pointer);

  // Copy source data into internal field
  auto source_view = pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
    reinterpret_cast<pcms::Real*>(source_data), source_size, 1);
  ctx->source_field.SetDOFHolderDataHost(source_view);

  // Apply conservative projection
  ctx->projection->Apply(ctx->source_field, ctx->target_field);

  // Copy result back to user buffer
  auto target_view = ctx->target_field.GetDOFHolderDataHost();
  auto flat = pcms::FlattenToRank1View(target_view);
  std::memcpy(target_data, flat.data_handle(),
              static_cast<std::size_t>(target_size) * sizeof(double));
}

void pcms_destroy_conservative_projection(
  PcmsConservativeProjectionHandle projection)
{
  if (projection.pointer != nullptr) {
    delete reinterpret_cast<ConservativeProjectionContext*>(
      projection.pointer);
  }
}

#else // !(PCMS_ENABLE_PETSC && PCMS_ENABLE_MESHFIELDS)

// Stub implementations when PETSc or MeshFields are unavailable

PcmsConservativeProjectionHandle pcms_create_conservative_projection(
  const char*, int, const char*, int)
{
  pcms::printError("Conservative projection requires PCMS_ENABLE_PETSC and "
                   "PCMS_ENABLE_MESHFIELDS\n");
  return {nullptr};
}

int pcms_conservative_projection_get_source_size(
  PcmsConservativeProjectionHandle)
{
  return 0;
}

int pcms_conservative_projection_get_target_size(
  PcmsConservativeProjectionHandle)
{
  return 0;
}

void pcms_conservative_projection_apply(PcmsConservativeProjectionHandle,
                                        void*, int, void*, int)
{
  pcms::printError("Conservative projection requires PCMS_ENABLE_PETSC and "
                   "PCMS_ENABLE_MESHFIELDS\n");
}

void pcms_destroy_conservative_projection(PcmsConservativeProjectionHandle) {}

#endif
