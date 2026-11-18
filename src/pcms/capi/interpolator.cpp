//
// Created by hasanm4 on 2/17/25.
//
#include <pcms/capi/kokkos.h>
#include <pcms/capi/interpolator.h>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/backends/point_cloud/geometry.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/field/point_evaluator.hpp>
#include <pcms/backends/omega_h/entity_layout.hpp>
#include <pcms/backends/omega_h/omega_h.hpp>
#include <pcms/utility/arrays.hpp>
#include <pcms/utility/entity_types.hpp>
#include <pcms/utility/mesh_geometry.hpp>
#include <pcms/utility/print.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <cmath>
#include <memory>

namespace pcms
{
using backend::omega_h::CreateGeometry;
using backend::omega_h::OmegaHEntityLayout;
using utility::get_entity_centroids;
} // namespace pcms

namespace
{

/// What a PcmsInterpolatorHandle points to: a function on the source sites,
/// which pcms_interpolate fills from its input, and its MLS evaluator at the
/// target points.
struct MLSInterpolation
{
  pcms::Function<pcms::Real> source;
  std::unique_ptr<pcms::PointEvaluator<pcms::Real>> evaluator;
};

/// The support search of the C API: `radius` is the squared starting radius
/// (what this API has always taken), and the radius always adapts.
pcms::NeighbourhoodOptions Supports(double radius, int min_req_supports)
{
  return {.radius = std::sqrt(radius), .min_supports = min_req_supports};
}

PcmsInterpolatorHandle MakeHandle(
  const pcms::FunctionSpace& space,
  std::shared_ptr<const pcms::BoundSource> bound)
{
  return {new MLSInterpolation{
    space.CreateFunction<pcms::Real>(),
    space.CreatePointEvaluator<pcms::Real>(std::move(bound))}};
}

} // namespace

//[[nodiscard]]
PcmsInterpolatorHandle pcms_create_interpolator(PcmsOmegaHMeshHandle oh_mesh,
                                                double radius)
{
  auto& mesh = *reinterpret_cast<Omega_h::Mesh*>(oh_mesh.mesh_handle);
  OMEGA_H_CHECK_PRINTF(
    mesh.dim() == 2, "Only 2D meshes are supported but found %d\n", mesh.dim());
  const auto mesh_geometry =
    pcms::CreateGeometry(mesh, pcms::csys::Cartesian::Create(2));
  auto space = pcms::CreatePolynomialReconstructionSpace(
    pcms::backend::point_cloud::CreatePointCloudGeometry(
      std::make_shared<pcms::OmegaHEntityLayout>(mesh_geometry, pcms::Face, 1),
      Supports(radius, 10)));
  // the vertices of the same mesh, whose supports follow from its adjacency
  const pcms::OmegaHEntityLayout vertices(mesh_geometry, pcms::Vertex, 1);
  return MakeHandle(
    *space, space->GetGeometry()->Bind(vertices.GetDOFHolderCoordinates(),
                                       pcms::DOFHolderOwners(vertices)));
}

PcmsInterpolatorHandle pcms_create_point_based_interpolator(
  void* source_points, int source_points_size, void* target_points,
  int target_points_size, double radius, int degree, int min_req_supports,
  double lambda, double decay_factor)
{
  const auto system = pcms::csys::Cartesian::Create(2);
  auto space = pcms::CreatePolynomialReconstructionSpace(
    pcms::backend::point_cloud::CreatePointCloudGeometry(
      pcms::CoordinateView<pcms::HostMemorySpace>(
        system, pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
                  static_cast<const pcms::Real*>(source_points),
                  source_points_size / 2, 2)),
      Supports(radius, min_req_supports)),
    {.degree = static_cast<unsigned>(degree),
     .lambda = lambda,
     .decay_factor = decay_factor});
  const Kokkos::View<const pcms::Real**, Kokkos::LayoutRight,
                     pcms::HostMemorySpace, Kokkos::MemoryUnmanaged>
    targets_host(static_cast<const pcms::Real*>(target_points),
                 target_points_size / 2, 2);
  Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> targets(
    "pcms_interpolator_targets", targets_host.extent(0), 2);
  pcms::DeepCopyMismatchLayouts(targets, targets_host);
  return MakeHandle(*space, space->GetGeometry()->Bind(
                              pcms::CoordinateView<pcms::DeviceMemorySpace>(
                                system, pcms::MakeConstRank2View(targets))));
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
  delete static_cast<MLSInterpolation*>(interpolator.pointer);
}

void pcms_interpolate(PcmsInterpolatorHandle interpolator, void* input,
                      int input_size, void* output, int output_size)
{
  auto& mls = *static_cast<MLSInterpolation*>(interpolator.pointer);
  const pcms::LO num_sources = mls.source.GetLayout().GetNumDofHolder();
  const pcms::LO num_targets = mls.evaluator->NumPoints();
  OMEGA_H_CHECK_PRINTF(
    input_size == num_sources,
    "Input array size does not match the source size %d != %d\n", input_size,
    num_sources);
  OMEGA_H_CHECK_PRINTF(
    output_size == num_targets,
    "Output array size does not match the target size %d != %d\n", output_size,
    num_targets);

  const pcms::ValueBasis& basis = mls.source.GetData().GetValueBasis();
  mls.source.SetDOFHolderDataHost(
    pcms::ValueView<const pcms::Real, pcms::HostMemorySpace>(
      basis, pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
               static_cast<const pcms::Real*>(input), input_size, 1)));
  Kokkos::View<pcms::Real**, pcms::DeviceMemorySpace> values(
    "pcms_interpolate_values", output_size, 1);
  mls.evaluator->Evaluate(mls.source,
                          pcms::ValueView<pcms::Real, pcms::DeviceMemorySpace>(
                            basis, pcms::detail::MakeRank2View(values)));
  Kokkos::View<pcms::Real**, Kokkos::LayoutRight, pcms::HostMemorySpace,
               Kokkos::MemoryUnmanaged>
    output_values(static_cast<pcms::Real*>(output), output_size, 1);
  pcms::DeepCopyMismatchLayouts(output_values, values);
}
