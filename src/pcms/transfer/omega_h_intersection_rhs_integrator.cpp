#include "pcms/transfer/omega_h_intersection_rhs_integrator.hpp"
#include "pcms/field/element_dispatch.h"
#include "pcms/transfer/omega_h_form_integrator_utils.hpp"
#include "pcms/transfer/petsc_utils.hpp"
#include "pcms/utility/assert.h"
#include <algorithm>
#include <Omega_h_for.hpp>
#include <Omega_h_shape.hpp>
#include <petscksp.h>

namespace pcms
{

namespace
{

// Results of the two-pass intersection kernel; the constructor moves these into
// the integrator's members.
struct Data
{
  Kokkos::View<Real**, DeviceMemorySpace> coords;
  Kokkos::View<PetscInt*, DeviceMemorySpace> node_gids;
  Kokkos::View<Real*, DeviceMemorySpace> coeffs;
  int ndof_per_elem = 0;
  PetscInt num_target_dofs = 0;
};

template <int TgtOrder>
Data BuildDataImpl(const OmegaHLagrangeLayout& source_layout,
                   const OmegaHLagrangeLayout& target_layout, int quad_order);

Data BuildData(const std::shared_ptr<const OmegaHLagrangeLayout>& source_layout,
               CoordinateSystem source_coordinate_system,
               const std::shared_ptr<const OmegaHLagrangeLayout>& target_layout,
               CoordinateSystem target_coordinate_system)
{
  detail::CheckOmegaHScalarLagrangeLayout(
    source_coordinate_system, source_layout, "OmegaHIntersectionRHSIntegrator",
    "source");
  detail::CheckOmegaHScalarLagrangeLayout(
    target_coordinate_system, target_layout, "OmegaHIntersectionRHSIntegrator",
    "target");

  // The integrand f_src * phi_target has polynomial degree source_order +
  // target_order on each intersection subtriangle; integrate it exactly (with a
  // 1-point floor so a P0->P0 pair still gets a valid rule).
  const int quad_order =
    std::max(1, source_layout->GetOrder() + target_layout->GetOrder());

  return detail::DispatchByOrder(target_layout->GetOrder(), [&](auto order_c) {
    constexpr int TgtOrder = decltype(order_c)::value;
    return BuildDataImpl<TgtOrder>(*source_layout, *target_layout, quad_order);
  });
}

template <int TgtOrder>
Data BuildDataImpl(const OmegaHLagrangeLayout& source_layout,
                   const OmegaHLagrangeLayout& target_layout, int quad_order)
{
  using Basis = detail::TargetTriBasis<TgtOrder>;
  constexpr int ndof = Basis::ndof;

  Omega_h::Mesh& source_mesh = source_layout.GetMesh();
  Omega_h::Mesh& target_mesh = target_layout.GetMesh();

  const auto intersections = intersectTargets(source_mesh, target_mesh);

  const auto& tgt_coords = target_mesh.coords();
  const auto& tgt_faces2nodes =
    target_mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  const auto& src_coords = source_mesh.coords();
  const auto& src_faces2nodes =
    source_mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;

  detail::IntegrationData ip_data(quad_order);
  const int npts = ip_data.size();
  auto bary_coords = ip_data.bary_coords; // device view
  auto weights = ip_data.weights;         // device view

  const auto global_to_local = target_layout.GetGlobalToLocalPermutation();

  const Omega_h::LOs tgt2src_offsets = intersections.tgt2src_offsets;
  const Omega_h::LOs tgt2src_indices = intersections.tgt2src_indices;
  const int nelems = target_mesh.nelems();

  // Pass 1: count integration points per target element.
  Omega_h::Write<Omega_h::LO> ip_counts(nelems, 0, "rhs_ip_counts");
  Kokkos::parallel_for(
    "rhs_count", nelems, KOKKOS_LAMBDA(int elm) {
      int count = 0;
      detail::ForEachIntersectionSubtriangle(
        elm, {tgt2src_offsets, tgt2src_indices}, tgt_coords, src_coords,
        tgt_faces2nodes, src_faces2nodes,
        [&](const Omega_h::Few<Omega_h::Vector<2>, 3>&,
            const r3d::Few<r3d::Vector<2>, 3>&,
            const r3d::Few<r3d::Vector<2>, 3>&, int,
            Omega_h::Real) { count += npts; });
      ip_counts[elm] = count;
    });
  Kokkos::fence();

  const auto ip_offsets = Omega_h::offset_scan(
    Omega_h::Read<Omega_h::LO>(ip_counts), "rhs_ip_offsets");
  const int num_pts = static_cast<int>(ip_offsets.last());

  // Pass 2: fill coords, node_gids, and coeffs on device. node_gids/coeffs hold
  // ndof (target DOFs per element) entries per integration point.
  Kokkos::View<Real**, DeviceMemorySpace> coords("rhs_coords", num_pts, 2);
  Kokkos::View<PetscInt*, DeviceMemorySpace> node_gids(
    "rhs_node_gids", static_cast<std::size_t>(num_pts) * ndof);
  Kokkos::View<Real*, DeviceMemorySpace> coeffs(
    "rhs_coeffs", static_cast<std::size_t>(num_pts) * ndof);

  Kokkos::parallel_for(
    "rhs_fill", nelems, KOKKOS_LAMBDA(int elm) {
      const auto tgt_verts = Omega_h::gather_verts<3>(tgt_faces2nodes, elm);
      const Omega_h::Matrix<2, 3> tgt_vert_mat =
        Omega_h::gather_vectors<3, 2>(tgt_coords, tgt_verts);
      Omega_h::Few<Omega_h::Vector<2>, 3> tgt_omh;
      for (int i = 0; i < 3; ++i)
        tgt_omh[i] = {tgt_vert_mat[i][0], tgt_vert_mat[i][1]};

      int ip_local = 0;
      const int offset = ip_offsets[elm];

      detail::ForEachIntersectionSubtriangle(
        elm, {tgt2src_offsets, tgt2src_indices}, tgt_coords, src_coords,
        tgt_faces2nodes, src_faces2nodes,
        [&](const Omega_h::Few<Omega_h::Vector<2>, 3>& tri,
            const r3d::Few<r3d::Vector<2>, 3>&,
            const r3d::Few<r3d::Vector<2>, 3>&, int, Omega_h::Real area) {
          for (int ip_idx = 0; ip_idx < npts; ++ip_idx) {
            const auto bary = bary_coords(ip_idx);
            const double w = weights(ip_idx);
            const auto pt = detail::GlobalFromBarycentric(bary, tri);

            Omega_h::Real basis[ndof];
            Basis::Values(pt, tgt_omh, basis);

            const int global_ip = offset + ip_local;
            coords(global_ip, 0) = pt[0];
            coords(global_ip, 1) = pt[1];
            for (int k = 0; k < ndof; ++k) {
              node_gids(global_ip * ndof + k) = static_cast<PetscInt>(
                Basis::Index(global_to_local, elm, tgt_verts, k));
              coeffs(global_ip * ndof + k) = basis[k] * w * 2.0 * area;
            }
            ++ip_local;
          }
        });
    });
  Kokkos::fence();

  Data d;
  d.coords = std::move(coords);
  d.node_gids = std::move(node_gids);
  d.coeffs = std::move(coeffs);
  d.ndof_per_elem = ndof;
  d.num_target_dofs =
    static_cast<PetscInt>(target_layout.GetNumOwnedDofHolder());
  return d;
}

} // namespace

// ---------------------------------------------------------------------------
// OmegaHIntersectionRHSIntegrator constructors
// ---------------------------------------------------------------------------

OmegaHIntersectionRHSIntegrator::OmegaHIntersectionRHSIntegrator(
  const FunctionSpace& source_space, const FunctionSpace& target_space)
  : OmegaHIntersectionRHSIntegrator(
      std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(
        source_space.GetLayout()),
      source_space.GetCoordinateSystem(),
      std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(
        target_space.GetLayout()),
      target_space.GetCoordinateSystem())
{
}

OmegaHIntersectionRHSIntegrator::OmegaHIntersectionRHSIntegrator(
  std::shared_ptr<const OmegaHLagrangeLayout> source_layout,
  CoordinateSystem source_coordinate_system,
  std::shared_ptr<const OmegaHLagrangeLayout> target_layout,
  CoordinateSystem target_coordinate_system)
{
  Data data = BuildData(source_layout, source_coordinate_system, target_layout,
                        target_coordinate_system);
  point_set_ = IntegrationPointSet<DeviceMemorySpace>(
    CoordinateSystem::Cartesian, std::move(data.coords));
  node_gids_ = std::move(data.node_gids);
  coeffs_ = std::move(data.coeffs);
  ndof_per_elem_ = data.ndof_per_elem;

  const PetscInt nnz = static_cast<PetscInt>(node_gids_.extent(0));
  PetscErrorCode ierr =
    createSeqVec(PETSC_COMM_WORLD, data.num_target_dofs, &vec_);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
  // VecSetPreallocationCOO takes the COO indices on the host
  // TODO: ask Todd/PETSc folks if there is a better way to do this for GPU
  // support
  auto node_gids_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, node_gids_);
  ierr = VecSetPreallocationCOO(vec_, nnz, node_gids_host.data());
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

OmegaHIntersectionRHSIntegrator::~OmegaHIntersectionRHSIntegrator()
{
  if (vec_) {
    VecDestroy(&vec_);
  }
}

// ---------------------------------------------------------------------------
// OmegaHIntersectionRHSIntegrator public interface
// ---------------------------------------------------------------------------

const IntegrationPointSet<DeviceMemorySpace>&
OmegaHIntersectionRHSIntegrator::GetIntegrationPoints() const noexcept
{
  return point_set_;
}

Vec OmegaHIntersectionRHSIntegrator::GetVector() const noexcept
{
  return vec_;
}

void OmegaHIntersectionRHSIntegrator::Assemble(
  Rank2View<const Real, DeviceMemorySpace> sampled_values)
{
  const int ndof = ndof_per_elem_;
  const std::size_t num_pts =
    static_cast<std::size_t>(node_gids_.extent(0)) / ndof;
  PCMS_ALWAYS_ASSERT(static_cast<std::size_t>(sampled_values.extent(0)) ==
                     num_pts);
  PCMS_ALWAYS_ASSERT(sampled_values.extent(1) >= 1);

  PetscErrorCode ierr = VecZeroEntries(vec_);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);

  auto sv = Kokkos::View<const Real**, Kokkos::LayoutRight, DeviceMemorySpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
    sampled_values.data_handle(), sampled_values.extent(0),
    sampled_values.extent(1));
  Kokkos::View<PetscScalar*, DeviceMemorySpace> coo_vals("rhs_coo_vals",
                                                         num_pts * ndof);
  auto coeffs = coeffs_;
  Kokkos::parallel_for(
    "rhs_coo_vals", static_cast<int>(num_pts), KOKKOS_LAMBDA(int i) {
      const PetscScalar f = static_cast<PetscScalar>(sv(i, 0));
      for (int j = 0; j < ndof; ++j) {
        coo_vals(i * ndof + j) =
          static_cast<PetscScalar>(coeffs(i * ndof + j)) * f;
      }
    });

  ierr = VecSetValuesCOO(vec_, coo_vals.data(), ADD_VALUES);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

// ---------------------------------------------------------------------------
// Builder
// ---------------------------------------------------------------------------

std::unique_ptr<LinearFormIntegrator> BuildOmegaHConservativeRHSIntegrator(
  const FunctionSpace& source_space, const FunctionSpace& target_space)
{
  return std::make_unique<OmegaHIntersectionRHSIntegrator>(source_space,
                                                           target_space);
}

} // namespace pcms
