#include "pcms/transfer/omega_h_intersection_rhs_integrator.hpp"
#include "pcms/field/element_dispatch.hpp"
#include "pcms/field/evaluator/omega_h_lagrange.hpp"
#include "pcms/field/function_space/lagrange.hpp"
#include "pcms/transfer/omega_h_form_integrator_utils.hpp"
#include "pcms/transfer/petsc_utils.hpp"
#include "pcms/utility/assert.hpp"
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
  std::shared_ptr<const OmegaHGeometry> source_geometry;
  Kokkos::View<Real**, DeviceMemorySpace> coords;
  Kokkos::View<LO*, DeviceMemorySpace> source_elements;
  Kokkos::View<PetscInt*, DeviceMemorySpace> node_gids;
  Kokkos::View<Real*, DeviceMemorySpace> coeffs;
  int ndof_per_elem = 0;
  PetscInt num_target_dofs = 0;
};

template <int Dim, int TgtOrder>
Data BuildDataImpl(const OmegaHLagrangeLayout& target_layout, int quad_order,
                   const OmegaHGeometry& source_geometry,
                   const OmegaHGeometry& target_geometry);

Data BuildData(const FunctionSpace& source_space,
               const FunctionSpace& target_space)
{
  // the intersection locates the target centroids with the source geometry
  const auto [source_layout, source_geometry] =
    detail::ResolveOmegaHScalarSpace(
      source_space, 0, 1, "OmegaHIntersectionRHSIntegrator", "source");
  const auto [target_layout, target_geometry] =
    detail::ResolveOmegaHScalarSpace(
      target_space, 0, 1, "OmegaHIntersectionRHSIntegrator", "target");

  const int dim = GetMesh(*target_geometry).dim();
  if (GetMesh(*source_geometry).dim() != dim) {
    throw pcms_error("OmegaHIntersectionRHSIntegrator: source and target mesh "
                     "dimensions differ");
  }

  // The integrand f_src * phi_target has polynomial degree source_order +
  // target_order on each intersection sub-simplex; integrate it exactly (with a
  // 1-point floor so a P0->P0 pair still gets a valid rule).
  const int quad_order =
    std::max(1, source_layout->GetOrder() + target_layout->GetOrder());

  Data data =
    detail::DispatchByOrder(target_layout->GetOrder(), [&](auto order_c) {
      constexpr int TgtOrder = decltype(order_c)::value;
      if (dim == 3) {
        return BuildDataImpl<3, TgtOrder>(*target_layout, quad_order,
                                          *source_geometry, *target_geometry);
      }
      return BuildDataImpl<2, TgtOrder>(*target_layout, quad_order,
                                        *source_geometry, *target_geometry);
    });
  data.source_geometry = source_geometry;
  return data;
}

template <int Dim, int TgtOrder>
Data BuildDataImpl(const OmegaHLagrangeLayout& target_layout, int quad_order,
                   const OmegaHGeometry& source_geometry,
                   const OmegaHGeometry& target_geometry)
{
  using Basis = detail::TargetSimplexBasis<Dim, TgtOrder>;
  constexpr int ndof = Basis::ndof;
  // Reference-to-physical Jacobian factor for a simplex: the reference simplex
  // measure is 1/Dim! (1/2 in 2D, 1/6 in 3D), so a physical sub-simplex of
  // measure `m` scales the reference quadrature weights by Dim! * m.
  constexpr Omega_h::Real ref_factor = (Dim == 3) ? 6.0 : 2.0;

  Omega_h::Mesh& source_mesh = GetMesh(source_geometry);
  Omega_h::Mesh& target_mesh = GetMesh(target_geometry);

  const auto intersections = intersectTargets(source_geometry, target_mesh);

  const auto& tgt_coords = target_mesh.coords();
  const auto& tgt_elems2nodes = target_mesh.ask_down(Dim, Omega_h::VERT).ab2b;
  const auto& src_coords = source_mesh.coords();
  const auto& src_elems2nodes = source_mesh.ask_down(Dim, Omega_h::VERT).ab2b;

  detail::IntegrationData<Dim> ip_data(quad_order);
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
      detail::ForEachIntersectionSubsimplex<Dim>(
        elm, {tgt2src_offsets, tgt2src_indices}, tgt_coords, src_coords,
        tgt_elems2nodes, src_elems2nodes,
        [&](const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>&, int,
            Omega_h::Real) { count += npts; });
      ip_counts[elm] = count;
    });
  Kokkos::fence();

  const auto ip_offsets = Omega_h::offset_scan(
    Omega_h::Read<Omega_h::LO>(ip_counts), "rhs_ip_offsets");
  const int num_pts = static_cast<int>(ip_offsets.last());

  // Pass 2: fill the coordinates and source element of each integration
  // point, node_gids, and coeffs on device. node_gids/coeffs hold ndof (target
  // DOFs per element) entries per integration point.
  Kokkos::View<Real**, DeviceMemorySpace> coords("rhs_coords", num_pts, Dim);
  Kokkos::View<LO*, DeviceMemorySpace> source_elements("rhs_source_elements",
                                                       num_pts);
  Kokkos::View<PetscInt*, DeviceMemorySpace> node_gids(
    "rhs_node_gids", static_cast<std::size_t>(num_pts) * ndof);
  Kokkos::View<Real*, DeviceMemorySpace> coeffs(
    "rhs_coeffs", static_cast<std::size_t>(num_pts) * ndof);

  Kokkos::parallel_for(
    "rhs_fill", nelems, KOKKOS_LAMBDA(int elm) {
      const auto tgt_verts =
        Omega_h::gather_verts<Dim + 1>(tgt_elems2nodes, elm);
      const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1> tgt_omh =
        Omega_h::gather_vectors<Dim + 1, Dim>(tgt_coords, tgt_verts);

      int ip_local = 0;
      const int offset = ip_offsets[elm];

      detail::ForEachIntersectionSubsimplex<Dim>(
        elm, {tgt2src_offsets, tgt2src_indices}, tgt_coords, src_coords,
        tgt_elems2nodes, src_elems2nodes,
        [&](const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>& sub, int src_elm,
            Omega_h::Real measure) {
          for (int ip_idx = 0; ip_idx < npts; ++ip_idx) {
            Omega_h::Vector<Dim + 1> bary;
            for (int d = 0; d < Dim + 1; ++d) {
              bary[d] = bary_coords(ip_idx, d);
            }
            const double w = weights(ip_idx);
            const auto pt = detail::GlobalFromBarycentric<Dim>(bary, sub);

            Omega_h::Real basis[ndof];
            Basis::Values(pt, tgt_omh, basis);

            const int global_ip = offset + ip_local;
            for (int d = 0; d < Dim; ++d) {
              coords(global_ip, d) = pt[d];
            }
            // the source element the point lies in, so the source is bound
            // to it without a search
            source_elements(global_ip) = src_elm;
            for (int k = 0; k < ndof; ++k) {
              node_gids(global_ip * ndof + k) = static_cast<PetscInt>(
                Basis::Index(global_to_local, elm, tgt_verts, k));
              coeffs(global_ip * ndof + k) =
                basis[k] * w * ref_factor * measure;
            }
            ++ip_local;
          }
        });
    });
  Kokkos::fence();

  Data d;
  d.coords = std::move(coords);
  d.source_elements = std::move(source_elements);
  d.node_gids = std::move(node_gids);
  d.coeffs = std::move(coeffs);
  d.ndof_per_elem = ndof;
  d.num_target_dofs =
    static_cast<PetscInt>(target_layout.GetNumDofHolderWithGhosts());
  return d;
}

} // namespace

// ---------------------------------------------------------------------------
// OmegaHIntersectionRHSIntegrator constructors
// ---------------------------------------------------------------------------

OmegaHIntersectionRHSIntegrator::OmegaHIntersectionRHSIntegrator(
  const FunctionSpace& source_space, const FunctionSpace& target_space)
{
  Data data = BuildData(source_space, target_space);
  source_geometry_ = std::move(data.source_geometry);
  coords_ = std::move(data.coords);
  source_elements_ = std::move(data.source_elements);
  node_gids_ = std::move(data.node_gids);
  coeffs_ = std::move(data.coeffs);
  ndof_per_elem_ = data.ndof_per_elem;

  const PetscInt nnz = static_cast<PetscInt>(node_gids_.extent(0));
  PetscErrorCode ierr =
    createSeqVec(PETSC_COMM_SELF, data.num_target_dofs, &vec_);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  // VecSetPreallocationCOO takes the COO indices on the host
  // TODO: ask Todd/PETSc folks if there is a better way to do this for GPU
  // support
  auto node_gids_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, node_gids_);
  ierr = VecSetPreallocationCOO(vec_, nnz, node_gids_host.data());
  CHKERRABORT(PETSC_COMM_SELF, ierr);
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

CoordinateView<DeviceMemorySpace>
OmegaHIntersectionRHSIntegrator::GetIntegrationPoints() const noexcept
{
  return CoordinateView<DeviceMemorySpace>(
    source_geometry_->GetCoordinateSystem(), MakeConstRank2View(coords_));
}

OwnerHints OmegaHIntersectionRHSIntegrator::IntegrationPointOwners() const
{
  return OwnerHints{source_geometry_, GetMesh(*source_geometry_).dim(),
                    source_elements_};
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
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  auto sv = sampled_values;
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
  CHKERRABORT(PETSC_COMM_SELF, ierr);
}

// ---------------------------------------------------------------------------
// Builder
// ---------------------------------------------------------------------------

std::unique_ptr<OmegaHIntersectionRHSIntegrator>
BuildOmegaHConservativeRHSIntegrator(const FunctionSpace& source_space,
                                     const FunctionSpace& target_space)
{
  return std::make_unique<OmegaHIntersectionRHSIntegrator>(source_space,
                                                           target_space);
}

} // namespace pcms
