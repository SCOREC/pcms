#include "pcms/transfer/omega_h_intersection_quadrature.hpp"
#include "pcms/field/element_dispatch.h"
#include "pcms/field/evaluation_request.h"
#include "pcms/transfer/omega_h_form_integrator_utils.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include <algorithm>
#include <Omega_h_for.hpp>
#include <Omega_h_shape.hpp>

namespace pcms
{

namespace
{

// Results of the two-pass intersection kernel; the constructor moves these into
// the integrator's members.
struct Data
{
  Kokkos::View<Real**, DeviceMemorySpace> coords;
  Kokkos::View<LO*, DeviceMemorySpace> source_elem_ids;
  Kokkos::View<PetscInt*, DeviceMemorySpace> node_gids;
  Kokkos::View<Real*, DeviceMemorySpace> coeffs;
  int ndof_per_elem = 0;
  PetscInt num_target_dofs = 0;
};

template <int Dim, int TgtOrder>
Data BuildDataImpl(const OmegaHLagrangeLayout& source_layout,
                   const OmegaHLagrangeLayout& target_layout, int quad_order,
                   const GridPointSearchVariant* source_search,
                   const IntersectionResults* precomputed);

Data BuildData(const std::shared_ptr<const OmegaHLagrangeLayout>& source_layout,
               CoordinateSystem source_coordinate_system,
               const std::shared_ptr<const OmegaHLagrangeLayout>& target_layout,
               CoordinateSystem target_coordinate_system,
               const GridPointSearchVariant* source_search,
               const MeshIntersection* intersection)
{
  detail::CheckOmegaHScalarLagrangeLayout(
    source_coordinate_system, source_layout, "OmegaHIntersectionQuadrature",
    "source");
  detail::CheckOmegaHScalarLagrangeLayout(
    target_coordinate_system, target_layout, "OmegaHIntersectionQuadrature",
    "target");

  const int dim = target_layout->GetMesh().dim();
  if (source_layout->GetMesh().dim() != dim) {
    throw pcms_error("OmegaHIntersectionQuadrature: source and target mesh "
                     "dimensions differ");
  }

  const IntersectionResults* precomputed = nullptr;
  if (intersection != nullptr) {
    if (!intersection->GetSourceDiscretization()->SameEntities(
          *source_layout->GetDiscretization()) ||
        !intersection->GetTargetDiscretization()->SameEntities(
          *target_layout->GetDiscretization())) {
      throw pcms_error("OmegaHIntersectionQuadrature: the supplied mesh "
                       "intersection is not between the source and target "
                       "meshes");
    }
    precomputed = &intersection->GetTargetToSource();
  }

  // The integrand f_src * phi_target has polynomial degree source_order +
  // target_order on each intersection sub-simplex; integrate it exactly (with a
  // 1-point floor so a P0->P0 pair still gets a valid rule).
  const int quad_order =
    std::max(1, source_layout->GetOrder() + target_layout->GetOrder());

  return detail::DispatchByOrder(target_layout->GetOrder(), [&](auto order_c) {
    constexpr int TgtOrder = decltype(order_c)::value;
    if (dim == 3) {
      return BuildDataImpl<3, TgtOrder>(*source_layout, *target_layout,
                                        quad_order, source_search, precomputed);
    }
    return BuildDataImpl<2, TgtOrder>(*source_layout, *target_layout,
                                      quad_order, source_search, precomputed);
  });
}

template <int Dim, int TgtOrder>
Data BuildDataImpl(const OmegaHLagrangeLayout& source_layout,
                   const OmegaHLagrangeLayout& target_layout, int quad_order,
                   const GridPointSearchVariant* source_search,
                   const IntersectionResults* precomputed)
{
  using Basis = detail::TargetSimplexBasis<Dim, TgtOrder>;
  constexpr int ndof = Basis::ndof;
  // Reference-to-physical Jacobian factor for a simplex: the reference simplex
  // measure is 1/Dim! (1/2 in 2D, 1/6 in 3D), so a physical sub-simplex of
  // measure `m` scales the reference quadrature weights by Dim! * m.
  constexpr Omega_h::Real ref_factor = (Dim == 3) ? 6.0 : 2.0;

  Omega_h::Mesh& source_mesh = source_layout.GetMesh();
  Omega_h::Mesh& target_mesh = target_layout.GetMesh();

  const IntersectionResults intersections =
    precomputed     ? *precomputed
    : source_search ? intersectTargets(source_mesh, target_mesh, *source_search)
                    : intersectTargets(source_mesh, target_mesh);

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

  // Pass 2: fill coords, node_gids, and coeffs on device. node_gids/coeffs hold
  // ndof (target DOFs per element) entries per integration point.
  Kokkos::View<Real**, DeviceMemorySpace> coords("rhs_coords", num_pts, Dim);
  Kokkos::View<LO*, DeviceMemorySpace> source_elem_ids("rhs_source_elems",
                                                       num_pts);
  Kokkos::View<PetscInt*, DeviceMemorySpace> node_gids(
    "rhs_node_gids", static_cast<std::size_t>(num_pts) * ndof);
  Kokkos::View<Real*, DeviceMemorySpace> coeffs(
    "rhs_coeffs", static_cast<std::size_t>(num_pts) * ndof);

  Kokkos::parallel_for(
    "rhs_fill", nelems, KOKKOS_LAMBDA(int elm) {
      const auto tgt_verts =
        Omega_h::gather_verts<Dim + 1>(tgt_elems2nodes, elm);
      const Omega_h::Matrix<Dim, Dim + 1> tgt_vert_mat =
        Omega_h::gather_vectors<Dim + 1, Dim>(tgt_coords, tgt_verts);
      Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1> tgt_omh;
      for (int i = 0; i < Dim + 1; ++i) {
        for (int d = 0; d < Dim; ++d) {
          tgt_omh[i][d] = tgt_vert_mat[i][d];
        }
      }

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
            source_elem_ids(global_ip) = src_elm;
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
  d.source_elem_ids = std::move(source_elem_ids);
  d.node_gids = std::move(node_gids);
  d.coeffs = std::move(coeffs);
  d.ndof_per_elem = ndof;
  d.num_target_dofs =
    static_cast<PetscInt>(target_layout.GetNumOwnedDofHolder());
  return d;
}

} // namespace

OmegaHIntersectionQuadrature::OmegaHIntersectionQuadrature(
  const FunctionSpace& source_space, const FunctionSpace& target_space,
  std::shared_ptr<const MeshIntersection> intersection)
  : OmegaHIntersectionQuadrature(
      std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(
        source_space.GetLayout()),
      source_space.GetCoordinateSystem(),
      std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(
        target_space.GetLayout()),
      target_space.GetCoordinateSystem(), SourceSearchFromSpace(source_space),
      intersection.get())
{
  source_evaluator_ =
    source_space.CreatePointEvaluator<Real>(EvaluationRequest::FromElements(
      GetIntegrationPoints(), GetSourceElementIds()));
}

OmegaHIntersectionQuadrature::OmegaHIntersectionQuadrature(
  std::shared_ptr<const OmegaHLagrangeLayout> source_layout,
  CoordinateSystem source_coordinate_system,
  std::shared_ptr<const OmegaHLagrangeLayout> target_layout,
  CoordinateSystem target_coordinate_system,
  const GridPointSearchVariant* source_search,
  const MeshIntersection* intersection)
  : source_layout_(std::move(source_layout)),
    target_layout_(std::move(target_layout))
{
  Data data =
    BuildData(source_layout_, source_coordinate_system, target_layout_,
              target_coordinate_system, source_search, intersection);
  coords_ = std::move(data.coords);
  source_elem_ids_ = std::move(data.source_elem_ids);
  target_dofs_ = std::move(data.node_gids);
  weights_ = std::move(data.coeffs);
  dofs_per_element_ = data.ndof_per_elem;
  num_target_dofs_ = data.num_target_dofs;
}

CoordinateView<DeviceMemorySpace>
OmegaHIntersectionQuadrature::GetIntegrationPoints() const noexcept
{
  return CoordinateView<DeviceMemorySpace>(CoordinateSystem::Cartesian,
                                           MakeConstRank2View(coords_));
}

LO OmegaHIntersectionQuadrature::GetNumPoints() const noexcept
{
  return static_cast<LO>(coords_.extent(0));
}

Kokkos::View<const LO*, DeviceMemorySpace>
OmegaHIntersectionQuadrature::GetSourceElementIds() const noexcept
{
  return source_elem_ids_;
}

Kokkos::View<const PetscInt*, DeviceMemorySpace>
OmegaHIntersectionQuadrature::GetTargetDofs() const noexcept
{
  return target_dofs_;
}

Kokkos::View<const Real*, DeviceMemorySpace>
OmegaHIntersectionQuadrature::GetWeights() const noexcept
{
  return weights_;
}

int OmegaHIntersectionQuadrature::GetDofsPerElement() const noexcept
{
  return dofs_per_element_;
}

PetscInt OmegaHIntersectionQuadrature::GetNumTargetDofs() const noexcept
{
  return num_target_dofs_;
}

const PointEvaluator<Real>* OmegaHIntersectionQuadrature::GetSourceEvaluator()
  const noexcept
{
  return source_evaluator_.get();
}

const std::shared_ptr<const OmegaHLagrangeLayout>&
OmegaHIntersectionQuadrature::GetSourceLayout() const noexcept
{
  return source_layout_;
}

const std::shared_ptr<const OmegaHLagrangeLayout>&
OmegaHIntersectionQuadrature::GetTargetLayout() const noexcept
{
  return target_layout_;
}

} // namespace pcms
