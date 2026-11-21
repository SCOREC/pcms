#include "pcms/transfer/omega_h_intersection_quadrature.hpp"
#include "pcms/geometry/detail/geometry.hpp"
#include "pcms/transfer/omega_h_form_integrator_utils.hpp"
#include "pcms/utility/arrays.hpp"
#include <Omega_h_for.hpp>
#include <Omega_h_shape.hpp>

namespace pcms::omega_h
{

namespace
{

template <int Dim>
Quadrature CreateIntersectionQuadratureImpl(
  std::shared_ptr<const OmegaHGeometry> source,
  std::shared_ptr<const OmegaHGeometry> target, int order)
{
  // the reference rule's weights sum to 1 / Dim!
  constexpr Real ReferenceFactor = (Dim == 3) ? 6.0 : 2.0;
  Omega_h::Mesh& target_mesh = target->GetMesh();
  const auto intersections = intersectTargets(*source, target_mesh);
  const auto tgt_coords = target_mesh.coords();
  const auto tgt_elems2nodes = target_mesh.ask_down(Dim, Omega_h::VERT).ab2b;
  const auto src_coords = source->GetMesh().coords();
  const auto src_elems2nodes =
    source->GetMesh().ask_down(Dim, Omega_h::VERT).ab2b;
  const transfer::detail::IntegrationData<Dim> rule(order);
  const int npts = rule.size();
  const auto bary_coords = rule.bary_coords;
  const auto rule_weights = rule.weights;
  const int nelems = target_mesh.nelems();

  Omega_h::Write<Omega_h::LO> counts(nelems, 0, "quadrature_counts");
  Kokkos::parallel_for(
    "CreateIntersectionQuadrature::count", nelems, KOKKOS_LAMBDA(int elm) {
      int count = 0;
      transfer::detail::ForEachIntersectionSubsimplex<Dim>(
        elm, intersections, tgt_coords, src_coords, tgt_elems2nodes,
        src_elems2nodes,
        [&](const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>&, int,
            Omega_h::Real) { count += npts; });
      counts[elm] = count;
    });
  const auto offsets =
    Omega_h::offset_scan(Omega_h::Read<Omega_h::LO>(counts), "offsets");
  const int num_points = offsets.last();

  Kokkos::View<Real**, DeviceMemorySpace> coords("quadrature_coords",
                                                 num_points, Dim);
  Kokkos::View<LO*, DeviceMemorySpace> target_elements("target_elements",
                                                       num_points);
  Kokkos::View<LO*, DeviceMemorySpace> source_elements("source_elements",
                                                       num_points);
  Kokkos::View<Real*, DeviceMemorySpace> weights("quadrature_weights",
                                                 num_points);
  Kokkos::parallel_for(
    "CreateIntersectionQuadrature::fill", nelems, KOKKOS_LAMBDA(int elm) {
      int i = offsets[elm];
      transfer::detail::ForEachIntersectionSubsimplex<Dim>(
        elm, intersections, tgt_coords, src_coords, tgt_elems2nodes,
        src_elems2nodes,
        [&](const Omega_h::Few<Omega_h::Vector<Dim>, Dim + 1>& sub, int src_elm,
            Omega_h::Real measure) {
          for (int q = 0; q < npts; ++q, ++i) {
            Omega_h::Vector<Dim + 1> bary;
            for (int k = 0; k <= Dim; ++k) {
              bary[k] = bary_coords(q, k);
            }
            const auto pt =
              transfer::detail::GlobalFromBarycentric<Dim>(bary, sub);
            for (int d = 0; d < Dim; ++d) {
              coords(i, d) = pt[d];
            }
            target_elements(i) = elm;
            source_elements(i) = src_elm;
            weights(i) = rule_weights(q) * ReferenceFactor * measure;
          }
        });
    });
  auto bound =
    target->Bind(CoordinateView<DeviceMemorySpace>(
                   target->GetCoordinateSystem(), MakeConstRank2View(coords)),
                 OwnerHints{target, Dim, target_elements});
  pcms::detail::ScaleByVolumeFactor(weights, *bound);
  return Quadrature{std::move(bound), weights,
                    OwnerHints{std::move(source), Dim, source_elements}};
}

} // namespace

Quadrature CreateIntersectionQuadrature(
  std::shared_ptr<const OmegaHGeometry> source,
  std::shared_ptr<const OmegaHGeometry> target, int order)
{
  // order 1 is the lowest rule; it is exact for order 0
  order = std::max(1, order);
  if (target->GetMesh().dim() == 3) {
    return CreateIntersectionQuadratureImpl<3>(std::move(source),
                                               std::move(target), order);
  }
  return CreateIntersectionQuadratureImpl<2>(std::move(source),
                                             std::move(target), order);
}

} // namespace pcms::omega_h
