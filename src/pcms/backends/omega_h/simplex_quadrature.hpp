#ifndef PCMS_BACKENDS_OMEGA_H_SIMPLEX_QUADRATURE_HPP
#define PCMS_BACKENDS_OMEGA_H_SIMPLEX_QUADRATURE_HPP

#include "pcms/configuration.hpp"
#ifdef PCMS_ENABLE_OMEGA_H_MESHFIELDS
#include "pcms/utility/assert.hpp"
#include <MeshField_Integrate.hpp>
#include <MeshField_Shape.hpp>
#include <Omega_h_defines.hpp>
#include <Kokkos_Core.hpp>
#include <string>

namespace pcms::omega_h::detail
{

template <int Dim>
struct SimplexTopology;
template <>
struct SimplexTopology<2>
{
  static constexpr MeshField::Mesh_Topology value = MeshField::Triangle;
};
template <>
struct SimplexTopology<3>
{
  static constexpr MeshField::Mesh_Topology value = MeshField::Tetrahedron;
};

/// @brief MeshFields' reference rule of a simplex on the device: Dim + 1
///        barycentric coordinates and a weight per point; the weights sum to
///        the reference measure 1 / Dim!.
template <int Dim>
struct IntegrationData
{
  Kokkos::View<Omega_h::Real * [Dim + 1]> bary_coords;
  Kokkos::View<Omega_h::Real*> weights;

  /// @param order the degree the rule integrates exactly
  /// @throws pcms_error if MeshFields has no rule of that order
  explicit IntegrationData(int order)
  {
    const auto* rule =
      MeshField::getIntegration<SimplexTopology<Dim>::value>()->getAccurate(
        order);
    if (rule == nullptr) {
      throw pcms_error("IntegrationData: MeshFields has no simplex rule of "
                       "order " +
                       std::to_string(order));
    }
    const auto points = rule->getPoints();
    bary_coords =
      Kokkos::View<Omega_h::Real * [Dim + 1]>("bary_coords", points.size());
    weights = Kokkos::View<Omega_h::Real*>("weights", points.size());
    auto bary_coords_host = Kokkos::create_mirror_view(bary_coords);
    auto weights_host = Kokkos::create_mirror_view(weights);
    for (std::size_t i = 0; i < points.size(); ++i) {
      // MeshFields stores the first Dim barycentric coordinates
      bary_coords_host(i, Dim) = 1.0;
      for (int d = 0; d < Dim; ++d) {
        bary_coords_host(i, d) = points[i].param[d];
        bary_coords_host(i, Dim) -= points[i].param[d];
      }
      weights_host(i) = points[i].weight;
    }
    Kokkos::deep_copy(bary_coords, bary_coords_host);
    Kokkos::deep_copy(weights, weights_host);
  }

  int size() const { return bary_coords.extent(0); }
};

} // namespace pcms::omega_h::detail

#endif // PCMS_ENABLE_OMEGA_H_MESHFIELDS
#endif // PCMS_BACKENDS_OMEGA_H_SIMPLEX_QUADRATURE_HPP
