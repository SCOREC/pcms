#include "pcms/backends/uniform_grid/uniform_grid_geometry.hpp"
#include "pcms/geometry/detail/geometry.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include "pcms/utility/detail/profile.hpp"
#include <variant>
#include <vector>

namespace pcms::uniform_grid
{

namespace
{

template <typename T>
using DeviceView = Kokkos::View<T, DeviceMemorySpace>;

/// @brief Fills the tensor rule of `nodes` and `rule_weights` on [-1, 1] over
///        each cell of the grid.
/// @param xi the points, cell by cell
/// @param weights their weights in the coordinate measure
template <unsigned Dim>
void FillCellRule(UniformGrid<Dim> grid, DeviceView<const Real*> nodes,
                  DeviceView<const Real*> rule_weights, DeviceView<Real**> xi,
                  DeviceView<Real*> weights)
{
  const int n = nodes.extent(0);
  int per_cell = 1;
  for (unsigned d = 0; d < Dim; ++d) {
    per_cell *= n;
  }
  pcms::detail::ForEach(
    "uniform_grid_cell_rule", static_cast<LO>(weights.extent(0)),
    KOKKOS_LAMBDA(const LO i) {
      const LO cell = i / per_cell;
      int digits = i % per_cell;
      const auto box = grid.GetCellBBOX(cell);
      Real weight = 1;
      for (unsigned d = 0; d < Dim; ++d) {
        const int k = digits % n;
        digits /= n;
        xi(i, d) = box.center[d] + box.half_width[d] * nodes(k);
        weight *= box.half_width[d] * rule_weights(k);
      }
      weights(i) = weight;
    });
}

} // namespace

template <unsigned Dim>
std::shared_ptr<BoundPoints> UniformGridGeometry<Dim>::BindImpl(
  const BindRequest& request) const
{
  PCMS_FUNCTION_TIMER;
  return std::make_shared<pcms::detail::UniformGridBinding<Dim>>(
    shared_from_this(), request, grid_);
}

template <unsigned Dim>
std::pair<Kokkos::View<Real**, DeviceMemorySpace>,
          Kokkos::View<Real*, DeviceMemorySpace>>
UniformGridGeometry<Dim>::CellRule(int order) const
{
  if (order < 1) {
    throw pcms_error("CreateQuadrature: the order must be positive");
  }
  std::vector<Real> nodes_host;
  std::vector<Real> weights_host;
  pcms::detail::GaussLegendre(order, nodes_host, weights_host);
  const auto per_axis = static_cast<int>(nodes_host.size());
  DeviceView<Real*> nodes("gauss_nodes", per_axis);
  DeviceView<Real*> rule("gauss_weights", per_axis);
  Kokkos::deep_copy(nodes,
                    Kokkos::View<const Real*, HostMemorySpace,
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                      nodes_host.data(), per_axis));
  Kokkos::deep_copy(rule, Kokkos::View<const Real*, HostMemorySpace,
                                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                            weights_host.data(), per_axis));
  size_t num_points = grid_.GetNumCells();
  for (unsigned d = 0; d < Dim; ++d) {
    num_points *= per_axis;
  }
  DeviceView<Real**> xi("quadrature_xi", num_points, Dim);
  DeviceView<Real*> weights("quadrature_weights", num_points);
  FillCellRule<Dim>(grid_, nodes, rule, xi, weights);
  return {xi, weights};
}

template <unsigned Dim>
Quadrature UniformGridGeometry<Dim>::CreateQuadrature(
  const pcms::detail::QuadratureRule& rule) const
{
  const auto* gauss = std::get_if<pcms::detail::GaussRule>(&rule);
  if (gauss == nullptr) {
    return Geometry::CreateQuadrature(rule);
  }
  auto [xi, weights] = CellRule(gauss->order);
  return pcms::detail::BindQuadrature(
    *this, CoordinateView<DeviceMemorySpace>(system_, MakeConstRank2View(xi)),
    weights);
}

template class UniformGridGeometry<1>;
template class UniformGridGeometry<2>;
template class UniformGridGeometry<3>;

template <unsigned Dim>
std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<Dim> grid, std::shared_ptr<const CoordinateSystem> system)
{
  return std::make_shared<UniformGridGeometry<Dim>>(grid, std::move(system));
}

template std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<1>, std::shared_ptr<const CoordinateSystem>);
template std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<2>, std::shared_ptr<const CoordinateSystem>);
template std::shared_ptr<const Geometry> CreateGeometry(
  UniformGrid<3>, std::shared_ptr<const CoordinateSystem>);

} // namespace pcms::uniform_grid
