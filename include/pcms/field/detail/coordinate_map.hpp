// Part of pcms/field/coordinate_map.hpp, which includes it after the
// declarations it needs.
#include "pcms/field/coordinate_map.hpp"
#ifndef PCMS_FIELD_DETAIL_COORDINATE_MAP_HPP
#define PCMS_FIELD_DETAIL_COORDINATE_MAP_HPP

#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include <Kokkos_Core.hpp>
#include <array>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

namespace pcms::detail
{

/// The largest dimension of the systems of a kernel map.
inline constexpr int MaxKernelMapDimension = 3;

/// @return {S, T, matches}: the source and target dimensions S, T <= 3 whose
///         typed arrays `Kernel` takes, and how many pairs it takes
template <typename Kernel, int Index = 0>
constexpr std::array<int, 3> KernelDimensions()
{
  if constexpr (Index == MaxKernelMapDimension * MaxKernelMapDimension) {
    return {0, 0, 0};
  } else {
    constexpr int S = Index / MaxKernelMapDimension + 1;
    constexpr int T = Index % MaxKernelMapDimension + 1;
    auto dims = KernelDimensions<Kernel, Index + 1>();
    if constexpr (std::is_invocable_v<
                    const Kernel&, const Kokkos::Array<Real, S>&,
                    Kokkos::Array<Real, T>&, Kokkos::Array<Real, T * S>&>) {
      dims = {S, T, dims[2] + 1};
    }
    return dims;
  }
}

/// A bound map that holds its mapped points and its Jacobian.
class BoundJacobianMap final : public BoundCoordinateMap
{
public:
  /// @param target the map's target system
  /// @param mapped [points][target dimension]
  /// @param jacobian [points][target dimension][source dimension]
  BoundJacobianMap(std::shared_ptr<const CoordinateSystem> target,
                   Kokkos::View<Real**, DeviceMemorySpace> mapped,
                   Kokkos::View<Real***, DeviceMemorySpace> jacobian)
    : BoundCoordinateMap(std::move(target), std::move(mapped)),
      jacobian_(std::move(jacobian))
  {
  }

  [[nodiscard]] View<3, const Real, DeviceMemorySpace> Jacobian() const override
  {
    return View<3, const Real, DeviceMemorySpace>(
      jacobian_.data(), jacobian_.extent(0), jacobian_.extent(1),
      jacobian_.extent(2));
  }

private:
  Kokkos::View<Real***, DeviceMemorySpace> jacobian_;
};

/// Runs `kernel` at every point of `in`: the images into `mapped` and the
/// Jacobians into `jacobian`. (A free function: CUDA does not allow device
/// lambdas in protected members.)
template <int S, int T, typename Kernel>
void EvaluateKernelMap(const Kernel& kernel,
                       Rank2View<const Real, DeviceMemorySpace> in,
                       Kokkos::View<Real**, DeviceMemorySpace> mapped,
                       Kokkos::View<Real***, DeviceMemorySpace> jacobian)
{
  ForEach(
    "KernelMap::Bind", static_cast<LO>(in.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      Kokkos::Array<Real, S> x = {};
      Kokkos::Array<Real, T> out = {};
      Kokkos::Array<Real, T * S> jac = {};
      for (int i = 0; i < S; ++i) {
        x[i] = in(p, i);
      }
      kernel(x, out, jac);
      for (int j = 0; j < T; ++j) {
        mapped(p, j) = out[j];
        for (int i = 0; i < S; ++i) {
          jacobian(p, j, i) = jac[S * j + i];
        }
      }
    });
}

/// The coordinate map of a kernel (CreateCoordinateMap) from S to T
/// coordinates.
template <typename Kernel, int S, int T>
class KernelMap final : public CoordinateMap
{
public:
  KernelMap(std::shared_ptr<const CoordinateSystem> source,
            std::shared_ptr<const CoordinateSystem> target, Kernel kernel)
    : CoordinateMap(std::move(source), std::move(target)),
      kernel_(std::move(kernel))
  {
    if (GetSourceCoordinateSystem()->Dimension() != S ||
        GetTargetCoordinateSystem()->Dimension() != T) {
      throw pcms_error("CreateCoordinateMap: the kernel maps " +
                       std::to_string(S) + " to " + std::to_string(T) +
                       " coordinates; the systems have other dimensions");
    }
  }

protected:
  [[nodiscard]] std::unique_ptr<BoundCoordinateMap> BindImpl(
    const CoordinateView<DeviceMemorySpace>& query_points) const override
  {
    const auto in = query_points.GetValues();
    const auto n = in.extent(0);
    Kokkos::View<Real**, DeviceMemorySpace> mapped(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "mapped_points"), n, T);
    Kokkos::View<Real***, DeviceMemorySpace> jacobian(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "map_jacobian"), n, T, S);
    EvaluateKernelMap<S, T>(kernel_, in, mapped, jacobian);
    return std::make_unique<BoundJacobianMap>(
      GetTargetCoordinateSystem(), std::move(mapped), std::move(jacobian));
  }

private:
  Kernel kernel_;
};

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_COORDINATE_MAP_HPP
