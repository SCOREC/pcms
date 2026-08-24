#ifndef PCMS_FIELD_BASIS_TRANSFORMATION_HPP
#define PCMS_FIELD_BASIS_TRANSFORMATION_HPP

#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/point_status.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include <Kokkos_Core.hpp>
#include <memory>
#include <type_traits>
#include <utility>

namespace pcms
{

// Re-expresses components from one ValueBasis to another at the point set the
// transformation was bound to.
class BoundBasisTransformation
{
public:
  [[nodiscard]] const ValueBasis& GetSourceBasis() const noexcept
  {
    return source_basis_;
  }
  [[nodiscard]] const ValueBasis& GetTargetBasis() const noexcept
  {
    return target_basis_;
  }

  [[nodiscard]] LO NumPoints() const noexcept { return num_points_; }

  [[nodiscard]] const PointStatusView& GetStatus() const noexcept
  {
    return status_;
  }

  void Apply(ValueView<const Real, DeviceMemorySpace> in,
             ValueView<Real, DeviceMemorySpace> out) const;

  virtual ~BoundBasisTransformation() = default;

  BoundBasisTransformation(const BoundBasisTransformation&) = delete;
  BoundBasisTransformation& operator=(const BoundBasisTransformation&) = delete;

protected:
  BoundBasisTransformation(ValueBasis source_basis, ValueBasis target_basis,
                           LO num_points);

  virtual void ApplyVectorImpl(
    Rank2View<const Real, DeviceMemorySpace> in,
    Rank2View<Real, DeviceMemorySpace> out) const = 0;

  PointStatusView status_; 
private:
  ValueBasis source_basis_;
  ValueBasis target_basis_;
  LO num_points_;
};

// basis-to-basis component transformation
class BasisTransformation
{
public:
  [[nodiscard]] virtual const ValueBasis& GetSourceBasis() const noexcept = 0;
  [[nodiscard]] virtual const ValueBasis& GetTargetBasis() const noexcept = 0;

  [[nodiscard]] virtual std::unique_ptr<BoundBasisTransformation> Bind(
    const CoordinateView<DeviceMemorySpace>& points) const = 0;

  virtual ~BasisTransformation() = default;
};

// basis bound to a set of coordinates
template <typename Basis>
class BoundBasis final : public BoundBasisTransformation
{
public:
  BoundBasis(Basis basis, const CoordinateView<DeviceMemorySpace>& points)
    : BoundBasisTransformation(basis.Source(), basis.Target(),
                               static_cast<LO>(points.GetValues().extent(0))),
      basis_(std::move(basis)),
      state_(basis_.Prepare(points, status_))
  {
  }

private:
  void ApplyVectorImpl(
    Rank2View<const Real, DeviceMemorySpace> in,
    Rank2View<Real, DeviceMemorySpace> out) const override
  {
    basis_.ApplyVector(state_, in, out, NumPoints());
  }

  Basis basis_;
  typename Basis::State state_;
};

// Basis prior to being bound to a set of coordinates
template <typename Basis>
class UnboundBasis final : public BasisTransformation
{
public:
  UnboundBasis() = default;

  template <typename Arg0, typename... Args,
            typename = std::enable_if_t<
              !std::is_same_v<std::decay_t<Arg0>, UnboundBasis> &&
              std::is_constructible_v<Basis, Arg0, Args...>>>
  explicit UnboundBasis(Arg0&& arg0, Args&&... args)
    : basis_(std::forward<Arg0>(arg0), std::forward<Args>(args)...)
  {
  }

  [[nodiscard]] const ValueBasis& GetSourceBasis() const noexcept override
  {
    return basis_.Source();
  }
  [[nodiscard]] const ValueBasis& GetTargetBasis() const noexcept override
  {
    return basis_.Target();
  }
  [[nodiscard]] std::unique_ptr<BoundBasisTransformation> Bind(
    const CoordinateView<DeviceMemorySpace>& points) const override
  {
    return std::make_unique<BoundBasis<Basis>>(basis_, points);
  }

private:
  Basis basis_;
};

namespace basis
{

struct CartesianToCylindrical
{
  struct State
  {
    Kokkos::View<Real*, DeviceMemorySpace> cos_theta;
    Kokkos::View<Real*, DeviceMemorySpace> sin_theta;
  };

  [[nodiscard]] const ValueBasis& Source() const noexcept;
  [[nodiscard]] const ValueBasis& Target() const noexcept;
  [[nodiscard]] State Prepare(const CoordinateView<DeviceMemorySpace>& points,
                              PointStatusView& status) const;
  void ApplyVector(const State& state,
                   Rank2View<const Real, DeviceMemorySpace> in,
                   Rank2View<Real, DeviceMemorySpace> out, LO n) const;
};

struct CylindricalToCartesian
{
  using State = CartesianToCylindrical::State;

  [[nodiscard]] const ValueBasis& Source() const noexcept;
  [[nodiscard]] const ValueBasis& Target() const noexcept;
  [[nodiscard]] State Prepare(const CoordinateView<DeviceMemorySpace>& points,
                              PointStatusView& status) const;
  void ApplyVector(const State& state,
                   Rank2View<const Real, DeviceMemorySpace> in,
                   Rank2View<Real, DeviceMemorySpace> out, LO n) const;
};

} // namespace basis

extern template class BoundBasis<basis::CartesianToCylindrical>;
extern template class UnboundBasis<basis::CartesianToCylindrical>;
extern template class BoundBasis<basis::CylindricalToCartesian>;
extern template class UnboundBasis<basis::CylindricalToCartesian>;

using CartesianToCylindricalBasis = UnboundBasis<basis::CartesianToCylindrical>;
using CylindricalToCartesianBasis = UnboundBasis<basis::CylindricalToCartesian>;

} // namespace pcms

#endif // PCMS_FIELD_BASIS_TRANSFORMATION_HPP
