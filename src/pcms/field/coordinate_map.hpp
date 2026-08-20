#ifndef PCMS_FIELD_COORDINATE_MAP_HPP
#define PCMS_FIELD_COORDINATE_MAP_HPP

#include "pcms/field/basis_transformation.hpp"
#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/point_status.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include <Kokkos_Core.hpp>
#include <memory>

namespace pcms
{

struct MappedPoints
{
  std::shared_ptr<const CoordinateSystem> system;
  Kokkos::View<Real**, DeviceMemorySpace> coords;
  PointStatusView status; 

  [[nodiscard]] CoordinateView<DeviceMemorySpace> View() const
  {
    return CoordinateView<DeviceMemorySpace>(system,
                                             MakeConstRank2View(coords));
  }

  [[nodiscard]] LO NumPoints() const noexcept
  {
    return static_cast<LO>(coords.extent(0));
  }
};

class CoordinateMap
{
public:
  [[nodiscard]] virtual std::shared_ptr<const CoordinateSystem>
  GetSourceCoordinateSystem() const noexcept = 0;
  [[nodiscard]] virtual std::shared_ptr<const CoordinateSystem>
  GetTargetCoordinateSystem() const noexcept = 0;

  [[nodiscard]] virtual MappedPoints Map(
    const CoordinateView<DeviceMemorySpace>& points) const = 0;

  [[nodiscard]] virtual bool SupportsBoundaryProjection() const noexcept
  {
    return false;
  }

  [[nodiscard]] virtual std::unique_ptr<BoundBasisTransformation>
  MakeBasisTransformation(
    const CoordinateView<DeviceMemorySpace>& source_points,
    const MappedPoints& mapped) const;

  virtual ~CoordinateMap() = default;

protected:
  void ValidateSourcePoints(const CoordinateView<DeviceMemorySpace>& points,
                            const char* who) const;
};

class CartesianToCylindrical final : public CoordinateMap
{
public:
  [[nodiscard]] std::shared_ptr<const CoordinateSystem>
  GetSourceCoordinateSystem() const noexcept override;
  [[nodiscard]] std::shared_ptr<const CoordinateSystem>
  GetTargetCoordinateSystem() const noexcept override;
  [[nodiscard]] MappedPoints Map(
    const CoordinateView<DeviceMemorySpace>& points) const override;
  [[nodiscard]] std::unique_ptr<BoundBasisTransformation>
  MakeBasisTransformation(
    const CoordinateView<DeviceMemorySpace>& source_points,
    const MappedPoints& mapped) const override;
};

class CylindricalToCartesian final : public CoordinateMap
{
public:
  [[nodiscard]] std::shared_ptr<const CoordinateSystem>
  GetSourceCoordinateSystem() const noexcept override;
  [[nodiscard]] std::shared_ptr<const CoordinateSystem>
  GetTargetCoordinateSystem() const noexcept override;
  [[nodiscard]] MappedPoints Map(
    const CoordinateView<DeviceMemorySpace>& points) const override;
  [[nodiscard]] std::unique_ptr<BoundBasisTransformation>
  MakeBasisTransformation(
    const CoordinateView<DeviceMemorySpace>& source_points,
    const MappedPoints& mapped) const override;
};

} // namespace pcms

#endif // PCMS_FIELD_COORDINATE_MAP_HPP
