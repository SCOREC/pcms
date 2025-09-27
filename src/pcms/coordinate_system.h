#ifndef PCMS_COORDINATE_SYSTEM_H
#define PCMS_COORDINATE_SYSTEM_H
#include "arrays.h"

namespace pcms
{

enum class CoordinateSystem
{
  Cartesian,
  Cylindrical,
  XGC,
  GNET,
  BEAMS3D
};

template <typename MemorySpace>
class CoordinateView
{
public:
  CoordinateView(CoordinateSystem cs,
                 Rank2View<const Real, MemorySpace> coords) noexcept
    : coordinate_system_(cs), coordinates_(coords)
  {
  }

  [[nodiscard]] CoordinateSystem GetCoordinateSystem() const noexcept
  {
    return coordinate_system_;
  }

  [[nodiscard]] Rank2View<const Real, MemorySpace> GetCoordinates()
    const noexcept
  {
    return coordinates_;
  }

  // would prefer if these operations were limited to use by
  // CoordinateTransformation as they are unsafe (i.e., you can break class
  // invariant) passkey pattern?
  [[nodiscard]] Rank2View<const Real, MemorySpace> GetCoordinates() noexcept
  {
    return coordinates_;
  }
  void SetCoordinateSystem(CoordinateSystem cs) noexcept
  {
    coordinate_system_ = cs;
  }

private:
  CoordinateSystem coordinate_system_;
  Rank2View<const Real, MemorySpace> coordinates_;
};

class CoordinateTransformation
{
public:
  virtual void Evaluate(const CoordinateView<HostMemorySpace> from,
                        CoordinateView<HostMemorySpace> to) = 0;
  // TODO
  // #IFDEF PCMS_HAS_DEVICE
  // virtual void Evaluate(const CoordinateView<DeviceMemorySpace> from,
  // CoordinateView<DeviceMemorySpace> to) = 0;

  virtual ~CoordinateTransformation() noexcept = default;
};

/*
class CartesianToCylindrical : public CoordinateTransformation {
public:
  void Evaluate(const CoordinateView<HostMemorySpace> from,
CoordinateView<HostMemorySpace> to) override; private: Rank2View<const Real,
HostMemorySpace> from_;
};

class CylindricalToCartesian : public CoordinateTransformation {
public:
  void Evaluate(const CoordinateView from, CoordinateView to) override;
private:
  Rank2View<const double> from_;
};

// this function creates a transformation object based on the input coordinate
systems
// trhows an exception if the transformation is not supported
std::unique_ptr<CoordinateTransformation>
CreateCoordinateTransformation(CoordinateSystem from, CoordinateSystem to);

// this is a convenience function that creates a transformation and applies it
void Transform(const CoordinateView from, CoordinateView to);
*/

} // namespace pcms

#endif // PCMS_COORDINATE_SYSTEM_H
