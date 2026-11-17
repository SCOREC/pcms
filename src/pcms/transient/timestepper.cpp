#include "pcms/transient/timestepper.hpp"

namespace pcms
{

FixedTimestepper::FixedTimestepper(Real dt) : dt_(dt) {}

Real FixedTimestepper::InitialStep() const
{
  return dt_;
}

Step FixedTimestepper::UpdateStep(Real, Real)
{
  return {true, dt_};
}

} // namespace pcms
