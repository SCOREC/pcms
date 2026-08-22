#include "pcms/transient/timestepper.hpp"

namespace pcms::transient
{

FixedTimestepper::FixedTimestepper(Real dt) : dt_(dt) {}

Real FixedTimestepper::InitialStep() const
{
  return dt_;
}

std::pair<bool, Real> FixedTimestepper::Update(Real, Real)
{
  return {true, dt_};
}

} //namespace pcms::transient
