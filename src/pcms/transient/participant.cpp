#include "pcms/transient/participant.hpp"

#include <utility>

namespace pcms::transient
{

InterfaceState::InterfaceState() = default;

InterfaceState::InterfaceState(std::size_t n, Real value) : dofs_(n, value) {}

InterfaceState::InterfaceState(std::vector<Real> dofs)
  : dofs_(std::move(dofs))
{
}

std::size_t InterfaceState::Size() const noexcept
{
  return dofs_.size();
}

Real& InterfaceState::operator[](std::size_t i) noexcept
{
  return dofs_[i];
}

Real InterfaceState::operator[](std::size_t i) const noexcept
{
  return dofs_[i];
}

std::span<const Real> InterfaceState::View() const noexcept
{
  return dofs_;
}

std::span<const Real> Participant::ReportQoI() const
{
  return {};
}

} // namespace pcms::transient
