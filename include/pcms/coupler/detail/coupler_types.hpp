#ifndef PCMS_COUPLER_DETAIL_COUPLER_TYPES_HPP
#define PCMS_COUPLER_DETAIL_COUPLER_TYPES_HPP

#include "pcms/coupler/coupler_types.hpp"

namespace pcms::detail
{

inline redev::Mode ToRedev(Mode mode) noexcept
{
  return mode == Mode::Deferred ? redev::Mode::Deferred
                                : redev::Mode::Synchronous;
}

} // namespace pcms::detail

#endif // PCMS_COUPLER_DETAIL_COUPLER_TYPES_HPP
