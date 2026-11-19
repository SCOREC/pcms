#ifndef PCMS_COUPLER_COUPLER_TYPES_HPP
#define PCMS_COUPLER_COUPLER_TYPES_HPP

#include <redev.h>

namespace pcms
{

using ProcessType = redev::ProcessType;
/// Whether a send or receive completes before returning (Synchronous) or
/// is completed by the communication phase (Deferred).
enum class Mode
{
  Deferred,
  Synchronous
};

} // namespace pcms

#include "pcms/coupler/detail/coupler_types.hpp"

#endif // PCMS_COUPLER_COUPLER_TYPES_HPP
