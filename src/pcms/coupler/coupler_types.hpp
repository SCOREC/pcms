#ifndef PCMS_COUPLER_TYPES_H_
#define PCMS_COUPLER_TYPES_H_

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

namespace detail
{
inline redev::Mode ToRedev(Mode mode) noexcept
{
  return mode == Mode::Deferred ? redev::Mode::Deferred
                                : redev::Mode::Synchronous;
}
} // namespace detail

} // namespace pcms

#endif // PCMS_COUPLER_TYPES_H_
