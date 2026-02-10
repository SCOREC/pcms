#include "pcms.h"
#include "pcms/types.h"

namespace pcms
{
static_assert(std::is_same_v<Real, redev::Real>,
              "pcms and redev real types must match");
static_assert(std::is_same_v<LO, redev::LO>,
              "pcms and redev LO types must match");
static_assert(std::is_same_v<GO, redev::GO>,
              "pcms and redev GO types must match");
#ifdef PCMS_ENABLE_OMEGA_H
static_assert(std::is_same_v<Real, Omega_h::Real>,
              "pcms and Omega_h real types must match");
static_assert(std::is_same_v<LO, Omega_h::LO>,
              "pcms and Omega_h LO types must match");
static_assert(std::is_same_v<GO, Omega_h::GO>,
              "pcms and Omega_h GO types must match");
#endif
} // namespace pcms
