#include "wdmcpl.h"
#include <redev.h>
#include "wdmcpl/coordinate_transform.h"
#include "wdmcpl/types.h"

namespace wdmcpl
{
static_assert(std::is_same_v<Real, redev::Real>,
              "wdmcpl and redev real types must match");
static_assert(std::is_same_v<LO, redev::LO>,
              "wdmcpl and redev LO types must match");
static_assert(std::is_same_v<GO, redev::GO>,
              "wdmcpl and redev GO types must match");
static_assert(std::is_same_v<Real, Omega_h::Real>,
              "wdmcpl and Omega_h real types must match");
static_assert(std::is_same_v<LO, Omega_h::LO>,
              "wdmcpl and Omega_h LO types must match");
static_assert(std::is_same_v<GO, Omega_h::GO>,
              "wdmcpl and Omega_h GO types must match");
} // namespace wdmcpl
