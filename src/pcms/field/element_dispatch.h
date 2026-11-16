#ifndef PCMS_FIELD_ELEMENT_DISPATCH_H
#define PCMS_FIELD_ELEMENT_DISPATCH_H

#include "pcms/utility/assert.h"
#include <string>
#include <type_traits>

namespace pcms::detail
{

// Runtime -> compile-time bridge for finite-element order.
//
// Turns a runtime Lagrange order into a compile-time std::integral_constant and
// invokes the generic callable `f` with it, so the callee can instantiate
// order-templated kernels (shape functions, DOF-per-element counts) without
// scattering `if (order == ...)` branches through device code. This is the same
// dispatch shape as MakeMeshFieldBackend's switch, factored out so both the
// transfer integrators and the MeshField evaluator backend can share it.
//
// Supported orders: 0 (one DOF per element) and 1 (one DOF per vertex). Add a
// case here when a higher-order element (P2, ...) lands; every consumer picks
// it up at once. `f` must return the same type for every supported order.
template <typename F>
auto DispatchByOrder(int order, F&& f)
{
  switch (order) {
    case 0: return f(std::integral_constant<int, 0>{});
    case 1: return f(std::integral_constant<int, 1>{});
    default:
      throw pcms_error("DispatchByOrder: unsupported element order " +
                       std::to_string(order));
  }
}

} // namespace pcms::detail

#endif // PCMS_FIELD_ELEMENT_DISPATCH_H
