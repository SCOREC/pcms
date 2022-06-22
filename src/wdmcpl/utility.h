#ifndef WDM_COUPLING_UTILITY_H
#define WDM_COUPLING_UTILITY_H
#include <tuple>
namespace wdmcpl {
template <typename Tuple, typename Func>
constexpr void tuple_for_each(Tuple&& tuple, Func&& f) {

  std::apply([&](auto&& ...x){(..., f(x));}, tuple);
}

// take a tuple of functor objects and call each
template <typename Tuple,typename...Args>
constexpr void tuple_call_each(Tuple&& tuple, Args&&... args) {

  std::apply([&](auto&& ...x){(..., x(args...));}, tuple);
}
}

#endif // WDM_COUPLING_UTILITY_H
