#ifndef PCMS_PYTHONAPI_NUMPY_ARRAY_TRANSFORM_HPP
#define PCMS_PYTHONAPI_NUMPY_ARRAY_TRANSFORM_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/types.hpp"
#include "pcms/pythonapi/numpy_view.hpp"

namespace py = pybind11;

// Mirror pybind11's built-in std::array type caster for Kokkos::Array so that
// def_readwrite / def_property members work with Python sequences.
namespace pybind11::detail
{
template <typename Type, size_t Size>
struct type_caster<Kokkos::Array<Type, Size>>
{
  using ArrayType = Kokkos::Array<Type, Size>;
  using value_conv = make_caster<Type>;

  bool load(handle src, bool convert)
  {
    if (!isinstance<sequence>(src))
      return false;
    auto seq = reinterpret_borrow<sequence>(src);
    if (seq.size() != Size)
      return false;
    size_t i = 0;
    for (auto it : seq) {
      value_conv conv;
      if (!conv.load(it, convert))
        return false;
      value[i++] = cast_op<Type>(conv);
    }
    return true;
  }

  template <typename T>
  static handle cast(T&& src, return_value_policy policy, handle parent)
  {
    list l(Size);
    size_t i = 0;
    for (auto& v : src)
      l[i++] = make_caster<Type>::cast(v, policy, parent);
    return l.release();
  }

  PYBIND11_TYPE_CASTER(ArrayType, const_name("List[") +
                                    make_caster<Type>::name + const_name("]"));
};
} // namespace pybind11::detail

namespace pcms
{

} // namespace pcms

#endif // PCMS_PYTHONAPI_NUMPY_ARRAY_TRANSFORM_HPP
