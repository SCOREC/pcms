#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <Omega_h_array.hpp>
#include "pcms/utility/arrays.h"
#include "pcms/utility/types.h"
#include "numpy_view.h"

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

// Omega_h arrays (Read/Write) have storage in the default (device) memory
// space.

// Copy a 1D numpy array into a new Omega_h::Read.
template <typename T>
Omega_h::Read<T> NumpyToOmegaHRead(const NumpyInput<T>& arr)
{
  return Omega_h::Read<T>(Omega_h::Write<T>(
    NumpyToKokkosView<1, T, DeviceMemorySpace>(arr, "numpy_to_omega_h")));
}

// Copy an Omega_h::Read into a new 1D numpy array.
template <typename T>
py::array_t<T> ToNumpy(Omega_h::Read<T> arr)
{
  return ToNumpy(Rank1View<const T, DeviceMemorySpace, Kokkos::layout_right>(
    arr.data(), arr.size()));
}

} // namespace pcms
