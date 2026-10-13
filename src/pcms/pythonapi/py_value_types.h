#ifndef PCMS_PYTHONAPI_PY_VALUE_TYPES_H
#define PCMS_PYTHONAPI_PY_VALUE_TYPES_H

// Python-side counterparts of the non-owning C++ view types. A C++
// CoordinateView or ValueView borrows its memory, so the Python objects pair
// the same tag (coordinate system or value basis) with a numpy array that owns
// or keeps alive that memory, and are converted to the C++ view at each call.

#include <memory>
#include <string>
#include <utility>

#include <Kokkos_Core.hpp>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "numpy_view.h"
#include "pcms/configuration.h"
#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"

namespace pcms
{

namespace py = pybind11;

// pybind11 holders cannot be shared_ptr<const T>. Python only reaches the
// const interface of these objects, so they are held as shared_ptr<T> and
// converted at the boundary.
template <typename T>
std::shared_ptr<T> ToPyHolder(const std::shared_ptr<const T>& p)
{
  return std::const_pointer_cast<T>(p);
}

// Coordinate system arguments default to None in Python, meaning Cartesian in
// the dimension of the data they are used with. C++ resolves that dimension
// from a deferred Cartesian system, which Python never sees.
inline std::shared_ptr<const CoordinateSystem> CartesianIfNone(
  std::shared_ptr<const CoordinateSystem> system)
{
  return system ? std::move(system) : csys::Cartesian::Deferred();
}

// Coordinates tagged with the coordinate system that gives them meaning.
// `values` is a read-only 2D float64 array (num_points × dim) that pcms owns
// or keeps alive.
struct PyCoordinateView
{
  std::shared_ptr<const CoordinateSystem> system;
  py::array_t<Real> values;
};

// Values tagged with the basis their components are expressed in. `values` is
// a 2D array (num_points × num_components): the caller's own array, or a
// read-only one returned by pcms. Calls that write it require it to be a
// writeable C-contiguous float64 array.
struct PyValueView
{
  ValueBasis basis;
  py::array values;
};

// Python view of C++ coordinates, sharing memory with them when they are host
// accessible; `owner` must keep that memory alive.
template <typename MemorySpace, typename Layout>
PyCoordinateView ToPyCoordinateView(
  const CoordinateView<MemorySpace, Layout>& view, py::handle owner)
{
  return PyCoordinateView{view.GetCoordinateSystem(),
                          ViewToNumpy(view.GetValues(), owner)};
}

// Python view of C++ values, sharing memory with them when they are host
// accessible; `owner` must keep that memory alive.
template <typename Element, typename MemorySpace, typename Layout>
PyValueView ToPyValueView(const ValueView<Element, MemorySpace, Layout>& view,
                          py::handle owner)
{
  return PyValueView{view.GetBasis(), ViewToNumpy(view.GetValues(), owner)};
}

// Device copy of the coordinates together with the CoordinateView over it.
// The view is valid while `storage` is alive.
struct DeviceCoordinates
{
  Kokkos::View<Real**, DeviceMemorySpace> storage;
  CoordinateView<DeviceMemorySpace> view;
};

inline DeviceCoordinates CopyCoordinatesToDevice(const PyCoordinateView& coords)
{
  auto storage = CopyNumpyToKokkosView<2, Real, DeviceMemorySpace>(
    NumpyInput<Real>::ensure(coords.values), "coordinates");
  return DeviceCoordinates{storage, CoordinateView<DeviceMemorySpace>(
                                      coords.system, MakeRank2View(storage))};
}

// Number of components of per-entity data shaped (n, num_components), where a
// 1D array of length n has one component. The call the data is passed to
// checks the extents against what it expects.
inline LO PerEntityComponents(const py::array& arr)
{
  if (arr.ndim() == 1) {
    return 1;
  }
  if (arr.ndim() == 2) {
    return static_cast<LO>(arr.shape(1));
  }
  throw std::invalid_argument(
    "expected an array of shape (n, num_components) or (n,), got " +
    detail::ShapeString(arr));
}

// Calls f with a device view over the caller's writeable output buffer and
// makes what f writes visible in that buffer. On host builds f writes the
// buffer in place; otherwise it writes a device copy that is copied back.
template <typename F>
void WithDeviceOutput(py::array& output, F&& f)
{
  auto output_host =
    NumpyOutputToHostComponentView<Real>(output, PerEntityComponents(output));
#ifdef PCMS_HAS_DISTINCT_DEVICE_MEMORY_SPACE
  // TODO: write in place on the device via DLPack.
  auto output_view = MakeUnmanagedKokkosView(output_host);
  auto output_device = Kokkos::View<Real**, DeviceMemorySpace>(
    "output_device", output_view.extent(0), output_view.extent(1));
  DeepCopyMismatchLayouts(output_device, output_view);
  std::forward<F>(f)(MakeRank2View(output_device));
  DeepCopyMismatchLayouts(output_view, output_device);
#else
  std::forward<F>(f)(output_host);
#endif
}

// Scalar values carry no basis, so for rank-0 values a plain array stands in
// for a ValueView. Values of higher rank must declare their basis.
inline void RequireScalarForPlainArray(const char* context,
                                       const ValueBasis& basis)
{
  if (basis.Rank() != 0) {
    throw py::type_error(std::string(context) + ": rank-" +
                         std::to_string(basis.Rank()) +
                         " values need a ValueView declaring their basis; a "
                         "plain array is accepted only for scalar values");
  }
}

} // namespace pcms

#endif // PCMS_PYTHONAPI_PY_VALUE_TYPES_H
