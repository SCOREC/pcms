#ifndef PCMS_PYTHONAPI_PY_VALUE_TYPES_H
#define PCMS_PYTHONAPI_PY_VALUE_TYPES_H

// Python-side counterparts of the non-owning C++ view types. A C++
// CoordinateView or ValueView borrows its memory, so the Python objects pair
// the same tag (coordinate system or value basis) with a numpy array that owns
// or keeps alive that memory, and are converted to the C++ view at each call.

#include <memory>
#include <optional>
#include <string>
#include <utility>

#include <Kokkos_Core.hpp>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "numpy_view.h"
#include "pcms/configuration.h"
#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"
#include "pcms/field/coordinate_view.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/types.h"

namespace pcms
{

namespace py = pybind11;

// Handle of an object that is immutable once built: Python holds it through
// the const pointer, so the const interface is the only one the bindings can
// reach, and the compiler enforces that. pybind11 cannot bind this directly:
// a bound class and its holder cannot be const, so the class is the handle.
template <typename T>
struct Immutable
{
  std::shared_ptr<const T> ptr;
};

// The Python class of a handle of T. Two handles to one object are two Python
// objects, so they compare equal, and hash alike, when they hold the same
// object (not when they are the same Python object, `is`).
template <typename T>
py::class_<Immutable<T>> ImmutableClass(py::handle scope, const char* name,
                                        const char* doc = nullptr)
{
  py::class_<Immutable<T>> cls(scope, name, doc);
  cls
    .def(
      "__eq__",
      [](const Immutable<T>& a, const Immutable<T>& b) {
        return a.ptr.get() == b.ptr.get();
      },
      py::is_operator())
    .def("__hash__", [](const Immutable<T>& self) {
      return std::hash<const T*>{}(self.ptr.get());
    });
  return cls;
}

// The Python classes of the concrete coordinate systems derive from
// CoordinateSystem, so Python sees the kind of any system pcms returns
// (isinstance(cs, pcms.csys.Cartesian)); their handles derive likewise.
struct CartesianHandle : Immutable<CoordinateSystem>
{};
struct CylindricalRThetaZHandle : Immutable<CoordinateSystem>
{};
struct CylindricalRZHandle : Immutable<CoordinateSystem>
{};

// The Python object of a coordinate system, of the class of its kind; None for
// null.
inline py::object ToPython(
  const std::shared_ptr<const CoordinateSystem>& system)
{
  if (system == nullptr) {
    return py::none();
  }
  if (dynamic_cast<const csys::Cartesian*>(system.get()) != nullptr) {
    return py::cast(CartesianHandle{{system}});
  }
  if (dynamic_cast<const csys::CylindricalRThetaZ*>(system.get()) != nullptr) {
    return py::cast(CylindricalRThetaZHandle{{system}});
  }
  if (dynamic_cast<const csys::CylindricalRZ*>(system.get()) != nullptr) {
    return py::cast(CylindricalRZHandle{{system}});
  }
  return py::cast(Immutable<CoordinateSystem>{system});
}

// Coordinate system arguments default to None in Python, meaning Cartesian in
// the dimension of the data they are used with. C++ resolves that dimension
// from a deferred Cartesian system, which Python never sees.
inline std::shared_ptr<const CoordinateSystem> CartesianIfNone(
  const std::optional<Immutable<CoordinateSystem>>& system)
{
  return system ? system->ptr : csys::Cartesian::Deferred();
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

/// @brief Read-only NumPy array over device point indices.
/// @param indices indices to expose
/// @param owner Python object that keeps `indices` alive
/// @return the array (empty if there are no indices)
inline py::array_t<LO> IndicesToNumpy(
  Kokkos::View<const LO*, DeviceMemorySpace> indices, py::object owner)
{
  if (indices.extent(0) == 0) {
    return py::array_t<LO>(0);
  }
  return ViewToNumpy(MakeRank1View(indices), owner);
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
