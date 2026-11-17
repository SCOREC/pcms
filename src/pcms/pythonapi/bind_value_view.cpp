#include <memory>
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pcms/field/coordinate_system.hpp"
#include "pcms/field/value_view.hpp"
#include "py_value_types.h"

namespace py = pybind11;

namespace pcms
{

namespace
{

VarianceSignature MakeVarianceSignature(const std::vector<Variance>& v)
{
  switch (v.size()) {
    case 0: return VarianceSignature{};
    case 1: return VarianceSignature{v[0]};
    case 2: return VarianceSignature{v[0], v[1]};
    case 3: return VarianceSignature{v[0], v[1], v[2]};
    case 4: return VarianceSignature{v[0], v[1], v[2], v[3]};
    default: break;
  }
  throw py::value_error("VarianceSignature: rank must be 0.." +
                        std::to_string(VarianceSignature::max_rank) + ", got " +
                        std::to_string(v.size()));
}

std::string Repr(const VarianceSignature& v)
{
  std::string s = "[";
  for (int i = 0; i < v.Rank(); ++i) {
    s += (i > 0 ? ", " : "");
    s += v[i] == Variance::Contravariant ? "Contravariant" : "Covariant";
  }
  return s + "]";
}

std::string Repr(const ValueBasis& b)
{
  if (b.Rank() == 0) {
    return "<pcms.ValueBasis scalar>";
  }
  return "<pcms.ValueBasis " +
         (b.system ? std::string(b.system->Kind()) : std::string("None")) +
         " " +
         (b.scaling == ComponentScaling::Physical ? "Physical" : "Natural") +
         " " + Repr(b.variance) + ">";
}

} // namespace

void bind_value_view_module(py::module& m)
{
  py::enum_<Variance>(m, "Variance",
                      "Index variance: Contravariant components are on the "
                      "basis dx/dq^i, Covariant on the dual grad q^i")
    .value("Contravariant", Variance::Contravariant)
    .value("Covariant", Variance::Covariant);

  py::enum_<ComponentScaling>(
    m, "ComponentScaling",
    "How components are scaled relative to the coordinate-induced basis")
    .value("Physical", ComponentScaling::Physical,
           "Components on the unit basis; all have the same units")
    .value("Natural", ComponentScaling::Natural,
           "Components on the unscaled coordinate basis");

  py::enum_<FieldValueType>(m, "FieldValueType")
    .value("Scalar", FieldValueType::Scalar)
    .value("Vector", FieldValueType::Vector)
    .value("Tensor", FieldValueType::Tensor);

  py::class_<VarianceSignature>(m, "VarianceSignature",
                                "The variance of each index of a value")
    .def(py::init<>(), "Rank 0 (scalar)")
    .def(py::init(&MakeVarianceSignature), py::arg("variance"),
         "From a list of Variance, one per index")
    .def("rank", &VarianceSignature::Rank)
    .def("__len__", &VarianceSignature::Rank)
    .def("__getitem__",
         [](const VarianceSignature& self, int i) {
           if (i < 0) {
             i += self.Rank();
           }
           if (i < 0 || i >= self.Rank()) {
             throw py::index_error("VarianceSignature index out of range");
           }
           return self[i];
         })
    .def(
      "__eq__",
      [](const VarianceSignature& a, const VarianceSignature& b) {
        return a == b;
      },
      py::is_operator())
    .def("__repr__", [](const VarianceSignature& self) {
      return "<pcms.VarianceSignature " + Repr(self) + ">";
    });

  auto values = m.def_submodule("values", "Common variance signatures");
  values.attr("Scalar") = values::Scalar;
  values.attr("Vector") = values::Vector;
  values.attr("Covector") = values::Covector;
  values.attr("Tensor") = values::Tensor;
  values.def("Of", &MakeVarianceSignature, py::arg("variance"),
             "Variance signature with the given variance per index");

  py::class_<ValueBasis>(
    m, "ValueBasis",
    "The basis a value's components are expressed in. The default is the "
    "scalar basis.")
    .def(py::init([](const std::shared_ptr<CoordinateSystem>& system,
                     ComponentScaling scaling, VarianceSignature variance) {
           return ValueBasis{system, scaling, variance};
         }),
         py::arg("system") = py::none(),
         py::arg("scaling") = ComponentScaling::Physical,
         py::arg("variance") = VarianceSignature{})
    .def_property(
      "system", [](const ValueBasis& self) { return ToPyHolder(self.system); },
      [](ValueBasis& self, const std::shared_ptr<CoordinateSystem>& system) {
        self.system = system;
      })
    .def_readwrite("scaling", &ValueBasis::scaling)
    .def_readwrite("variance", &ValueBasis::variance)
    .def("rank", &ValueBasis::Rank)
    .def(
      "__eq__",
      [](const ValueBasis& a, const ValueBasis& b) {
        return SameValueBasis(a, b);
      },
      py::is_operator(),
      "Whether the two bases describe the same components (SameValueBasis)")
    .def("__repr__", [](const ValueBasis& self) { return Repr(self); });

  m.def("same_value_basis", &SameValueBasis, py::arg("a"), py::arg("b"));

  py::class_<PyValueView>(m, "ValueView",
                          "Values (num_points × num_components) tagged with "
                          "the basis their components are in")
    .def(py::init([](const ValueBasis& basis, const py::object& values) {
           py::array array = py::array::ensure(values);
           if (!array) {
             throw py::type_error("ValueView: values must be array-like");
           }
           if (array.ndim() != 2) {
             throw py::value_error("ValueView: expected a 2D array, got a " +
                                   std::to_string(array.ndim()) + "D array");
           }
           detail::ValidateValueSemantics(basis,
                                          static_cast<int>(array.shape(1)));
           return PyValueView{basis, array};
         }),
         py::arg("basis"), py::arg("values"),
         "Tag a 2D array with a basis. The array is not copied: a ValueView "
         "passed as an output is written in place.")
    .def("get_basis", [](const PyValueView& self) { return self.basis; })
    .def("get_values", [](const PyValueView& self) { return self.values; })
    .def("rank", [](const PyValueView& self) { return self.basis.Rank(); })
    .def("get_value_type",
         [](const PyValueView& self) {
           return ValueTypeOfRank(self.basis.Rank());
         })
    .def(
      "extent",
      [](const PyValueView& self, int r) {
        if (r < 0 || r > 1) {
          throw py::index_error("ValueView.extent: rank must be 0 or 1");
        }
        return self.values.shape(r);
      },
      py::arg("r"))
    .def("size", [](const PyValueView& self) { return self.values.size(); })
    .def(
      "__array__",
      [](const PyValueView& self, const py::object& dtype,
         const py::object& copy) {
        auto np = py::module_::import("numpy");
        if (!copy.is_none() && copy.cast<bool>()) {
          return np.attr("array")(self.values, dtype);
        }
        return np.attr("asarray")(self.values, dtype);
      },
      py::arg("dtype") = py::none(), py::arg("copy") = py::none(),
      "The values, so numpy functions accept a ValueView")
    .def("__repr__", [](const PyValueView& self) {
      return "<pcms.ValueView " + std::to_string(self.values.shape(0)) + "x" +
             std::to_string(self.values.shape(1)) + " " + Repr(self.basis) +
             ">";
    });
}

} // namespace pcms
