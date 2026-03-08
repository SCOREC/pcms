#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Omega_h_graph.hpp>
#include <pcms/transfer/mls_interpolation.hpp>
#include <pcms/transfer/adj_search.hpp>
#include "numpy_array_transform.h"

namespace py = pybind11;

namespace pcms
{

void bind_mls_interpolation_module(py::module& m)
{
  py::class_<SupportResults>(m, "SupportResults")
    .def(py::init<>())
    .def(py::init([](py::array_t<Omega_h::LO> supports_ptr,
                     py::array_t<Omega_h::LO> supports_idx,
                     py::array_t<Omega_h::Real> radii2) {
           auto supports_ptr_read =
             numpy_to_omega_h_read<Omega_h::LO>(supports_ptr);
           auto supports_idx_read =
             numpy_to_omega_h_read<Omega_h::LO>(supports_idx);
           auto radii2_write = numpy_to_omega_h_write<Omega_h::Real>(radii2);
           return SupportResults{supports_ptr_read, supports_idx_read,
                                 radii2_write};
         }),
         py::arg("supports_ptr"), py::arg("supports_idx"), py::arg("radii2"))
    .def_property(
      "supports_ptr",
      [](const SupportResults& self) {
        return omega_h_read_to_numpy(self.supports_ptr);
      },
      [](SupportResults& self, py::array_t<Omega_h::LO> value) {
        self.supports_ptr = numpy_to_omega_h_read<Omega_h::LO>(value);
      })
    .def_property(
      "supports_idx",
      [](const SupportResults& self) {
        return omega_h_read_to_numpy(self.supports_idx);
      },
      [](SupportResults& self, py::array_t<Omega_h::LO> value) {
        self.supports_idx = numpy_to_omega_h_read<Omega_h::LO>(value);
      })
    .def_property(
      "radii2",
      [](const SupportResults& self) {
        return omega_h_read_to_numpy(Omega_h::read(self.radii2));
      },
      [](SupportResults& self, py::array_t<Omega_h::Real> value) {
        self.radii2 = numpy_to_omega_h_write<Omega_h::Real>(value);
      });

  py::enum_<pcms::RadialBasisFunction>(m, "RadialBasisFunction")
    .value("RBF_GAUSSIAN", pcms::RadialBasisFunction::RBF_GAUSSIAN)
    .value("RBF_C4", pcms::RadialBasisFunction::RBF_C4)
    .value("RBF_CONST", pcms::RadialBasisFunction::RBF_CONST)
    .value("NO_OP", pcms::RadialBasisFunction::NO_OP)
    .export_values();

  m.def(
    "mls_interpolation",
    [](py::array_t<Omega_h::Real> source_values,
       py::array_t<Omega_h::Real> source_coordinates,
       py::array_t<Omega_h::Real> target_coordinates,
       const SupportResults& support, Omega_h::LO dim, Omega_h::LO degree,
       pcms::RadialBasisFunction bf, double lambda_reg, double tol,
       double decay_factor) {
      auto source_values_read =
        numpy_to_omega_h_read<Omega_h::Real>(source_values);
      auto source_coordinates_read =
        numpy_to_omega_h_read<Omega_h::Real>(source_coordinates);
      auto target_coordinates_read =
        numpy_to_omega_h_read<Omega_h::Real>(target_coordinates);

      auto interpolated_values = pcms::mls_interpolation(
        source_values_read, source_coordinates_read, target_coordinates_read,
        support, dim, degree, bf, lambda_reg, tol, decay_factor);

      return omega_h_read_to_numpy(Omega_h::Reals(interpolated_values));
    },
    py::arg("source_values"), py::arg("source_coordinates"),
    py::arg("target_coordinates"), py::arg("support"), py::arg("dim"),
    py::arg("degree"),
    py::arg("radial_basis_function") = pcms::RadialBasisFunction::NO_OP,
    py::arg("lambda_reg") = 0.0, py::arg("tol") = 1e-6,
    py::arg("decay_factor") = 5.0);
}

} // namespace pcms
