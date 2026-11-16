#include <pybind11/pybind11.h>
#include <pcms/field/evaluator/mls_interpolation.hpp>

namespace py = pybind11;

namespace pcms
{

void bind_mls_interpolation_module(py::module& m)
{
  // RadialBasisFunction is part of the MLSOptions API used by
  // PolynomialReconstructionFunctionSpace.
  py::enum_<pcms::RadialBasisFunction>(m, "RadialBasisFunction")
    .value("RBF_GAUSSIAN", pcms::RadialBasisFunction::RBF_GAUSSIAN)
    .value("RBF_C4", pcms::RadialBasisFunction::RBF_C4)
    .value("RBF_CONST", pcms::RadialBasisFunction::RBF_CONST)
    .value("NO_OP", pcms::RadialBasisFunction::NO_OP)
    .export_values();
}

} // namespace pcms
