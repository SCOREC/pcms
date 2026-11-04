#include "pcms/field/function_spaces.hpp"
#include "pcms/configuration.hpp"
#include "pcms/field/function_space/polynomial_reconstruction.hpp"
#include "pcms/field/layout/empty.hpp"
#include "pcms/utility/assert.hpp"
#ifdef PCMS_ENABLE_OMEGA_H
#include "pcms/backends/omega_h/discretization.hpp"
#include "pcms/backends/uniform_grid/discretization.hpp"
#include "pcms/field/function_space/lagrange.hpp"
#include "pcms/field/function_space/spline.hpp"
#endif

#include <string>
#include <typeinfo>

namespace pcms
{
namespace
{
[[noreturn]] void Unsupported(const char* factory, const Discretization& d)
{
  throw pcms_error(std::string(factory) +
                   ": unsupported discretization type " + typeid(d).name());
}

void RequireDiscretization(const char* factory,
                           const std::shared_ptr<const Discretization>& d)
{
  if (d == nullptr) {
    throw pcms_error(std::string(factory) + ": discretization is null");
  }
}
} // namespace

std::shared_ptr<FunctionSpace> CreateLagrangeSpace(
  std::shared_ptr<const Discretization> discretization, int order,
  int num_components, std::shared_ptr<const CoordinateSystem> coordinate_system,
  LagrangeOptions options)
{
  RequireDiscretization("CreateLagrangeSpace", discretization);
#ifdef PCMS_ENABLE_OMEGA_H
  using backend::omega_h::OmegaHDiscretization;
  using backend::uniform_grid::UniformGridDiscretization;
  if (const auto* oh =
        dynamic_cast<const OmegaHDiscretization*>(discretization.get())) {
    auto impl = pcms::field::LagrangeFunctionSpace::DefaultBackend;
    if (options.implementation == LagrangeImplementation::MeshFields) {
      impl = pcms::field::LagrangeFunctionSpace::Backend::MeshFields;
    } else if (options.implementation == LagrangeImplementation::OmegaH) {
      impl = pcms::field::LagrangeFunctionSpace::Backend::OmegaH;
    }
    if (oh->GetOwnedMask().has_value()) {
      return pcms::field::LagrangeFunctionSpace::FromMesh(
        oh->GetMesh(), order, num_components, std::move(coordinate_system),
        *oh->GetOwnedMask(), oh->GetGlobalIdName(), impl,
        std::move(options.layout_name));
    }
    return pcms::field::LagrangeFunctionSpace::FromMesh(
      oh->GetMesh(), order, num_components, std::move(coordinate_system),
      oh->GetGlobalIdName(), impl, std::move(options.layout_name));
  }
  if (options.implementation != LagrangeImplementation::Default) {
    throw pcms_error("CreateLagrangeSpace: LagrangeOptions::implementation "
                     "only applies to Omega_h discretizations");
  }
  if (const auto* ug = dynamic_cast<const UniformGridDiscretization<2>*>(
        discretization.get())) {
    return pcms::field::LagrangeFunctionSpace::FromUniformGrid(
      ug->GetGrid(), num_components, std::move(coordinate_system), order,
      std::move(options.layout_name));
  }
  if (const auto* ug = dynamic_cast<const UniformGridDiscretization<3>*>(
        discretization.get())) {
    return pcms::field::LagrangeFunctionSpace::FromUniformGrid(
      ug->GetGrid(), num_components, std::move(coordinate_system), order,
      std::move(options.layout_name));
  }
#else
  (void)order;
  (void)num_components;
  (void)coordinate_system;
  (void)options;
#endif
  Unsupported("CreateLagrangeSpace", *discretization);
}

std::shared_ptr<FunctionSpace> CreateSplineSpace(
  std::shared_ptr<const Discretization> discretization,
  std::shared_ptr<const CoordinateSystem> coordinate_system)
{
  RequireDiscretization("CreateSplineSpace", discretization);
#ifdef PCMS_ENABLE_OMEGA_H
  if (const auto* ug =
        dynamic_cast<const backend::uniform_grid::UniformGridDiscretization<2>*>(
          discretization.get())) {
    return pcms::field::SplineFunctionSpace::FromUniformGrid(ug->GetGrid(),
                                                std::move(coordinate_system));
  }
#else
  (void)coordinate_system;
#endif
  Unsupported("CreateSplineSpace", *discretization);
}

std::shared_ptr<FunctionSpace> CreatePolynomialReconstructionSpace(
  CoordinateView<HostMemorySpace> points, MLSOptions options)
{
  return pcms::field::PolynomialReconstructionFunctionSpace::Create(std::move(points),
                                                       options);
}

std::shared_ptr<FunctionSpace> CreatePolynomialReconstructionSpace(
  std::shared_ptr<const Discretization> discretization, int source_entity_dim,
  std::shared_ptr<const CoordinateSystem> coordinate_system,
  MLSOptions options)
{
  RequireDiscretization("CreatePolynomialReconstructionSpace", discretization);
#ifdef PCMS_ENABLE_OMEGA_H
  if (const auto* oh =
        dynamic_cast<const backend::omega_h::OmegaHDiscretization*>(
          discretization.get())) {
    return pcms::field::PolynomialReconstructionFunctionSpace::FromMesh(
      oh->GetMesh(), source_entity_dim, std::move(coordinate_system), options);
  }
#else
  (void)source_entity_dim;
  (void)coordinate_system;
  (void)options;
#endif
  Unsupported("CreatePolynomialReconstructionSpace", *discretization);
}

std::shared_ptr<const FieldLayout> CreateEmptyLayout(
  std::shared_ptr<const CoordinateSystem> coordinate_system, int num_components)
{
  return std::make_shared<pcms::field::EmptyFieldLayout>(std::move(coordinate_system),
                                            num_components);
}

} // namespace pcms
