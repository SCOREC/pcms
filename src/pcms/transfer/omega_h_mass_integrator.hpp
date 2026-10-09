#ifndef PCMS_TRANSFER_OMEGA_H_MASS_INTEGRATOR_HPP
#define PCMS_TRANSFER_OMEGA_H_MASS_INTEGRATOR_HPP

#include "pcms/field/function_space.h"
#include "pcms/field/layout/omega_h_lagrange.h"
#include "pcms/transfer/bilinear_form_integrator.hpp"
#include "pcms/transfer/mass_matrix_type.hpp"
#include <memory>
#include <vector>

namespace pcms
{

// Mass-matrix bilinear form integrator for Lagrange spaces on Omega_h 2D
// simplex meshes. Supports order-0 (diagonal, M_ee = area(e)) and order-1
// (consistent 3x3 element blocks) target spaces. With
// MassMatrixType::Lumped the order-1 matrix is row-sum lumped into a
// diagonal matrix (order-0 is already diagonal, so lumping is a no-op).
//
// The matrix is assembled once at construction into a PETSc sparse AIJ matrix.
// Sparsity is derived directly from the element node GIDs — no separate COO
// utility is needed. The integrator owns the matrix; callers receive a
// non-owning handle via GetMatrix(). PETSc reference-counts the matrix, so a
// KSP that calls KSPSetOperators(GetMatrix()) safely outlives this object.
class OmegaHMassIntegrator : public BilinearFormIntegrator
{
public:
  explicit OmegaHMassIntegrator(
    const FunctionSpace& target_space,
    MassMatrixType mass_type = MassMatrixType::Consistent);
  OmegaHMassIntegrator(
    std::shared_ptr<const OmegaHLagrangeLayout> target_layout,
    std::shared_ptr<const CoordinateSystem> coordinate_system,
    MassMatrixType mass_type = MassMatrixType::Consistent);
  ~OmegaHMassIntegrator();

  Mat GetMatrix() const noexcept override;
  bool IsDiagonal() const noexcept override;

private:
  Mat mat_ = nullptr;
  bool diagonal_ = false;
};

// Builds an OmegaHMassIntegrator for the given target space.
// target_space must use OmegaHLagrangeLayout, scalar, Cartesian, order-0 or
// order-1.
std::unique_ptr<BilinearFormIntegrator> BuildOmegaHMassIntegrator(
  const FunctionSpace& target_space,
  MassMatrixType mass_type = MassMatrixType::Consistent);

} // namespace pcms

#endif // PCMS_TRANSFER_OMEGA_H_MASS_INTEGRATOR_HPP
