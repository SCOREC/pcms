#ifndef PCMS_TRANSFER_MASS_SMOOTHER_HPP
#define PCMS_TRANSFER_MASS_SMOOTHER_HPP

#include <memory>
#include <vector>

#include <petscvec.h>

#include "pcms/field/field_layout.h"
#include "pcms/field/field_operator.hpp"
#include "pcms/field/function_space.h"
#include "pcms/transfer/bilinear_form_integrator.hpp"
#include "pcms/transfer/mass_matrix_type.hpp"
#include "pcms/utility/types.h"

namespace pcms
{

/// Conservative, bounds-preserving smoother: one Apply performs one sweep
/// u <- inv(M_L) M u with M the consistent mass matrix of the layout and M_L
/// its row sums (Farrell et al., CMAME 198 (2009), Eqs. 37-38, applied to
/// the whole field). Call Apply repeatedly for more sweeps.
///
/// Construction throws pcms_error if the matrix is diagonal, has a
/// non-positive row sum, or has a negative entry, since each breaks the
/// smoothing/convexity guarantees. Serial (PETSC_COMM_SELF) only.
class MassSmoother : public FieldOperator<Real>
{
public:
  /// @param layout the layout the smoothed fields live on
  /// @param mass consistent mass matrix of that layout, rows indexed by the
  ///        layout's permuted index; ownership is taken
  MassSmoother(std::shared_ptr<const FieldLayout> layout,
               std::unique_ptr<BilinearFormIntegrator> mass);
  MassSmoother(const FunctionSpace& space,
               std::unique_ptr<BilinearFormIntegrator> mass);
  ~MassSmoother() override;

  MassSmoother(const MassSmoother&) = delete;
  MassSmoother& operator=(const MassSmoother&) = delete;

  /// One sweep, in place. Each component is smoothed with the same matrix.
  /// Throws pcms_error if the field is not on the construction layout.
  void Apply(Field<Real>& field) const override;

private:
  std::shared_ptr<const FieldLayout> layout_;
  std::unique_ptr<BilinearFormIntegrator> mass_;
  Vec lumped_mass_ = nullptr;
  Vec q_ = nullptr;
  Vec tmp_ = nullptr;
  mutable std::vector<Real> host_values_;
};

/// Consistent mass integrator for the space, dispatched on its layout type.
/// Throws pcms_error for layouts without a mass integrator.
std::unique_ptr<BilinearFormIntegrator> BuildMassIntegrator(
  const FunctionSpace& space,
  MassMatrixType mass_type = MassMatrixType::Consistent);

/// MassSmoother over BuildMassIntegrator(space).
std::unique_ptr<MassSmoother> BuildMassSmoother(const FunctionSpace& space);

} // namespace pcms

#endif // PCMS_TRANSFER_MASS_SMOOTHER_HPP
