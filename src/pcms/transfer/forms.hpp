#ifndef PCMS_TRANSFER_FORMS_HPP
#define PCMS_TRANSFER_FORMS_HPP

#include "pcms/field/coordinate_map.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/geometry/quadrature.hpp"
#include "pcms/utility/types.hpp"
#include <Kokkos_Core.hpp>
#include <memory>

namespace pcms::transfer::detail
{

/// @brief A Real function of `space` whose basis is a scalar, or an array of
///        its components if it has several.
inline Function<Real> ScratchFunction(const FunctionSpace& space)
{
  ValueBasis basis;
  if (space.GetLayout()->GetNumComponents() > 1) {
    basis.kind = ValueKind::Array;
  }
  return space.CreateFunction<Real>("", std::move(basis));
}

// The internal forms of the Galerkin projections, written on the evaluation of
// a space at the points of a quadrature (E) and the quadrature weights (W);
// they know nothing of the kind of the space. Their interfaces are meant to
// be published later.

/// The linear form v -> integral of f v over a quadrature, for samples f at
/// its points: the functional E^T W f of a space.
class LinearForm
{
public:
  /// @param quadrature points bound on the space's geometry, one weight each
  /// @throws pcms_error if the quadrature is not on the space's geometry
  LinearForm(std::shared_ptr<const FunctionSpace> space, Quadrature quadrature);

  /// functional = E^T W samples (the evaluator must support the transpose).
  /// `samples` holds one row per point and one
  /// column per component, tagged with the basis of the functional's data;
  /// the functional is a field of the space's layout (overwritten).
  void Assemble(ValueView<const Real, DeviceMemorySpace> samples,
                Field<Real>& functional) const;

  /// W samples, in a scratch the form owns (valid until the next call)
  [[nodiscard]] ValueView<const Real, DeviceMemorySpace> Weight(
    ValueView<const Real, DeviceMemorySpace> samples) const;

  [[nodiscard]] const std::shared_ptr<const FunctionSpace>& Space()
    const noexcept
  {
    return space_;
  }
  [[nodiscard]] const Quadrature& GetQuadrature() const noexcept
  {
    return quadrature_;
  }
  /// E: the space at the quadrature's points
  [[nodiscard]] const PointEvaluator<Real>& Evaluator() const noexcept
  {
    return *evaluator_;
  }

private:
  std::shared_ptr<const FunctionSpace> space_;
  Quadrature quadrature_;
  std::unique_ptr<PointEvaluator<Real>> evaluator_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> weighted_; // [points][comps]
};

/// The mass form (u, v) -> integral of u v over a quadrature: E^T W E.
class BilinearForm
{
public:
  /// @throws pcms_error as LinearForm, or if the evaluator of the space does
  ///         not support the transpose
  BilinearForm(std::shared_ptr<const FunctionSpace> space,
               Quadrature quadrature);

  /// functional = E^T W E u (overwritten)
  void Apply(const Field<Real>& u, Field<Real>& functional) const;

  /// @return whether the evaluator provides the matrix E (the diagonal and
  ///         the COO assembly need it)
  [[nodiscard]] bool SupportsMatrix() const noexcept
  {
    return form_.Evaluator().SupportsMatrix();
  }

  /// diagonal = the diagonal of E^T W E as a field of the space's layout.
  /// Single rank: contributions to ghost holders are not reduced.
  /// @throws pcms_error unless SupportsMatrix() or if the layout has ghosts
  void AssembleDiagonal(Field<Real>& diagonal) const;

  /// lumped = E^T W E 1, the row sums
  void AssembleLumped(Field<Real>& lumped) const;

  /// E^T W E as coordinate entries with global DOF-holder ids; entries with
  /// equal indices are to be summed
  struct Coo
  {
    Kokkos::View<GO*, DeviceMemorySpace> rows;
    Kokkos::View<GO*, DeviceMemorySpace> columns;
    Kokkos::View<Real*, DeviceMemorySpace> values;
  };
  /// @throws pcms_error unless SupportsMatrix()
  [[nodiscard]] Coo AssembleCoo() const;

  [[nodiscard]] const LinearForm& Form() const noexcept { return form_; }

private:
  LinearForm form_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> samples_; // [points][comps]
};

} // namespace pcms::transfer::detail

#endif // PCMS_TRANSFER_FORMS_HPP
