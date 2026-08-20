#ifndef PCMS_TRANSFER_TRANSFORMED_TRANSFER_OPERATOR_HPP
#define PCMS_TRANSFER_TRANSFORMED_TRANSFER_OPERATOR_HPP

#include "pcms/field/basis_transformation.hpp"
#include "pcms/field/field.h"
#include "pcms/field/field_data.h"
#include "pcms/field/function_space.h"
#include "pcms/field/out_of_bounds_policy.h"
#include "pcms/transfer/transfer_operator.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include "pcms/utility/profile.h"
#include <Kokkos_Core.hpp>
#include <memory>
#include <type_traits>
#include <utility>

namespace pcms
{

namespace detail
{

// Restores rows whose pre-transformation components all equal the FILL
// sentinel (NaN sentinels included), so a basis rotation never mixes sentinel
// values into real components. Free function because CUDA extended lambdas
// may not be defined inside private member functions.
//
// NOTE: detection compares against the sentinel because evaluators expose no
// out-of-bounds mask — an in-bounds row whose components all exactly equal
// the sentinel would be misclassified. Obtaining a real mask is future work
// (the Map-time PointStatus of coordinate maps covers map-domain
// violations, not evaluator-mesh ones).
inline void RestoreFillRows(
  const Kokkos::View<Real**, DeviceMemorySpace>& pre_transform,
  Rank2View<Real, DeviceMemorySpace> values, Real fill_value)
{
  const auto num_cols = static_cast<int>(pre_transform.extent(1));
  Kokkos::parallel_for(
    "transformed_transfer_restore_fill_rows",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(pre_transform.extent(0))),
    KOKKOS_LAMBDA(const LO i) {
      bool filled = true;
      for (int c = 0; c < num_cols; ++c) {
        const Real value = pre_transform(i, c);
        filled =
          filled && ((value == fill_value) ||
                     (Kokkos::isnan(value) && Kokkos::isnan(fill_value)));
      }
      if (filled) {
        for (int c = 0; c < num_cols; ++c) {
          values(i, c) = fill_value;
        }
      }
    });
}

} // namespace detail

// Composes a value-basis transformation around a transfer operator: run the
// wrapped transfer, then re-express the transferred components in the
// target's declared basis at the target DOF-holder points. Both spaces share
// one coordinate system; only the basis the components are expressed in
// changes (borrowed-basis storage over an ordinary transfer).
//
// The wrapped operator is named at construction and built here, never handed
// in already built: this operator binds the transformation at the target's
// DOF-holder points and then trusts the wrapped operator's output rows to be
// exactly those holders in that order. Owning construction is what makes that
// correspondence hold; a caller supplying a finished operator would have to
// know the point set to keep it, and nothing could check that they had.
//
// The wrapped operator's uncommitted apply is what lets this wrapper rotate
// without a staging round trip.
//
// Transfer BETWEEN coordinate systems is not composed here. Until a coupling
// needs it as an operator, it is built by hand: map the query points with a
// CoordinateMap, create the source's point evaluator on the mapped points,
// and rotate the evaluated components with a transformation bound to the
// same points (see the manual-composition test).
//
// The caller supplies an unbound rule naming the mathematics; this operator
// binds it at the points it actually writes. Value semantics come entirely
// from the two fields at Apply time: the basis actually written must match
// the target's declared basis.
//
// The wrapped operator type is named per constructor call with
// std::in_place_type (explicit template arguments on constructors are not
// expressible in C++), so a TransformedTransferOperator<T> is one concrete
// type whichever operator it wraps. Wrappable operators (Interpolator) are
// TransferOperator<T>s constructible as InnerOp(src, tgt, policy, args...).
template <typename T>
class TransformedTransferOperator : public TransferOperator<T>
{
public:
  /// Builds an InnerOp over (source_space, target_space, policy,
  /// inner_args...) and binds value_transformation at the target layout's
  /// DOF-holder coordinates.
  template <typename InnerOp, typename... Args>
  TransformedTransferOperator(
    std::in_place_type_t<InnerOp>, const FunctionSpace& source_space,
    const FunctionSpace& target_space,
    std::shared_ptr<const BasisTransformation> value_transformation,
    OutOfBoundsPolicy policy = {}, Args&&... inner_args)
    : policy_(policy),
      inner_(std::make_unique<InnerOp>(source_space, target_space, policy,
                                       std::forward<Args>(inner_args)...))
  {
    static_assert(std::is_base_of_v<TransferOperator<T>, InnerOp>,
                  "InnerOp must be a TransferOperator<T>");
    if (value_transformation == nullptr) {
      throw pcms_error(
        "TransformedTransferOperator: value_transformation must not be null "
        "(use the inner operator directly when nothing is transformed)");
    }
    value_transformation_ = value_transformation->Bind(
        target_space.GetLayout()->GetDOFHolderCoordinates());
  }

  void Apply(const Field<T>& source, Field<T>& target) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto& sd = source.GetData();
    const auto& td = target.GetData();
    if (sd.GetValueType() != td.GetValueType()) {
      throw pcms_error(
        "TransformedTransferOperator: source and target value types differ");
    }
    const bool needs_transformation =
      !SameValueBasis(sd.GetValueBasis(), td.GetValueBasis());
    if (!needs_transformation) {
      // Bases agree: the inner committed path, zero wrapper overhead.
      inner_->Apply(source, target);
      return;
    }
    if constexpr (!std::is_same_v<T, Real>) {
      throw pcms_error(
        "TransformedTransferOperator: basis changes require T == Real");
    } else {
      if (!SameValueBasis(sd.GetValueBasis(),
                          value_transformation_->GetSourceBasis()) ||
          !SameValueBasis(td.GetValueBasis(),
                          value_transformation_->GetTargetBasis())) {
        throw pcms_error(
          "TransformedTransferOperator: the value transformation does not "
          "bridge the source's stored basis to the target's declared basis");
      }
      // The bridge check pins the transformation's target basis to the
      // target's declared basis, so the keyed apply (transfer + rotate + FILL
      // restore) produces exactly what the target declares; commit it.
      const auto n = static_cast<size_t>(value_transformation_->NumPoints());
      const auto width =
        static_cast<size_t>(target.GetLayout().GetNumComponents());
      if (rotated_scratch_.extent(0) != n ||
          rotated_scratch_.extent(1) != width) {
        rotated_scratch_ = Kokkos::View<Real**, DeviceMemorySpace>(
          "transformed_transfer_rotated", n, width);
      }
      Apply(this->MakeTransferKey(), source, MakeRank2View(rotated_scratch_));
      target.SetDOFHolderDataUnchecked(MakeConstRank2View(rotated_scratch_));
    }
  }

  // Uncommitted apply (keyed): with no target declaration to consult, the
  // output basis is the value transformation's target basis when the source's
  // stored basis differs from it, else the source's stored basis. See
  // TransferOperator for the contract.
  void Apply(TransferKey key, const Field<T>& source,
             Rank2View<T, DeviceMemorySpace> out) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto& sd = source.GetData();
    const bool needs_transformation = !SameValueBasis(
      sd.GetValueBasis(), value_transformation_->GetTargetBasis());
    if (!needs_transformation) {
      inner_->Apply(key, source, out);
      return;
    }
    if constexpr (!std::is_same_v<T, Real>) {
      throw pcms_error(
        "TransformedTransferOperator: basis changes require T == Real");
    } else {
      if (!SameValueBasis(sd.GetValueBasis(),
                          value_transformation_->GetSourceBasis())) {
        throw pcms_error(
          "TransformedTransferOperator: the value transformation does not "
          "start at the source's stored basis");
      }
      if (transfer_scratch_.extent(0) != out.extent(0) ||
          transfer_scratch_.extent(1) != out.extent(1)) {
        transfer_scratch_ = Kokkos::View<Real**, DeviceMemorySpace>(
          "transformed_transfer_scratch", out.extent(0), out.extent(1));
      }
      // The inner operator validates the buffer extents against the target
      // layout.
      inner_->Apply(key, source, MakeRank2View(transfer_scratch_));
      value_transformation_->Apply(
        ValueView<const Real, DeviceMemorySpace>(
          sd.GetValueBasis(), MakeConstRank2View(transfer_scratch_)),
        ValueView<Real, DeviceMemorySpace>(
          value_transformation_->GetTargetBasis(), out));
      if (policy_.mode == OutOfBoundsMode::FILL) {
        detail::RestoreFillRows(transfer_scratch_, out, policy_.fill_value);
      }
    }
  }

private:
  OutOfBoundsPolicy policy_;
  // Per-call scratch (repo operator idiom): concurrent Apply calls on one
  // operator are not thread-safe.
  mutable Kokkos::View<Real**, DeviceMemorySpace> transfer_scratch_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> rotated_scratch_;
  std::unique_ptr<TransferOperator<T>> inner_;
  std::unique_ptr<BoundBasisTransformation> value_transformation_;
};

} // namespace pcms

#endif // PCMS_TRANSFER_TRANSFORMED_TRANSFER_OPERATOR_HPP
