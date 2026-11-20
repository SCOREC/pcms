#ifndef PCMS_FIELD_FUNCTION_SPACE_SPARSE_ROWS_HPP
#define PCMS_FIELD_FUNCTION_SPACE_SPARSE_ROWS_HPP

// Private: spaces whose evaluation at a set of bound points is a sparse
// matrix of DOF weights built when the points are bound (the 1D periodic
// spaces and the product spaces). One evaluator serves them all: Evaluate is
// the product of the matrix with the DOFs, EvaluateTranspose its transpose,
// GetMatrix the matrix itself.

#include "pcms/field/detail/field_ops.hpp"
#include "pcms/field/detail/simple_field_data.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/utility/arrays.hpp"
#include <Kokkos_Core.hpp>
#include <functional>
#include <optional>
#include <memory>
#include <set>

namespace pcms::detail
{

/// The rows of the evaluation matrix of a set of bound points. Row
/// p * num_out + c is the weights of the DOFs that component c of the output at
/// point p combines (num_out is 1 for Value and the dimension for Gradient).
struct SparseRows
{
  Kokkos::View<LO*, DeviceMemorySpace> offsets; // [rows + 1]
  Kokkos::View<GO*, DeviceMemorySpace> columns; // global ids of the DOF holders
  Kokkos::View<LO*, DeviceMemorySpace> holders; // local holders of the columns
  Kokkos::View<Real*, DeviceMemorySpace> values;
  int num_out = 1;
};

/// out(p, c) = sum of row (p, c) times the DOFs. (Free functions: CUDA does not
/// allow device lambdas in the virtual functions of the evaluators.)
inline void ApplySparseRows(const SparseRows& rows,
                            Rank2View<const Real, DeviceMemorySpace> dofs,
                            Rank2View<Real, DeviceMemorySpace> out)
{
  const auto offsets = rows.offsets;
  const auto holders = rows.holders;
  const auto values = rows.values;
  const int num_out = rows.num_out;
  ForEach(
    "ApplySparseRows", static_cast<LO>(out.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      for (int c = 0; c < num_out; ++c) {
        Real sum = 0;
        for (LO j = offsets(p * num_out + c); j < offsets(p * num_out + c + 1);
             ++j) {
          sum += values(j) * dofs(holders(j), 0);
        }
        out(p, c) = sum;
      }
    });
}

/// total(holder) += sum over the entries of the holder times in(p, 0).
inline void ApplySparseRowsTranspose(
  const SparseRows& rows, Rank2View<const Real, DeviceMemorySpace> in,
  Kokkos::View<Real**, DeviceMemorySpace> total)
{
  const auto offsets = rows.offsets;
  const auto holders = rows.holders;
  const auto values = rows.values;
  Kokkos::parallel_for(
    "ApplySparseRowsTranspose",
    Kokkos::RangePolicy<DefaultExecutionSpace>(0,
                                               static_cast<LO>(in.extent(0))),
    KOKKOS_LAMBDA(const LO p) {
      for (LO j = offsets(p); j < offsets(p + 1); ++j) {
        Kokkos::atomic_add(&total(holders(j), 0), values(j) * in(p, 0));
      }
    });
}

/// Evaluates through SparseRows; Value and Gradient of scalar fields.
class SparseRowEvaluator final : public PointEvaluator<Real>
{
public:
  SparseRowEvaluator(std::shared_ptr<const BoundPoints> bound,
                     std::shared_ptr<const FieldLayout> layout,
                     DifferentialOperator op, SparseRows rows)
    : PointEvaluator<Real>(std::move(bound), op, layout),
      layout_(std::move(layout)),
      rows_(std::move(rows))
  {
  }

  [[nodiscard]] bool SupportsTranspose() const noexcept override
  {
    return true;
  }
  [[nodiscard]] bool SupportsMatrix() const noexcept override { return true; }

  void EvaluateImpl(const Field<Real>& field,
                    ValueView<Real, DeviceMemorySpace> values) const override
  {
    ApplySparseRows(rows_, field.GetDOFHolderData().GetValues(),
                    values.GetValues());
  }

  void EvaluateTransposeImpl(Rank2View<const Real, DeviceMemorySpace> values,
                             Field<Real>& functional) const override
  {
    Kokkos::View<Real**, DeviceMemorySpace> total(
      "sparse_rows_transpose", layout_->GetNumDofHolder(), 1);
    ApplySparseRowsTranspose(rows_, values, total);
    AddToField(functional, total);
  }

  [[nodiscard]] PointMatrix<Real> GetMatrixImpl() const override
  {
    return {rows_.offsets, rows_.columns, rows_.values};
  }

private:
  std::shared_ptr<const FieldLayout> layout_;
  SparseRows rows_;
};

/// The evaluator source of a space whose evaluators are SparseRowEvaluators
/// of rows a function makes from the bound points.
class SparseRowEvaluatorFactory final : public PointEvaluatorFactory
{
public:
  using MakeRows = std::function<SparseRows(
    const std::shared_ptr<const BoundPoints>&, DifferentialOperator)>;

  SparseRowEvaluatorFactory(std::shared_ptr<const Geometry> geometry,
                            std::shared_ptr<const FieldLayout> layout,
                            std::set<DifferentialOperator> supported,
                            MakeRows make_rows)
    : geometry_(std::move(geometry)),
      layout_(std::move(layout)),
      supported_(std::move(supported)),
      make_rows_(std::move(make_rows))
  {
  }

  [[nodiscard]] std::shared_ptr<const Geometry> GetGeometry()
    const noexcept override
  {
    return geometry_;
  }
  [[nodiscard]] bool Supports(DifferentialOperator op) const noexcept override
  {
    return supported_.count(op) != 0;
  }

protected:
  PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, std::shared_ptr<const BoundPoints> bound,
    DifferentialOperator op) const override
  {
    if (value_type != Type::Real) {
      throw pcms_error("CreatePointEvaluator: this space has Real fields only");
    }
    auto rows = make_rows_(bound, op);
    return std::unique_ptr<PointEvaluator<Real>>(
      std::make_unique<SparseRowEvaluator>(std::move(bound), layout_, op,
                                           std::move(rows)));
  }

private:
  std::shared_ptr<const Geometry> geometry_;
  std::shared_ptr<const FieldLayout> layout_;
  std::set<DifferentialOperator> supported_;
  MakeRows make_rows_;
};

/// A space with a natural quadrature a function makes.
class QuadratureSpace final : public FunctionSpace
{
public:
  QuadratureSpace(std::shared_ptr<const FieldLayout> layout,
                  std::shared_ptr<const PointEvaluatorFactory> evaluator,
                  std::shared_ptr<const Geometry> geometry,
                  std::function<Quadrature()> make_quadrature)
    : FunctionSpace(layout, std::move(evaluator), std::move(geometry)),
      layout_(std::move(layout)),
      make_quadrature_(std::move(make_quadrature))
  {
  }

  [[nodiscard]] Quadrature CreateQuadrature() const override
  {
    return make_quadrature_();
  }

protected:
  FieldDataVariant CreateFieldData(Type type, ValueBasis basis) const override
  {
    if (type != Type::Real) {
      throw pcms_error("this space has Real fields only");
    }
    return std::make_unique<SimpleFieldData<Real>>(layout_, std::move(basis));
  }

private:
  std::shared_ptr<const FieldLayout> layout_;
  std::function<Quadrature()> make_quadrature_;
};

} // namespace pcms::detail

#endif // PCMS_FIELD_FUNCTION_SPACE_SPARSE_ROWS_HPP
