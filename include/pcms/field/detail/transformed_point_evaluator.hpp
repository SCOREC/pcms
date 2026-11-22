#ifndef PCMS_FIELD_DETAIL_TRANSFORMED_POINT_EVALUATOR_HPP
#define PCMS_FIELD_DETAIL_TRANSFORMED_POINT_EVALUATOR_HPP

#include "pcms/field/field.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/point_evaluator.hpp"
#include "pcms/field/value_transform.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/geometry/bound_points.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/memory_spaces.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include "pcms/utility/detail/profile.hpp"
#include <functional>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

namespace pcms::detail
{

/// The value transforms, in order, that re-express values in a basis.
using ValueTransforms =
  std::function<std::vector<const ValueTransform*>(const ValueBasis&)>;

/// Point evaluator that re-expresses the values of an inner evaluator through
/// the value transforms of their basis (those of a bind with a coordinate
/// map). Values without transforms are written by the inner evaluator
/// directly; others are transformed in place, through a buffer only where a
/// transform changes the number of components.
/// PointEvaluatorFactory::CreatePointEvaluator creates it.
template <typename T>
class TransformedPointEvaluator final : public PointEvaluator<T>
{
public:
  /// @param bound BoundPoints of the query points
  /// @param inner evaluator created from `bound` that writes the values
  ///        before the transforms; not null
  /// @param transforms the transforms of the values of `inner`
  TransformedPointEvaluator(std::shared_ptr<const BoundPoints> bound,
                            std::unique_ptr<PointEvaluator<T>> inner,
                            ValueTransforms transforms)
    : PointEvaluator<T>(std::move(bound), inner->Operator(),
                        inner->AcceptedLayout()),
      inner_(std::move(inner)),
      transforms_(std::move(transforms))
  {
  }

  // Public because CUDA forbids extended lambdas in non-public member
  // functions. User should call Evaluate
  // TODO: add passkey or other mechanism to avoid user calling this
  void EvaluateImpl(const Field<T>& field,
                    ValueView<T, DeviceMemorySpace> values) const override
  {
    PCMS_FUNCTION_TIMER;
    const auto [bases, transforms] = Chain(field.GetData().GetValueBasis());
    if (transforms.empty()) {
      inner_->Evaluate(field, values);
      return;
    }
    if constexpr (std::is_same_v<T, Real>) {
      const auto out = values.GetValues();
      const int dout = bases.back().system->Dimension();
      const std::size_t blocks = out.extent(1) / dout;
      // the output itself for as many components as it has
      const auto storage = [&](const ValueBasis& basis) {
        const int dim = basis.system->Dimension();
        return dim == dout ? out : Buffer(out.extent(0), blocks * dim);
      };
      inner_->Evaluate(field, ValueView<Real, DeviceMemorySpace>(
                                bases.front(), storage(bases.front())));
      for (std::size_t k = 1; k < bases.size(); ++k) {
        CheckOneBufferShape(bases[k - 1], bases[k], dout);
        transforms[k - 1]->Apply(
          ValueView<const Real, DeviceMemorySpace>(bases[k - 1],
                                                   storage(bases[k - 1])),
          ValueView<Real, DeviceMemorySpace>(bases[k], storage(bases[k])));
      }
    } else {
      throw pcms_error("TransformedPointEvaluator::Evaluate: basis "
                       "transformations require T == Real");
    }
  }

  /// @return whether the inner evaluator supports the transpose
  [[nodiscard]] bool SupportsTranspose() const noexcept override
  {
    return inner_->SupportsTranspose();
  }

  // Public for the same reason as EvaluateImpl.
  void EvaluateTransposeImpl(Rank2View<const T, DeviceMemorySpace> values,
                             Field<T>& functional) const override
  {
    PCMS_FUNCTION_TIMER;
    // the inner evaluator may count the points this evaluator fills as found
    Kokkos::View<T**, DeviceMemorySpace> masked;
    bool writable = false;
    if (const auto filled = this->FilledPoints(); filled.extent(0) > 0) {
      masked = Kokkos::View<T**, DeviceMemorySpace>(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "transpose_masked"),
        values.extent(0), values.extent(1));
      CopyDeviceRank2ViewToDeviceView(masked, values);
      FillRows(MakeRank2View(masked), filled, T{0});
      values = MakeConstRank2View(masked);
      writable = true;
    }
    const auto [bases, transforms] =
      Chain(functional.GetData().GetValueBasis());
    if (transforms.empty()) {
      inner_->EvaluateTranspose(
        ValueView<const T, DeviceMemorySpace>(bases.front(), values),
        functional);
      return;
    }
    if constexpr (std::is_same_v<T, Real>) {
      // E^T T_1^T ... T_n^T
      const int dout = bases.back().system->Dimension();
      const std::size_t blocks = values.extent(1) / dout;
      Rank2View<const Real, DeviceMemorySpace> in = values;
      for (std::size_t k = bases.size() - 1; k > 0; --k) {
        CheckOneBufferShape(bases[k - 1], bases[k], dout);
        const int dim = bases[k - 1].system->Dimension();
        if (dim == dout && !writable) {
          masked = Kokkos::View<T**, DeviceMemorySpace>(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, "transpose_values"),
            values.extent(0), values.extent(1));
          writable = true;
        }
        const auto out = dim == dout ? MakeRank2View(masked)
                                     : Buffer(values.extent(0), blocks * dim);
        transforms[k - 1]->ApplyTranspose(
          ValueView<const Real, DeviceMemorySpace>(bases[k], in),
          ValueView<Real, DeviceMemorySpace>(bases[k - 1], out));
        in = out;
      }
      inner_->EvaluateTranspose(
        ValueView<const Real, DeviceMemorySpace>(bases.front(), in),
        functional);
    } else {
      throw pcms_error("TransformedPointEvaluator::EvaluateTranspose: basis "
                       "transformations require T == Real");
    }
  }

  [[nodiscard]] ValueBasis OutputBasis(const ValueBasis& stored) const override
  {
    return Chain(stored).first.back();
  }

private:
  // the buffer of the values whose dimension is not the output's,
  // [points][columns], reallocated on a change of shape
  Rank2View<Real, DeviceMemorySpace> Buffer(std::size_t points,
                                            std::size_t columns) const
  {
    if (buffer_.extent(0) != points || buffer_.extent(1) != columns) {
      buffer_ = Kokkos::View<Real**, DeviceMemorySpace>(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "transformed_values"),
        points, columns);
    }
    return MakeRank2View(buffer_);
  }

  // one buffer holds the values of one dimension other than the output's
  static void CheckOneBufferShape(const ValueBasis& a, const ValueBasis& b,
                                  int dout)
  {
    const int da = a.system->Dimension();
    const int db = b.system->Dimension();
    PCMS_ALWAYS_ASSERT(da == dout || db == dout || da == db);
  }

  // the inner evaluator's basis for `stored`, then the output basis of each
  // of its transforms; and the transforms
  [[nodiscard]] std::pair<std::vector<ValueBasis>,
                          std::vector<const ValueTransform*>>
  Chain(const ValueBasis& stored) const
  {
    std::vector<ValueBasis> bases{inner_->OutputBasis(stored)};
    auto transforms = transforms_(bases.front());
    for (const auto* transform : transforms) {
      bases.push_back(transform->OutputBasis(bases.back()));
    }
    return {std::move(bases), std::move(transforms)};
  }

  std::unique_ptr<PointEvaluator<T>> inner_;
  ValueTransforms transforms_;
  mutable Kokkos::View<Real**, DeviceMemorySpace> buffer_;
};

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_TRANSFORMED_POINT_EVALUATOR_HPP
