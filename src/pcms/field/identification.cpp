#include "pcms/field/identification.hpp"
#include "pcms/field/function_space.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/memory_spaces.hpp"
#include <algorithm>

namespace pcms
{

namespace
{

/// Copies the columns `axes` of `points` (a free function: CUDA does not allow
/// device lambdas in members).
Kokkos::View<Real**, DeviceMemorySpace> Columns(
  Rank2View<const Real, DeviceMemorySpace> points,
  Kokkos::View<const int*, DeviceMemorySpace> axes)
{
  Kokkos::View<Real**, DeviceMemorySpace> out("twist_shift_inputs",
                                              points.extent(0), axes.extent(0));
  const auto count = static_cast<int>(axes.extent(0));
  pcms::detail::ForEach(
    "TwistShiftInputs", static_cast<LO>(points.extent(0)),
    KOKKOS_LAMBDA(const LO p) {
      for (int k = 0; k < count; ++k) {
        out(p, k) = points(p, axes(k));
      }
    });
  return out;
}

void AddShift(Rank2View<Real, DeviceMemorySpace> points,
              Rank1View<const Real, DeviceMemorySpace> wraps, int axis,
              Kokkos::View<const Real**, DeviceMemorySpace> shift)
{
  pcms::detail::ForEach(
    "TwistShift", static_cast<LO>(points.extent(0)),
    KOKKOS_LAMBDA(const LO p) { points(p, axis) += wraps(p) * shift(p, 0); });
}

/// @return the number of the points `missing` that have wraps
LO CountWrapped(Kokkos::View<const LO*, DeviceMemorySpace> missing,
                Kokkos::View<const Real*, DeviceMemorySpace> wraps)
{
  LO count = 0;
  Kokkos::parallel_reduce(
    "TwistShiftWrappedMissing",
    Kokkos::RangePolicy<DeviceMemorySpace::execution_space>(
      0, static_cast<LO>(missing.extent(0))),
    KOKKOS_LAMBDA(const LO i, LO& sum) { sum += wraps(missing(i)) != 0; },
    count);
  return count;
}

class TwistShiftMap final : public IdentificationMap
{
public:
  TwistShiftMap(int shift_axis, std::shared_ptr<const Function<Real>> shift,
                const std::vector<int>& input_axes)
    : shift_axis_(shift_axis),
      shift_(std::move(shift)),
      input_axes_("twist_shift_axes", input_axes.size()),
      max_axis_(shift_axis)
  {
    const auto host = Kokkos::create_mirror_view(input_axes_);
    for (std::size_t k = 0; k < input_axes.size(); ++k) {
      host(k) = input_axes[k];
      max_axis_ = std::max(max_axis_, input_axes[k]);
    }
    Kokkos::deep_copy(input_axes_, host);
  }

  void Apply(Rank2View<Real, DeviceMemorySpace> points,
             Rank1View<const Real, DeviceMemorySpace> wraps) const override
  {
    if (max_axis_ >= static_cast<int>(points.extent(1))) {
      throw pcms_error("twist-shift: an axis of the map is not an axis of the "
                       "coordinate system");
    }
    const auto inputs = Columns(points, input_axes_);
    const auto geometry = shift_->GetSpace().GetGeometry();
    // S is needed only at the points that wrapped
    const auto bound = geometry->Bind(
      CoordinateView<DeviceMemorySpace>(geometry->GetCoordinateSystem(),
                                        MakeConstRank2View(inputs)),
      OutOfBoundsPolicy{OutOfBoundsMode::FILL});
    if (CountWrapped(bound->MissingPoints(), wraps) > 0) {
      throw pcms_error("twist-shift: a point crosses the identification "
                       "outside the domain of the shift");
    }
    const auto evaluator = shift_->GetSpace().CreatePointEvaluator<Real>(bound);
    Kokkos::View<Real**, DeviceMemorySpace> shift("twist_shift",
                                                  points.extent(0), 1);
    evaluator->Evaluate(
      *shift_, ValueView<Real, DeviceMemorySpace>(
                 evaluator->OutputBasis(shift_->GetData().GetValueBasis()),
                 MakeRank2View(shift)));
    AddShift(points, wraps, shift_axis_, shift);
  }

private:
  int shift_axis_;
  std::shared_ptr<const Function<Real>> shift_;
  Kokkos::View<int*, DeviceMemorySpace> input_axes_;
  int max_axis_;
};

} // namespace

std::shared_ptr<const IdentificationMap> CreateTwistShiftMap(
  int shift_axis, std::shared_ptr<const Function<Real>> shift,
  std::vector<int> input_axes)
{
  const bool negative =
    shift_axis < 0 ||
    std::ranges::any_of(input_axes, [](int a) { return a < 0; });
  if (negative || shift == nullptr ||
      shift->GetSpace().GetGeometry()->GetCoordinateSystem()->Dimension() !=
        static_cast<int>(input_axes.size())) {
    throw pcms_error("CreateTwistShiftMap: needs axes that are not negative "
                     "and a shift Function on a geometry of one dimension per "
                     "input axis");
  }
  return std::make_shared<const TwistShiftMap>(shift_axis, std::move(shift),
                                               input_axes);
}

} // namespace pcms
