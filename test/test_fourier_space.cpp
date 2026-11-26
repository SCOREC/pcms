#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Kokkos_Core.hpp>
#include <pcms/field/field_ops.hpp>
#include <pcms/field/function_space.hpp>
#include <pcms/field/detail/simple_field_data.hpp>
#include <pcms/geometry/geometry.hpp>
#include <pcms/transfer/projection.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/transfer/forms.hpp>
#include <pcms/transfer/linear_solver.hpp>
#include <cmath>
#include <numbers>
#include <vector>

// A real Fourier space on the periodic interval [0, 2 pi): the modes 1,
// cos(kx), sin(kx) for k = 1..K, an equispaced trapezoid rule as its
// quadrature, and an evaluator with Evaluate and EvaluateTranspose only (no
// matrix). It is here to show that the projection machinery composes with a
// space that is not an Omega_h one; it is a test fixture, not a space PCMS
// provides.

namespace
{

using pcms::DeviceMemorySpace;
using pcms::LO;
using pcms::Real;

constexpr Real TwoPi = 2 * std::numbers::pi_v<Real>;

// mode m: 0 -> 1, 2k - 1 -> cos(kx), 2k -> sin(kx)
KOKKOS_INLINE_FUNCTION Real Mode(int m, Real x)
{
  if (m == 0) {
    return 1.0;
  }
  const int k = (m + 1) / 2;
  return m % 2 == 1 ? Kokkos::cos(k * x) : Kokkos::sin(k * x);
}

class FourierLayout final : public pcms::FieldLayout
{
public:
  FourierLayout(std::shared_ptr<const pcms::Geometry> geometry, int harmonics)
    : size_(2 * harmonics + 1),
      owned_("owned", size_),
      gids_("gids", size_),
      zeros_("zeros", size_),
      coords_("coords", size_, 1)
  {
    for (int i = 0; i < size_; ++i) {
      owned_(i) = true;
      gids_(i) = i;
    }
    SetGeometry(std::move(geometry));
    SetCoordinateSystem(pcms::csys::CreateCartesian(1));
  }
  int GetNumComponents() const override { return 1; }
  LO GetNumDofHolder() const override { return size_; }
  pcms::GO GetNumGlobalDofHolder() const override { return size_; }
  bool IsDistributed() const override { return false; }
  pcms::Rank1View<const bool, pcms::HostMemorySpace> GetOwnedHost()
    const override
  {
    return pcms::MakeRank1View(owned_);
  }
  pcms::GlobalIDView<pcms::HostMemorySpace> GetGidsHost() const override
  {
    return pcms::GlobalIDView<pcms::HostMemorySpace>(gids_.data(), size_);
  }
  pcms::EntOffsetsArray GetEntOffsets() const override
  {
    const auto n = static_cast<std::size_t>(size_);
    return {0, n, n, n, n};
  }
  pcms::CoordinateView<DeviceMemorySpace> GetDOFHolderCoordinates()
    const override
  {
    return pcms::CoordinateView<DeviceMemorySpace>(
      GetCoordinateSystem(), pcms::MakeConstRank2View(coords_));
  }
  int GetDimension() const override { return 1; }
  pcms::Rank1View<const LO, pcms::HostMemorySpace>
  GetDOFHolderClassificationDimensionsHost() const override
  {
    return pcms::MakeRank1View(zeros_);
  }
  pcms::Rank1View<const LO, pcms::HostMemorySpace>
  GetDOFHolderClassificationIdsHost() const override
  {
    return pcms::MakeRank1View(zeros_);
  }

private:
  int size_;
  Kokkos::View<bool*, pcms::HostMemorySpace> owned_;
  Kokkos::View<pcms::GO*, pcms::HostMemorySpace> gids_;
  Kokkos::View<LO*, pcms::HostMemorySpace> zeros_;
  Kokkos::View<Real**, DeviceMemorySpace> coords_;
};

class PeriodicBound final : public pcms::BoundPoints
{
public:
  PeriodicBound(std::shared_ptr<const pcms::Geometry> geometry,
                const pcms::BindRequest& request)
    : BoundPoints(std::move(geometry), request)
  {
  }

protected:
  void RebindImpl() override {}
};

// The periodic interval: every point is found. Its quadrature of `order` is
// the trapezoid rule with order + 1 points, which integrates the trigonometric
// polynomials up to degree `order` exactly.
class PeriodicGeometry final : public pcms::Geometry
{
public:
  PeriodicGeometry() : system_(pcms::csys::CreateCartesian(1)) {}
  const std::shared_ptr<const pcms::CoordinateSystem>& GetCoordinateSystem()
    const noexcept override
  {
    return system_;
  }

protected:
  pcms::Quadrature CreateQuadrature(
    const pcms::detail::QuadratureRule& rule) const override
  {
    const int n = std::get<pcms::detail::GaussRule>(rule).order + 1;
    Kokkos::View<Real**, DeviceMemorySpace> x("trapezoid_points", n, 1);
    Kokkos::View<Real*, DeviceMemorySpace> w("trapezoid_weights", n);
    Kokkos::parallel_for(
      "trapezoid", Kokkos::RangePolicy<pcms::DefaultExecutionSpace>(0, n),
      KOKKOS_LAMBDA(const int i) {
        x(i, 0) = TwoPi * i / n;
        w(i) = TwoPi / n;
      });
    return pcms::Quadrature{Bind(pcms::CoordinateView<DeviceMemorySpace>(
                              system_, pcms::MakeConstRank2View(x))),
                            w,
                            {}};
  }

  std::shared_ptr<pcms::BoundPoints> BindImpl(
    const pcms::BindRequest& request) const override
  {
    return std::make_shared<PeriodicBound>(shared_from_this(), request);
  }

private:
  std::shared_ptr<const pcms::CoordinateSystem> system_;
};

class FourierEvaluator final : public pcms::PointEvaluator<Real>
{
public:
  FourierEvaluator(std::shared_ptr<const pcms::BoundPoints> bound,
                   std::shared_ptr<const pcms::FieldLayout> layout)
    : PointEvaluator(std::move(bound), pcms::DifferentialOperator::Value,
                     layout),
      layout_(std::move(layout))
  {
  }
  bool SupportsTranspose() const noexcept override { return true; }

  void EvaluateImpl(
    const pcms::Field<Real>& field,
    pcms::ValueView<Real, DeviceMemorySpace> values) const override
  {
    const auto out = values.GetValues();
    const auto dofs = field.GetDOFHolderData().GetValues();
    const auto points = Bound().LocatedPoints().GetValues();
    const int modes = layout_->GetNumDofHolder();
    Kokkos::parallel_for(
      "fourier_evaluate",
      Kokkos::RangePolicy<pcms::DefaultExecutionSpace>(
        0, static_cast<LO>(points.extent(0))),
      KOKKOS_LAMBDA(const LO p) {
        Real sum = 0;
        for (int m = 0; m < modes; ++m) {
          sum += Mode(m, points(p, 0)) * dofs(m, 0);
        }
        out(p, 0) = sum;
      });
  }

  void EvaluateTransposeImpl(
    pcms::Rank2View<const Real, DeviceMemorySpace> values,
    pcms::Field<Real>& functional) const override
  {
    const auto current = functional.GetDOFHolderData().GetValues();
    const auto points = Bound().LocatedPoints().GetValues();
    const int modes = layout_->GetNumDofHolder();
    const auto num_points = static_cast<LO>(points.extent(0));
    Kokkos::View<Real**, DeviceMemorySpace> sum("fourier_transpose", modes, 1);
    Kokkos::parallel_for(
      "fourier_transpose",
      Kokkos::RangePolicy<pcms::DefaultExecutionSpace>(0, modes),
      KOKKOS_LAMBDA(const int m) {
        Real total = current(m, 0);
        for (LO p = 0; p < num_points; ++p) {
          total += Mode(m, points(p, 0)) * values(p, 0);
        }
        sum(m, 0) = total;
      });
    functional.SetDOFHolderDataUnchecked(pcms::MakeConstRank2View(sum));
  }

private:
  std::shared_ptr<const pcms::FieldLayout> layout_;
};

class FourierEvaluatorFactory final : public pcms::PointEvaluatorFactory
{
public:
  FourierEvaluatorFactory(std::shared_ptr<const pcms::Geometry> geometry,
                          std::shared_ptr<const pcms::FieldLayout> layout)
    : geometry_(std::move(geometry)), layout_(std::move(layout))
  {
  }
  std::shared_ptr<const pcms::Geometry> GetGeometry() const noexcept override
  {
    return geometry_;
  }
  bool Supports(pcms::DifferentialOperator op) const noexcept override
  {
    return op == pcms::DifferentialOperator::Value;
  }

protected:
  pcms::PointEvaluatorVariant CreatePointEvaluatorImpl(
    pcms::Type, std::shared_ptr<const pcms::BoundPoints> bound,
    pcms::DifferentialOperator) const override
  {
    return std::unique_ptr<pcms::PointEvaluator<Real>>(
      std::make_unique<FourierEvaluator>(std::move(bound), layout_));
  }

private:
  std::shared_ptr<const pcms::Geometry> geometry_;
  std::shared_ptr<const pcms::FieldLayout> layout_;
};

// The Fourier space: its natural quadrature integrates the product of two of
// its functions, the trigonometric polynomials of degree 2N.
class FourierSpace final : public pcms::FunctionSpace
{
public:
  FourierSpace(std::shared_ptr<const pcms::FieldLayout> layout,
               std::shared_ptr<const pcms::PointEvaluatorFactory> evaluator,
               std::shared_ptr<const pcms::Geometry> geometry, int harmonics)
    : pcms::FunctionSpace(std::move(layout), std::move(evaluator), geometry),
      geometry_(std::move(geometry)),
      harmonics_(harmonics)
  {
  }

  [[nodiscard]] pcms::Quadrature CreateQuadrature() const override
  {
    return pcms::CreateGaussQuadrature(geometry_, 2 * harmonics_);
  }

private:
  std::shared_ptr<const pcms::Geometry> geometry_;
  int harmonics_;
};

std::shared_ptr<const pcms::FunctionSpace> MakeFourierSpace(
  const std::shared_ptr<const pcms::Geometry>& geometry, int harmonics)
{
  const auto layout = std::make_shared<FourierLayout>(geometry, harmonics);
  return std::make_shared<FourierSpace>(
    layout, std::make_shared<FourierEvaluatorFactory>(geometry, layout),
    geometry, harmonics);
}

void SetCoefficients(pcms::Field<Real>& field, const std::vector<Real>& c)
{
  Kokkos::View<Real**, pcms::HostMemorySpace> values("c", c.size(), 1);
  for (std::size_t i = 0; i < c.size(); ++i) {
    values(i, 0) = c[i];
  }
  field.SetDOFHolderDataHost(pcms::ValueView<const Real, pcms::HostMemorySpace>(
    field.GetData().GetValueBasis(), pcms::MakeConstRank2View(values)));
}

} // namespace

TEST_CASE("Fourier space: the mass matrix is diagonal and CG takes two "
          "iterations",
          "[fourier][projection]")
{
  const auto geometry = std::make_shared<PeriodicGeometry>();
  const int harmonics = 4;
  const auto space = MakeFourierSpace(geometry, harmonics);
  const pcms::transfer::detail::BilinearForm mass(space);
  REQUIRE_FALSE(mass.SupportsMatrix());
  const int n = 2 * harmonics + 1;
  for (int j = 0; j < n; ++j) {
    auto e = space->CreateFunction<Real>();
    std::vector<Real> unit(n, 0.0);
    unit[j] = 1.0;
    SetCoefficients(e, unit);
    auto column = space->CreateFunction<Real>();
    mass.Apply(e, column);
    const auto values = column.GetDOFHolderDataHost().GetValues();
    for (int i = 0; i < n; ++i) {
      const Real expected = i != j ? 0.0 : (j == 0 ? TwoPi : TwoPi / 2);
      CHECK(values(i, 0) == Catch::Approx(expected).margin(1e-13));
    }
  }
  // two distinct eigenvalues (2 pi and pi): CG converges in two iterations
  pcms::transfer::detail::LinearSolver solver(mass, pcms::LinearSolveOptions{});
  auto b = space->CreateFunction<Real>();
  SetCoefficients(b, {1, 2, -3, 4, 5, -6, 7, 8, -9});
  auto x = space->CreateFunction<Real>();
  const int iterations = solver.Solve(b, x);
  CHECK(iterations <= 2);
  const auto values = x.GetDOFHolderDataHost().GetValues();
  CHECK(values(0, 0) == Catch::Approx(1 / TwoPi).epsilon(1e-12));
  CHECK(values(1, 0) == Catch::Approx(2 / (TwoPi / 2)).epsilon(1e-12));
  std::printf("[fourier] mass solve: %d CG iterations\n", iterations);
}

TEST_CASE("Fourier space: projection onto a space with more or fewer modes "
          "gives the analytic coefficients",
          "[fourier][projection]")
{
  const auto geometry = std::make_shared<PeriodicGeometry>();
  const int source_harmonics = 2;
  const auto source_space = MakeFourierSpace(geometry, source_harmonics);
  const std::vector<Real> source_coefficients = {0.5, 2.0, -1.0, 0.25, 3.0};
  auto source = source_space->CreateFunction<Real>();
  SetCoefficients(source, source_coefficients);
  for (const int target_harmonics : {1, 2, 4}) {
    CAPTURE(target_harmonics);
    const auto target_space = MakeFourierSpace(geometry, target_harmonics);
    pcms::ProjectionOptions options;
    // the source times the target: trigonometric degree up to the two
    // harmonics
    options.quadrature = pcms::CreateGaussQuadrature(
      geometry, source_harmonics + target_harmonics);
    auto target = target_space->CreateFunction<Real>();
    pcms::CreateProjection(source_space, target_space, options)
      ->Apply(source, target);
    const auto values = target.GetDOFHolderDataHost().GetValues();
    REQUIRE(values.extent(0) == 2 * target_harmonics + 1);
    for (std::size_t m = 0; m < values.extent(0); ++m) {
      // the modes are orthogonal: the projection truncates or pads with zeros
      const Real expected =
        m < source_coefficients.size() ? source_coefficients[m] : 0.0;
      CHECK(values(m, 0) == Catch::Approx(expected).margin(1e-12));
    }
  }
  // the assembled method needs the evaluation matrix, which this evaluator
  // does not provide
  const auto target_space = MakeFourierSpace(geometry, 2);
  pcms::ProjectionOptions assembled;
  assembled.solve.method = pcms::LinearSolveOptions::Method::Assembled;
  CHECK_THROWS(pcms::CreateProjection(source_space, target_space, assembled));
}
