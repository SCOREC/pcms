#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <pcms/backends/uniform_grid/geometry.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/configuration.hpp>
#include <pcms/field/detail/simple_field_data.hpp>
#include <pcms/field/field_metadata.hpp>
#include <pcms/field/field_ops.hpp>
#include <pcms/transfer/copy.hpp>
#include <pcms/transfer/interpolator.hpp>
#include <pcms/transfer/l2_inner_product.hpp>
#include <pcms/transfer/projection.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <pcms/utility/coordinate_systems/cylindrical.hpp>
#include <pcms/utility/coordinate_systems/parameter.hpp>

using Catch::Matchers::ContainsSubstring;
using pcms::Duality;
using pcms::FieldMetadata;
using pcms::Real;
using pcms::Units;

namespace
{

// P1 on `cells` cells of [0, 1]: cells + 1 DOFs
std::shared_ptr<const pcms::FunctionSpace> Space(int cells)
{
  pcms::UniformGrid<1> grid;
  grid.bot_left = {0.0};
  grid.edge_length = {1.0};
  grid.divisions = {cells};
  return pcms::CreateLagrangeSpace(
    pcms::uniform_grid::CreateGeometry(grid, pcms::csys::CreateCartesian(1)), 1,
    1);
}

void SetDofs(pcms::Field<Real>& field, const std::vector<Real>& c)
{
  Kokkos::View<Real**, pcms::HostMemorySpace> values("c", c.size(), 1);
  for (std::size_t i = 0; i < c.size(); ++i) {
    values(i, 0) = c[i];
  }
  field.SetDOFHolderDataHost(pcms::ValueView<const Real, pcms::HostMemorySpace>(
    field.GetData().GetValueBasis(), pcms::MakeConstRank2View(values)));
}

// Every transfer between periodic Lagrange spaces: a copy, an interpolator and
// the projections with the matrix-free solve
struct Transfers
{
  std::shared_ptr<const pcms::FunctionSpace> source_space = Space(8);
  std::shared_ptr<const pcms::FunctionSpace> target_space = Space(16);
  std::shared_ptr<const pcms::FunctionSpace> same_space = Space(8);
  std::vector<
    std::pair<const char*, std::unique_ptr<pcms::TransferOperator<Real>>>>
    operators;
  Transfers()
  {
    operators.emplace_back("copy",
                           pcms::CreateCopy<Real>(source_space, source_space));
    operators.emplace_back("interpolator", pcms::CreateInterpolator<Real>(
                                             source_space, target_space));
    operators.emplace_back("projection",
                           pcms::CreateProjection(source_space, target_space));
    pcms::ProjectionOptions lumped;
    lumped.mass_matrix_type = pcms::MassMatrixType::Lumped;
    operators.emplace_back(
      "lumped projection",
      pcms::CreateProjection(source_space, target_space, lumped));
    pcms::ProjectionOptions control_variate;
    control_variate.control_variate = true;
    operators.emplace_back(
      "control variate projection",
      pcms::CreateProjection(source_space, target_space, control_variate));
  }
  // the space of the target of operator `index`
  const pcms::FunctionSpace& TargetOf(std::size_t index) const
  {
    return index == 0 ? *source_space : *target_space;
  }
};

// the text of a field's units, "none" when they are unknown
std::string UnitsOf(const pcms::Field<Real>& field)
{
  const auto& units = field.GetData().GetValueBasis().units;
  return units ? units->ToString() : "none";
}

// A FieldFactory (no evaluation, no space) with no DOF holders
class TestFieldFactory : public pcms::FieldFactory
{
public:
  TestFieldFactory()
    : layout_(pcms::CreateEmptyLayout(pcms::csys::CreateCartesian(1), 1))
  {
  }

  [[nodiscard]] std::shared_ptr<const pcms::FieldLayout> GetLayout()
    const noexcept override
  {
    return layout_;
  }

protected:
  [[nodiscard]] pcms::FieldVariant CreateFieldImpl(
    pcms::Type /*storage_type*/, pcms::ValueBasis basis) const override
  {
    return WrapField<Real>(
      layout_,
      std::make_unique<pcms::detail::SimpleFieldData<Real>>(layout_, basis));
  }

  [[nodiscard]] pcms::FieldVariant CreateFieldImpl(
    pcms::FieldDataVariant data) const override
  {
    return WrapField<Real>(
      layout_,
      std::move(std::get<std::unique_ptr<pcms::FieldData<Real>>>(data)));
  }

private:
  std::shared_ptr<const pcms::FieldLayout> layout_;
};

// units checks only run with a units library
constexpr bool UnitsEnabled =
#ifdef PCMS_ENABLE_UNITS
  true;
#else
  false;
#endif

} // namespace

TEST_CASE("fields carry metadata, which travels with them", "[metadata]")
{
  const auto space = Space(8);
  auto plain = space->CreateFunction<Real>("b");
  REQUIRE(plain.GetName() == "b");
  REQUIRE_FALSE(plain.GetData().GetValueBasis().units.has_value());
  REQUIRE_FALSE(plain.GetMetadata().duality.has_value());

  auto field = space->CreateFunction<Real>(
    FieldMetadata{.name = "n", .duality = Duality::Primal},
    {.units = Units::Parse("m^-3")});
  REQUIRE(field.GetName() == "n");
  REQUIRE(UnitsOf(field) == Units::Parse("m^-3").ToString());
  REQUIRE(field.GetMetadata().duality == Duality::Primal);

  const auto moved = std::move(field);
  REQUIRE(moved.GetName() == "n");
  REQUIRE(moved.GetData().GetValueBasis().units.has_value());
  REQUIRE(moved.GetMetadata().duality == Duality::Primal);
}

TEST_CASE("creation takes metadata and a value basis", "[metadata]")
{
  const auto space = Space(8);

  SECTION("metadata with the default scalar basis")
  {
    const auto f = space->CreateFunction<Real>(
      FieldMetadata{.name = "p", .duality = Duality::Dual});
    REQUIRE(f.GetName() == "p");
    REQUIRE(f.GetMetadata().duality == Duality::Dual);
    REQUIRE(f.GetData().GetValueBasis().Rank() == 0);
    REQUIRE(f.GetData().GetValueBasis().scaling ==
            pcms::ComponentScaling::Physical);
  }
  SECTION("metadata with a vector basis: the space's Cartesian system, "
          "physical components")
  {
    const auto f = space->CreateFunction<Real>(
      FieldMetadata{.name = "v"},
      {.variance = pcms::values::Vector, .units = Units::Parse("m/s")});
    REQUIRE(f.GetName() == "v");
    const auto& basis = f.GetData().GetValueBasis();
    REQUIRE(basis.units.has_value());
    REQUIRE(basis.Rank() == 1);
    REQUIRE(pcms::detail::SameCoordinateSystem(basis.system,
                                               space->GetCoordinateSystem()));
    REQUIRE(basis.scaling == pcms::ComponentScaling::Physical);
  }
  SECTION("a vector on a system that is not orthonormal states its scaling")
  {
    for (const auto& system : {pcms::csys::CreateCylindricalRThetaZ(),
                               pcms::csys::CreateParameterSystem("xi", 1)}) {
      CAPTURE(system->Kind());
      REQUIRE_THROWS_WITH(
        space->CreateFunction<Real>(
          "v", {.system = system, .variance = pcms::values::Vector}),
        ContainsSubstring("scaling must be stated"));
    }
  }
  SECTION("the name shorthand takes the basis too")
  {
    const auto f = space->CreateFunction<Real>(
      "v", {.scaling = pcms::ComponentScaling::Natural,
            .variance = pcms::values::Vector});
    REQUIRE(f.GetName() == "v");
    REQUIRE_FALSE(f.GetData().GetValueBasis().units.has_value());
    REQUIRE(f.GetData().GetValueBasis().scaling ==
            pcms::ComponentScaling::Natural);
  }
  SECTION("a name alone and no arguments still work")
  {
    REQUIRE(space->CreateFunction<Real>("u").GetName() == "u");
    REQUIRE(space->CreateFunction<Real>().GetName().empty());
  }
  SECTION("wrapped data takes metadata")
  {
    auto data = std::make_unique<pcms::detail::SimpleFieldData<Real>>(
      space->GetLayout(), pcms::ValueBasis{});
    const auto f = space->CreateFunction<Real>(
      FieldMetadata{.name = "w", .duality = Duality::Primal}, std::move(data));
    REQUIRE(f.GetName() == "w");
    REQUIRE(f.GetMetadata().duality == Duality::Primal);
  }
}

TEST_CASE("a FieldFactory takes metadata and a value basis", "[metadata]")
{
  const TestFieldFactory factory;

  const auto f = factory.CreateField<Real>(
    FieldMetadata{.name = "q", .duality = Duality::Dual},
    {.units = Units::Parse("C")});
  REQUIRE(f.GetName() == "q");
  REQUIRE(f.GetData().GetValueBasis().units.has_value());
  REQUIRE(f.GetMetadata().duality == Duality::Dual);

  const auto vector =
    factory.CreateField<Real>("e", {.variance = pcms::values::Vector});
  REQUIRE(vector.GetName() == "e");
  REQUIRE(vector.GetData().GetValueBasis().Rank() == 1);
  REQUIRE(vector.GetData().GetValueBasis().scaling ==
          pcms::ComponentScaling::Physical);

  const auto named = factory.CreateField<Real>("s");
  REQUIRE(named.GetName() == "s");
  REQUIRE_FALSE(named.GetMetadata().duality.has_value());

  auto data = std::make_unique<pcms::detail::SimpleFieldData<Real>>(
    factory.GetLayout(), pcms::ValueBasis{});
  const auto wrapped = factory.CreateField<Real>(
    FieldMetadata{.name = "w", .duality = Duality::Primal}, std::move(data));
  REQUIRE(wrapped.GetMetadata().duality == Duality::Primal);
}

TEST_CASE("transfers check units", "[metadata]")
{
  Transfers transfers;
  for (std::size_t i = 0; i < transfers.operators.size(); ++i) {
    const auto& [name, op] = transfers.operators[i];
    CAPTURE(name);
    const auto& target_space = transfers.TargetOf(i);
    const auto make = [&](const char* source_units, const char* target_units) {
      pcms::ValueBasis s;
      pcms::ValueBasis t;
      if (source_units != nullptr) {
        s.units = Units::Parse(source_units);
      }
      if (target_units != nullptr) {
        t.units = Units::Parse(target_units);
      }
      auto source = transfers.source_space->CreateFunction<Real>("source", s);
      SetDofs(source, std::vector<Real>(9, 1.0));
      return std::make_pair(std::move(source),
                            target_space.CreateFunction<Real>("target", t));
    };

    SECTION("unknown units are not checked")
    {
      const std::vector<std::pair<const char*, const char*>> cases{
        {nullptr, nullptr}, {"T", nullptr}, {nullptr, "T"}};
      for (const auto& [s, t] : cases) {
        auto [source, target] = make(s, t);
        REQUIRE_NOTHROW(op->Apply(source, target));
      }
    }
    SECTION("the same units pass")
    {
      auto [source, target] = make("T", "kg/(A*s^2)");
      REQUIRE_NOTHROW(op->Apply(source, target));
    }
    if constexpr (UnitsEnabled) {
      SECTION("different dimensions throw")
      {
        auto [source, target] = make("T", "m");
        REQUIRE_THROWS_WITH(op->Apply(source, target),
                            ContainsSubstring("units [T] and [m]") &&
                              ContainsSubstring("different dimensions"));
      }
      SECTION("a different scale throws, and nothing is converted")
      {
        auto [source, target] = make("gauss", "T");
        REQUIRE_THROWS_WITH(
          op->Apply(source, target),
          ContainsSubstring("differ in scale or offset") &&
            ContainsSubstring("PCMS does not convert units, so convert one of "
                              "the fields"));
      }
    } else {
      SECTION("without a units library nothing is checked")
      {
        auto [source, target] = make("T", "m");
        REQUIRE_NOTHROW(op->Apply(source, target));
      }
    }
  }
}

TEST_CASE("transfers check the duality", "[metadata]")
{
  Transfers transfers;
  const auto make = [&](std::optional<Duality> source_duality,
                        std::optional<Duality> target_duality) {
    auto source = transfers.source_space->CreateFunction<Real>(
      FieldMetadata{.name = "source", .duality = source_duality});
    SetDofs(source, std::vector<Real>(9, 1.0));
    return std::make_pair(
      std::move(source),
      transfers.target_space->CreateFunction<Real>(
        FieldMetadata{.name = "target", .duality = target_duality}));
  };
  const auto& op = transfers.operators[1].second; // the interpolator
  REQUIRE(op->SupportsTranspose());

  SECTION("Apply refuses a dual side, in either direction")
  {
    auto [a, b] = make(Duality::Dual, std::nullopt);
    REQUIRE_THROWS_WITH(op->Apply(a, b), ContainsSubstring("'source' is dual"));
    auto [c, d] = make(std::nullopt, Duality::Dual);
    REQUIRE_THROWS_WITH(op->Apply(c, d), ContainsSubstring("'target' is dual"));
  }
  SECTION("Apply takes primal and unknown dualities")
  {
    auto [a, b] = make(Duality::Primal, Duality::Primal);
    REQUIRE_NOTHROW(op->Apply(a, b));
    auto [c, d] = make(std::nullopt, std::nullopt);
    REQUIRE_NOTHROW(op->Apply(c, d));
    auto [e, f] = make(Duality::Primal, std::nullopt);
    REQUIRE_NOTHROW(op->Apply(e, f));
  }
  SECTION("ApplyTranspose refuses a primal side, in either direction")
  {
    // loads are on the target space, the functional on the source space
    const auto loads = [&](std::optional<Duality> e) {
      auto f = transfers.target_space->CreateFunction<Real>(
        FieldMetadata{.name = "loads", .duality = e});
      SetDofs(f, std::vector<Real>(17, 1.0));
      return f;
    };
    const auto functional = [&](std::optional<Duality> e) {
      return transfers.source_space->CreateFunction<Real>(
        FieldMetadata{.name = "functional", .duality = e});
    };
    auto l1 = loads(Duality::Primal);
    auto f1 = functional(std::nullopt);
    REQUIRE_THROWS_WITH(op->ApplyTranspose(l1, f1),
                        ContainsSubstring("'loads' is primal"));
    auto l2 = loads(std::nullopt);
    auto f2 = functional(Duality::Primal);
    REQUIRE_THROWS_WITH(op->ApplyTranspose(l2, f2),
                        ContainsSubstring("'functional' is primal"));
    auto l3 = loads(Duality::Dual);
    auto f3 = functional(Duality::Dual);
    REQUIRE_NOTHROW(op->ApplyTranspose(l3, f3));
    auto l4 = loads(std::nullopt);
    auto f4 = functional(std::nullopt);
    REQUIRE_NOTHROW(op->ApplyTranspose(l4, f4));
  }
}

TEST_CASE("a transfer does not change the target's metadata", "[metadata]")
{
  Transfers transfers;
  const auto& op = transfers.operators[1].second;
  auto source = transfers.source_space->CreateFunction<Real>(
    FieldMetadata{.name = "source", .duality = Duality::Primal},
    {.units = Units::Parse("m^-3")});
  SetDofs(source, std::vector<Real>(9, 1.0));
  auto target = transfers.target_space->CreateFunction<Real>("target");
  op->Apply(source, target);
  REQUIRE(target.GetName() == "target");
  REQUIRE_FALSE(target.GetData().GetValueBasis().units.has_value());
  REQUIRE_FALSE(target.GetMetadata().duality.has_value());
}

TEST_CASE("ApplyTranspose checks units", "[metadata]")
{
  Transfers transfers;
  // the interpolator and the projections of the periodic Lagrange spaces
  // support the transpose; the copy does not
  for (std::size_t i = 1; i < 4; ++i) {
    const auto& [name, op] = transfers.operators[i];
    CAPTURE(name);
    REQUIRE(op->SupportsTranspose());
    const auto run = [&](const char* loads_units,
                         const char* functional_units) {
      pcms::ValueBasis l;
      pcms::ValueBasis f;
      if (loads_units != nullptr) {
        l.units = Units::Parse(loads_units);
      }
      if (functional_units != nullptr) {
        f.units = Units::Parse(functional_units);
      }
      auto loads = transfers.target_space->CreateFunction<Real>("loads", l);
      auto functional =
        transfers.source_space->CreateFunction<Real>("functional", f);
      SetDofs(loads, std::vector<Real>(17, 1.0));
      op->ApplyTranspose(loads, functional);
    };

    REQUIRE_NOTHROW(run(nullptr, nullptr));
    REQUIRE_NOTHROW(run("N", nullptr));
    REQUIRE_NOTHROW(run("N", "kg*m/s^2"));
    if constexpr (UnitsEnabled) {
      REQUIRE_THROWS_WITH(run("N", "m"),
                          ContainsSubstring("units [N] and [m]"));
      REQUIRE_THROWS_WITH(run("N", "dyn"),
                          ContainsSubstring("differ in scale or offset"));
    } else {
      REQUIRE_NOTHROW(run("N", "m"));
    }
  }
}

TEST_CASE("a transfer without a transpose still reports it", "[metadata]")
{
  Transfers transfers;
  const auto& copy = transfers.operators[0].second;
  REQUIRE_FALSE(copy->SupportsTranspose());
  auto a = transfers.source_space->CreateFunction<Real>("a");
  auto b = transfers.source_space->CreateFunction<Real>("b");
  REQUIRE_THROWS_WITH(copy->ApplyTranspose(a, b),
                      ContainsSubstring("not supported"));
}

TEST_CASE("a transfer refuses a source and a target that share their data",
          "[metadata]")
{
  Transfers transfers;
  auto f = transfers.source_space->CreateFunction<Real>("f");
  // the interpolator's transpose and the projection
  REQUIRE_THROWS_WITH(transfers.operators[1].second->ApplyTranspose(f, f),
                      ContainsSubstring("share their data"));
  REQUIRE_THROWS_WITH(transfers.operators[2].second->Apply(f, f),
                      ContainsSubstring("share their data"));
}

TEST_CASE("L2Integral::Distance checks the units of both fields; Inner does "
          "not",
          "[metadata]")
{
  const auto space = Space(8);
  const auto other = Space(16);
  auto u = space->CreateFunction<Real>("u", {.units = Units::Parse("m")});
  auto v = other->CreateFunction<Real>("v", {.units = Units::Parse("cm")});
  auto w = other->CreateFunction<Real>("w", {.units = Units::Parse("m")});
  const pcms::L2Integral integral(space, other);
  REQUIRE_NOTHROW(integral.Distance(u, w));
  REQUIRE_NOTHROW(integral.Inner(u, v));
  if constexpr (UnitsEnabled) {
    REQUIRE_THROWS_WITH(integral.Distance(u, v),
                        ContainsSubstring("L2Integral::Distance"));
  } else {
    REQUIRE_NOTHROW(integral.Distance(u, v));
  }
}

TEST_CASE("Axpy checks the value bases", "[metadata]")
{
  const auto space = Space(8);
  auto x = space->CreateFunction<Real>("x", {.units = Units::Parse("m")});
  auto same = space->CreateFunction<Real>("y", {.units = Units::Parse("m")});
  auto unset = space->CreateFunction<Real>("z");
  auto vector =
    space->CreateFunction<Real>("v", {.variance = pcms::values::Vector});
  REQUIRE_NOTHROW(pcms::Axpy(1.0, x, same));
  REQUIRE_NOTHROW(pcms::Axpy(1.0, x, unset));
  REQUIRE_THROWS_WITH(pcms::Axpy(1.0, x, vector),
                      ContainsSubstring("value bases differ"));
  if constexpr (UnitsEnabled) {
    auto kelvin =
      space->CreateFunction<Real>("k", {.units = Units::Parse("K")});
    REQUIRE_THROWS_WITH(pcms::Axpy(1.0, x, kelvin),
                        ContainsSubstring("different dimensions"));
  }
}
