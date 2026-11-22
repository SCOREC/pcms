#include "periodic_test_utils.hpp"
#include <catch2/matchers/catch_matchers_string.hpp>
#include <pcms/configuration.hpp>
#include <pcms/field/detail/simple_field_data.hpp>
#include <pcms/field/field_metadata.hpp>
#include <pcms/transfer/copy.hpp>
#include <pcms/transfer/interpolator.hpp>
#include <pcms/transfer/projection.hpp>
#include <pcms/utility/coordinate_systems/cartesian.hpp>

using Catch::Matchers::ContainsSubstring;
using pcms::Extent;
using pcms::FieldMetadata;
using pcms::Real;
using pcms::Units;

namespace
{

std::shared_ptr<const pcms::FunctionSpace> Space(int cells)
{
  return pcms::CreateLagrangeSpace(pcms::test::PeriodicGrid(cells, 1.0, 0.0), 1,
                                   1);
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
  const auto& units = field.GetMetadata().units;
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
  REQUIRE_FALSE(plain.GetMetadata().units.has_value());
  REQUIRE_FALSE(plain.GetMetadata().extent.has_value());

  auto field = space->CreateFunction<Real>(FieldMetadata{
    .name = "n", .units = Units::Parse("m^-3"), .extent = Extent::Intensive});
  REQUIRE(field.GetName() == "n");
  REQUIRE(UnitsOf(field) == Units::Parse("m^-3").ToString());
  REQUIRE(field.GetMetadata().extent == Extent::Intensive);

  const auto moved = std::move(field);
  REQUIRE(moved.GetName() == "n");
  REQUIRE(moved.GetMetadata().units.has_value());
  REQUIRE(moved.GetMetadata().extent == Extent::Intensive);
}

TEST_CASE("creation takes metadata and a value basis", "[metadata]")
{
  const auto space = Space(8);

  SECTION("metadata with the default scalar basis")
  {
    const auto f = space->CreateFunction<Real>(
      FieldMetadata{.name = "p", .extent = Extent::Extensive});
    REQUIRE(f.GetName() == "p");
    REQUIRE(f.GetMetadata().extent == Extent::Extensive);
    REQUIRE(f.GetData().GetValueBasis().Rank() == 0);
  }
  SECTION("metadata with a vector basis: system and scaling are derived")
  {
    const auto f = space->CreateFunction<Real>(
      FieldMetadata{.name = "v", .units = Units::Parse("m/s")},
      {.variance = pcms::values::Vector});
    REQUIRE(f.GetName() == "v");
    REQUIRE(f.GetMetadata().units.has_value());
    const auto& basis = f.GetData().GetValueBasis();
    REQUIRE(basis.Rank() == 1);
    REQUIRE(pcms::detail::SameCoordinateSystem(basis.system,
                                               space->GetCoordinateSystem()));
    REQUIRE(basis.scaling == pcms::ComponentScaling::Physical);
  }
  SECTION("the name shorthand takes the basis too")
  {
    const auto f = space->CreateFunction<Real>(
      "v", {.scaling = pcms::ComponentScaling::Natural,
            .variance = pcms::values::Vector});
    REQUIRE(f.GetName() == "v");
    REQUIRE_FALSE(f.GetMetadata().units.has_value());
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
      FieldMetadata{.name = "w", .extent = Extent::Intensive}, std::move(data));
    REQUIRE(f.GetName() == "w");
    REQUIRE(f.GetMetadata().extent == Extent::Intensive);
  }
}

TEST_CASE("a FieldFactory takes metadata and a value basis", "[metadata]")
{
  const TestFieldFactory factory;

  const auto f = factory.CreateField<Real>(FieldMetadata{
    .name = "q", .units = Units::Parse("C"), .extent = Extent::Extensive});
  REQUIRE(f.GetName() == "q");
  REQUIRE(f.GetMetadata().units.has_value());
  REQUIRE(f.GetMetadata().extent == Extent::Extensive);

  const auto vector =
    factory.CreateField<Real>("e", {.variance = pcms::values::Vector});
  REQUIRE(vector.GetName() == "e");
  REQUIRE(vector.GetData().GetValueBasis().Rank() == 1);
  REQUIRE(vector.GetData().GetValueBasis().scaling ==
          pcms::ComponentScaling::Physical);

  const auto named = factory.CreateField<Real>("s");
  REQUIRE(named.GetName() == "s");
  REQUIRE_FALSE(named.GetMetadata().extent.has_value());

  auto data = std::make_unique<pcms::detail::SimpleFieldData<Real>>(
    factory.GetLayout(), pcms::ValueBasis{});
  const auto wrapped = factory.CreateField<Real>(
    FieldMetadata{.name = "w", .extent = Extent::Intensive}, std::move(data));
  REQUIRE(wrapped.GetMetadata().extent == Extent::Intensive);
}

TEST_CASE("transfers check units", "[metadata]")
{
  Transfers transfers;
  for (std::size_t i = 0; i < transfers.operators.size(); ++i) {
    const auto& [name, op] = transfers.operators[i];
    CAPTURE(name);
    const auto& target_space = transfers.TargetOf(i);
    const auto make = [&](const char* source_units, const char* target_units) {
      FieldMetadata s{.name = "source"};
      FieldMetadata t{.name = "target"};
      if (source_units != nullptr) {
        s.units = Units::Parse(source_units);
      }
      if (target_units != nullptr) {
        t.units = Units::Parse(target_units);
      }
      auto source = transfers.source_space->CreateFunction<Real>(s);
      pcms::test::SetDofs(source, std::vector<Real>(8, 1.0));
      return std::make_pair(std::move(source),
                            target_space.CreateFunction<Real>(t));
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
                            ContainsSubstring("transfer from 'source' [T] to "
                                              "'target' [m]") &&
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

TEST_CASE("transfers check the extent", "[metadata]")
{
  Transfers transfers;
  const auto make = [&](std::optional<Extent> source_extent,
                        std::optional<Extent> target_extent) {
    auto source = transfers.source_space->CreateFunction<Real>(
      FieldMetadata{.name = "source", .extent = source_extent});
    pcms::test::SetDofs(source, std::vector<Real>(8, 1.0));
    return std::make_pair(
      std::move(source),
      transfers.target_space->CreateFunction<Real>(
        FieldMetadata{.name = "target", .extent = target_extent}));
  };
  const auto& op = transfers.operators[1].second; // the interpolator
  REQUIRE(op->SupportsTranspose());

  SECTION("Apply refuses an extensive side, in either direction")
  {
    auto [a, b] = make(Extent::Extensive, std::nullopt);
    REQUIRE_THROWS_WITH(op->Apply(a, b),
                        ContainsSubstring("'source' is extensive"));
    auto [c, d] = make(std::nullopt, Extent::Extensive);
    REQUIRE_THROWS_WITH(op->Apply(c, d),
                        ContainsSubstring("'target' is extensive"));
  }
  SECTION("Apply takes intensive and unknown extents")
  {
    auto [a, b] = make(Extent::Intensive, Extent::Intensive);
    REQUIRE_NOTHROW(op->Apply(a, b));
    auto [c, d] = make(std::nullopt, std::nullopt);
    REQUIRE_NOTHROW(op->Apply(c, d));
    auto [e, f] = make(Extent::Intensive, std::nullopt);
    REQUIRE_NOTHROW(op->Apply(e, f));
  }
  SECTION("ApplyTranspose refuses an intensive side, in either direction")
  {
    // loads are on the target space, the functional on the source space
    const auto loads = [&](std::optional<Extent> e) {
      auto f = transfers.target_space->CreateFunction<Real>(
        FieldMetadata{.name = "loads", .extent = e});
      pcms::test::SetDofs(f, std::vector<Real>(16, 1.0));
      return f;
    };
    const auto functional = [&](std::optional<Extent> e) {
      return transfers.source_space->CreateFunction<Real>(
        FieldMetadata{.name = "functional", .extent = e});
    };
    auto l1 = loads(Extent::Intensive);
    auto f1 = functional(std::nullopt);
    REQUIRE_THROWS_WITH(op->ApplyTranspose(l1, f1),
                        ContainsSubstring("'loads' is intensive"));
    auto l2 = loads(std::nullopt);
    auto f2 = functional(Extent::Intensive);
    REQUIRE_THROWS_WITH(op->ApplyTranspose(l2, f2),
                        ContainsSubstring("'functional' is intensive"));
    auto l3 = loads(Extent::Extensive);
    auto f3 = functional(Extent::Extensive);
    REQUIRE_NOTHROW(op->ApplyTranspose(l3, f3));
    auto l4 = loads(std::nullopt);
    auto f4 = functional(std::nullopt);
    REQUIRE_NOTHROW(op->ApplyTranspose(l4, f4));
  }
}

TEST_CASE("a transfer teaches the target what it does not know", "[metadata]")
{
  Transfers transfers;
  const auto& op = transfers.operators[1].second;
  auto source = transfers.source_space->CreateFunction<Real>(
    FieldMetadata{.name = "source",
                  .units = Units::Parse("m^-3"),
                  .extent = Extent::Intensive});
  pcms::test::SetDofs(source, std::vector<Real>(8, 1.0));

  SECTION("unknown units and extent are learned")
  {
    auto target = transfers.target_space->CreateFunction<Real>("target");
    op->Apply(source, target);
    REQUIRE(target.GetName() == "target");
    REQUIRE(target.GetMetadata().units.has_value());
    REQUIRE(UnitsOf(target) == UnitsOf(source));
    REQUIRE(target.GetMetadata().extent == Extent::Intensive);
  }
  SECTION("known values are never overwritten")
  {
    // the same units, written differently: the target keeps its own
    auto target = transfers.target_space->CreateFunction<Real>(
      FieldMetadata{.name = "target", .units = Units::Parse("1/m^3")});
    op->Apply(source, target);
    REQUIRE(UnitsOf(target) == Units::Parse("1/m^3").ToString());
    // the extent is learned alone
    REQUIRE(target.GetMetadata().extent == Extent::Intensive);
    if constexpr (!UnitsEnabled) {
      // without a units library nothing is compared, so units that differ
      // pass, and the known ones stay
      auto metre = transfers.target_space->CreateFunction<Real>(
        FieldMetadata{.name = "metre", .units = Units::Parse("m")});
      op->Apply(source, metre);
      REQUIRE(UnitsOf(metre) == "m");
    }
  }
  SECTION("what a failed transfer was refused for is not learned")
  {
    auto bad = transfers.source_space->CreateFunction<Real>(
      FieldMetadata{.name = "bad",
                    .units = Units::Parse("m^-3"),
                    .extent = Extent::Extensive});
    auto target = transfers.target_space->CreateFunction<Real>("target");
    REQUIRE_THROWS(op->Apply(bad, target));
    REQUIRE_FALSE(target.GetMetadata().units.has_value());
    REQUIRE_FALSE(target.GetMetadata().extent.has_value());
  }
  SECTION("ApplyTranspose learns too")
  {
    auto loads = transfers.target_space->CreateFunction<Real>(
      FieldMetadata{.name = "loads",
                    .units = Units::Parse("N"),
                    .extent = Extent::Extensive});
    pcms::test::SetDofs(loads, std::vector<Real>(16, 1.0));
    auto functional =
      transfers.source_space->CreateFunction<Real>("functional");
    op->ApplyTranspose(loads, functional);
    REQUIRE(functional.GetMetadata().units.has_value());
    REQUIRE(functional.GetMetadata().extent == Extent::Extensive);
  }
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
      FieldMetadata l{.name = "loads"};
      FieldMetadata f{.name = "functional"};
      if (loads_units != nullptr) {
        l.units = Units::Parse(loads_units);
      }
      if (functional_units != nullptr) {
        f.units = Units::Parse(functional_units);
      }
      auto loads = transfers.target_space->CreateFunction<Real>(l);
      auto functional = transfers.source_space->CreateFunction<Real>(f);
      pcms::test::SetDofs(loads, std::vector<Real>(16, 1.0));
      op->ApplyTranspose(loads, functional);
    };

    REQUIRE_NOTHROW(run(nullptr, nullptr));
    REQUIRE_NOTHROW(run("N", nullptr));
    REQUIRE_NOTHROW(run("N", "kg*m/s^2"));
    if constexpr (UnitsEnabled) {
      REQUIRE_THROWS_WITH(run("N", "m"),
                          ContainsSubstring("transposed transfer from 'loads' "
                                            "[N] to 'functional' [m]"));
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
