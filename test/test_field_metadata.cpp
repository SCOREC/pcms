#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <pcms/field/field_metadata.hpp>
#include <pcms/field/function_spaces.hpp>
#include <pcms/transfer/copy.hpp>
#include <pcms/transfer/interpolator.hpp>
#include <pcms/transfer/projection.hpp>
#include "field_test_utils.hpp"
#include <optional>
#include <utility>

using Catch::Matchers::ContainsSubstring;
using pcms::Duality;
using pcms::FieldMetadata;
using pcms::Real;

namespace
{

// A copy, an interpolator and a projection between P1 spaces of two meshes of
// the unit square
struct Transfers
{
  explicit Transfers(Omega_h::Library& lib)
    : source_mesh(pcms::test::BuildSquareMesh(lib, 4)),
      target_mesh(pcms::test::BuildSquareMesh(lib, 8)),
      source_space(pcms::test::MakeP1Space(source_mesh)),
      target_space(pcms::test::MakeP1Space(target_mesh))
  {
  }

  [[nodiscard]] std::unique_ptr<pcms::TransferOperator<Real>> Copy() const
  {
    return pcms::CreateCopy<Real>(source_space, source_space);
  }

  [[nodiscard]] std::unique_ptr<pcms::TransferOperator<Real>> Interpolator()
    const
  {
    return pcms::CreateInterpolator<Real>(source_space, target_space);
  }

  [[nodiscard]] std::unique_ptr<pcms::TransferOperator<Real>> Projection() const
  {
    return pcms::CreateProjection(source_space, target_space);
  }

  Omega_h::Mesh source_mesh;
  Omega_h::Mesh target_mesh;
  std::shared_ptr<const pcms::FunctionSpace> source_space;
  std::shared_ptr<const pcms::FunctionSpace> target_space;
};

void SetOnes(pcms::Field<Real>& field)
{
  pcms::test::SetField(field, OMEGA_H_LAMBDA(Real, Real) { return 1.0; });
}

} // namespace

TEST_CASE("fields carry metadata, which travels with them", "[metadata]")
{
  Omega_h::Library lib;
  const Transfers transfers(lib);
  auto plain = transfers.source_space->CreateFunction<Real>("b");
  REQUIRE(plain.GetName() == "b");
  REQUIRE_FALSE(plain.GetMetadata().duality.has_value());

  auto field = transfers.source_space->CreateFunction<Real>(
    FieldMetadata{.name = "n", .duality = Duality::Primal});
  REQUIRE(field.GetName() == "n");
  REQUIRE(field.GetMetadata().duality == Duality::Primal);

  const auto moved = std::move(field);
  REQUIRE(moved.GetName() == "n");
  REQUIRE(moved.GetMetadata().duality == Duality::Primal);
}

TEST_CASE("transfers check the duality", "[metadata]")
{
  Omega_h::Library lib;
  const Transfers transfers(lib);
  const auto make = [&](std::optional<Duality> source_duality,
                        std::optional<Duality> target_duality) {
    auto source = transfers.source_space->CreateFunction<Real>(
      FieldMetadata{.name = "source", .duality = source_duality});
    SetOnes(source);
    return std::make_pair(
      std::move(source),
      transfers.target_space->CreateFunction<Real>(
        FieldMetadata{.name = "target", .duality = target_duality}));
  };
  const auto op = transfers.Interpolator();
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
      SetOnes(f);
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
  Omega_h::Library lib;
  const Transfers transfers(lib);
  const auto op = transfers.Interpolator();
  auto source = transfers.source_space->CreateFunction<Real>(
    FieldMetadata{.name = "source", .duality = Duality::Primal});
  SetOnes(source);
  auto target = transfers.target_space->CreateFunction<Real>("target");
  op->Apply(source, target);
  REQUIRE(target.GetName() == "target");
  REQUIRE_FALSE(target.GetMetadata().duality.has_value());
}

TEST_CASE("a transfer without a transpose still reports it", "[metadata]")
{
  Omega_h::Library lib;
  const Transfers transfers(lib);
  const auto copy = transfers.Copy();
  REQUIRE_FALSE(copy->SupportsTranspose());
  auto a = transfers.source_space->CreateFunction<Real>("a");
  auto b = transfers.source_space->CreateFunction<Real>("b");
  REQUIRE_THROWS_WITH(copy->ApplyTranspose(a, b),
                      ContainsSubstring("not supported"));
}

TEST_CASE("a transfer refuses a source and a target that share their data",
          "[metadata]")
{
  Omega_h::Library lib;
  const Transfers transfers(lib);
  auto f = transfers.source_space->CreateFunction<Real>("f");
  REQUIRE_THROWS_WITH(transfers.Interpolator()->ApplyTranspose(f, f),
                      ContainsSubstring("share their data"));
  REQUIRE_THROWS_WITH(transfers.Projection()->Apply(f, f),
                      ContainsSubstring("share their data"));
}
