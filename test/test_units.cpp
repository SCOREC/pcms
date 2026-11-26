#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "pcms/configuration.hpp"
#include "pcms/field/value_view.hpp"
#include "pcms/field/units.hpp"
#include "pcms/utility/assert.hpp"

using Catch::Matchers::ContainsSubstring;
using pcms::Units;
using pcms::UnitSyntax;

namespace
{
bool Same(const Units& a, const Units& b)
{
  return pcms::detail::SameUnits(a, b);
}
bool SameDimensions(const Units& a, const Units& b)
{
  return pcms::detail::SameDimensions(a, b);
}
} // namespace

TEST_CASE("Units keep what they are given, whatever the build")
{
  REQUIRE_THROWS_AS(Units::Parse(""), pcms::pcms_error);
  REQUIRE_THROWS_AS(Units::Parse("  "), pcms::pcms_error);
  REQUIRE_NOTHROW(Units::Parse("T"));
  REQUIRE_NOTHROW(Units::Parse(" T "));
}

#ifndef PCMS_ENABLE_UNITS
TEST_CASE("Units without a units library are their text")
{
  REQUIRE(Units::Parse("  kg/(A*s^2) ").ToString() == "kg/(A*s^2)");
  REQUIRE(Units::Parse("not a unit").ToString() == "not a unit");
  REQUIRE(Units::Parse("m2 s-1", UnitSyntax::CF).ToString() == "m2 s-1");
  REQUIRE(Units::Parse("mixed", UnitSyntax::IMAS).ToString() == "mixed");
  // nothing is compared
  REQUIRE(Same(Units::Parse("T"), Units::Parse("m")));
  REQUIRE(
    pcms::detail::SameValueBasis(pcms::ValueBasis{.units = Units::Parse("T")},
                                 pcms::ValueBasis{.units = Units::Parse("m")}));
}
#else
TEST_CASE("Units parse SI, IMAS and fusion notations")
{
  REQUIRE(Same(Units::Parse("T"), Units::Parse("kg/(A*s^2)")));
  REQUIRE(Same(Units::Parse("Tesla"), Units::Parse("T")));
  REQUIRE(Same(Units::Parse("m/s"), Units::Parse("m.s^-1")));
  REQUIRE(Same(Units::Parse("A.m^-2"), Units::Parse("A/m^2")));
  REQUIRE(Same(Units::Parse("m^-3"), Units::Parse("1/m^3")));
  REQUIRE(SameDimensions(Units::Parse("  eV "), Units::Parse("J")));
  REQUIRE_FALSE(Same(Units::Parse("eV"), Units::Parse("J")));
}

TEST_CASE("Units compare dimensions, scale and offset")
{
  const auto tesla = Units::Parse("T");
  const auto gauss = Units::Parse("gauss");
  REQUIRE(SameDimensions(tesla, gauss));
  REQUIRE_FALSE(Same(tesla, gauss));
  REQUIRE_FALSE(SameDimensions(tesla, Units::Parse("m")));
  REQUIRE(SameDimensions(Units::Parse("degC"), Units::Parse("K")));
  REQUIRE_FALSE(Same(Units::Parse("degC"), Units::Parse("K")));
}

TEST_CASE("Units round-trip through ToString")
{
  for (const char* text : {"T", "gauss", "m/s", "A/m^2", "eV", "1/m^3", "Pa"}) {
    const auto units = Units::Parse(text);
    REQUIRE(Same(Units::Parse(units.ToString()), units));
  }
}

TEST_CASE("Units parse CF strings")
{
  REQUIRE(Same(Units::Parse("m2 s-1", UnitSyntax::CF), Units::Parse("m^2/s")));
  REQUIRE(
    Same(Units::Parse("J kg-1 K-1", UnitSyntax::CF), Units::Parse("J/(kg*K)")));
  REQUIRE(Same(Units::Parse("kg m-2 s-1", UnitSyntax::CF),
               Units::Parse("kg/(m^2*s)")));
  // Without the CF syntax LLNL/units reads this as 1/mol^2.
  REQUIRE(Same(Units::Parse("mol mol-1", UnitSyntax::CF), Units::Parse("1")));
  REQUIRE(Same(Units::Parse("degrees_Celsius", UnitSyntax::CF),
               Units::Parse("degC")));
  REQUIRE_THROWS_AS(Units::Parse("days since 2000-01-01", UnitSyntax::CF),
                    pcms::pcms_error);
}

TEST_CASE("Units parse IMAS strings")
{
  REQUIRE(Same(Units::Parse("-", UnitSyntax::IMAS), Units::Parse("1")));
  REQUIRE(
    Same(Units::Parse("W.m^-3", UnitSyntax::IMAS), Units::Parse("W/m^3")));
  REQUIRE_THROWS_WITH(Units::Parse("mixed", UnitSyntax::IMAS),
                      ContainsSubstring("leave the units unspecified"));
  REQUIRE_THROWS_AS(Units::Parse("as_parent", UnitSyntax::IMAS),
                    pcms::pcms_error);
}

TEST_CASE("Units reject what is not a unit")
{
  REQUIRE_THROWS_AS(Units::Parse(""), pcms::pcms_error);
  REQUIRE_THROWS_AS(Units::Parse("  "), pcms::pcms_error);
  REQUIRE_THROWS_AS(Units::Parse("-"), pcms::pcms_error);
  REQUIRE_THROWS_WITH(Units::Parse("UNITLESS"),
                      ContainsSubstring("\"UNITLESS\" is not a unit"));
}

TEST_CASE("value bases compare units")
{
  const auto basis = [](const char* units) {
    pcms::ValueBasis b;
    if (units != nullptr) {
      b.units = Units::Parse(units);
    }
    return b;
  };
  const auto tesla = basis("T");

  SECTION("the same or unset units match")
  {
    REQUIRE(pcms::detail::SameValueBasis(tesla, tesla));
    REQUIRE(pcms::detail::SameValueBasis(tesla, basis("kg/(A*s^2)")));
    REQUIRE(pcms::detail::SameValueBasis(tesla, basis(nullptr)));
    REQUIRE(pcms::detail::SameValueBasis(basis(nullptr), tesla));
    REQUIRE(pcms::detail::UnitsMismatch(tesla, basis(nullptr)).empty());
  }
  SECTION("different dimensions")
  {
    REQUIRE_FALSE(pcms::detail::SameValueBasis(tesla, basis("m")));
    REQUIRE_THAT(pcms::detail::UnitsMismatch(tesla, basis("m")),
                 ContainsSubstring("units [T] and [m]") &&
                   ContainsSubstring("different dimensions"));
  }
  SECTION("a different scale")
  {
    REQUIRE_FALSE(pcms::detail::SameValueBasis(basis("gauss"), tesla));
    REQUIRE_THAT(pcms::detail::UnitsMismatch(basis("gauss"), tesla),
                 ContainsSubstring("differ in scale or offset"));
  }
}
#endif // PCMS_ENABLE_UNITS
