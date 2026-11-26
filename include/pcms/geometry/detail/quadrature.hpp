// Part of pcms/geometry/quadrature.hpp, which includes it after the
// declarations it needs.
#include "pcms/geometry/quadrature.hpp"
#ifndef PCMS_GEOMETRY_DETAIL_QUADRATURE_HPP
#define PCMS_GEOMETRY_DETAIL_QUADRATURE_HPP

#include <cstdint>
#include <span>
#include <variant>

namespace pcms::detail
{

struct GaussRule
{
  int order;
};

struct MonteCarloRule
{
  int samples_per_cell;
  std::uint64_t seed;
};

/// @brief The rule a quadrature factory asks a geometry for.
using QuadratureRule = std::variant<GaussRule, MonteCarloRule>;

/// @brief The quadrature of `geometry` by `rule`, for the factories and the
///        kinds built on other geometries.
[[nodiscard]] Quadrature CreateQuadrature(const Geometry& geometry,
                                          const QuadratureRule& rule);

/// @brief Replaces each value by its sum over the ranks the entities of
///        `geometry` are partitioned over, e.g. a quadrature sum. Collective.
void SumOverRanks(const Geometry& geometry, std::span<Real> values);

/// @brief Replaces each value by its maximum over the ranks the entities of
///        `geometry` are partitioned over. Collective.
void MaxOverRanks(const Geometry& geometry, std::span<Real> values);

} // namespace pcms::detail

#endif // PCMS_GEOMETRY_DETAIL_QUADRATURE_HPP
