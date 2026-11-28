#ifndef PCMS_FIELD_DETAIL_FIELD_OPS_HPP
#define PCMS_FIELD_DETAIL_FIELD_OPS_HPP

#include "pcms/field/field.hpp"
#include "pcms/utility/assert.hpp"
#include <cmath>
#include <span>
#include <string>
#include <typeinfo>

namespace pcms::detail
{

// The coefficient algebra of the solvers and projections over Fields: each
// function checks its operands, then calls the FieldData coefficient
// operation of the same name. They act on the owned coefficients; those that
// write synchronize the field, and the reductions are collective over the
// ranks the data is partitioned over.

template <typename T>
void CheckSameLayout(const Field<T>& a, const Field<T>& b, const char* what)
{
  if (&a.GetLayout() != &b.GetLayout()) {
    throw pcms_error(std::string(what) +
                     ": the fields must have the same layout");
  }
}

/// @throws pcms_error unless `a` and `b` have the same layout, value basis and
///         kind of field data
template <typename T>
void CheckOperands(const Field<T>& a, const Field<T>& b, const char* what)
{
  CheckSameLayout(a, b, what);
  if (!SameValueBasis(a.GetData().GetValueBasis(),
                      b.GetData().GetValueBasis())) {
    throw pcms_error(std::string(what) + ": the fields' value bases differ");
  }
  if (typeid(a.GetData()) != typeid(b.GetData())) {
    throw pcms_error(std::string(what) +
                     ": the fields' data are of different kinds");
  }
}

/// @brief z = a x + b y
template <typename T>
void LinearSum(T a, const Field<T>& x, T b, const Field<T>& y, Field<T>& z)
{
  CheckOperands(x, z, "LinearSum");
  CheckOperands(y, z, "LinearSum");
  LinearSum(a, x.GetData(), b, y.GetData(), z.GetData());
  z.Synchronize();
}

/// @brief z = c x
template <typename T>
void Scale(T c, const Field<T>& x, Field<T>& z)
{
  CheckOperands(x, z, "Scale");
  Scale(c, x.GetData(), z.GetData());
  z.Synchronize();
}

/// @brief z = c
template <typename T>
void Const(T c, Field<T>& z)
{
  Const(c, z.GetData());
  z.Synchronize();
}

/// @brief z = x / y, entry by entry
template <typename T>
void Div(const Field<T>& x, const Field<T>& y, Field<T>& z)
{
  CheckOperands(x, z, "Div");
  CheckOperands(y, z, "Div");
  Div(x.GetData(), y.GetData(), z.GetData());
  z.Synchronize();
}

/// @return the sum over the owned coefficients of every rank of x_i y_i
template <typename T>
[[nodiscard]] T DotProd(const Field<T>& x, const Field<T>& y)
{
  CheckOperands(x, y, "DotProd");
  T sum = DotProdLocal(x.GetData(), y.GetData());
  SumOverRanks(x.GetData(), std::span<T>(&sum, 1));
  return sum;
}

/// @return the smallest owned coefficient over every rank
template <typename T>
[[nodiscard]] T Min(const Field<T>& x)
{
  T negated = -MinLocal(x.GetData());
  MaxOverRanks(x.GetData(), std::span<T>(&negated, 1));
  return -negated;
}

/// @return the largest owned coefficient magnitude over every rank
template <typename T>
[[nodiscard]] T MaxNorm(const Field<T>& x)
{
  T max = MaxNormLocal(x.GetData());
  MaxOverRanks(x.GetData(), std::span<T>(&max, 1));
  return max;
}

/// @return the Euclidean norm of the coefficients
template <typename T>
[[nodiscard]] T Norm2(const Field<T>& x)
{
  return std::sqrt(DotProd(x, x));
}

/// @brief The pairing of a functional (e.g. the moments E^T W f of a function
///        f) with a function: l(u) = sum_i b_i u_i. Basis-independent: the
///        functional is dual to the coefficients.
template <typename T>
[[nodiscard]] T Apply(const Field<T>& functional, const Function<T>& u)
{
  return DotProd(functional, static_cast<const Field<T>&>(u));
}

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_FIELD_OPS_HPP
