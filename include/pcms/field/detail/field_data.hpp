// Part of pcms/field/field_data.hpp, which includes it after declaring
// FieldData and before defining it.
#include "pcms/field/field_data.hpp"
#ifndef PCMS_FIELD_DETAIL_FIELD_DATA_HPP
#define PCMS_FIELD_DETAIL_FIELD_DATA_HPP

#include <span>

namespace pcms::detail
{

// The FieldData coefficient operations, for the coefficient algebra of
// pcms/field/detail/field_ops.hpp, which checks the operands.

/// @brief z = a x + b y
template <typename T>
void LinearSum(T a, const FieldData<T>& x, T b, const FieldData<T>& y,
               FieldData<T>& z)
{
  z.LinearSum(a, x, b, y);
}

/// @brief z = c x
template <typename T>
void Scale(T c, const FieldData<T>& x, FieldData<T>& z)
{
  z.Scale(c, x);
}

/// @brief z = c
template <typename T>
void Const(T c, FieldData<T>& z)
{
  z.Const(c);
}

/// @brief z = x / y, entry by entry
template <typename T>
void Div(const FieldData<T>& x, const FieldData<T>& y, FieldData<T>& z)
{
  z.Div(x, y);
}

/// @return the sum of x_i y_i over the coefficients this rank owns
template <typename T>
[[nodiscard]] T DotProdLocal(const FieldData<T>& x, const FieldData<T>& y)
{
  return x.DotProdLocal(y);
}

/// @return the smallest coefficient this rank owns
template <typename T>
[[nodiscard]] T MinLocal(const FieldData<T>& x)
{
  return x.MinLocal();
}

/// @return the largest coefficient magnitude this rank owns
template <typename T>
[[nodiscard]] T MaxNormLocal(const FieldData<T>& x)
{
  return x.MaxNormLocal();
}

/// @brief Replaces each value by its sum over the ranks `data` is
///        partitioned over. Collective.
template <typename T>
void SumOverRanks(const FieldData<T>& data, std::span<T> values)
{
  data.SumOverRanks(values);
}

/// @brief Replaces each value by its maximum over the ranks `data` is
///        partitioned over. Collective.
template <typename T>
void MaxOverRanks(const FieldData<T>& data, std::span<T> values)
{
  data.MaxOverRanks(values);
}

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_FIELD_DATA_HPP
