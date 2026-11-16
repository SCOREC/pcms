#ifndef PCMS_MPI_TYPE_H
#define PCMS_MPI_TYPE_H

#include "pcms/utility/types.h"
#include <complex>
#include <mpi.h>

namespace pcms
{

template <class T>
[[nodiscard]] constexpr MPI_Datatype GetMPIType(T) noexcept
{
  if constexpr (std::is_same_v<T, char>) {
    return MPI_CHAR;
  } else if constexpr (std::is_same_v<T, signed short int>) {
    return MPI_SHORT;
  } else if constexpr (std::is_same_v<T, signed int>) {
    return MPI_INT;
  } else if constexpr (std::is_same_v<T, signed long>) {
    return MPI_LONG;
  } else if constexpr (std::is_same_v<T, signed long long>) {
    return MPI_LONG_LONG;
  } else if constexpr (std::is_same_v<T, signed char>) {
    return MPI_SIGNED_CHAR;
  } else if constexpr (std::is_same_v<T, unsigned char>) {
    return MPI_UNSIGNED_CHAR;
  } else if constexpr (std::is_same_v<T, unsigned short>) {
    return MPI_UNSIGNED_SHORT;
  } else if constexpr (std::is_same_v<T, unsigned int>) {
    return MPI_UNSIGNED;
  } else if constexpr (std::is_same_v<T, unsigned long>) {
    return MPI_UNSIGNED_LONG;
  } else if constexpr (std::is_same_v<T, unsigned long long>) {
    return MPI_UNSIGNED_LONG_LONG;
  } else if constexpr (std::is_same_v<T, float>) {
    return MPI_FLOAT;
  } else if constexpr (std::is_same_v<T, double>) {
    return MPI_DOUBLE;
  } else if constexpr (std::is_same_v<T, long double>) {
    return MPI_LONG_DOUBLE;
  } else if constexpr (std::is_same_v<T, wchar_t>) {
    return MPI_WCHAR;
  } else if constexpr (std::is_same_v<T, int8_t>) {
    return MPI_INT8_T;
  } else if constexpr (std::is_same_v<T, int16_t>) {
    return MPI_INT16_T;
  } else if constexpr (std::is_same_v<T, int32_t>) {
    return MPI_INT32_T;
  } else if constexpr (std::is_same_v<T, int64_t>) {
    return MPI_INT64_T;
  } else if constexpr (std::is_same_v<T, uint8_t>) {
    return MPI_UINT8_T;
  } else if constexpr (std::is_same_v<T, uint16_t>) {
    return MPI_UINT16_T;
  } else if constexpr (std::is_same_v<T, uint32_t>) {
    return MPI_UINT32_T;
  } else if constexpr (std::is_same_v<T, uint64_t>) {
    return MPI_UINT64_T;
  } else if constexpr (std::is_same_v<T, bool>) {
    return MPI_CXX_BOOL;
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    return MPI_CXX_FLOAT_COMPLEX;
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    return MPI_CXX_DOUBLE_COMPLEX;
  } else if constexpr (std::is_same_v<T, std::complex<long double>>) {
    return MPI_CXX_LONG_DOUBLE_COMPLEX;
  } else {
    static_assert(detail::dependent_always_false<T>::value,
                  "type has unknown map to MPI_Datatype");
    return {};
  }
}

} // namespace pcms

#endif // PCMS_MPI_TYPE_H
