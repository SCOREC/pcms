#ifndef PCMS_UTILITY_DETAIL_ASSERT_HPP
#define PCMS_UTILITY_DETAIL_ASSERT_HPP

namespace pcms::detail
{

// from scorec/core/pcu_fail.h
void Pcms_Assert_Fail(const char* msg) __attribute__((noreturn));

/// Prints "<cond> failed at <file> + <line> on rank <rank>[: <msg>]" and
/// aborts; the rank is the one in MPI_COMM_WORLD, or -1 if MPI is not
/// initialized.
[[noreturn]] void AssertFail(const char* cond, const char* file, int line,
                             const char* msg);
/// As AssertFail, with the rank given.
[[noreturn]] void AssertFailOnRank(const char* cond, const char* file, int line,
                                   int rank, const char* msg);

} // namespace pcms::detail

#endif // PCMS_UTILITY_DETAIL_ASSERT_HPP
