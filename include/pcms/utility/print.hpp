#ifndef PCMS_UTILITY_PRINT_HPP
#define PCMS_UTILITY_PRINT_HPP
#include <cstdio>

namespace pcms
{

/// Redirects PCMS's informational output (stdout by default).
void setStdout(FILE* out);
/// Redirects PCMS's error output (stderr by default).
void setStderr(FILE* err);

} // namespace pcms

#endif // PCMS_UTILITY_PRINT_HPP
