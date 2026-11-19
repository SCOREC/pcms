#include "pcms/utility/print.hpp"
#include "pcms/utility/detail/print.hpp"

#include <cstdio>
#include <cassert>

namespace pcms
{

FILE* pcms_stdout = stdout;
FILE* pcms_stderr = stderr;

namespace detail
{
FILE* getStdout()
{
  return pcms_stdout;
}
FILE* getStderr()
{
  return pcms_stderr;
}
} // namespace detail

void setStdout(FILE* out)
{
  assert(out != nullptr);
  pcms_stdout = out;
}

void setStderr(FILE* err)
{
  assert(err != nullptr);
  pcms_stderr = err;
}
} // namespace pcms
