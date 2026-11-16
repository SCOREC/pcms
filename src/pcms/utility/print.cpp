#include "pcms/utility/print.h"

#include <cstdio>
#include <cassert>

namespace pcms
{

FILE* pcms_stdout = stdout;
FILE* pcms_stderr = stderr;

FILE* getStdout()
{
  return pcms_stdout;
}
FILE* getStderr()
{
  return pcms_stderr;
}

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
