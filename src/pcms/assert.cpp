#include "pcms/assert.h"
#include "pcms/print.h"
#include <cstdio>
#include <cstdlib>
namespace pcms
{
void Pcms_Assert_Fail(const char* msg)
{
  printError("%s", msg);
  abort();
}
} // namespace pcms
