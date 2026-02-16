#include "pcms/utility/print.h"

namespace pcms
{
void Pcms_Assert_Fail(const char* msg)
{
  printError("%s", msg);
  abort();
}
} // namespace pcms
