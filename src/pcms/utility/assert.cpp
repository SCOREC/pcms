#include "pcms/utility/print.hpp"

namespace pcms
{
void Pcms_Assert_Fail(const char* msg)
{
  printError("%s", msg);
  abort();
}
} // namespace pcms
