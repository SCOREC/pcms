#include "wdmcpl/assert.h"
#include <cstdio>
#include <cstdlib>
namespace wdmcpl {
void Wdmcpl_Assert_Fail(const char* msg) {
  fprintf(stderr, "%s", msg);
  abort();
}
}
