#include "pcms/assert.h"
#include <cstdio>
#include <cstdlib>
namespace pcms {
void Pcms_Assert_Fail(const char* msg) {
  fprintf(stderr, "%s", msg);
  abort();
}
}
