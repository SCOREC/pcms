#include "coupling.h"

namespace coupler {


void close_engines(adios2::Engine engine[], const int num) {
  for(int i = 0; i < num; i++) {
    engine[i].Close();
  }
}

}
