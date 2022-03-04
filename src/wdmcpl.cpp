#include "wdmcpl.h"
#include <redev.h>

namespace wdmcpl {
  void init(int argc, char** argv, MPI_Comm comm) {
    redev::RCBPtn ptn;
    auto isRendezvous = true;
    auto noParticipant = true;
    redev::Redev(comm,ptn,isRendezvous,noParticipant);
  }
}
