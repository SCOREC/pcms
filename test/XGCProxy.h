#ifndef WDM_COUPLING_XGCPROXY_H
#define WDM_COUPLING_XGCPROXY_H
#include <iostream>
#include "wdmcpl/wdmcpl_static.h"
class XGCProxy
{
public:
  void operator()(wdmcpl::Run, int timestep, int rkstage) {
      std::cout<<"Running XGC PIC proxy\n";
  }
  void operator()(wdmcpl::NativeToInternal){std::cout<<"XGC To Internal\n";}
  void operator()(wdmcpl::InternalToNative){std::cout<<"XGC From Internal\n";}
};

#endif // WDM_COUPLING_XGCPROXY_H
