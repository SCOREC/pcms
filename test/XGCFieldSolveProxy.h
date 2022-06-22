#ifndef WDM_COUPLING_XGCFIELDSOLVEPROXY_H
#define WDM_COUPLING_XGCFIELDSOLVEPROXY_H
#include <iostream>
#include "wdmcpl/wdmcpl_static.h"
class XGCFieldSolveProxy
{
public:
  void operator()(wdmcpl::Run, int timestep, int rkstage)
  {
      std::cout << "Running XGC FieldSolveProxy\n";
  }
  void operator()(wdmcpl::NativeToInternal){std::cout<<"FieldSolve To Internal\n";}
  void operator()(wdmcpl::InternalToNative){std::cout<<"FieldSolve To Native\n";}
};

#endif // WDM_COUPLING_XGCFIELDSOLVEPROXY_H
