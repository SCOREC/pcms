#ifndef WDM_COUPLING_GENEPROXY_H
#define WDM_COUPLING_GENEPROXY_H
#include <iostream>
#include "wdmcpl/wdmcpl_static.h"

class GeneProxy
{
public:
  void operator()(wdmcpl::Run, int timestep, int rkstage)
  {
      std::cout<<"Running GENE PROXY\n";
  }
  void operator()(wdmcpl::NativeToInternal){
    std::cout<<"GENE TO INTERNAL\n";
  }
  void operator()(wdmcpl::InternalToNative){
    std::cout<<"GENE FROM INTERNAL\n";
  }

private:

};

#endif // WDM_COUPLING_GENEPROXY_H
