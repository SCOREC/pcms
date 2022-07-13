#ifndef WDM_COUPLING_XGCPROXY_H
#define WDM_COUPLING_XGCPROXY_H
#include <Omega_h_mesh.hpp>

namespace wdmcpl
{
class XGCProxy
{
public:
  XGCProxy(Omega_h::Mesh mesh);
  void poisson_solve_axisymmetric(){};
  void poisson_solve_non_axisymmetric(){};
  void push(){};

private:
  Omega_h::Mesh mesh_;
  // in real xgc field data is not stored in omega_h mesh, and is replicated
  // here by storing field data in arrays
  Omega_h::Reals potential_;
  Omega_h::Reals density_;
};
} // namespace wdmcpl

#endif // WDM_COUPLING_XGCPROXY_H
