#include "XGCProxy.h"
namespace wdmcpl {

XGCProxy::XGCProxy(Omega_h::Mesh mesh) : mesh_(std::move(mesh)) {
  mesh_.add_tag<Omega_h::Real>(0,"density",1);
  mesh_.add_tag<Omega_h::Real>(0,"electrostatic potential",1);
}
}
