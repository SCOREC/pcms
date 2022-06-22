//#include "GeneProxy.h"
//#include "XGCProxy.h"
//#include "XGCFieldSolveProxy.h"
#include "wdmcpl/wdmcpl_static.h"
#include "wdmcpl/interpolation.h"
#include <catch2/catch.hpp>

#include <Omega_h_library.hpp>
#include <Omega_h_array_ops.hpp>
#include <Omega_h_comm.hpp>
#include <Omega_h_mesh.hpp>

using namespace wdmcpl;

template <typename MeshT, typename... CoupledFields>
class GenericProxyApp
{
public:
  static_assert((... && CoupledFieldC<CoupledFields>));
  constexpr GenericProxyApp(std::string_view name, MeshT mesh,
                            CoupledFields... coupled_fields)
    : name_(name), mesh_(mesh), coupled_fields_(std::move(coupled_fields)...)
  {}
  void operator()(wdmcpl::Run, int timestep, int rkstage)
  {
    std::cout << "Running " << name_ << " Proxy\n";
  }
  void operator()(wdmcpl::NativeToInternal)
  {
    std::cout << name_ << " NativeToInternal\n";
    tuple_call_each(coupled_fields_, wdmcpl::NativeToInternal{});
  }
  void operator()(wdmcpl::InternalToNative)
  {
    std::cout << name_ << " InternalToNative\n";
    tuple_call_each(coupled_fields_, wdmcpl::InternalToNative{});
  }

private:
  std::string name_;
  MeshT mesh_;
  std::tuple<CoupledFields...> coupled_fields_;
};

struct OverlapCombiner
{
  // TODO generic OverlapCombiner. Current thoughts are that combiner can take a
  // (compile time) list of indexes that correspond to the fields in each solver
  // in the list of solvers. Negative index means that field is not involved in
  // the combination. The operator also needs to take a variadic list of
  // coupling codes. So, in the case of Gene/XGC if both fields are stored with
  // the same index. Then, we would have <0,0> (Gene/XGC) Density overlaps and
  // <1,1> (Gene/XGC) current density overlaps. Or, we assume ordering of fields
  // is same in Application functor.
  //
  // NOTE: Combiner operates on the internal data structures, so we want to
  // avoid making the "user" i.e. someone adding a new coupling code from
  // writing this hence why the generic version is needed
  OverlapCombiner(double c1, double c2)
  {
    inside_coordinate_ = c1;
    outside_coordinate_ = c2;
    if (c2 < c1) {
      std::swap(inside_coordinate_, outside_coordinate_);
    }
  }
  // loop over relevant fields and combine them
  template <typename T1, typename T2>
  void operator()(T1&, T2&)
  {
    std::cout << "Overlap Combiner (" << inside_coordinate_ << ","
              << outside_coordinate_ << ") (needs impl)\n";
  };

private:
  double inside_coordinate_;
  double outside_coordinate_;
};

struct None
{
  template <typename... Args>
  constexpr void operator()(Args...)
  {}
};

// can we write the field transfer to only need 2
// 1) structured mesh using coordinate transforms for application specific stuff
// 2) unstructured mesh using coordinate transforms for application specific
// stuff what are the steps of the field transfer? from native, to internal loop
// over each of the nodes in the internal mesh to get the positions in new
// coordinate system loop over each
template <int order = 1>
struct LagrangeFieldTransfer
{
  void operator()()
  {
    // if mesh doesn't have
    // for point in original mesh
    double* old_coords;
    double* old_data;
    double new_coord;
    auto new_value =
      lagrange_interpolation<order>(old_coords, old_data, new_coord);
  }
  template <typename FromField, typename ToField>
  void operator()(FromField, ToField)
  {}
};
auto gene_proxy(const Omega_h::Mesh* internal_mesh)
{
  std::vector<double> gene_mesh;
  Field gene_native_density{
    GENEFieldFollowing{},
    std::move(gene_mesh),
    std::vector<double>{},
  };
  Field gene_internal_density{Cylindrical{}, internal_mesh,
                              std::vector<double>{}};
  CoupledField gene_density{"Density", std::move(gene_native_density),
                            std::move(gene_internal_density),
                            LagrangeFieldTransfer{}, LagrangeFieldTransfer{}};
  return GenericProxyApp{"GENE", std::move(gene_density)};
}
auto xgc_proxy(const Omega_h::Mesh* internal_mesh)
{
  Field xgc_native_density{Cylindrical{}, internal_mesh, std::vector<double>{}};
  CoupledField xgc_density{"Density", xgc_native_density, None{}, None{},
                           None{}};
  return GenericProxyApp{"XGC", std::move(xgc_density)};
}
auto xgc_field_proxy(const Omega_h::Mesh* internal_mesh)
{
  Field xgc_native_density{Cylindrical{}, internal_mesh, std::vector<double>{}};
  CoupledField xgc_density{"Density", std::move(xgc_native_density), None{}, None{},
                          None{}};
  return GenericProxyApp{"XGC Field", xgc_density};
}

auto proxy_coupler()
{
  // get basic omega_h mesh
  auto lib = Omega_h::Library();
  auto internal_mesh = std::make_unique<Omega_h::Mesh>();

  auto GeneGProxy = gene_proxy(internal_mesh.get());
  auto XGCGProxy = xgc_proxy(internal_mesh.get());
  auto XGCGFieldSolveProxy = xgc_field_proxy(internal_mesh.get());

  ApplicationGroup pic_solve_group(OverlapCombiner{0.6, 0.8}, GeneGProxy,
                                   XGCGProxy);
  ApplicationGroup field_solve_group(None{}, XGCGFieldSolveProxy);
  // charge density & E-Field Coupling
  Coupling gene_xgc_coupling(pic_solve_group, field_solve_group,
                             std::move(internal_mesh));

  return Coupler{std::move(gene_xgc_coupling)};
}

TEST_CASE("run proxy coupler")
{
  auto coupler = proxy_coupler();
  constexpr int nsteps = 1;
  constexpr int rksteps = 1;
  coupler.Run(nsteps, rksteps);

  // for vertex in target_mesh.vertices:
  //    result = source_mesh.evaluate(source_field, method, vertex.coordinate)
  //    target_mesh.set_field(target_field, vertex)
}