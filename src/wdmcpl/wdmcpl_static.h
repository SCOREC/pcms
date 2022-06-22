#ifndef WDM_COUPLING_WDMCPL_STATIC_H
#define WDM_COUPLING_WDMCPL_STATIC_H
#include "wdmcpl/coordinate_transform.h"
#include "wdmcpl/utility.h"
#include <iostream>
namespace wdmcpl
{
// Tag types
struct Run
{};
struct NativeToInternal
{};
struct InternalToNative
{};
// Concepts
template <typename T>
concept Application = std::invocable<T, Run, int, int> &&
  std::invocable<T, NativeToInternal> && std::invocable<T, InternalToNative>;

template<typename T>
concept CoupledFieldC = std::invocable<T, NativeToInternal> && std::invocable<T, InternalToNative>;

template <typename T, typename... Args>
concept Combiner = std::invocable<T, Args&...> &&(...&& Application<Args>);

template <typename FieldT, typename T>
class FieldSender
{
public:
  void operator()() {}
};
template <typename FieldT, typename T>
class FieldReciever
{
public:
  void operator()() {}
};

template <typename FieldFrom, typename FieldTo>
struct FieldTransfer
{
  template <typename Method> // Method can be enum selector rather than structs
  // to put everything in one place
  void operator()(){
    // not general, but for fusion codes we have a
    // two step proces:
    // step 1 "interpolate" to the poloidal plane
    // step 2 "interpolate" to native positions on
    // field_to_ mesh
  };
  const FieldFrom& field_from_;
  FieldTo& field_to_;
};

template <typename NativeField, typename InternalField, typename NTIFieldTransfer, typename ITNFieldTransfer>
struct CoupledField
{
  CoupledField(std::string_view name, NativeField native_field,
               InternalField internal_field,
               NTIFieldTransfer native_to_internal,
               ITNFieldTransfer internal_to_native)
    : name_(name),
      native_to_internal_(std::move(native_to_internal)),
      internal_to_native_(std::move(internal_to_native)),
      native_field_(std::move(native_field)),
      internal_field_(std::move(internal_field))
  {}
  void operator()(NativeToInternal){
    std::cout<<"\tN2I "<<name_<<"\n";
    native_to_internal_(native_field_,internal_field_);
  }
  void operator()(InternalToNative){
    std::cout<<"\tI2N "<<name_<<"\n";
    internal_to_native_(internal_field_,native_field_);
  }
  std::string name_;
  NTIFieldTransfer native_to_internal_;
  ITNFieldTransfer internal_to_native_;
  NativeField native_field_;
  InternalField internal_field_;
};

template <typename CoordinateSystemT, typename FieldStorageT, typename MeshStorageT>
struct Field
{
  Field(CoordinateSystemT, MeshStorageT mesh, FieldStorageT field)
    : mesh_(std::move(mesh)), field_(std::move(field))
  {}
  // coordinate system that field will use
  //using CoordinateSystem = CoordinateSystemT;
  CoordinateSystemT coordinate_system_;
  FieldStorageT field_; // field
  MeshStorageT mesh_; // mesh
};
// application concept needs to own fields
template <typename CombinerT, typename... Applications>
class ApplicationGroup
{
  // static assert here gives much better error message than constraining the
  // class
  static_assert((... && Application<Applications>),
                "All applications must abide by application concept.");
  static_assert(Combiner<CombinerT, Applications...>,
                "Combiner must abide by combiner concept.");

public:
  ApplicationGroup(CombinerT combiner, Applications... applications)
    : applications_(std::move(applications)...), combiner_(std::move(combiner))
  {}

  ////
  template <typename... Args>
  void Run(Args&&... args)
  {
    this->Solve(std::forward<Args>(args)...);
    NativeToInternal(std::forward<Args>(args)...);
    Combine();
  }

  template <typename... Args>
  constexpr void Solve(Args&&... args)
  {
    tuple_call_each(applications_, wdmcpl::Run{}, std::forward<Args>(args)...);
  }
  template <typename... Args>
  constexpr void NativeToInternal(Args&&... args)
  {
    tuple_call_each(applications_, wdmcpl::NativeToInternal{});
  }
  template <typename... Args>
  constexpr void InternalToNative(Args&&... args)
  {
    tuple_call_each(applications_, wdmcpl::InternalToNative{});
  }
  // constexpr void Combine() { std::cout<<"Need to impliment combiner!\n";}
  constexpr auto Combine() { return std::apply(combiner_, applications_); }

private:
  std::tuple<Applications...> applications_;
  CombinerT combiner_;
};

template <typename ApplicationGroupA, typename ApplicationGroupB, typename MeshT>
struct Coupling
{
  Coupling(ApplicationGroupA sga, ApplicationGroupB sgb, MeshT internal_mesh)
    : application_group_a_{std::move(sga)}, application_group_b_{std::move(sgb)}, internal_mesh_{std::move(internal_mesh)}
  {}
  // step the coupling object
  template <typename... Args>
  void operator()(Args&&... args)
  {
    application_group_a_.Run(std::forward<Args>(args)...);
    application_group_b_.InternalToNative(std::forward<Args>(args)...);
    application_group_b_.Run(std::forward<Args>(args)...);
    application_group_a_.InternalToNative(std::forward<Args>(args)...);
  };
  ApplicationGroupA application_group_a_;
  ApplicationGroupB application_group_b_;
  MeshT internal_mesh_;
};

template <typename... Couplings>
class Coupler
{
public:
  Coupler(Couplings... cpl) : couplings_(std::move(cpl)...) {}
  void Run(int nsteps, int nrksteps)
  {
    for (int step = 0; step < nsteps; ++step) {
      for (int rkstage = 0; rkstage < nrksteps; ++rkstage) {
        tuple_call_each(couplings_, step, rkstage);
      }
    }
  }

private:
  std::tuple<Couplings...> couplings_;
};
} // namespace wdmcpl
#endif // WDM_COUPLING_WDMCPL_STATIC_H
