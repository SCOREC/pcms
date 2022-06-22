#ifndef WDM_COUPLING_WDMCPL_DYNAMIC_H
#define WDM_COUPLING_WDMCPL_DYNAMIC_H
#include <vector>
#include <memory>
#include <Omega_h_mesh.hpp>
#include <iostream>
namespace wdmcpl {

class FieldTransfer
{};
class Field {};

struct ICoupledField {

  virtual void NativeToInternal() = 0;
  virtual void InternalToNative() = 0;
};

template <typename NTIFieldTransfer, typename ITNFieldTransfer>
class CoupledField : public ICoupledField
{
  CoupledField(std::string_view name, std::unique_ptr<Field> native_field,
               std::unique_ptr<Field> internal_field,
               NTIFieldTransfer native_to_internal,
               ITNFieldTransfer internal_to_native)
    : name_(name),
      native_to_internal_(std::move(native_to_internal)),
      internal_to_native_(std::move(internal_to_native)),
      native_field_(std::move(native_field)),
      internal_field_(std::move(internal_field))
  {}
public:
  void NativeToInternal() {
    std::cout<<"\tN2I "<<name_<<"\n";
    native_to_internal_(native_field_,internal_field_);
  }
  void InternalToNative() {
    std::cout<<"\tI2N "<<name_<<"\n";
    internal_to_native_(internal_field_,native_field_);
  }
private:
  std::string name_;
  NTIFieldTransfer native_to_internal_;
  ITNFieldTransfer internal_to_native_;
  std::unique_ptr<Field> native_field_;
  std::unique_ptr<Field> internal_field_;
};

/*
class CoupledField
{
public:

  virtual void NativeToInternal() = 0;
  virtual void InternalToNative() = 0;
private:
  std::unique_ptr<Field> native_field_;
  std::unique_ptr<Field> internal_field_;
};*/

class Application{
public:
  virtual void Run(int step, int rkstep) = 0;
  void NativeToInternal() {
    for(auto& coupled_field: coupled_fields_) {
      coupled_field->NativeToInternal();
    }
  }
  void InternalToNative() {
    for(auto& coupled_field: coupled_fields_) {
      coupled_field->InternalToNative();
    }
  }
private:
  std::vector<std::unique_ptr<ICoupledField>> coupled_fields_;
};

class ApplicationGroup{
    public:
      void Run(int step, int rkstep) {
        for (auto& app : applications_) {
          app->Run(step,rkstep);
        }
        NativeToInternal();
        Combine();
      }
    private:
      void NativeToInternal(){
      }
      void Combine() {}
      std::vector<std::unique_ptr<Application>> applications_;
      std::unique_ptr<Omega_h::Mesh> internal_mesh_;

  };

  class Coupling{
    public:
      void Run(int step, int rkstep) {
      }

    private:
      std::unique_ptr<ApplicationGroup> application_group_a_;
      std::unique_ptr<ApplicationGroup> application_group_b_;

  };
  class Coupler {
    public:
      void Run(int nsteps, int nrksteps){
        for(int step=0; step<nsteps; ++step)
        {
          for (int rkstep=0; rkstep<nrksteps; ++rkstep) {
            for(const auto& coupling : couplings_) {
              coupling->Run(step,rkstep);
            }
          }
        }
      }
    private:
      std::vector<std::unique_ptr<Coupling>> couplings_;
  };

}

#endif // WDM_COUPLING_WDMCPL_DYNAMIC_H
