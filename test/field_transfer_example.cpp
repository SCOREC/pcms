#include <Omega_h_build.hpp>
#include <Omega_h_vtk.hpp>
#include <wdmcpl/transfer_field.h>
#include <wdmcpl/omega_h_field.h>

using wdmcpl::Real;
using Field = wdmcpl::OmegaHField<Real,Real>;
using wdmcpl::set;
using wdmcpl::get_nodal_coordinates;
using wdmcpl::get_nodal_data;
using wdmcpl::make_array_view;
using wdmcpl::interpolate_field;
using wdmcpl::copy_field;

struct MeanCombiner {
  void operator()(std::vector<std::reference_wrapper<Field>> fields, Field& combined) const {
    auto field_size = combined.Size();
    Omega_h::Write<Real> combined_array(field_size);
    for( auto& field : fields) {
      assert(field.get().Size() == field_size);
      auto field_array = get_nodal_data(field.get());
      Omega_h::parallel_for(field_size, OMEGA_H_LAMBDA(int i){
                                         combined_array[i] += field_array[i];
                                        });
    }
    auto num_fields = fields.size();
    Omega_h::parallel_for(field_size, OMEGA_H_LAMBDA(int i){
      combined_array[i] /= num_fields;
    });
    set(combined, make_array_view(Omega_h::Read(combined_array)));
  }
};

int main(int argc, char** argv) {

  auto lib = Omega_h::Library{&argc, &argv};
  auto world = lib.world();
  auto internal_mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 40, 40, 0, false);
  auto app_mesh = Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);
  assert(internal_mesh.dim() == 2 && app_mesh.dim() == 2);
  Field app_a_field("app_a_field",app_mesh);
  Field app_b_field("app_b_field",app_mesh);
  Field internal_app_a_field("internal_app_a_field",internal_mesh);
  Field internal_app_b_field("internal_app_b_field",internal_mesh);
  Field internal_combined("internal_combined",internal_mesh);
  auto app_nodal_coords = get_nodal_coordinates(app_a_field);
  Omega_h::Write<double> values_a(app_a_field.Size());
  Omega_h::Write<double> values_b(app_b_field.Size());
  assert(app_a_field.Size() == app_b_field.Size());
  assert(app_a_field.Size() == app_nodal_coords.size()/2);
  for(int i=0; i<app_nodal_coords.size()/2; ++i) {
    auto x = app_nodal_coords[i*2]-0.5;
    auto y = app_nodal_coords[i*2+1]-0.5;
    auto r = sqrt(x*x+y*y);
    r = r>1E-3?r:1E-3;
    std::cout<<r<<"\n";
    //values[i] = 1/r;
    values_b[i] = sin(x)*sin(y);
    values_a[i] = cos(x)*cos(y);
  }
  //auto read_values = Omega_h::Read(values);
  set(app_a_field, make_array_view(Omega_h::Read(values_a)));
  set(app_b_field, make_array_view(Omega_h::Read(values_b)));
  //copy_field(app_a_field, app_b_field);

  //set(app_b_field, make_array_view(read_values));
  interpolate_field(app_a_field,internal_app_a_field, wdmcpl::Lagrange<1>{});
  interpolate_field(app_b_field,internal_app_b_field, wdmcpl::NearestNeighbor{});

  // combine after interpolation
  MeanCombiner combiner{};
  std::vector<std::reference_wrapper<Field>> internal_fields{internal_app_a_field, internal_app_b_field};
  combiner(internal_fields, internal_combined);

  Omega_h::vtk::write_parallel("internal_mesh.vtk", &internal_mesh, internal_mesh.dim());
  Omega_h::vtk::write_parallel("app_mesh.vtk", &app_mesh, app_mesh.dim());
  return 0;
}