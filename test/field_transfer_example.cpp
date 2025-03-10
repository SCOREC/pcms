#include <Omega_h_build.hpp>
#include <Omega_h_vtk.hpp>
#include <pcms/transfer_field.h>
#include "pcms/adapter/omega_h/omega_h_field.h"
// for transfer operation dummy test
#include <pcms.h>
#include <chrono>

using pcms::Real;
using OHField = pcms::OmegaHField<Real, Real>;
using OHShim = pcms::OmegaHFieldAdapter<Real, Real>;
using pcms::copy_field;
using pcms::get_nodal_coordinates;
using pcms::get_nodal_data;
using pcms::interpolate_field;
using pcms::make_array_view;
using pcms::set_nodal_data;

inline constexpr int num_trials = 1000;

struct MeanCombiner
{
  void operator()(const std::vector<std::reference_wrapper<OHField>>& fields,
                  OHField& combined) const
  {
    auto field_size = combined.Size();
    Omega_h::Write<Real> combined_array(field_size);
    for (auto& field : fields) {
      assert(field.get().Size() == field_size);
      auto field_array = get_nodal_data(field.get());
      Omega_h::parallel_for(
        field_size,
        OMEGA_H_LAMBDA(int i) { combined_array[i] += field_array[i]; });
    }
    auto num_fields = fields.size();
    Omega_h::parallel_for(
      field_size, OMEGA_H_LAMBDA(int i) { combined_array[i] = combined_array[i] / num_fields; });
    set_nodal_data(combined, make_array_view(Omega_h::Read(combined_array)));
  }
};
void SetApplicationFields(const OHField& app_a_field,
                          const OHField& app_b_field);

using pcms::CoupledField;
using pcms::ProcessType;
using pcms::FieldCommunicator;
using pcms::InternalField;

void test_standalone(Omega_h::Mesh& internal_mesh, Omega_h::Mesh& app_mesh)
{
  OHField app_a_field("app_a_field", app_mesh);
  OHField app_b_field("app_b_field", app_mesh);
  SetApplicationFields(app_a_field, app_b_field);

  OHField internal_app_a_field("internal_app_a_field", internal_mesh);
  OHField internal_app_b_field("internal_app_b_field", internal_mesh);
  OHField internal_combined("internal_combined", internal_mesh);

  for (int i = 0; i < num_trials; ++i) {
    interpolate_field(app_a_field, internal_app_a_field, pcms::Lagrange<1>{});
    interpolate_field(app_b_field, internal_app_b_field,
                      pcms::NearestNeighbor{});

    // combine after interpolation
    MeanCombiner combiner{};
    std::vector<std::reference_wrapper<OHField>> internal_fields{
      internal_app_a_field, internal_app_b_field};

    combiner(internal_fields, internal_combined);
  }
}
void SetApplicationFields(const OHField& app_a_field,
                          const OHField& app_b_field)
{
  auto app_nodal_coords = get_nodal_coordinates(app_a_field);
  Omega_h::Write<double> values_a(app_a_field.Size());
  Omega_h::Write<double> values_b(app_b_field.Size());
  assert(app_a_field.Size() == app_b_field.Size());
  assert(app_a_field.Size() == app_nodal_coords.size() / 2);
  for (int i = 0; i < app_nodal_coords.size() / 2; ++i) {
    auto x = (app_nodal_coords[i * 2] - 0.5) * 3.1415 * 2;
    auto y = (app_nodal_coords[i * 2 + 1] - 0.5) * 3.1415 * 2;
    values_b[i] = sin(x) * sin(y);
    values_a[i] = cos(x) * cos(y);
  }
  set_nodal_data(app_a_field, make_array_view(Omega_h::Read(values_a)));
  set_nodal_data(app_b_field, make_array_view(Omega_h::Read(values_b)));
}

int main(int argc, char** argv)
{

  auto lib = Omega_h::Library{&argc, &argv};
  auto world = lib.world();
  auto internal_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 40, 40, 0, false);
  auto app_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);

  assert(internal_mesh.dim() == 2 && app_mesh.dim() == 2);
  auto point1 = std::chrono::steady_clock::now();
  test_standalone(internal_mesh, app_mesh);
  auto point2 = std::chrono::steady_clock::now();

  Omega_h::vtk::write_parallel("internal_mesh.vtk", &internal_mesh,
                               internal_mesh.dim());
  Omega_h::vtk::write_parallel("app_mesh.vtk", &app_mesh, app_mesh.dim());
  auto point3 = std::chrono::steady_clock::now();
  std::cout << "Test Standalone: "
            << std::chrono::duration<double>(point2 - point1).count() << "\n";
  std::cout << "Write Files: "
            << std::chrono::duration<double>(point3 - point2).count() << "\n";

  return 0;
}
