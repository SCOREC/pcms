#include <Omega_h_build.hpp>
#include <Omega_h_vtk.hpp>
#include <wdmcpl/transfer_field.h>
#include <wdmcpl/omega_h_field.h>
// for transfer operation dummy test
#include <wdmcpl.h>
#include <chrono>

using wdmcpl::Real;
using OHField = wdmcpl::OmegaHField<Real, Real>;
using OHShim = wdmcpl::OmegaHFieldShim<Real, Real>;
using wdmcpl::copy_field;
using wdmcpl::get_nodal_coordinates;
using wdmcpl::get_nodal_data;
using wdmcpl::interpolate_field;
using wdmcpl::make_array_view;
using wdmcpl::set;


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
      field_size, OMEGA_H_LAMBDA(int i) { combined_array[i] /= num_fields; });
    set(combined, make_array_view(Omega_h::Read(combined_array)));
  }
};
void SetApplicationFields(const OHField& app_a_field,
                          const OHField& app_b_field);

using wdmcpl::CoupledField;
using wdmcpl::FieldEvaluationMethod;
using wdmcpl::FieldTransferMethod;
using wdmcpl::GatherOperation;
using wdmcpl::TransferOptions;
using wdmcpl::detail::InternalField;

void test_gather_operation(Omega_h::Mesh& internal_mesh,
                           Omega_h::Mesh& app_mesh)
{
  OHField app_a_field("gather_app_a", app_mesh);
  OHField app_b_field("gather_app_b", app_mesh);
  SetApplicationFields(app_a_field, app_b_field);
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  redev::Redev redev{mpi_comm, redev::ClassPtn{}};
  InternalField combined_field(OHField("gather_combined", internal_mesh));
  TransferOptions transfer_lag1{FieldTransferMethod::Interpolate,
                                   FieldEvaluationMethod::Lagrange1};
  TransferOptions transfer_nn{FieldTransferMethod::Interpolate,
                                   FieldEvaluationMethod::NearestNeighbor};
  std::vector<CoupledField> coupled_fields;
  coupled_fields.reserve(2);
  coupled_fields.emplace_back("gather_app_a",
                              OHShim("gather_app_a", app_mesh),
                              wdmcpl::FieldCommunicator(),
                              internal_mesh,
                              transfer_lag1, transfer_lag1);
  coupled_fields.emplace_back("gather_app_b",
                              OHShim("gather_app_b", app_mesh),
                              wdmcpl::FieldCommunicator(),
                              internal_mesh,
                              transfer_nn, transfer_lag1);
  GatherOperation gather(
    coupled_fields, combined_field,
    [](nonstd::span<const std::reference_wrapper<InternalField>> fields,
       InternalField& combined_field) {
      std::vector<std::reference_wrapper<OHField>> typed_fields;
      typed_fields.reserve(fields.size());

      std::transform(fields.begin(), fields.end(),
                     std::back_inserter(typed_fields),
                     [](const std::reference_wrapper<InternalField>& internal) {
                       return std::ref(std::get<OHField>(internal.get()));
                     });
      auto& typed_combined_field = std::get<OHField>(combined_field);
      std::invoke(MeanCombiner{}, typed_fields, typed_combined_field);
    });
  for(int i=0; i<100; ++i) {
    gather.Run();
  }
}

void test_standalone(Omega_h::Mesh& internal_mesh, Omega_h::Mesh& app_mesh)
{
  OHField app_a_field("app_a_field", app_mesh);
  OHField app_b_field("app_b_field", app_mesh);
  SetApplicationFields(app_a_field, app_b_field);

  OHField internal_app_a_field("internal_app_a_field", internal_mesh);
  OHField internal_app_b_field("internal_app_b_field", internal_mesh);
  OHField internal_combined("internal_combined", internal_mesh);
  // copy_field(app_a_field, app_b_field);

  for(int i=0; i<100; ++i) {
  // set(app_b_field, make_array_view(read_values));
  interpolate_field(app_a_field, internal_app_a_field, wdmcpl::Lagrange<1>{});
  interpolate_field(app_b_field, internal_app_b_field,
                    wdmcpl::NearestNeighbor{});

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
  // auto read_values = Omega_h::Read(values);
  set(app_a_field, make_array_view(Omega_h::Read(values_a)));
  set(app_b_field, make_array_view(Omega_h::Read(values_b)));
}

int main(int argc, char** argv)
{

  auto lib = Omega_h::Library{&argc, &argv};
  auto world = lib.world();
  auto internal_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 40, 40, 0, false);
  auto app_mesh =
    Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1, 1, 1, 10, 10, 0, false);

  // we can try this when search structure is working better
  // auto internal_mesh = Omega_h::Mesh{&lib};
  // Omega_h::binary::read("./d3d-full_9k_sfc.osh",lib.world(), &internal_mesh);
  // auto app_mesh = Omega_h::Mesh{&lib};
  // Omega_h::binary::read("./d3d-full_9k_sfc.osh", lib.world(), &app_mesh);
  assert(internal_mesh.dim() == 2 && app_mesh.dim() == 2);
  auto point1 = std::chrono::steady_clock::now();
  test_standalone(internal_mesh, app_mesh);
  auto point2 = std::chrono::steady_clock::now();
  test_gather_operation(internal_mesh, app_mesh);
  auto point3 = std::chrono::steady_clock::now();

  Omega_h::vtk::write_parallel("internal_mesh.vtk", &internal_mesh,
                               internal_mesh.dim());
  Omega_h::vtk::write_parallel("app_mesh.vtk", &app_mesh, app_mesh.dim());
  auto point4 = std::chrono::steady_clock::now();
  std::cout<<"Test Standalone: "<<std::chrono::duration<double>(point2-point1).count()<<"\n";
  std::cout<<"Test Gather: "<<std::chrono::duration<double>(point3-point2).count()<<"\n";
  std::cout<<"Write Files: "<<std::chrono::duration<double>(point4-point3).count()<<"\n";


  return 0;
}