#include "pcms/transfer/load_vector_integrator.hpp"
#include "pcms/utility/assert.h"

namespace pcms
{
// TODO:: create a function that can sample sobol sequences instead of
// reading from file
Kokkos::View<MeshField::Real* [3]> read_sobol_barycentric_samples_from_file(
  std::string file_path)
{
  std::ifstream infile(file_path);
  if (!infile) {
    throw std::runtime_error("Could not open sobol sample file");
  }

  std::vector<MeshField::Real> buffer;
  size_t nrows = 0;
  const size_t ncols = 3;
  std::string line;

  // skip header row
  std::getline(infile, line);

  // Read the file and store values in a flat vector
  while (std::getline(infile, line)) {
    if (line.empty())
      continue;

    std::istringstream iss(line);
    MeshField::Real b0, b1, b2;

    if (!(iss >> b0 >> b1 >> b2)) {
      throw std::runtime_error("Invalid barycentric row in file: " + file_path);
    }

    buffer.push_back(b0);
    buffer.push_back(b1);
    buffer.push_back(b2);
    ++nrows;
  }

  Kokkos::View<MeshField::Real* [3]> device_samples("device_data", nrows);

  auto host_samples = Kokkos::create_mirror_view(device_samples);

  // Fill host view from buffer
  for (size_t i = 0; i < nrows; ++i) {
    for (size_t j = 0; j < ncols; ++j) {
      host_samples(i, j) = buffer[i * ncols + j];
    }
  }

  Kokkos::deep_copy(device_samples, host_samples);

  return device_samples;
}

Kokkos::View<MeshField::Real* [3]> generate_uniform_random_barycentric_coords(
  const int npoints_each_tri)
{

  Kokkos::Random_XorShift64_Pool<> rand_pool(1200);
  Kokkos::View<MeshField::Real* [3]> barycentric_coords(
    "uniform random barycentric coordinates", npoints_each_tri);
  Omega_h::parallel_for(
    npoints_each_tri, OMEGA_H_LAMBDA(int i) {
      // taken from  Shape Distributions (ACM Transactions on Graphics, Vol.
      // 21, No. 4, October 2002.) page 814 Eq 1
      auto rand_gen = rand_pool.get_state();
      Omega_h::Real r1 = rand_gen.drand();
      Omega_h::Real r2 = rand_gen.drand();
      Omega_h::Real sqrt_r1 = Kokkos::sqrt(r1);
      barycentric_coords(i, 0) = 1.0 - sqrt_r1;
      barycentric_coords(i, 1) = sqrt_r1 * (1.0 - r2);
      barycentric_coords(i, 2) = sqrt_r1 * r2;
      rand_pool.free_state(rand_gen);
    });

  return barycentric_coords;
}

Kokkos::View<pcms::Real* [2]> global_coords_from_ref_barycentric_coords(
  Omega_h::Mesh& target_mesh,
  const Kokkos::View<MeshField::Real* [3]>& ref_barycentric_coords)
{

  int nelems = target_mesh.nelems();

  const auto& faces2nodes =
    target_mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  const auto& coordinates = target_mesh.coords();
  const int npoints_each_tri = ref_barycentric_coords.extent(0);

  Kokkos::View<pcms::Real* [2]> global_coords(
    "stores global coordinates of  reference samples in each element",
    nelems * npoints_each_tri);

  Omega_h::parallel_for(
    nelems, OMEGA_H_LAMBDA(int id) {
      const auto elem_verts = Omega_h::gather_verts<3>(faces2nodes, id);
      const Omega_h::Few<Omega_h::Vector<2>, 3> verts_coords =
        Omega_h::gather_vectors<3, 2>(coordinates, elem_verts);

      int base_idx = id * npoints_each_tri;
      for (int i = 0; i < npoints_each_tri; ++i) {

        global_coords(base_idx + i, 0) =
          verts_coords[0][0] * ref_barycentric_coords(i, 0) +
          verts_coords[1][0] * ref_barycentric_coords(i, 1) +
          verts_coords[2][0] * ref_barycentric_coords(i, 2);

        global_coords(base_idx + i, 1) =
          verts_coords[0][1] * ref_barycentric_coords(i, 0) +
          verts_coords[1][1] * ref_barycentric_coords(i, 1) +
          verts_coords[2][1] * ref_barycentric_coords(i, 2);
      }
    });

  return global_coords;
}

Kokkos::View<pcms::GridPointSearch2D::Result*> localize_points_in_mesh(
  Omega_h::Mesh& mesh, const Kokkos::View<pcms::Real* [2]>& points)
{

  PCMS_ALWAYS_ASSERT(mesh.dim() == 2);
  const auto mesh_bbox = Omega_h::get_bounding_box<2>(&mesh);

  const auto npoints = points.extent(0);

  Kokkos::Array<Omega_h::Real, 2> edge_length;

  edge_length = {mesh_bbox.max[0] - mesh_bbox.min[0],
                 mesh_bbox.max[1] - mesh_bbox.min[1]};

  int nx = edge_length[0] * Kokkos::sqrt(mesh.nelems());
  int ny = edge_length[1] * Kokkos::sqrt(mesh.nelems());

  pcms::GridPointSearch2D search_element(mesh, nx, ny);
  auto results = search_element(points);

  return results;
}

Omega_h::Reals evaluate_field_from_point_localization(
  Omega_h::Mesh& mesh, const Omega_h::Reals& nodal_field_values,
  const Kokkos::View<pcms::GridPointSearch2D::Result*>& results)
{

  const auto& npoints = results.extent(0);

  const auto& faces2nodes = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  const auto& coordinates = mesh.coords();

  Omega_h::Write<Omega_h::Real> field_values_at_points(
    npoints, 0.0, "stores field values at the given points");

  Omega_h::parallel_for(
    npoints, OMEGA_H_LAMBDA(int id) {
      field_values_at_points[id] = 0.0;

      const int tri = results[id].element_id;
      if (tri < 0)
        return;
      const auto el_verts = Omega_h::gather_verts<3>(faces2nodes, tri);

      for (int j = 0; j < 3; ++j) {
        const auto node_id = el_verts[j];
        field_values_at_points[id] +=
          nodal_field_values[node_id] * results[id].parametric_coords[j];
      }
    });
  Kokkos::fence();
  return Omega_h::read(field_values_at_points);
}

KOKKOS_INLINE_FUNCTION
Omega_h::Vector<3> compute_elm_load_vector_mc(
  const Kokkos::View<MeshField::Real* [3]>& ref_barycentric_coords,
  const Omega_h::Reals& field_values_at_points,
  const Omega_h::Real& elmVolume, const int elm)
{

  const int npoints_each_tri = ref_barycentric_coords.extent(0);
  Omega_h::Real sum0 = 0.0;
  Omega_h::Real sum1 = 0.0;
  Omega_h::Real sum2 = 0.0;
  const int base = elm * npoints_each_tri;
  for (int i = 0; i < npoints_each_tri; ++i) {
    const auto f = field_values_at_points[base + i];
    sum0 += ref_barycentric_coords(i, 0) * f;
    sum1 += ref_barycentric_coords(i, 1) * f;
    sum2 += ref_barycentric_coords(i, 2) * f;
  }
  return elmVolume / npoints_each_tri * Omega_h::Vector<3>{sum0, sum1, sum2};
}

Kokkos::View<MeshField::Real*> buildLoadVectorMC(
  Omega_h::Mesh& target_mesh, const Omega_h::Reals& field_values_at_points,
  const int npoints_each_tri, SamplingMethod method,
  const std::string& sobol_filename)
{

  Kokkos::View<MeshField::Real* [3]> ref_barycentric_coords;
  if (method == SamplingMethod::SOBOL) {
    if (sobol_filename.empty()) {
                throw std::runtime_error("Could not open sobol sample file");
    }
    ref_barycentric_coords =
      read_sobol_barycentric_samples_from_file(sobol_filename);
  } else {
    ref_barycentric_coords =
      generate_uniform_random_barycentric_coords(npoints_each_tri);
  }

  if (ref_barycentric_coords.extent(0) !=
      static_cast<std::size_t>(npoints_each_tri)) {
    throw std::runtime_error(
      "Sobol sample count does not match npoints_each_tri.");
  }

  int subVectorSize = 3;

  Omega_h::Reals elementsArea;
  elementsArea = Omega_h::measure_elements_real(&target_mesh);

  Kokkos::View<MeshField::Real*> elmLoadVector("elmLoadVector",
                                               target_mesh.nelems() * 3);

  Kokkos::parallel_for(
    "eval load vector using MC", target_mesh.nelems(),
    KOKKOS_LAMBDA(const int elm) {
      auto result = compute_elm_load_vector_mc(
        ref_barycentric_coords, field_values_at_points, elementsArea[elm], elm);

      for (int i = 0; i < 3; ++i) {
        elmLoadVector(elm * subVectorSize + i) = result[i];
      }
    });

  return elmLoadVector;
}
} // namespace pcms
