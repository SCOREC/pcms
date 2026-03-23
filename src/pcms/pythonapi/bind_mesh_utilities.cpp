#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Omega_h_mesh.hpp>
#include <Omega_h_array_ops.hpp>
#include <pcms/interpolator/mls_interpolation.hpp>
#include <pcms/interpolator/adj_search.hpp>
#include "numpy_array_transform.h"

namespace py = pybind11;

namespace pcms
{

// Helper function to compute entity centroids
py::array_t<double> compute_entity_centroids(Omega_h::Mesh& mesh,
                                             Omega_h::Int entity_dim)
{
  const auto dim = mesh.dim();
  const auto nents = mesh.nents(entity_dim);
  const auto ent2verts = mesh.ask_down(entity_dim, Omega_h::VERT).ab2b;
  const auto coords = mesh.coords();
  Omega_h::Write<Omega_h::Real> centroids(dim * nents, 0.0, "entity centroids");

  if (dim == 2 && entity_dim == 2) {
    Kokkos::parallel_for(
      "entity_centroids_tri2d", nents, OMEGA_H_LAMBDA(const Omega_h::LO id) {
        const auto current_el_verts = Omega_h::gather_verts<3>(ent2verts, id);
        const Omega_h::Few<Omega_h::Vector<2>, 3> current_el_vert_coords =
          Omega_h::gather_vectors<3, 2>(coords, current_el_verts);
        auto centroid = Omega_h::average(current_el_vert_coords);
        int index = 2 * id;
        centroids[index] = centroid[0];
        centroids[index + 1] = centroid[1];
      });
  } else {
    throw std::runtime_error(
      "Centroid computation not implemented for this dimension combination");
  }

  return omega_h_read_to_numpy(Omega_h::read(centroids));
}

// Helper function to determine patch size based on interpolation degree
Omega_h::LO patch_size_for_degree(Omega_h::LO interp_degree)
{
  if (interp_degree == 1)
    return 4;
  if (interp_degree == 2)
    return 10;
  if (interp_degree == 3)
    return 12;
  return 12;
}

// Map entity field to vertices using simple averaging
py::array_t<double> map_entity_field_to_vertices_average(
  Omega_h::Mesh& mesh, py::array_t<double> entity_field_values,
  Omega_h::Int entity_dim)
{

  const Omega_h::LO nverts = mesh.nverts();
  const Omega_h::LO nentities = mesh.nents(entity_dim);

  // Get entity to vertices adjacency
  auto entity_to_verts_omega_h = mesh.ask_verts_of(entity_dim);
  const Omega_h::LO verts_per_entity =
    entity_to_verts_omega_h.size() / nentities;

  // Convert numpy array to Omega_h array
  auto entity_field_omega_h =
    numpy_to_omega_h_read<Omega_h::Real>(entity_field_values);

  // Initialize vertex field and count arrays
  Omega_h::Write<Omega_h::Real> vertex_field_values(nverts, 0.0);
  Omega_h::Write<Omega_h::Real> vertex_counts(nverts, 0.0);

  // Parallel accumulation
  Kokkos::parallel_for(
    "accumulate_to_vertices", nentities,
    OMEGA_H_LAMBDA(const Omega_h::LO entity_id) {
      const Omega_h::Real entity_value = entity_field_omega_h[entity_id];
      for (Omega_h::LO v_local = 0; v_local < verts_per_entity; ++v_local) {
        const Omega_h::LO vert_id =
          entity_to_verts_omega_h[entity_id * verts_per_entity + v_local];
        Kokkos::atomic_add(&vertex_field_values[vert_id], entity_value);
        Kokkos::atomic_add(&vertex_counts[vert_id], 1.0);
      }
    });

  // Parallel averaging
  Kokkos::parallel_for(
    "average_vertices", nverts, OMEGA_H_LAMBDA(const Omega_h::LO vert_id) {
      if (vertex_counts[vert_id] > 0.0) {
        vertex_field_values[vert_id] /= vertex_counts[vert_id];
      }
    });

  return omega_h_read_to_numpy(Omega_h::read(vertex_field_values));
}

// Map entity field to vertices using MLS interpolation
py::array_t<double> map_entity_field_to_vertices_mls(
  Omega_h::Mesh& mesh, py::array_t<double> entity_field_values,
  Omega_h::Int entity_dim, Omega_h::LO interp_degree, Omega_h::LO min_support,
  double lambda_reg)
{

  const Omega_h::Int dim = mesh.dim();

  // Compute source coordinates (entity centroids) and flatten
  auto centroids_flat = compute_entity_centroids(mesh, entity_dim);
  auto source_coordinates =
    numpy_to_omega_h_read<Omega_h::Real>(centroids_flat);

  // Get target coordinates (vertex coordinates)
  auto target_coordinates = mesh.coords();

  // Determine minimum support size if not provided
  if (min_support < 0) {
    min_support = patch_size_for_degree(interp_degree);
  }

  // Get patches (vertex-centered neighborhoods)
  auto patches = mesh.get_vtx_patches(min_support, entity_dim);
  auto supports_ptr = patches.a2ab;
  auto supports_idx = patches.ab2b;

  // Create radii array (all ones for now, as in Python version)
  const Omega_h::LO num_vertices = mesh.nverts();
  Omega_h::Write<Omega_h::Real> radii2_write(num_vertices, 1.0);

  // Create SupportResults
  SupportResults support{supports_ptr, supports_idx, radii2_write};

  // Convert entity field values to Omega_h array
  auto entity_field_omega_h =
    numpy_to_omega_h_read<Omega_h::Real>(entity_field_values);

  // Call MLS interpolation
  auto interpolated_values = pcms::mls_interpolation(
    entity_field_omega_h, source_coordinates, target_coordinates, support, dim,
    interp_degree, pcms::RadialBasisFunction::NO_OP, lambda_reg,
    1e-6, // tol
    5.0   // decay_factor
  );

  return omega_h_read_to_numpy(Omega_h::Reals(interpolated_values));
}

void bind_mesh_utilities_module(py::module& m)
{
  m.def("compute_entity_centroids", &compute_entity_centroids, py::arg("mesh"),
        py::arg("entity_dim"),
        "Compute centroids of mesh entities (faces or elements). "
        "Returns flattened array of shape (nentities * dim,).");

  // Backward compatibility alias
  m.def("compute_entity_centroids_few", &compute_entity_centroids,
        py::arg("mesh"), py::arg("entity_dim"),
        "Compute entity centroids using gather_verts/gather_vectors/average. "
        "Returns flattened array of shape (nentities * dim,). "
        "Note: This is an alias for compute_entity_centroids.");

  m.def(
    "patch_size_for_degree", &patch_size_for_degree, py::arg("interp_degree"),
    "Get recommended minimum patch size for given MLS interpolation degree. "
    "Returns 4 for degree 1, 10 for degree 2, and 12 for degree 3.");

  m.def(
    "map_entity_field_to_vertices_average",
    &map_entity_field_to_vertices_average, py::arg("mesh"),
    py::arg("entity_field_values"), py::arg("entity_dim"),
    "Map field values defined on mesh entities (faces/elements) to vertices "
    "by simple averaging. Each vertex value is the average of all adjacent "
    "entity values.\n\n"
    "Parameters:\n"
    "  mesh: OmegaHMesh object\n"
    "  entity_field_values: numpy array of field values at entities (shape: "
    "nentities,)\n"
    "  entity_dim: dimension of entities (2 for faces, 3 for elements)\n\n"
    "Returns:\n"
    "  numpy array of field values at vertices (shape: nverts,)");

  m.def(
    "map_entity_field_to_vertices_mls", &map_entity_field_to_vertices_mls,
    py::arg("mesh"), py::arg("entity_field_values"), py::arg("entity_dim"),
    py::arg("interp_degree") = 1, py::arg("min_support") = -1,
    py::arg("lambda_reg") = 1e-5,
    "Map field values defined on mesh entities (faces/elements) to vertices "
    "using Moving Least Squares (MLS) interpolation. This provides smoother "
    "and more accurate interpolation than simple averaging.\n\n"
    "Parameters:\n"
    "  mesh: OmegaHMesh object\n"
    "  entity_field_values: numpy array of field values at entities (shape: "
    "nentities,)\n"
    "  entity_dim: dimension of entities (2 for faces, 3 for elements)\n"
    "  interp_degree: polynomial degree for MLS (default: 1)\n"
    "  min_support: minimum patch size, auto-selected if < 0 (default: -1)\n"
    "  lambda_reg: regularization parameter (default: 1e-5)\n\n"
    "Returns:\n"
    "  numpy array of field values at vertices (shape: nverts,)");
}

} // namespace pcms
