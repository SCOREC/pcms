#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Omega_h_mesh.hpp>
#include <Omega_h_array_ops.hpp>
#include <pcms/field/evaluator/mls_interpolation.hpp>
#include <pcms/localization/adj_search.hpp>
#include <pcms/utility/mesh_geometry.h>
#include "numpy_array_transform.h"

namespace py = pybind11;

namespace pcms
{

namespace
{

// Number of components of per-entity values shaped (nentities,) or
// (nentities, ncomps). A 1D array has one component.
LO EntityComponents(const py::array& values, LO nentities)
{
  if (values.ndim() == 1 && values.shape(0) == nentities) {
    return 1;
  }
  if (values.ndim() == 2 && values.shape(0) == nentities &&
      values.shape(1) > 0) {
    return static_cast<LO>(values.shape(1));
  }
  throw std::invalid_argument("entity_field_values must have shape (" +
                              std::to_string(nentities) + ",) or (" +
                              std::to_string(nentities) + ", ncomps), got " +
                              detail::ShapeString(values));
}

// Vertex values with the same rank as the entity values they came from.
py::array ShapeLike(py::array_t<double> vertex_values,
                    const py::array& entity_values, LO ncomps)
{
  if (entity_values.ndim() == 1) {
    return vertex_values;
  }
  return vertex_values.reshape(
    std::vector<py::ssize_t>{vertex_values.size() / ncomps, ncomps});
}

} // namespace

// Helper function to compute entity centroids
py::array compute_entity_centroids(Omega_h::Mesh& mesh, Omega_h::Int entity_dim)
{
  const py::ssize_t dim = mesh.dim();
  return ToNumpy(pcms::get_entity_centroids(mesh, entity_dim))
    .reshape(std::vector<py::ssize_t>{mesh.nents(entity_dim), dim});
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
py::array map_entity_field_to_vertices_average(
  Omega_h::Mesh& mesh, const NumpyInput<double>& entity_field_values,
  Omega_h::Int entity_dim)
{

  const Omega_h::LO nverts = mesh.nverts();
  const Omega_h::LO nentities = mesh.nents(entity_dim);
  const Omega_h::LO ncomps = EntityComponents(entity_field_values, nentities);

  // Get entity to vertices adjacency
  auto entity_to_verts_omega_h = mesh.ask_verts_of(entity_dim);
  const Omega_h::LO verts_per_entity =
    entity_to_verts_omega_h.size() / nentities;

  // Convert numpy array to Omega_h array
  auto entity_field_omega_h = NumpyToOmegaHRead(entity_field_values);

  // Initialize vertex field and count arrays
  Omega_h::Write<Omega_h::Real> vertex_field_values(nverts * ncomps, 0.0);
  Omega_h::Write<Omega_h::Real> vertex_counts(nverts, 0.0);

  // Parallel accumulation
  Kokkos::parallel_for(
    "accumulate_to_vertices", nentities,
    OMEGA_H_LAMBDA(const Omega_h::LO entity_id) {
      for (Omega_h::LO v_local = 0; v_local < verts_per_entity; ++v_local) {
        const Omega_h::LO vert_id =
          entity_to_verts_omega_h[entity_id * verts_per_entity + v_local];
        for (Omega_h::LO c = 0; c < ncomps; ++c) {
          Kokkos::atomic_add(&vertex_field_values[vert_id * ncomps + c],
                             entity_field_omega_h[entity_id * ncomps + c]);
        }
        Kokkos::atomic_add(&vertex_counts[vert_id], 1.0);
      }
    });

  // Parallel averaging
  Kokkos::parallel_for(
    "average_vertices", nverts, OMEGA_H_LAMBDA(const Omega_h::LO vert_id) {
      if (vertex_counts[vert_id] > 0.0) {
        for (Omega_h::LO c = 0; c < ncomps; ++c) {
          vertex_field_values[vert_id * ncomps + c] /= vertex_counts[vert_id];
        }
      }
    });

  return ShapeLike(ToNumpy(Omega_h::read(vertex_field_values)),
                   entity_field_values, ncomps);
}

// Map entity field to vertices using MLS interpolation
py::array map_entity_field_to_vertices_mls(
  Omega_h::Mesh& mesh, const NumpyInput<double>& entity_field_values,
  Omega_h::Int entity_dim, Omega_h::LO interp_degree, Omega_h::LO min_support,
  double lambda_reg)
{

  const Omega_h::Int dim = mesh.dim();
  const Omega_h::LO nentities = mesh.nents(entity_dim);
  const Omega_h::LO ncomps = EntityComponents(entity_field_values, nentities);

  // Source coordinates are the entity centroids
  auto source_coordinates = pcms::get_entity_centroids(mesh, entity_dim);

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
  auto entity_field_omega_h = NumpyToOmegaHRead(entity_field_values);

  // mls_interpolation handles one component at a time
  Omega_h::Write<Omega_h::Real> vertex_field_values(num_vertices * ncomps);
  for (Omega_h::LO c = 0; c < ncomps; ++c) {
    Omega_h::Write<Omega_h::Real> component_values(nentities);
    Kokkos::parallel_for(
      "extract_component", nentities, OMEGA_H_LAMBDA(const Omega_h::LO i) {
        component_values[i] = entity_field_omega_h[i * ncomps + c];
      });

    auto interpolated_values = pcms::mls_interpolation(
      Omega_h::read(component_values), source_coordinates, target_coordinates,
      support, dim, interp_degree, pcms::RadialBasisFunction::NO_OP, lambda_reg,
      1e-6, // tol
      5.0   // decay_factor
    );

    Kokkos::parallel_for(
      "insert_component", num_vertices, OMEGA_H_LAMBDA(const Omega_h::LO v) {
        vertex_field_values[v * ncomps + c] = interpolated_values[v];
      });
  }

  return ShapeLike(ToNumpy(Omega_h::read(vertex_field_values)),
                   entity_field_values, ncomps);
}

void bind_mesh_utilities_module(py::module& m)
{
  m.def("compute_entity_centroids", &compute_entity_centroids, py::arg("mesh"),
        py::arg("entity_dim"),
        "Compute centroids of mesh entities (faces or elements). "
        "Returns an array of shape (nentities, dim).");

  // Backward compatibility alias
  m.def("compute_entity_centroids_few", &compute_entity_centroids,
        py::arg("mesh"), py::arg("entity_dim"),
        "Compute entity centroids using gather_verts/gather_vectors/average. "
        "Returns an array of shape (nentities, dim). "
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
    "(nentities,) or (nentities, ncomps))\n"
    "  entity_dim: dimension of entities (2 for faces, 3 for elements)\n\n"
    "Returns:\n"
    "  numpy array of field values at vertices (shape: (nverts,) or "
    "(nverts, ncomps), matching entity_field_values)");

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
    "(nentities,) or (nentities, ncomps))\n"
    "  entity_dim: dimension of entities (2 for faces, 3 for elements)\n"
    "  interp_degree: polynomial degree for MLS (default: 1)\n"
    "  min_support: minimum patch size, auto-selected if < 0 (default: -1)\n"
    "  lambda_reg: regularization parameter (default: 1e-5)\n\n"
    "Returns:\n"
    "  numpy array of field values at vertices (shape: (nverts,) or "
    "(nverts, ncomps), matching entity_field_values)");
}

} // namespace pcms
