//
// Created by Fuad Hasan on 11/11/25.
// This test verifies interpolation on an LTX reactor mesh
// The XGC mesh is a point cloud and the DEGAS2 mesh is in Omega_h format
// From XGC to DEGAS2 Centroids and Vice versa with two different data
//

#include <catch2/catch_session.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_vtk.hpp>
#include <pcms/interpolator/interpolation_base.h>
#include <pcms/print.h>

#include <vector>
#include <fstream>

std::vector<double> read_xgc_mesh_nodes(std::string filename);
void read_data(Omega_h::HostWrite<Omega_h::Real> density,
               Omega_h::HostWrite<Omega_h::Real> temperature,
               std::string filename);
void write_xgc_mesh_as_vtu(
  const std::string& output_filename, const std::vector<double>& node_coords,
  const std::string& connectivity_file,
  const std::vector<Omega_h::HostWrite<Omega_h::Real>>& node_data,
  const std::vector<std::string>& data_names);
double compute_l2_norm(const Omega_h::HostWrite<Omega_h::Real>& data,
                       const Omega_h::HostWrite<Omega_h::Real>& reference);

// FIXME: What's the way to avoid global variables?
std::string degas2_mesh_filename = "";
std::string ltx_mesh_base_filename = "";
std::string data_root_dir = "";

int main(const int argc, char* argv[])
{
  Catch::Session session;
  pcms::printInfo("Starting interpolation test on LTX mesh...");
  using namespace Catch::Clara;

  const auto cli =
    session.cli() |
    Opt(degas2_mesh_filename, "degas2_mesh_filename")["--degas2_mesh"](
      "Degas2 mesh file in Omega_h binary format") |
    Opt(ltx_mesh_base_filename, "ltx_mesh_base_filename")["--ltx_mesh"](
      "LTX mesh file in XGC mesh format (needs both .node and .ele)") |
    Opt(data_root_dir, "data_root_dir")["--data_root"](
      "Root directory for data files. It needs to these 8 files:"
      "degas2_data_0.original.txt degas2_data_0.txt "
      "degas2_interpolated_data_0.original.txt degas2_interpolated_data_0.txt "
      "xgc_data_0.original.txt xgc_data_0.txt "
      "xgc_interpolated_data_0.original.txt xgc_interpolated_data_0.txt");

  session.cli(cli);
  int returnCode = session.applyCommandLine(argc, argv);
  if (returnCode != 0)
    return returnCode;

  return session.run(argc, argv);
}

TEST_CASE("Test Interpolation on LTX Mesh", "[interpolation]")
{
  // ---------------------------- Loading Mesh ------------------- //
  auto lib = Omega_h::Library{};
  Omega_h::Mesh degas2_mesh(&lib);
  Omega_h::binary::read(degas2_mesh_filename, lib.world(), &degas2_mesh);

  // --------------------- Initialize Interpolators -------------- //
  const int degas2_num_elems = degas2_mesh.nelems();
  const auto degas2_mesh_centroids_host =
    Omega_h::HostRead(getCentroids(degas2_mesh));
  printf("[INFO] Degas2 Mesh loaded from %s with %d elements\n",
         degas2_mesh_filename.c_str(), degas2_num_elems);
  const auto degas2_mesh_centroids_view =
    pcms::Rank1View<const double, pcms::HostMemorySpace>(
      degas2_mesh_centroids_host.data(), degas2_mesh_centroids_host.size());

  auto xgc_mesh_points = read_xgc_mesh_nodes(ltx_mesh_base_filename + ".node");
  const int xgc_num_nodes = xgc_mesh_points.size() / 2;
  printf("[INFO] XGC Mesh loaded from %s with %d points\n",
         ltx_mesh_base_filename.c_str(), xgc_num_nodes);
  const auto xgc_mesh_points_view =
    pcms::Rank1View<const double, pcms::HostMemorySpace>(
      xgc_mesh_points.data(), xgc_mesh_points.size());

  auto xgc_to_degas2_interpolator =
    MLSPointCloudInterpolation(xgc_mesh_points_view, degas2_mesh_centroids_view,
                               2, 0.00001, 10, 1, true, 0.0, 50.0);
  auto degas2_to_xgc_interpolator =
    MLSPointCloudInterpolation(degas2_mesh_centroids_view, xgc_mesh_points_view,
                               2, 0.1, 10, 1, true, 0.0, 50.0);
  printf("[INFO] Interpolators initialized.\n");

  // ---------------------- Load Data ---------------------- //
  Omega_h::HostWrite<Omega_h::Real> density_at_xgc_nodes(xgc_num_nodes);
  Omega_h::HostWrite<Omega_h::Real> temp_at_xgc_nodes(xgc_num_nodes);
  Omega_h::HostWrite<Omega_h::Real> density_at_degas2_centroids(
    degas2_num_elems);
  Omega_h::HostWrite<Omega_h::Real> temp_at_degas2_centroids(degas2_num_elems);

  read_data(density_at_xgc_nodes, temp_at_xgc_nodes,
            data_root_dir + "/xgc_data_0.original.txt");
  printf("[INFO] XGC node data loaded: %s\n",
         (data_root_dir + "/xgc_data_0.original.txt").c_str());
  read_data(density_at_degas2_centroids, temp_at_degas2_centroids,
            data_root_dir + "/degas2_data_0.txt");
  printf("[INFO] Degas2 centroid data loaded: %s\n",
         (data_root_dir + "/degas2_data_0.txt").c_str());

  // ------------------ First Interpolation ------------------ //
  const auto density_at_xgc_nodes_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(density_at_xgc_nodes.data(),
                                                   density_at_xgc_nodes.size());
  const auto temp_at_xgc_nodes_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(temp_at_xgc_nodes.data(),
                                                   temp_at_xgc_nodes.size());

  const auto density_at_degas2_centroids_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(
      density_at_degas2_centroids.data(), density_at_degas2_centroids.size());
  const auto temp_at_degas2_centroids_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(
      temp_at_degas2_centroids.data(), temp_at_degas2_centroids.size());

  Omega_h::HostWrite<Omega_h::Real> interpolated_xgc_density(degas2_num_elems);
  Omega_h::HostWrite<Omega_h::Real> interpolated_xgc_temp(degas2_num_elems);
  Omega_h::HostWrite<Omega_h::Real> interpolated_degas2_density(xgc_num_nodes);
  Omega_h::HostWrite<Omega_h::Real> interpolated_degas2_temp(xgc_num_nodes);

  const auto interpolated_xgc_density_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(
      interpolated_xgc_density.data(), interpolated_xgc_density.size());
  const auto interpolated_xgc_temp_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(
      interpolated_xgc_temp.data(), interpolated_xgc_temp.size());
  const auto interpolated_degas2_density_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(
      interpolated_degas2_density.data(), interpolated_degas2_density.size());
  const auto interpolated_degas2_temp_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(
      interpolated_degas2_temp.data(), interpolated_degas2_temp.size());

  xgc_to_degas2_interpolator.eval(density_at_xgc_nodes_view,
                                  interpolated_xgc_density_view);
  xgc_to_degas2_interpolator.eval(temp_at_xgc_nodes_view,
                                  interpolated_xgc_temp_view);
  printf("[INFO] Interpolated data from XGC nodes to Degas2 centroids.\n");

  degas2_to_xgc_interpolator.eval(density_at_degas2_centroids_view,
                                  interpolated_degas2_density_view);
  degas2_to_xgc_interpolator.eval(temp_at_degas2_centroids_view,
                                  interpolated_degas2_temp_view);
  printf("[INFO] Interpolated data from Degas2 centroids to XGC nodes");

  // ------------------ Interpolation Back ------------------ //
  Omega_h::HostWrite<Omega_h::Real> interpolated_back_density_at_xgc_nodes(
    xgc_num_nodes);
  Omega_h::HostWrite<Omega_h::Real> interpolated_back_temp_at_xgc_nodes(
    xgc_num_nodes);
  Omega_h::HostWrite<Omega_h::Real>
    interpolated_back_density_at_degas2_centroids(degas2_num_elems);
  Omega_h::HostWrite<Omega_h::Real> interpolated_back_temp_at_degas2_centroids(
    degas2_num_elems);

  const auto interpolated_back_density_at_xgc_nodes_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(
      interpolated_back_density_at_xgc_nodes.data(),
      interpolated_back_density_at_xgc_nodes.size());
  const auto interpolated_back_temp_at_xgc_nodes_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(
      interpolated_back_temp_at_xgc_nodes.data(),
      interpolated_back_temp_at_xgc_nodes.size());
  const auto interpolated_back_density_at_degas2_centroids_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(
      interpolated_back_density_at_degas2_centroids.data(),
      interpolated_back_density_at_degas2_centroids.size());
  const auto interpolated_back_temp_at_degas2_centroids_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(
      interpolated_back_temp_at_degas2_centroids.data(),
      interpolated_back_temp_at_degas2_centroids.size());

  degas2_to_xgc_interpolator.eval(interpolated_xgc_density_view,
                                  interpolated_back_density_at_xgc_nodes_view);
  degas2_to_xgc_interpolator.eval(interpolated_xgc_temp_view,
                                  interpolated_back_temp_at_xgc_nodes_view);
  printf("[INFO] Interpolated back data from Degas2 centroids to XGC nodes.\n");
  xgc_to_degas2_interpolator.eval(
    interpolated_degas2_density_view,
    interpolated_back_density_at_degas2_centroids_view);
  xgc_to_degas2_interpolator.eval(
    interpolated_degas2_temp_view,
    interpolated_back_temp_at_degas2_centroids_view);
  printf("[INFO] Interpolated back data from XGC nodes to Degas2 centroids.\n");

  // ------------------ Write Output VTU Files ------------------ //
  degas2_mesh.add_tag(Omega_h::FACE, "density_at_degas2_centroids", 1,
                      Omega_h::Reals(density_at_degas2_centroids));
  degas2_mesh.add_tag(Omega_h::FACE, "temperature_at_degas2_centroids", 1,
                      Omega_h::Reals(temp_at_degas2_centroids));
  degas2_mesh.add_tag(
    Omega_h::FACE, "interpolated_back_density_at_degas2_centroids", 1,
    Omega_h::Reals(interpolated_back_density_at_degas2_centroids));
  degas2_mesh.add_tag(
    Omega_h::FACE, "interpolated_back_temperature_at_degas2_centroids", 1,
    Omega_h::Reals(interpolated_back_temp_at_degas2_centroids));
  degas2_mesh.add_tag(Omega_h::FACE, "interpolated_xgc_density", 1,
                      Omega_h::Reals(interpolated_xgc_density));
  degas2_mesh.add_tag(Omega_h::FACE, "interpolated_xgc_temperature", 1,
                      Omega_h::Reals(interpolated_xgc_temp));

  std::string output_vtu_filename = "degas2_mesh.vtu";
  Omega_h::vtk::write_vtu(output_vtu_filename, &degas2_mesh);
  printf("[INFO] Wrote Degas2 mesh with all data in %s.\n",
         output_vtu_filename.c_str());

  std::string output_xgc_vtu_filename = "xgc_mesh.vtu";
  std::vector xgc_node_data_arrays = {density_at_xgc_nodes,
                                      temp_at_xgc_nodes,
                                      interpolated_back_density_at_xgc_nodes,
                                      interpolated_back_temp_at_xgc_nodes,
                                      interpolated_degas2_density,
                                      interpolated_degas2_temp};
  std::vector<std::string> xgc_node_data_names = {
    "density_at_xgc_nodes",
    "temperature_at_xgc_nodes",
    "interpolated_back_density_at_xgc_nodes",
    "interpolated_back_temperature_at_xgc_nodes",
    "interpolated_degas2_density",
    "interpolated_degas2_temperature"};

  write_xgc_mesh_as_vtu(output_xgc_vtu_filename, xgc_mesh_points,
                        ltx_mesh_base_filename + ".ele", xgc_node_data_arrays,
                        xgc_node_data_names);
  printf("[INFO] Wrote XGC mesh with all data in %s.\n",
         output_xgc_vtu_filename.c_str());

  // ------------------ Verification ------------------ //
  double tol = 10.0 / 100.0; // 10 percent tolerance
  for (int i = 0; i < xgc_num_nodes; ++i) {
    CHECK_THAT(interpolated_back_density_at_xgc_nodes[i],
               Catch::Matchers::WithinRel(density_at_xgc_nodes[i], tol) ||
                 Catch::Matchers::WithinAbs(density_at_xgc_nodes[i], tol));

    CHECK_THAT(interpolated_back_temp_at_xgc_nodes[i],
               Catch::Matchers::WithinRel(temp_at_xgc_nodes[i], tol) ||
                 Catch::Matchers::WithinAbs(temp_at_xgc_nodes[i], tol));
  }

  for (int i = 0; i < degas2_num_elems; ++i) {
    CHECK_THAT(
      interpolated_back_density_at_degas2_centroids[i],
      Catch::Matchers::WithinRel(density_at_degas2_centroids[i], tol) ||
        Catch::Matchers::WithinAbs(density_at_degas2_centroids[i], tol));

    CHECK_THAT(interpolated_back_temp_at_degas2_centroids[i],
               Catch::Matchers::WithinRel(temp_at_degas2_centroids[i], tol) ||
                 Catch::Matchers::WithinAbs(temp_at_degas2_centroids[i], tol));
  }

  double l2_norm_density_xgc = compute_l2_norm(
    interpolated_back_density_at_xgc_nodes, density_at_xgc_nodes);
  double l2_norm_temp_xgc =
    compute_l2_norm(interpolated_back_temp_at_xgc_nodes, temp_at_xgc_nodes);
  double l2_norm_density_degas2 = compute_l2_norm(
    interpolated_back_density_at_degas2_centroids, density_at_degas2_centroids);
  double l2_norm_temp_degas2 = compute_l2_norm(
    interpolated_back_temp_at_degas2_centroids, temp_at_degas2_centroids);

  printf("[INFO] L2 Norms of errors after interpolation back:\n");
  printf("       Density\t(XGC -> Degas2 --> XGC):\t%e\n", l2_norm_density_xgc);
  printf("       Temperature\t(XGC -> Degas2 --> XGC):\t%e\n",
         l2_norm_temp_xgc);
  printf("       Density\t(Degas2 -> XGC --> Degas2):\t%e\n",
         l2_norm_density_degas2);
  printf("       Temperature\t(Degas2 -> XGC --> Degas2):\t%e\n",
         l2_norm_temp_degas2);

  // require that l2 norms are within 0.02
  double l2_tol = 0.02;
  CHECK(l2_norm_density_xgc < l2_tol);
  CHECK(l2_norm_temp_xgc < l2_tol);
  CHECK(l2_norm_density_degas2 < l2_tol);
  CHECK(l2_norm_temp_degas2 < l2_tol);
}

std::vector<double> read_xgc_mesh_nodes(std::string filename)
{
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open file: " + filename);
  }

  std::string line;
  std::getline(file, line); // Read the first line
  std::istringstream header_stream(line);
  int num_nodes, dim, dummy1, dummy2;
  header_stream >> num_nodes >> dim >> dummy1 >> dummy2;

  assert(dim == 2 && "Expected 2D coordinates in the file");

  std::vector<double> nodes;
  nodes.reserve(num_nodes * dim);

  while (std::getline(file, line)) {
    std::istringstream line_stream(line);
    int node_id;
    double x, y;
    int boundary_attr; // not used
    line_stream >> node_id >> x >> y >> boundary_attr;

    nodes.push_back(x);
    nodes.push_back(y);
  }

  file.close();
  return nodes;
}

void read_data(Omega_h::HostWrite<Omega_h::Real> density,
               Omega_h::HostWrite<Omega_h::Real> temperature,
               std::string filename)
{
  assert(density.size() == temperature.size() &&
         "Density and Temperature vectors must be of the same size");
  assert(density.size() > 0 && "Data vectors must not be empty");

  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open file: " + filename);
  }

  int data_count = 0;

  int index;
  double dens, temp;
  while (file >> index >> dens >> temp) {

    OMEGA_H_CHECK_PRINTF((index > 0 && index <= density.size()),
                         "Index %d out of bounds (1 to %d)\n", index,
                         density.size());

    density[index - 1] = dens;
    temperature[index - 1] = temp;
    data_count++;
  }

  file.close();
  OMEGA_H_CHECK_PRINTF(data_count == density.size(),
                       "Expected %d data points, but read %d\n", density.size(),
                       data_count);
}

std::vector<std::array<int, 3>> read_xgc_mesh_triangle_connectivity(
  const std::string& filename)
{
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("Unable to open file: " + filename);
  }

  std::string line;
  std::getline(file, line); // Read the first line
  std::istringstream header_stream(line);
  int num_elements, dummy1, dummy2;
  header_stream >> num_elements >> dummy1 >> dummy2;

  std::vector<std::array<int, 3>> connectivity;
  connectivity.reserve(num_elements);

  while (std::getline(file, line)) {
    std::istringstream line_stream(line);
    int element_id, v1, v2, v3;
    line_stream >> element_id >> v1 >> v2 >> v3;

    // Convert 1-based indices to 0-based indices
    connectivity.push_back({v1 - 1, v2 - 1, v3 - 1});
  }

  file.close();
  return connectivity;
}

void write_xgc_mesh_as_vtu(
  const std::string& output_filename, const std::vector<double>& node_coords,
  const std::string& connectivity_file,
  const std::vector<Omega_h::HostWrite<Omega_h::Real>>& node_data,
  const std::vector<std::string>& data_names)
{
  const size_t n_points = node_coords.size() / 2;
  assert(node_coords.size() % 2 == 0);
  assert(node_data.size() == data_names.size());
  for (const auto& data : node_data) {
    assert(data.size() == n_points);
  }

  std::vector<std::array<int, 3>> connectivity =
    read_xgc_mesh_triangle_connectivity(connectivity_file);

  const size_t n_cells = connectivity.size();

  std::ofstream file(output_filename);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + output_filename);
  }

  // Write VTU header
  file << "<?xml version=\"1.0\"?>\n";
  file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
          "byte_order=\"LittleEndian\">\n";
  file << "  <UnstructuredGrid>\n";
  file << "    <Piece NumberOfPoints=\"" << n_points << "\" NumberOfCells=\""
       << n_cells << "\">\n";

  // Write points (coordinates) - VTK requires 3D coordinates, so we add 0.0 for
  // Z
  file << "      <Points>\n";
  file << "        <DataArray type=\"Float64\" Name=\"Points\" "
          "NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (size_t i = 0; i < n_points; ++i) {
    file << "          " << node_coords[2 * i] << " " << node_coords[2 * i + 1]
         << " 0.0\n";
  }
  file << "        </DataArray>\n";
  file << "      </Points>\n";

  // Write point data
  file << "      <PointData>\n";
  for (size_t d = 0; d < node_data.size(); ++d) {
    const auto& data_name = data_names[d];
    const auto& data_values = node_data[d];

    file << "        <DataArray type=\"Float64\" Name=\"" << data_name
         << "\" format=\"ascii\">\n";
    for (size_t i = 0; i < data_values.size(); ++i) {
      file << "          " << data_values[i] << "\n";
    }
    file << "        </DataArray>\n";
  }
  file << "      </PointData>\n";

  // Write cells (connectivity)
  file << "      <Cells>\n";

  // Write connectivity
  file << "        <DataArray type=\"Int32\" Name=\"connectivity\" "
          "format=\"ascii\">\n";
  for (const auto& triangle : connectivity) {
    file << "          " << triangle[0] << " " << triangle[1] << " "
         << triangle[2] << "\n";
  }
  file << "        </DataArray>\n";

  // Write offsets
  file
    << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for (size_t i = 1; i <= n_cells; ++i) {
    file << "          " << i * 3 << "\n";
  }
  file << "        </DataArray>\n";

  // Write cell types (5 for triangles)
  file
    << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (size_t i = 0; i < n_cells; ++i) {
    file << "          5\n";
  }
  file << "        </DataArray>\n";

  file << "      </Cells>\n";

  // Write footer
  file << "    </Piece>\n";
  file << "  </UnstructuredGrid>\n";
  file << "</VTKFile>\n";

  file.close();
}

double compute_l2_norm(const Omega_h::HostWrite<Omega_h::Real>& data,
                       const Omega_h::HostWrite<Omega_h::Real>& reference)
{
  assert(data.size() == reference.size() &&
         "Data and reference must be of the same size");

  double sum_sq_diff = 0.0;
  double sum_sq_ref = 0.0;

  for (size_t i = 0; i < data.size(); ++i) {
    const double diff = data[i] - reference[i];
    sum_sq_diff += diff * diff;
    sum_sq_ref += reference[i] * reference[i];
  }

  double l2_norm = std::sqrt(sum_sq_diff);
  double l2_ref = std::sqrt(sum_sq_ref);

  if (l2_ref > 0.0) {
    return l2_norm / l2_ref; // Return relative L2 norm
  } else {
    return l2_norm; // If reference is zero, return absolute L2 norm
  }
}
