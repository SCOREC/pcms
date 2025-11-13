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

// FIXME: What's the way to avoid global variables?
std::string degas2_mesh_filename = "";
std::string ltx_mesh_base_filename = "";
std::string ltx_node_data_filename = "";
std::string degas2_centroid_data_filename = "";

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
    Opt(ltx_node_data_filename,
        "xgc_node_data_filename")["--xgc_node_data"]("XGC node data file") |
    Opt(degas2_centroid_data_filename,
        "degas2_centroid_data_filename")["--degas2_centroid_data"](
      "Degas2 mesh centroid data file");

  session.cli(cli);
  int returnCode = session.applyCommandLine(argc, argv);
  if (returnCode != 0)
    return returnCode;

  return session.run(argc, argv);
}

TEST_CASE("Test Interpolation on LTX Mesh", "[interpolation]")
{
  auto lib = Omega_h::Library{};
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(degas2_mesh_filename, lib.world(), &mesh);
  const int degas2_num_elems = mesh.nelems();
  const auto degas2_mesh_centroids_host = Omega_h::HostRead(getCentroids(mesh));
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
                               2, 0.1, 10, 1, true, 0.0, 5.0);
  auto degas2_to_xgc_interpolator =
    MLSPointCloudInterpolation(degas2_mesh_centroids_view, xgc_mesh_points_view,
                               2, 0.1, 10, 1, true, 0.0, 5.0);
  printf("[INFO] Interpolators initialized.\n");

  // Read Fields from Files
  Omega_h::HostWrite<Omega_h::Real> density_at_xgc_nodes(xgc_num_nodes);
  Omega_h::HostWrite<Omega_h::Real> temp_at_xgc_nodes(xgc_num_nodes);
  Omega_h::HostWrite<Omega_h::Real> density_at_degas2_centroids(
    degas2_num_elems);
  Omega_h::HostWrite<Omega_h::Real> temp_at_degas2_centroids(degas2_num_elems);

  read_data(density_at_xgc_nodes, temp_at_xgc_nodes, ltx_node_data_filename);
  printf("[INFO] Data files loaded: %s and %s\n",
         ltx_node_data_filename.c_str(), degas2_centroid_data_filename.c_str());

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

  xgc_to_degas2_interpolator.eval(density_at_xgc_nodes_view,
                                  density_at_degas2_centroids_view);
  xgc_to_degas2_interpolator.eval(temp_at_xgc_nodes_view,
                                  temp_at_degas2_centroids_view);
  printf("[INFO] Interpolated data from XGC nodes to Degas2 centroids.\n");

  Omega_h::HostWrite<Omega_h::Real> density_at_xgc_nodes_interpolated_back(
    xgc_num_nodes);
  Omega_h::HostWrite<Omega_h::Real> temp_at_xgc_nodes_interpolated_back(
    xgc_num_nodes);
  const auto density_at_xgc_nodes_interpolated_back_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(
      density_at_xgc_nodes_interpolated_back.data(),
      density_at_xgc_nodes_interpolated_back.size());
  const auto temp_at_xgc_nodes_interpolated_back_view =
    pcms::Rank1View<double, pcms::HostMemorySpace>(
      temp_at_xgc_nodes_interpolated_back.data(),
      temp_at_xgc_nodes_interpolated_back.size());

  degas2_to_xgc_interpolator.eval(density_at_degas2_centroids_view,
                                  density_at_xgc_nodes_interpolated_back_view);
  degas2_to_xgc_interpolator.eval(temp_at_degas2_centroids_view,
                                  temp_at_xgc_nodes_interpolated_back_view);
  printf("[INFO] Interpolated data back from Degas2 centroids to XGC nodes.\n");

  // write to VTK in Degas2 mesh
  auto density_at_node_centroids_read =
    Omega_h::Read<Omega_h::Real>(density_at_degas2_centroids);
  auto temp_at_node_centroids_read =
    Omega_h::Read<Omega_h::Real>(temp_at_degas2_centroids);
  mesh.add_tag(Omega_h::FACE, "interpolated_density", 1,
               density_at_node_centroids_read);
  mesh.add_tag(Omega_h::FACE, "interpolated_temperature", 1,
               temp_at_node_centroids_read);

  std::string output_vtu_filename = "degas2_mesh.vtu";
  Omega_h::vtk::write_vtu(output_vtu_filename, &mesh);
  printf("[INFO] Wrote Degas2 mesh with interpolated data in %s.\n",
         output_vtu_filename.c_str());

  std::string output_xgc_vtu_filename = "xgc_mesh.vtu";
  write_xgc_mesh_as_vtu(
    output_xgc_vtu_filename, xgc_mesh_points, ltx_mesh_base_filename + ".ele",
    {density_at_xgc_nodes, temp_at_xgc_nodes,
     density_at_xgc_nodes_interpolated_back,
     temp_at_xgc_nodes_interpolated_back},
    {"original_density", "original_temperature", "interpolated_back_density",
     "interpolated_back_temperature"});
  printf(
    "[INFO] Wrote XGC mesh with original and interpolated back data in %s.\n",
    output_xgc_vtu_filename.c_str());

  // Compare original and interpolated back data at XGC nodes
  for (int i = 0; i < xgc_num_nodes; ++i) {
    double tol_percent = 10.0;
    REQUIRE_THAT(
      density_at_xgc_nodes_interpolated_back[i],
      Catch::Matchers::WithinRel(density_at_xgc_nodes[i], tol_percent / 100.0));

    REQUIRE_THAT(
      temp_at_xgc_nodes_interpolated_back[i],
      Catch::Matchers::WithinRel(temp_at_xgc_nodes[i], tol_percent / 100.0));
  }
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
    double x, y, z;
    line_stream >> node_id >> x >> y >> z;

    if (std::abs(z) > 1e-9 || std::isnan(z)) {
      // the ltx mesh has some 1.0 in the node file, which is fine
      if (!((std::abs(z) - 1.0) < 1e-9)) {
        printf("[WARNING] Non-zero Z coordinate found (%f) for node ID %d. "
               "Expected Z=0 for 2D mesh.\n",
               z, node_id);
      }
    }

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
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream line_stream(line);
    int index;
    double dens, temp;

    // Read index, density, and temperature from the line
    line_stream >> index >> dens >> temp;

    assert(index > 0 && index <= density.size() &&
           "Index out of bounds in data file");

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
