#include "pcms/utility/eqdsk.h"
#include "pcms/utility/assert.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cctype>

namespace pcms
{

// ============================================================================
// File readers for EQDSKData
// ============================================================================

void read_geqdsk(const std::string& filename, EQDSKData& eqdsk)
{
  constexpr double psi_factor = 1.0; // Default psi factor

  std::ifstream fin(filename);

  if (!fin) {
    throw std::runtime_error("Cannot open G-EQDSK file: " + filename);
  }

  // ----------------------------------------
  // Read header
  // ----------------------------------------
  std::string first_line;
  std::getline(fin, first_line);

  if (first_line.length() < 60) {
    first_line.resize(60, ' ');
  }

  int mw = std::stoi(first_line.substr(52, 4));
  int mh = std::stoi(first_line.substr(56, 4));

  // ----------------------------------------
  // Read geometry
  // ----------------------------------------
  double rdim, zdim, rzero, r_min, zmid;
  fin >> rdim >> zdim >> rzero >> r_min >> zmid;

  double ssimag, ssibry, dummy;
  fin >> dummy >> dummy >> ssimag >> ssibry >> dummy;
  fin >> dummy >> dummy >> dummy >> dummy >> dummy;
  fin >> dummy >> dummy >> dummy >> dummy >> dummy;

  // ----------------------------------------
  // Read fpol (poloidal current function)
  // ----------------------------------------
  std::vector<double> fpol(mw);
  for (int i = 0; i < mw; i++) {
    fin >> fpol[i];
  }

  // Skip pres, ffprim, pprime arrays (3 * mw values)
  for (int i = 0; i < 3 * mw; i++) {
    fin >> dummy;
  }

  // ----------------------------------------
  // Read psirz
  // ----------------------------------------
  std::vector<std::vector<double>> psirz(mw, std::vector<double>(mh));
  for (int j = 0; j < mh; j++) {
    for (int i = 0; i < mw; i++) {
      fin >> psirz[i][j];
    }
  }

  // Skip qpsi array
  for (int i = 0; i < mw; i++) {
    fin >> dummy;
  }

  // Skip boundary and limiter info
  int nbdry, nlim;
  fin >> nbdry >> nlim;
  for (int i = 0; i < nbdry + nlim; i++) {
    fin >> dummy >> dummy;
  }

  fin.close();

  // ----------------------------------------
  // Build EQDSKData
  // ----------------------------------------
  ssimag *= psi_factor;
  ssibry *= psi_factor;
  double z_min = zmid - zdim / 2.0;

  // Initialize grid structure
  // grid.divisions stores cell counts, which is one less than point counts
  // mw = NW = number of R grid points, mh = NH = number of Z grid points
  eqdsk.grid.divisions[0] = static_cast<LO>(mw - 1);
  eqdsk.grid.divisions[1] = static_cast<LO>(mh - 1);
  eqdsk.grid.edge_length[0] = static_cast<Real>(rdim);
  eqdsk.grid.edge_length[1] = static_cast<Real>(zdim);
  eqdsk.grid.bot_left[0] = static_cast<Real>(r_min);
  eqdsk.grid.bot_left[1] = static_cast<Real>(z_min);

  // Allocate and copy PSIZR
  eqdsk.PSIZR = Kokkos::View<Real*, DeviceMemorySpace>(
    "PSIZR", static_cast<size_t>(mw) * mh);
  auto psizr_host = Kokkos::create_mirror_view(eqdsk.PSIZR);
  for (int j = 0; j < mh; ++j) {
    for (int i = 0; i < mw; ++i) {
      psizr_host(i + j * mw) = static_cast<Real>(psirz[i][j] * psi_factor);
    }
  }
  Kokkos::deep_copy(eqdsk.PSIZR, psizr_host);

  // Build psi_grid
  eqdsk.psi_grid = Kokkos::View<Real*, DeviceMemorySpace>("psi_grid", mw);
  auto psi_grid_host = Kokkos::create_mirror_view(eqdsk.psi_grid);
  for (int i = 0; i < mw; ++i) {
    psi_grid_host(i) = static_cast<Real>(
      double(i) / double(mw - 1) * (ssibry - ssimag) + ssimag);
  }
  Kokkos::deep_copy(eqdsk.psi_grid, psi_grid_host);

  // Allocate and copy poloidal current function
  eqdsk.I_psi = Kokkos::View<Real*, DeviceMemorySpace>("I_psi", mw);
  auto I_psi_host = Kokkos::create_mirror_view(eqdsk.I_psi);
  for (int i = 0; i < mw; ++i) {
    I_psi_host(i) = static_cast<Real>(fpol[i]);
  }
  Kokkos::deep_copy(eqdsk.I_psi, I_psi_host);
}

void read_eqd(const std::string& filename, EQDSKData& eqdsk)
{
  std::ifstream fin(filename);

  if (!fin) {
    throw std::runtime_error("Cannot open EQD file: " + filename);
  }

  // Read header
  std::string eq_header;
  std::getline(fin, eq_header);

  // Read dimensions
  int mw, mh, mpsi;
  fin >> mw >> mh >> mpsi;

  // Read geometry
  double r_min, r_max, z_min, z_max;
  fin >> r_min >> r_max >> z_min >> z_max;

  double rmaxis, zmaxis, eq_axis_b;
  fin >> rmaxis >> zmaxis >> eq_axis_b;

  // Skip dummy variables
  double dummy;
  fin >> dummy >> dummy >> dummy;

  // Read psi_grid
  std::vector<double> psi_grid_vec(mpsi);
  for (int i = 0; i < mpsi; i++) {
    fin >> psi_grid_vec[i];
  }

  std::cout << "Axis (R,Z,B) = " << rmaxis << " " << zmaxis << " " << eq_axis_b
            << "\\n";

  // Read I(psi)
  std::vector<double> eq_I_vec(mpsi);
  for (int i = 0; i < mpsi; i++) {
    fin >> eq_I_vec[i];
  }

  // Read psi(R,Z)
  std::vector<std::vector<double>> psirz(mw, std::vector<double>(mh));
  for (int j = 0; j < mh; j++) {
    for (int i = 0; i < mw; i++) {
      fin >> psirz[i][j];
    }
  }

  // Verify end flag
  int end_flag;
  fin >> end_flag;
  if (end_flag != -1) {
    throw std::runtime_error("Wrong EQD file format");
  }

  fin.close();

  // ----------------------------------------
  // Build EQDSKData
  // ----------------------------------------

  // Initialize grid structure
  // grid.divisions stores cell counts, which is one less than point counts
  eqdsk.grid.divisions[0] = static_cast<LO>(mw - 1);
  eqdsk.grid.divisions[1] = static_cast<LO>(mh - 1);
  eqdsk.grid.edge_length[0] = static_cast<Real>(r_max - r_min);
  eqdsk.grid.edge_length[1] = static_cast<Real>(z_max - z_min);
  eqdsk.grid.bot_left[0] = static_cast<Real>(r_min);
  eqdsk.grid.bot_left[1] = static_cast<Real>(z_min);

  // Allocate and copy PSIZR
  eqdsk.PSIZR = Kokkos::View<Real*, DeviceMemorySpace>(
    "PSIZR", static_cast<size_t>(mw) * mh);
  auto psizr_host = Kokkos::create_mirror_view(eqdsk.PSIZR);
  for (int j = 0; j < mh; ++j) {
    for (int i = 0; i < mw; ++i) {
      psizr_host(i + j * mw) = static_cast<Real>(psirz[i][j]);
    }
  }
  Kokkos::deep_copy(eqdsk.PSIZR, psizr_host);

  // Allocate and copy psi_grid
  eqdsk.psi_grid = Kokkos::View<Real*, DeviceMemorySpace>("psi_grid", mpsi);
  auto psi_grid_host = Kokkos::create_mirror_view(eqdsk.psi_grid);
  for (int i = 0; i < mpsi; ++i) {
    psi_grid_host(i) = static_cast<Real>(psi_grid_vec[i]);
  }
  Kokkos::deep_copy(eqdsk.psi_grid, psi_grid_host);

  // Allocate and copy poloidal current function
  eqdsk.I_psi = Kokkos::View<Real*, DeviceMemorySpace>("I_psi", mpsi);
  auto I_psi_host = Kokkos::create_mirror_view(eqdsk.I_psi);
  for (int i = 0; i < mpsi; ++i) {
    I_psi_host(i) = static_cast<Real>(eq_I_vec[i]);
  }
  Kokkos::deep_copy(eqdsk.I_psi, I_psi_host);

  std::cout << "EQD file loaded successfully\n";
}

// ============================================================================
// EQDSKData Constructor
// ============================================================================
EQDSKData::EQDSKData(const std::string& filename)
{
  // Determine file type by extension
  std::string extension;
  size_t dot_pos = filename.find_last_of('.');
  if (dot_pos != std::string::npos) {
    extension = filename.substr(dot_pos);
    // Convert to lowercase for case-insensitive comparison
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) { return std::tolower(c); });
  }

  // Check for EQD format extensions
  if (extension == ".eqd") {
    read_eqd(filename, *this);
  }
  // Default to G-EQDSK format for .geq, .geqdsk, .g, or any other extension
  else {
    read_geqdsk(filename, *this);
  }
}

// Manual constructor for programmatic construction
EQDSKData::EQDSKData(const Uniform2DGrid& grid_in,
                     const std::vector<Real>& psirz,
                     const std::vector<Real>& psi_grid_vals,
                     const std::vector<Real>& I_psi_vals)
  : grid(grid_in)
{
  // grid.divisions stores cell counts, so number of grid points = divisions + 1
  const int mw = grid.divisions[0] + 1;
  const int mh = grid.divisions[1] + 1;

  // Validate input dimensions
  if (psirz.size() != static_cast<size_t>(mw) * mh) {
    throw std::runtime_error("psirz size mismatch: expected " +
                             std::to_string(mw * mh) + " but got " +
                             std::to_string(psirz.size()));
  }

  if (psi_grid_vals.size() != I_psi_vals.size()) {
    throw std::runtime_error("psi_grid and I_psi size mismatch: " +
                             std::to_string(psi_grid_vals.size()) + " vs " +
                             std::to_string(I_psi_vals.size()));
  }

  const int npsi = static_cast<int>(psi_grid_vals.size());

  // Allocate and copy PSIZR
  PSIZR = Kokkos::View<Real*, DeviceMemorySpace>("PSIZR",
                                                 static_cast<size_t>(mw) * mh);
  auto psizr_host = Kokkos::create_mirror_view(PSIZR);
  for (int i = 0; i < mw * mh; ++i) {
    psizr_host(i) = psirz[i];
  }
  Kokkos::deep_copy(PSIZR, psizr_host);

  // Allocate and copy psi_grid
  psi_grid = Kokkos::View<Real*, DeviceMemorySpace>("psi_grid", npsi);
  auto psi_grid_host = Kokkos::create_mirror_view(psi_grid);
  for (int i = 0; i < npsi; ++i) {
    psi_grid_host(i) = psi_grid_vals[i];
  }
  Kokkos::deep_copy(psi_grid, psi_grid_host);

  // Allocate and copy I_psi
  I_psi = Kokkos::View<Real*, DeviceMemorySpace>("I_psi", npsi);
  auto I_psi_host = Kokkos::create_mirror_view(I_psi);
  for (int i = 0; i < npsi; ++i) {
    I_psi_host(i) = I_psi_vals[i];
  }
  Kokkos::deep_copy(I_psi, I_psi_host);
}

} // namespace pcms
