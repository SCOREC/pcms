#ifndef PCMS_UTILITY_EQDSK_H
#define PCMS_UTILITY_EQDSK_H

#include "pcms/utility/types.h"
#include "pcms/utility/uniform_grid.h"
#include "pcms/utility/memory_spaces.h"
#include <vector>
#include <string>
#include <array>

namespace pcms
{

/**
 * @brief Data structure for EQDSK equilibrium data (simplified for field
 * evaluation)
 *
 * Contains the essential data needed for poloidal flux field evaluation:
 * - Uniform 2D grid structure (R-Z computational domain)
 * - PSIZR: Poloidal flux values on grid points
 * - psi_grid: Normalized poloidal flux grid values
 * - I_psi: Poloidal current function at psi_grid points (F(psi) for G-EQDSK,
 *   I(psi) for EQD format)
 *
 * All grid data is stored in column-major order (R varies fastest).
 */
struct EQDSKData
{
  // Grid structure for the R-Z computational domain
  // grid.divisions[0] = number of R cells (NW - 1)
  // grid.divisions[1] = number of Z cells (NH - 1)
  // grid.edge_length[0] = RDIM (horizontal dimension in meter)
  // grid.edge_length[1] = ZDIM (vertical dimension in meter)
  // grid.bot_left[0] = RLEFT (minimum R in meter)
  // grid.bot_left[1] = Z_min (minimum Z in meter)
  Uniform2DGrid grid;

  // PSIZR: Poloidal flux in Weber/rad on the rectangular grid points
  Kokkos::View<Real*, DeviceMemorySpace> PSIZR;

  // Psi grid: Normalized poloidal flux grid values
  Kokkos::View<Real*, DeviceMemorySpace> psi_grid;

  // I_psi: Poloidal current function I(psi) at psi_grid points
  Kokkos::View<Real*, DeviceMemorySpace> I_psi;

  // Constructors
  // Constructor from file
  EQDSKData(const std::string& filename);

  // Manual constructor for programmatic construction
  // grid: Uniform2DGrid structure defining the R-Z computational domain
  //   grid.divisions = number of cells in R and Z
  // psirz: 2D flux array (flattened in column-major order, size =
  //   (grid.divisions[0]+1)*(grid.divisions[1]+1))
  // psi_grid_vals: psi grid values (size = npsi)
  // I_psi_vals: poloidal current function values (size = npsi)
  EQDSKData(const Uniform2DGrid& grid, const std::vector<Real>& psirz,
            const std::vector<Real>& psi_grid_vals,
            const std::vector<Real>& I_psi_vals);

  // Accessor methods for backward compatibility
  [[nodiscard]] int GetNW() const noexcept { return grid.divisions[0] + 1; }
  [[nodiscard]] int GetNH() const noexcept { return grid.divisions[1] + 1; }
  [[nodiscard]] Real GetRDIM() const noexcept { return grid.edge_length[0]; }
  [[nodiscard]] Real GetZDIM() const noexcept { return grid.edge_length[1]; }
  [[nodiscard]] Real GetRLEFT() const noexcept { return grid.bot_left[0]; }
  [[nodiscard]] Real GetZMID() const noexcept
  {
    return grid.bot_left[1] + grid.edge_length[1] / 2.0;
  }

  // Grid bounds
  [[nodiscard]] Real GetZMin() const noexcept { return grid.bot_left[1]; }

  [[nodiscard]] Real GetZMax() const noexcept
  {
    return grid.bot_left[1] + grid.edge_length[1];
  }

  [[nodiscard]] Real GetRMin() const noexcept { return grid.bot_left[0]; }

  [[nodiscard]] Real GetRMax() const noexcept
  {
    return grid.bot_left[0] + grid.edge_length[0];
  }

  [[nodiscard]] Real GetDeltaR() const noexcept
  {
    const int num_R_cells = grid.divisions[0];
    return (num_R_cells > 0)
             ? grid.edge_length[0] / static_cast<Real>(num_R_cells)
             : 0.0;
  }

  [[nodiscard]] Real GetDeltaZ() const noexcept
  {
    const int num_Z_cells = grid.divisions[1];
    return (num_Z_cells > 0)
             ? grid.edge_length[1] / static_cast<Real>(num_Z_cells)
             : 0.0;
  }

  [[nodiscard]] int GetPsiIndex(int i_R, int i_Z) const noexcept
  {
    return i_R + i_Z * (grid.divisions[0] + 1);
  }

  [[nodiscard]] Real GetPsi(int i_R, int i_Z) const
  {
    return PSIZR[GetPsiIndex(i_R, i_Z)];
  }

  [[nodiscard]] Real GetR(int i_R) const noexcept
  {
    return grid.bot_left[0] + i_R * GetDeltaR();
  }

  [[nodiscard]] Real GetZ(int i_Z) const noexcept
  {
    return grid.bot_left[1] + i_Z * GetDeltaZ();
  }

  // Accessor for poloidal current data size
  [[nodiscard]] int GetNPsi() const noexcept
  {
    return static_cast<int>(psi_grid.extent(0));
  }

  // Accessor for psi grid value
  [[nodiscard]] Real GetPsiGrid(int i) const { return psi_grid[i]; }

  // Accessor for poloidal current at psi grid point
  [[nodiscard]] Real GetI(int i) const { return I_psi[i]; }
};

} // namespace pcms

#endif // PCMS_UTILITY_EQDSK_H
