#include <pcms/gnet/gnet_diagnostics.hpp>
#include <Omega_h_for.hpp>

namespace pcms::gnet
{

namespace
{
constexpr int dim = MAX_GNET_DIM;
}

Histogram2D calculate_velocity_distribution(
  const Omega_h::Reals& f5d, const Omega_h::Reals& sampled_points,
  const pcms::UniformGrid<MAX_GNET_DIM>& grid5d)
{

  int num_points = f5d.size();
  double umin = grid5d.bot_left[3];
  double vmin = grid5d.bot_left[4];
  double nbins_u = grid5d.divisions[3];
  double nbins_v = grid5d.divisions[4];

  double umax = umin + grid5d.edge_length[3];
  double vmax = vmin + grid5d.edge_length[4];
  double du = (umax - umin) / nbins_u;
  double dv = (vmax - vmin) / nbins_v;

  Histogram2D hist("velocity hist", grid5d.divisions[3], grid5d.divisions[4]);
  Kokkos::deep_copy(hist, 0.0);

  Omega_h::parallel_for(
    num_points, OMEGA_H_LAMBDA(const int i) {
      double u = sampled_points[i * dim + 3];
      double v = sampled_points[i * dim + 4];
      double fval = f5d[i];

      // Bin indices
      int uidx = static_cast<int>((u - umin) / du);
      int vidx = static_cast<int>((v - vmin) / dv);

      // check valid bin
      if (uidx >= 0 && uidx < nbins_u && vidx >= 0 && vidx < nbins_v) {
        Kokkos::atomic_add(&hist(uidx, vidx), fval);
      }
    });

  return hist;
}

Histogram2D calculate_2d_distribution(
  const Omega_h::Reals& f5d, const Omega_h::Reals& sampled_points,
  const pcms::UniformGrid<MAX_GNET_DIM>& grid5d)
{

  int num_points = f5d.size();
  double rho_min = grid5d.bot_left[0];
  double theta_min = grid5d.bot_left[1];
  double nbins_rho = grid5d.divisions[0];
  double nbins_theta = grid5d.divisions[1];

  double rho_max = rho_min + grid5d.edge_length[0];
  double theta_max = theta_min + grid5d.edge_length[1];
  double drho = (rho_max - rho_min) / nbins_rho;
  double dtheta = (theta_max - theta_min) / nbins_theta;

  Histogram2D hist("velocity hist", grid5d.divisions[0], grid5d.divisions[1]);
  Kokkos::deep_copy(hist, 0.0);

  Omega_h::parallel_for(
    num_points, OMEGA_H_LAMBDA(const int i) {
      double rho = sampled_points[i * dim + 0];
      double theta = sampled_points[i * dim + 1];
      double fval = f5d[i];

      // Bin indices
      int rho_idx = static_cast<int>((rho - rho_min) / drho);
      int theta_idx = static_cast<int>((theta - theta_min) / dtheta);

      // check valid bin
      if (rho_idx >= 0 && rho_idx < nbins_rho && theta_idx >= 0 &&
          theta_idx < nbins_theta) {
        Kokkos::atomic_add(&hist(rho_idx, theta_idx), fval);
      }
    });

  return hist;
}

Kokkos::View<double**, Kokkos::LayoutRight> get2DArray(
  const Omega_h::Reals& f5d, const Omega_h::Reals& sampled_points)
{

  int num_points = f5d.size();
  Kokkos::View<double**, Kokkos::LayoutRight> array(
    "stores the coordinates and values for GTC marker particles", num_points,
    dim + 1);
  Kokkos::deep_copy(array, 0.0);

  Omega_h::parallel_for(
    num_points, OMEGA_H_LAMBDA(const int i) {
      int base_idx_spts = i * dim;
      for (int d = 0; d < dim; ++d) {
        array(i, d) = sampled_points[base_idx_spts + d];
      }
      array(i, dim) = f5d[i];
    });

  return array;
}

} // namespace pcms::gnet
