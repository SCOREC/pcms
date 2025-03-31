#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/interpolator/mls_interpolation.hpp>
#include <pcms/interpolator/mls_interpolation_impl.hpp>
#include <pcms/interpolator/pcms_interpolator_view_utils.hpp>
#include "KokkosBatched_SVD_Decl.hpp"
#include "KokkosBatched_SVD_Serial_Impl.hpp"
#include <KokkosBlas2_serial_gemv_impl.hpp>
#include <KokkosBatched_Gemm_Decl.hpp>
#include <pcms/interpolator/pcms_interpolator_aliases.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_library.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <vector>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <adios2.h>
#include <mpi.h>
#include <fstream>
#include <string>

using namespace Omega_h;
using namespace pcms;

std::vector<double> readTeData(const std::string& filePath, MPI_Comm comm = MPI_COMM_WORLD) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    std::vector<double> Te_data;

    // Ensure MPI is initialized
    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized) {
        MPI_Init(nullptr, nullptr);
    }

    // Initialize ADIOS2
    adios2::ADIOS adios(comm);
    adios2::IO io = adios.DeclareIO("BoutIO");

    // Ensure the engine is set to BP format
    if (!io.InConfigFile()) {
        io.SetEngine("BP");
    }

    // Open the file in STREAMING mode
    adios2::Engine reader = io.Open(filePath, adios2::Mode::Read, comm);
    if (!reader) {
        std::cerr << "Error: Failed to open file " << filePath << std::endl;
        return Te_data;
    }

    // BeginStep (needed for streaming mode)
    adios2::StepStatus status = reader.BeginStep();
    if (status != adios2::StepStatus::OK) {
        std::cerr << "Error: Failed to start step in ADIOS2 reader.\n";
        reader.Close();
        return Te_data;
    }

    // Inquire the variable "Te"
    auto var_Te = io.InquireVariable<double>("Te");
    if (!var_Te) {
        std::cerr << "Error: Could not find variable 'Te' in the BOUT file.\n";
        reader.EndStep();
        reader.Close();
        return Te_data;
    }

    // Get the shape of the dataset
    auto shape = var_Te.Shape();

    // Debugging: Print detected shape
    if (rank == 0) {  // Only print from one rank to avoid excessive output
        std::cout << "Detected shape of 'Te': ";
        for (auto s : shape) {
            std::cout << s << " ";
        }
        std::cout << std::endl;
    }

    // Ensure shape is valid (must be 3D)
    if (shape.size() != 3) {
        std::cerr << "Error: Unexpected shape for 'Te'. Expected 3D array (Nx, Ny, Nz).\n";
        reader.EndStep();
        reader.Close();
        return Te_data;
    }

    size_t Nx = shape[0];   // 20 (number of x grid points)
    size_t Ny = shape[1];   // 10 (number of y grid points)
    size_t Nz = shape[2];   // 20 (number of z grid points)

    // Distribute workload across MPI ranks
    size_t N_per_rank = Nx / size;
    size_t start_x = rank * N_per_rank;
    size_t count_x = (rank == size - 1) ? Nx - rank * N_per_rank : N_per_rank;

    // Resize output vector to hold (count_x * Ny * Nz) values for this rank
    Te_data.resize(count_x * Ny * Nz);

    // Set selection for distributed x-range
    var_Te.SetSelection({{start_x, 0, 0}, {count_x, Ny, Nz}});

    // Read data into the vector
    reader.Get(var_Te, Te_data.data(), adios2::Mode::Sync);

    // EndStep must be called before closing
    reader.EndStep();
    reader.Close();

    // // Print first 10 values for debugging (only from rank 0)
    // if (rank == 0) {
    //     std::cout << "Te values (first 10 entries):\n";
    //     for (size_t i = 0; i < std::min(size_t(10), Te_data.size()); ++i) {
    //         std::cout << "Te[" << i << "] = " << Te_data[i] << "\n";
    //     }
    // }

    return Te_data;
}

std::vector<std::pair<std::vector<double>, double>> generateSourcePointsWithTe(const std::vector<double>& Te_data) {
    const int radial = 20;      // R direction points (Nx)
    const int poloidal = 20;    // Z direction points (Nz)
    const int toroidal = 10;    // θ direction points (Ny)

    double R_min = 1.4, R_max = 2.6;
    // double Z_min = -0.6, Z_max = 0.6;
    double Z_min = 0.0, Z_max = 1.2;
    double dy = 0.18; 
    double dR = (R_max - R_min) / (radial - 1);
    double dZ = (Z_max - Z_min) / (poloidal - 1);

    std::vector<std::pair<std::vector<double>, double>> source_points_with_Te;

    for (int i = 0; i < radial; ++i) {
        double R = R_min + i * dR;
        for (int j = 0; j < poloidal; ++j) {
            double Z = Z_min + j * dZ;
            for (int k = 0; k < toroidal; ++k) {
                double theta = k * dy;

                // Convert to Cartesian coordinates
                double x = R * cos(theta);
                double y = R * sin(theta);
                double z = Z;

                // Compute correct Te index (column-major: k + Ny * (j + Nz * i))
                size_t index = k + toroidal * (j + poloidal * i);

                // Bounds check to prevent errors
                if (index >= Te_data.size()) {
                    std::cerr << "Index out of bounds: " << index << " for Te_data size " << Te_data.size() << std::endl;
                    continue;
                }

                double Te_value = Te_data[index];

                // Store (x, y, z) with corresponding Te
                source_points_with_Te.push_back({{x, y, z}, Te_value});
            }
        }
    }
    return source_points_with_Te;
}

// n² search to find neighbors (Brute force search)
inline SupportResults findNeighbors(
    const std::vector<std::vector<double>>& host_source_points,
    const std::vector<std::vector<double>>& host_target_points,
    double cutoffDistance) {

    SupportResults results;
    int numTargets = host_target_points.size();
    int numSources = host_source_points.size();
    double cutoffDistanceSq = cutoffDistance * cutoffDistance;

    if (numSources == 0 || numTargets == 0) {
        std::cerr << "Warning: No source or target points provided!\n";
        return results;
    }

    // Allocate storage
    Write<LO> nSupports(numTargets, 0, "number of supports for each target");
    Write<Real> radii2(numTargets, cutoffDistanceSq, "squared radii of supports");

    // Create Kokkos Views for device storage
    Kokkos::View<double**> source_points("source_points", numSources, 3);
    Kokkos::View<double**> target_points("target_points", numTargets, 3);

    // Create host mirrors
    auto host_source_view = Kokkos::create_mirror_view(source_points);
    auto host_target_view = Kokkos::create_mirror_view(target_points);

    // Copy data from std::vector to host mirror views
    for (int i = 0; i < numSources; ++i) {
        for (int d = 0; d < 3; ++d) {  // 3D coordinates
            host_source_view(i, d) = host_source_points[i][d];
        }
    }

    for (int i = 0; i < numTargets; ++i) {
        for (int d = 0; d < 3; ++d) {
            host_target_view(i, d) = host_target_points[i][d];
        }
    }

    // Copy from host mirror to device views
    Kokkos::deep_copy(source_points, host_source_view);
    Kokkos::deep_copy(target_points, host_target_view);

    // Count neighbors 
    Kokkos::parallel_for("count neighbors", numTargets, KOKKOS_LAMBDA(int i) {
        int count = 0;
        for (int j = 0; j < numSources; ++j) {
            double dx = target_points(i, 0) - source_points(j, 0);
            double dy = target_points(i, 1) - source_points(j, 1);
            double dz = target_points(i, 2) - source_points(j, 2);
            double dist_sq = dx * dx + dy * dy + dz * dz;

            if (dist_sq <= cutoffDistanceSq) {
                count++;  
            }
        }
        nSupports[i] = count;
    });

    // Compute prefix sum (scan index)
    Write<LO> supports_ptr(numTargets + 1, 0, "scan index");
    LO total_supports = 0;
    Kokkos::parallel_scan(
        "compute scan index", numTargets,
        KOKKOS_LAMBDA(int j, int& update, bool final) {
            update += nSupports[j];
            if (final) {
                supports_ptr[j + 1] = update;
            }
        },
        total_supports);

    // Allocate storage for actual neighbor indices and distances
    Write<LO> supports_idx(total_supports, 0, "neighbor indices");
    Write<Real> radii_vals(total_supports, 0, "neighbor squared distances");

    // Fill neighbor indices 
    Kokkos::parallel_for("fill neighbor indices", numTargets, KOKKOS_LAMBDA(int i) {
        int offset = supports_ptr[i];  
        int count = 0;
        for (int j = 0; j < numSources; ++j) {
            double dx = target_points(i, 0) - source_points(j, 0);
            double dy = target_points(i, 1) - source_points(j, 1);
            double dz = target_points(i, 2) - source_points(j, 2);
            double dist_sq = dx * dx + dy * dy + dz * dz;

            if (dist_sq <= cutoffDistanceSq) {
                supports_idx[offset + count] = j;
                radii_vals[offset + count] = dist_sq;
                count++;
            }
        }
    });

    // Assign results
    results.supports_ptr = supports_ptr;
    results.supports_idx = supports_idx;
    results.radii2 = radii_vals;

    return results;
}

inline void test_interpolation_point_to_mesh(
    Mesh& mesh, double cutoffDistance,
    Reals source_values, Write<Real>& target_values,
    Reals source_coordinates, Reals target_coordinates,
    const SupportResults& support)
{
  using ExecSpace     = Kokkos::DefaultExecutionSpace;
  using TeamPolicy    = Kokkos::TeamPolicy<ExecSpace>;
  using MemberType    = TeamPolicy::member_type;
  using ScratchSpace  = typename MemberType::scratch_memory_space;

  using MatView = Kokkos::View<double**, Kokkos::LayoutRight, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using VecView = Kokkos::View<double*, ScratchSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

  const int dim = mesh.dim();
  const size_t num_targets = target_values.size();
  constexpr int basis_size = 10; // Quadratic: [1, x, y, z, x², y², z², xy, xz, yz]

  auto d_supports_ptr = support.supports_ptr;
  auto d_supports_idx = support.supports_idx;
  auto d_radii2       = support.radii2;

  Write<Real> temp_values(num_targets, "svd interpolated values");

  // Find max support
  LO max_support = 0;
  auto host_ptr = HostRead<LO>(d_supports_ptr);
  for (size_t i = 0; i < num_targets; ++i)
    max_support = std::max(max_support, host_ptr[i + 1] - host_ptr[i]);

  const int bytes_per_double = sizeof(double);
  const int scratch_bytes = bytes_per_double * (
    max_support * basis_size +   // A
    max_support +                // b
    basis_size +                 // x
    max_support +                // weights
    basis_size +                // sigma
    basis_size * basis_size +   // Vt
    max_support * max_support + // U
    max_support +               // work
    basis_size                  // Ut_b
  );

  TeamPolicy policy(num_targets, Kokkos::AUTO);
  policy = policy.set_scratch_size(1, Kokkos::PerTeam(scratch_bytes));

  Kokkos::parallel_for("MLS-SVD Interpolation", policy, KOKKOS_LAMBDA(const MemberType& team) {
    const int i = team.league_rank();
    const LO start = d_supports_ptr[i];
    const LO end   = d_supports_ptr[i + 1];
    const int nsupports = end - start;

    if (nsupports < basis_size) return;

    MatView A(team.team_scratch(1), nsupports, basis_size);
    VecView b(team.team_scratch(1), nsupports);
    VecView x(team.team_scratch(1), basis_size);
    VecView weights(team.team_scratch(1), nsupports);
    MatView U(team.team_scratch(1), nsupports, nsupports);
    MatView Vt(team.team_scratch(1), basis_size, basis_size);
    VecView sigma(team.team_scratch(1), basis_size);
    VecView work(team.team_scratch(1), nsupports);
    VecView Ut_b(team.team_scratch(1), basis_size);

    const double h = 0.5 * cutoffDistance;

    for (int j = 0; j < nsupports; ++j) {
      const int idx = d_supports_idx[start + j];
      double x_s = source_coordinates[idx * dim + 0];
      double y_s = source_coordinates[idx * dim + 1];
      double z_s = source_coordinates[idx * dim + 2];
      double r2  = d_radii2[start + j];
      double w   = exp(-r2 / (h * h));
      weights(j) = w;

      A(j, 0) = w * 1.0;
      A(j, 1) = w * x_s;
      A(j, 2) = w * y_s;
      A(j, 3) = w * z_s;
      A(j, 4) = w * x_s * x_s;
      A(j, 5) = w * y_s * y_s;
      A(j, 6) = w * z_s * z_s;
      A(j, 7) = w * x_s * y_s;
      A(j, 8) = w * x_s * z_s;
      A(j, 9) = w * y_s * z_s;

      b(j) = w * source_values[idx];
    }

    // SVD solve
    if (team.team_rank() == 0) {
      KokkosBatched::SerialSVD::invoke(
        KokkosBatched::SVD_USV_Tag(), A, U, sigma, Vt, work);
    }
    team.team_barrier();

    for (int k = 0; k < basis_size; ++k) {
      Ut_b(k) = 0.0;
      for (int j = 0; j < nsupports; ++j)
        Ut_b(k) += U(j, k) * b(j);
    }

    for (int irow = 0; irow < basis_size; ++irow) {
      double val = 0.0;
      for (int j = 0; j < basis_size; ++j)
        if (sigma(j) > 1e-10)
          val += Vt(j, irow) * (Ut_b(j) / sigma(j));
      x(irow) = val;
    }

    double xt = target_coordinates[i * dim + 0];
    double yt = target_coordinates[i * dim + 1];
    double zt = target_coordinates[i * dim + 2];

    Real interpolated =
      x(0) + x(1)*xt + x(2)*yt + x(3)*zt +
      x(4)*xt*xt + x(5)*yt*yt + x(6)*zt*zt +
      x(7)*xt*yt + x(8)*xt*zt + x(9)*yt*zt;

    temp_values[i] = interpolated;
  });

  target_values = temp_values;
}

    TEST_CASE("test point to mesh mls svd") {

    std::string filePath = "/lore/elahis/pcmsrelated/BOUT.dmp.bp";
    // Read the first time step (index 0)
    std::vector<double> Te_data = readTeData(filePath);

    auto lib = Library{};
    auto world = lib.world();

    auto mesh = build_box(world, OMEGA_H_SIMPLEX, 2.6, 2.6, 1.2, 14, 14, 20, false);

    // // Debug
    // auto host_coords = HostRead<Real>(mesh.coords());

    // double min_zz = std::numeric_limits<double>::max();
    // double max_zz = std::numeric_limits<double>::lowest();

    // for (int i = 0; i < mesh.nverts(); ++i) {
    //     double z = host_coords[i * 3 + 2];
    //     min_zz = std::min(min_zz, z);
    //     max_zz = std::max(max_zz, z);
    // }
    // std::cout << "Z range after shifting: [" << min_zz << ", " << max_zz << "]\n";

    Real cutoffDistance = 0.45;
    cutoffDistance = cutoffDistance * cutoffDistance; 

    const int dim = mesh.dim();
    const auto& target_coordinates = mesh.coords();
    const auto ntargets = mesh.nverts();

    auto source_data = generateSourcePointsWithTe(Te_data);
    size_t numSources = source_data.size();

    // // Debugging Print first few entries for debugging
    // std::cout << "First few source points (x, y, z, Te):" << std::endl;
    // for (size_t i = 0; i < std::min(numSources, size_t(10)); ++i) {
    //     std::cout << "Point " << i << ": ("
    //             << source_data[i].first[0] << ", "
    //             << source_data[i].first[1] << ", "
    //             << source_data[i].first[2] << ") "
    //             << " Te = " << source_data[i].second << std::endl;
    // }

    // Convert mesh coordinates to host format
    std::vector<std::vector<double>> host_target_data;
    auto host_target_coords = HostRead<Real>(target_coordinates);

    // Filtering target points
    for (size_t i = 0; i < ntargets; ++i) {
        double x = host_target_coords[i * 3];
        double y = host_target_coords[i * 3 + 1];
        double z = host_target_coords[i * 3 + 2];

        double r = std::sqrt(x * x + y * y);
        double theta = std::atan2(y, x); 

        // Normalize theta to [0, 2π] if needed
        if (theta < 0) theta += 2 * M_PI;

        if (r >= 1.4 && r <= 2.6 &&
            z >= 0.0 && z <= 1.2 &&
            theta >= 0.0 && theta <= 1.62) {
            host_target_data.push_back({x, y, z});
        }
    }

    // Update `ntargets` to reflect the filtered number of points
    const auto ntargets_filtered = host_target_data.size();

    // Define new variable for filtered target coordinates
    Write<Real> target_coordinates_filtered(ntargets_filtered * 3, 0.0, "target coordinates filtered");

    // Create a Kokkos View for host_target_data
    Kokkos::View<double**> host_target_data_view("host_target_data_view", ntargets_filtered, 3);

    // Copy host_target_data into the Kokkos View
    auto host_mirror = Kokkos::create_mirror_view(host_target_data_view);
    for (size_t i = 0; i < ntargets_filtered; ++i) {
        for (int d = 0; d < 3; ++d) {
            host_mirror(i, d) = host_target_data[i][d];
        }
    }
    Kokkos::deep_copy(host_target_data_view, host_mirror);

    // Fill target_coordinates_filtered using Kokkos parallel_for
    Kokkos::parallel_for("fill target_coordinates_filtered", ntargets_filtered, KOKKOS_LAMBDA(int i) {
        target_coordinates_filtered[i * 3] = host_target_data_view(i, 0);
        target_coordinates_filtered[i * 3 + 1] = host_target_data_view(i, 1);
        target_coordinates_filtered[i * 3 + 2] = host_target_data_view(i, 2);
    });

    // // Debugging
    // auto host_target_coords_filtered = HostRead<Real>(target_coordinates_filtered);

    // std::cout << "Filtered target coordinates (first 10):" << std::endl;
    // for (size_t i = 0; i < std::min(ntargets_filtered, size_t(10)); ++i) {
    //     std::cout << "Target " << i << ": ("
    //             << host_target_coords_filtered[i * 3] << ", "
    //             << host_target_coords_filtered[i * 3 + 1] << ", "
    //             << host_target_coords_filtered[i * 3 + 2] << ")" << std::endl;
    // }
    // // Debugging: Check min/max of r, theta, z
    // double min_r = std::numeric_limits<double>::max();
    // double max_r = std::numeric_limits<double>::lowest();
    // double min_theta = std::numeric_limits<double>::max();
    // double max_theta = std::numeric_limits<double>::lowest();
    // double min_z = std::numeric_limits<double>::max();
    // double max_z = std::numeric_limits<double>::lowest();
    // auto host_target_coords_filtered = HostRead<Real>(target_coordinates_filtered);
    // for (size_t i = 0; i < ntargets_filtered; ++i) {
    //     double x = host_target_coords_filtered[i * 3];
    //     double y = host_target_coords_filtered[i * 3 + 1];
    //     double z = host_target_coords_filtered[i * 3 + 2];
    //     double r = std::sqrt(x * x + y * y);
    //     double theta = std::atan2(y, x);
    //     if (theta < 0) theta += 2 * M_PI;
    //     min_r = std::min(min_r, r);
    //     max_r = std::max(max_r, r);
    //     min_theta = std::min(min_theta, theta);
    //     max_theta = std::max(max_theta, theta);
    //     min_z = std::min(min_z, z);
    //     max_z = std::max(max_z, z);
    // }
    // std::cout << "Filtered Target Range:\n";
    // std::cout << "  r     ∈ [" << min_r << ", " << max_r << "]\n";
    // std::cout << "  θ     ∈ [" << min_theta << ", " << max_theta << "] (rad)\n";
    // std::cout << "  z     ∈ [" << min_z << ", " << max_z << "]\n";

    // Extract {x, y, z} from source_data before calling findNeighbors
    std::vector<std::vector<double>> source_coordinates;
    for (const auto& p : source_data) {
        source_coordinates.push_back({p.first[0], p.first[1], p.first[2]}); 
    }

    // Now call findNeighbors with the correct data format
    SupportResults support = findNeighbors(source_coordinates, host_target_data, cutoffDistance);

    auto host_supports_ptr = HostRead<LO>(support.supports_ptr);
    auto host_supports_idx = HostRead<LO>(support.supports_idx);
    auto host_radii2 = HostRead<Real>(support.radii2);

    // // Debugging
    // std::cout << "First 10 targets' neighbor counts:" << std::endl;
    // for (size_t i = 0; i < std::min(ntargets_filtered, size_t(100)); ++i) {
    //     int num_neighbors = host_supports_ptr[i + 1] - host_supports_ptr[i];
    //     std::cout << "Target " << i << " has " << num_neighbors << " neighbors." << std::endl;
    // }

    // Ensure at least one neighbor is found
    CHECK(host_supports_ptr[ntargets_filtered] > 0);

    // Ensure supports_ptr size is correct
    CHECK(support.supports_ptr.size() == ntargets_filtered + 1);

    // Ensure supports_ptr is non-decreasing
    for (size_t i = 0; i < ntargets_filtered; ++i) {
        CHECK(host_supports_ptr[i] <= host_supports_ptr[i + 1]);
    }

    // Ensure first and last values are correct
    CHECK(host_supports_ptr[0] == 0);
    CHECK(host_supports_ptr[ntargets_filtered] == support.supports_idx.size());

    // Check that each target has at least one neighbor
    for (size_t i = 0; i < ntargets_filtered; ++i) {
        int num_neighbors = host_supports_ptr[i + 1] - host_supports_ptr[i];
        CHECK(num_neighbors >= 0);
    }

    // Check that found neighbors are actually within the cutoff distance
    const auto host_source_data = source_data;

    for (size_t i = 0; i < ntargets_filtered; ++i) {
        int start_idx = host_supports_ptr[i];
        int end_idx = host_supports_ptr[i + 1];

        for (int j = start_idx; j < end_idx; ++j) {
            int neighbor_index = host_supports_idx[j];

            // Extract only x, y, z from source data
            double dx = host_target_data[i][0] - host_source_data[neighbor_index].first[0];
            double dy = host_target_data[i][1] - host_source_data[neighbor_index].first[1];
            double dz = host_target_data[i][2] - host_source_data[neighbor_index].first[2];

            double dist_sq = dx * dx + dy * dy + dz * dz;
            CHECK(dist_sq <= cutoffDistance);
        }
    }

    // Check that every neighbor is counted only once per target
    for (size_t i = 0; i < ntargets_filtered; ++i) {
        int start_idx = host_supports_ptr[i];
        int end_idx = host_supports_ptr[i + 1];

        std::unordered_set<int> unique_neighbors;

        for (int j = start_idx; j < end_idx; ++j) {
            int neighbor_index = host_supports_idx[j];
            CHECK(unique_neighbors.find(neighbor_index) == unique_neighbors.end());
            unique_neighbors.insert(neighbor_index);
        }
    }

    // Convert source points to a GPU-compatible Kokkos View
    Kokkos::View<double**> source_view("source_view", numSources, 4);
    auto host_source_view = Kokkos::create_mirror_view(source_view);

    for (size_t i = 0; i < numSources; ++i) {
        for (int d = 0; d < 3; ++d) {
            host_source_view(i, d) = source_data[i].first[d];  // Extract {x, y, z}
        }
        host_source_view(i, 3) = source_data[i].second;  // Store Te value as the 4th column
    }

    // Copy data from host to device
    Kokkos::deep_copy(source_view, host_source_view);

    // Define a new variable for flattened source coordinates + Te
    Write<Real> source_data_flat(4 * numSources, 0.0, "source data flat");

    // Convert source points to flat `Reals` array for MLS
    Kokkos::parallel_for("convert source_data", numSources, KOKKOS_LAMBDA(const size_t i) {
        source_data_flat[i * 4] = source_view(i, 0);      // x
        source_data_flat[i * 4 + 1] = source_view(i, 1);  // y
        source_data_flat[i * 4 + 2] = source_view(i, 2);  // z
        source_data_flat[i * 4 + 3] = source_view(i, 3);  // Te
    });


    // Convert target points to a Kokkos View
    Kokkos::View<double**> target_view("target_view", ntargets_filtered, 3);
    auto host_target_view = Kokkos::create_mirror_view(target_view);

    for (size_t i = 0; i < ntargets_filtered; ++i) {
        for (int d = 0; d < dim; ++d) {
            host_target_view(i, d) = host_target_data[i][d];
        }
    }

    // Copy from host to device
    Kokkos::deep_copy(target_view, host_target_view);

    // SECTION("Linear MLS interpolation Te source to target") {
    //     int degree = 1;
    //     LO min_num_supports = 10;

    //     Write<Real> source_values(numSources, 0.0, "source Te values");
    //     Write<Real> target_values(ntargets_filtered, 0.0, "exact Te values");

    //     // Assign Te values from source points
    //     Kokkos::parallel_for("assign source Te values", numSources, KOKKOS_LAMBDA(const size_t i) {
    //         source_values[i] = source_view(i, 3);  // Te is stored in column 3
    //     });

    //     // Initialize target_values to 0 (they will be interpolated)
    //     Kokkos::parallel_for("initialize target Te values", ntargets_filtered, KOKKOS_LAMBDA(int i) {
    //         target_values[i] = 0.0;  // Placeholder, will be interpolated
    //     });

    //     // Flatten source coordinates
    //     Write<Real> source_coordinates_flat(numSources * 3, 0.0, "source coordinates");
    //     Kokkos::parallel_for("flatten source coordinates", numSources, KOKKOS_LAMBDA(const size_t i) {
    //         source_coordinates_flat[i * 3] = source_view(i, 0);  // x
    //         source_coordinates_flat[i * 3 + 1] = source_view(i, 1);  // y
    //         source_coordinates_flat[i * 3 + 2] = source_view(i, 2);  // z
    //     });

    //     // Flatten target coordinates
    //     Write<Real> target_coordinates_flat(ntargets_filtered * 3, 0.0, "target coordinates filtered");
    //     Kokkos::parallel_for("flatten target coordinates", ntargets_filtered, KOKKOS_LAMBDA(const size_t i) {
    //         target_coordinates_flat[i * 3] = target_view(i, 0);
    //         target_coordinates_flat[i * 3 + 1] = target_view(i, 1);
    //         target_coordinates_flat[i * 3 + 2] = target_view(i, 2);
    //     });

    //     test_interpolation_point_to_mesh(
    //         mesh, cutoffDistance, degree, min_num_supports,
    //         Reals(source_values), target_values,
    //         Reals(source_coordinates_flat), Reals(target_coordinates_flat)
    //     );

    //     // // Debug: Print first 10 target values
    //     // auto host_target_values_debug = HostRead<Real>(target_values);
    //     // std::cout << "\nTarget Te values after interpolation (first 10 targets)11111:\n";
    //     // for (size_t i = 0; i < std::min(ntargets_filtered, size_t(10)); ++i) {
    //     //     std::cout << "Target " << i << " | Te = " << host_target_values_debug[i] << "\n";
    //     // }
    // }

    SECTION("Two-way MLS interpolation: Source → Target → Source") {

        // First Interpolation: Source → Target
        Write<Real> source_values(numSources, 0.0, "source Te values");
        Write<Real> target_values(ntargets_filtered, 0.0, "interpolated Te at targets");

        // Assign `Te` values from source points
        Kokkos::parallel_for("assign source Te values", numSources, KOKKOS_LAMBDA(const size_t i) {
            source_values[i] = source_view(i, 3);  // Te is in column 3
        });

        // Flatten source and target coordinates
        Write<Real> source_coordinates_flat(numSources * 3, 0.0, "source coordinates");
        Write<Real> target_coordinates_flat(ntargets_filtered * 3, 0.0, "target coordinates filtered");

        Kokkos::parallel_for("flatten source coordinates", numSources, KOKKOS_LAMBDA(const size_t i) {
            source_coordinates_flat[i * 3] = source_view(i, 0);
            source_coordinates_flat[i * 3 + 1] = source_view(i, 1);
            source_coordinates_flat[i * 3 + 2] = source_view(i, 2);
        });

        Kokkos::parallel_for("flatten target coordinates", ntargets_filtered, KOKKOS_LAMBDA(const size_t i) {
            target_coordinates_flat[i * 3] = target_view(i, 0);
            target_coordinates_flat[i * 3 + 1] = target_view(i, 1);
            target_coordinates_flat[i * 3 + 2] = target_view(i, 2);
        });

        // Debug
        auto host_source_values_debug = HostRead<Real>(source_values);
        std::cout << "\nSource Te values before first interpolation (first 10 sources):\n";
        for (size_t i = 0; i < std::min(numSources, size_t(10)); ++i) {
            std::cout << "Source " << i << " | Te = " << host_source_values_debug[i] << "\n";
        }
        // Debug
        std::cout << "\nChecking neighbors per target (first 10 targets):\n";
        for (size_t i = 0; i < std::min(ntargets_filtered, size_t(10)); ++i) {
            int num_neighbors = host_supports_ptr[i + 1] - host_supports_ptr[i];
            std::cout << "Target " << i << " has " << num_neighbors << " neighbors.\n";
        }

        // Perform first MLS interpolation: `source → target`
        test_interpolation_point_to_mesh(
            mesh, cutoffDistance, Reals(source_values), target_values,
            Reals(source_coordinates_flat), Reals(target_coordinates_flat),
            support
        );

        // Debug: Print first 10 target values after first interpolation
        auto host_target_values_debug = HostRead<Real>(target_values);
        std::cout << "\nTarget Te values after first interpolation (first 10 targets):\n";
        for (size_t i = 0; i < std::min(ntargets_filtered, size_t(10)); ++i) {
            std::cout << "Target " << i << " | Te = " << host_target_values_debug[i] << "\n";
        }

        // Define new source points (previous target points) and their Te values
        Write<Real> new_source_values(ntargets_filtered, 0.0, "new source Te values");

        // Use Reals as an intermediate buffer before copying
        Reals temp_target_values = Reals(target_values);
        Kokkos::parallel_for("copy target_values -> new_source_values", ntargets_filtered, KOKKOS_LAMBDA(int i) {
            new_source_values[i] = temp_target_values[i]; 
        });

        // Define new source coordinates
        Write<Real> new_source_coordinates(ntargets_filtered * 3, 0.0, "new source coordinates");
        Reals temp_target_coordinates = Reals(target_coordinates_flat);
        Kokkos::parallel_for("copy target_coordinates -> new_source_coordinates", ntargets_filtered, KOKKOS_LAMBDA(int i) {
            new_source_coordinates[i * 3] = temp_target_coordinates[i * 3];
            new_source_coordinates[i * 3 + 1] = temp_target_coordinates[i * 3 + 1];
            new_source_coordinates[i * 3 + 2] = temp_target_coordinates[i * 3 + 2];
        });

        // Define new target points (previous source points), initialize Te to zero
        Write<Real> recovered_source_values(numSources, 0.0, "recovered Te at sources");
        Write<Real> new_target_coordinates(numSources * 3, 0.0, "new target coordinates");

        Reals temp_source_coordinates = Reals(source_coordinates_flat);
        Kokkos::parallel_for("copy source_coordinates -> new_target_coordinates", numSources, KOKKOS_LAMBDA(int i) {
            new_target_coordinates[i * 3] = temp_source_coordinates[i * 3];
            new_target_coordinates[i * 3 + 1] = temp_source_coordinates[i * 3 + 1];
            new_target_coordinates[i * 3 + 2] = temp_source_coordinates[i * 3 + 2];
        });

        // Convert `Reals` to host `std::vector<std::vector<double>>` format
        std::vector<std::vector<double>> host_new_source_data;
        std::vector<std::vector<double>> host_new_target_data;

        auto host_new_source_coords = HostRead<Real>(new_source_coordinates);
        auto host_new_target_coords = HostRead<Real>(new_target_coordinates);

        // Fill new source data
        for (size_t i = 0; i < ntargets_filtered; ++i) {
            host_new_source_data.push_back({
                host_new_source_coords[i * 3],
                host_new_source_coords[i * 3 + 1],
                host_new_source_coords[i * 3 + 2]
            });
        }

        // Fill new target data
        for (size_t i = 0; i < numSources; ++i) {
            host_new_target_data.push_back({
                host_new_target_coords[i * 3],
                host_new_target_coords[i * 3 + 1],
                host_new_target_coords[i * 3 + 2]
            });
        }

        // // Debug: Print Te values before second interpolation
        // std::cout << "\nTarget Te values before second interpolation (first 10 targets):\n";
        // for (size_t i = 0; i < std::min(ntargets_filtered, size_t(10)); ++i) {
        //     std::cout << "Target " << i << " | Te = " << host_target_values_debug[i] << "\n";
        // }

        Real second_cutoffDistance = cutoffDistance * 1.5;

        // Find new neighbors (target → source)
        SupportResults support_target_to_source = findNeighbors(
            host_new_source_data, 
            host_new_target_data,
            second_cutoffDistance
        );

        // Ensure neighbors are found
        CHECK(support_target_to_source.supports_idx.size() > 0);

        // Debug: Check neighbors per new target (original sources)
        auto host_supports_ptr_new = HostRead<LO>(support_target_to_source.supports_ptr);
        std::cout << "\nChecking neighbors per source (first 10 sources):\n";
        for (size_t i = 0; i < std::min(numSources, size_t(10)); ++i) {
            int num_neighbors = host_supports_ptr_new[i + 1] - host_supports_ptr_new[i];
            std::cout << "Source " << i << " has " << num_neighbors << " neighbors.\n";
        } 

        // Perform second MLS interpolation: `target → source`
        test_interpolation_point_to_mesh(
            mesh, cutoffDistance, Reals(new_source_values), recovered_source_values,
            Reals(new_source_coordinates), Reals(new_target_coordinates),
            support_target_to_source
        );

        // Compare original vs recovered values
        auto host_original_source_values = HostRead<Real>(source_values);
        auto host_recovered_source_values = HostRead<Real>(recovered_source_values);
        std::cout << "\nComparing Original vs. Recovered Te values (first 10 sources):\n";
        for (size_t i = 0; i < std::min(numSources, size_t(10)); ++i) {
            std::cout << "Source " << i
                    << " | Original Te = " << host_original_source_values[i]
                    << " | Recovered Te = " << host_recovered_source_values[i] << "\n";
        }

        // Read source data from host memory
        auto host_source_coords = HostRead<Real>(source_coordinates_flat);
        auto host_original_Te = HostRead<Real>(source_values);
        auto host_recovered_Te = HostRead<Real>(recovered_source_values);

        // Read target data from host memory
        auto host_target_coords = HostRead<Real>(target_coordinates_flat);
        auto host_target_Te = HostRead<Real>(target_values);  // Target values after first interpolation

        // Save source data to source_points.txt
        std::ofstream sourceFile("/lore/elahis/pcmsrelated/source_points.txt");
        sourceFile << "X_source Y_source Z_source Te_Original Te_Recovered\n";
        for (size_t i = 0; i < numSources; ++i) {
            double x = host_source_coords[i * 3];
            double y = host_source_coords[i * 3 + 1];
            double z = host_source_coords[i * 3 + 2];
            double te_orig = host_original_Te[i];
            double te_recovered = host_recovered_Te[i];
            sourceFile << x << " " << y << " " << z << " " << te_orig << " " << te_recovered << "\n";
        }
        sourceFile.close();

        // Save target data to target_points.txt
        std::ofstream targetFile("/lore/elahis/pcmsrelated/target_points.txt");
        targetFile << "X_target Y_target Z_target Te_Target_Interpolated\n";
        for (size_t i = 0; i < ntargets_filtered; ++i) {
            double x = host_target_coords[i * 3];
            double y = host_target_coords[i * 3 + 1];
            double z = host_target_coords[i * 3 + 2];
            double te = host_target_Te[i];
            targetFile << x << " " << y << " " << z << " " << te << "\n";
        }
        targetFile.close();

        //Debug
        double min_orig = std::numeric_limits<double>::max();
        double max_orig = std::numeric_limits<double>::lowest();
        double min_recv = std::numeric_limits<double>::max();
        double max_recv = std::numeric_limits<double>::lowest();
        for (size_t i = 0; i < numSources; ++i) {
            double te_orig = host_original_Te[i];
            double te_recv = host_recovered_Te[i];
            if (te_orig < min_orig) min_orig = te_orig;
            if (te_orig > max_orig) max_orig = te_orig;
            if (te_recv < min_recv) min_recv = te_recv;
            if (te_recv > max_recv) max_recv = te_recv;
        }
        double min_target = std::numeric_limits<double>::max();
        double max_target = std::numeric_limits<double>::lowest();
        for (size_t i = 0; i < ntargets_filtered; ++i) {
            double te = host_target_Te[i];
            if (te < min_target) min_target = te;
            if (te > max_target) max_target = te;
        }
        std::cout << "==== Temperature Range Summary ====\n";
        std::cout << "Te_Original  : min = " << min_orig << ", max = " << max_orig << "\n";
        std::cout << "Te_Recovered : min = " << min_recv << ", max = " << max_recv << "\n";
        std::cout << "Te_Target    : min = " << min_target << ", max = " << max_target << "\n";

        std::cout << "Saved source data to source_points.txt\n";
        std::cout << "Saved target data to target_points.txt\n";
    }
}