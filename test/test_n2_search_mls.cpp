#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <pcms/interpolator/mls_interpolation.hpp>
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

using namespace Omega_h;
using namespace pcms;

// Convert cylindrical coordinates (R, Z, θ) to Cartesian (x, y, z)
std::vector<std::vector<double>> generateSourcePoints() {
    const int radial = 20, toroidal = 10, poloidal = 20;
    double R_min = 1.4, R_max = 2.6, Z_min = -0.6, Z_max = 0.6;
    double dy = 0.18;
    double dR = (R_max - R_min) / (radial - 1);
    double dZ = (Z_max - Z_min) / (poloidal - 1);

    std::vector<std::vector<double>> source_points;
    for (int i = 0; i < radial; ++i) {
        double R = R_min + i * dR;
        for (int j = 0; j < poloidal; ++j) {
            double Z_val = Z_min + j * dZ;
            for (int k = 0; k < toroidal; ++k) {
                double theta = k * dy;
                source_points.push_back({R * cos(theta), R * sin(theta), Z_val});
            }
        }
    }
    return source_points;
}

KOKKOS_INLINE_FUNCTION
double evaluatePolynomial(const Coord& p, int degree) {
  auto x = p.x;
  auto y = p.y;
  auto z = p.z; 

  if (degree == 0) {
    return 3.0;
  } else if (degree == 1) {
    return x + y + z; 
  } else if (degree == 2) {
    return (x * x) + (y * y) + (z * z); 
  } else if (degree == 3) {
    return (x * x * x) + (y * y * y) + (z * z * z); 
  } else {
    printf("No polynomials with degree = %d\n", degree);
    return 0.0;
  }
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

inline void test_interpolation_point_to_mesh(Mesh& mesh, Real cutoffDistance, int degree, LO min_num_supports,
          Reals source_values, Reals exact_target_values,
          Reals source_coordinates, Reals target_coordinates)
{
  int dim = mesh.dim();
  assert(dim==3);
  Real tolerance = 0.05;

  std::vector<RadialBasisFunction> rbf_types = {
    RadialBasisFunction::RBF_GAUSSIAN, RadialBasisFunction::RBF_C4,
    RadialBasisFunction::RBF_CONST
  };

  // Ensure `source_coordinates` and `target_coordinates` are valid
  REQUIRE(source_coordinates.size() % dim == 0);
  REQUIRE(target_coordinates.size() % dim == 0);

  size_t num_sources = source_coordinates.size() / dim;
  size_t num_targets = target_coordinates.size() / dim;

  // Convert `Reals` to host std::vector<std::vector<double>> format
  std::vector<std::vector<double>> host_source_data(num_sources, std::vector<double>(dim));
  std::vector<std::vector<double>> host_target_data(num_targets, std::vector<double>(dim));

  // Read `Reals` into `HostRead<Real>` once to avoid redundant memory copies
  auto host_source_coords = HostRead<Real>(source_coordinates);
  auto host_target_coords = HostRead<Real>(target_coordinates);

  for (size_t i = 0; i < num_sources; ++i) {
      for (int d = 0; d < dim; ++d) {
          host_source_data[i][d] = host_source_coords[i * dim + d];
      }
  }
  for (size_t i = 0; i < num_targets; ++i) {
      for (int d = 0; d < dim; ++d) {
          host_target_data[i][d] = host_target_coords[i * dim + d];
      }
  }

  // Perform neighbor search
  SupportResults support = findNeighbors(host_source_data, host_target_data, cutoffDistance);

  // Ensure neighbor search results are valid
  CHECK(support.supports_ptr.size() == num_targets + 1);
  CHECK(support.supports_idx.size() > 0);  // Ensure we found at least one neighbor

  for (const auto& rbf : rbf_types) {
    auto approx_target_values =
      mls_interpolation(source_values, source_coordinates, target_coordinates,
                        support, dim, degree, support.radii2, rbf);

    // Read exact and interpolated values from device memory
    auto host_approx_target_values = HostRead<Real>(approx_target_values);
    auto host_exact_target_values = HostRead<Real>(exact_target_values);

    int m = exact_target_values.size();
    int n = approx_target_values.size();

    REQUIRE(m == n);

    // // Debugging
    // for (size_t i = 0; i < 10; ++i) { // Check first 10 values
    //     std::cout << "Target " << i << " exact: " << host_exact_target_values[i] 
    //             << ", approx: " << host_approx_target_values[i] << "\n";
    // }

    // double min_x = std::numeric_limits<double>::max();
    // double max_x = std::numeric_limits<double>::lowest();
    // double min_y = min_x, min_z = min_x;
    // double max_y = max_x, max_z = max_x;

    // for (size_t i = 0; i < num_targets; ++i) {
    //     min_x = std::min(min_x, host_target_data[i][0]);
    //     min_y = std::min(min_y, host_target_data[i][1]);
    //     min_z = std::min(min_z, host_target_data[i][2]);
    //     max_x = std::max(max_x, host_target_data[i][0]);
    //     max_y = std::max(max_y, host_target_data[i][1]);
    //     max_z = std::max(max_z, host_target_data[i][2]);
    // }

    // std::cout << "Filtered target bounding box: "
    //         << "[" << min_x << ", " << max_x << "] x "
    //         << "[" << min_y << ", " << max_y << "] x "
    //         << "[" << min_z << ", " << max_z << "]\n";

    auto host_supports_ptr = HostRead<LO>(support.supports_ptr);

    for (size_t i = 0; i < m; ++i) {
        int num_neighbors = host_supports_ptr[i + 1] - host_supports_ptr[i];
        if (num_neighbors >= min_num_supports) {
            CHECK_THAT(
                host_exact_target_values[i],
                Catch::Matchers::WithinAbs(host_approx_target_values[i], tolerance));
        }
    }
  }
}

TEST_CASE("point to mesh mls interp brute force search") {
    auto lib = Library{};
    auto world = lib.world();

    auto mesh = build_box(world, OMEGA_H_SIMPLEX, 5.2, 5.2, 1.2, 28, 28, 14, false);

    Real cutoffDistance = 0.2;
    cutoffDistance = cutoffDistance * cutoffDistance;  // Squared for efficiency

    const int dim = mesh.dim();
    const auto& target_coordinates = mesh.coords();
    const auto ntargets = mesh.nverts();

    // Generate source points
    auto source_data = generateSourcePoints();
    size_t numSources = source_data.size();

    // // Convert mesh coordinates to host format
    // std::vector<std::vector<double>> host_target_data(ntargets, std::vector<double>(dim));
    // auto host_target_coords = HostRead<Real>(target_coordinates);

    // for (size_t i = 0; i < ntargets; ++i) {
    //     for (int d = 0; d < dim; ++d) {
    //         host_target_data[i][d] = host_target_coords[i * dim + d];
    //     }
    // }

    // Convert mesh coordinates to host format
    std::vector<std::vector<double>> host_target_data;
    auto host_target_coords = HostRead<Real>(target_coordinates);

    // Compute bounding box of source points
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    double min_z = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = std::numeric_limits<double>::lowest();
    double max_z = std::numeric_limits<double>::lowest();

    for (const auto& point : source_data) {
        min_x = std::min(min_x, point[0]);
        min_y = std::min(min_y, point[1]);
        min_z = std::min(min_z, point[2]);
        max_x = std::max(max_x, point[0]);
        max_y = std::max(max_y, point[1]);
        max_z = std::max(max_z, point[2]);
    }

    // Now filter target points using this bounding box
    for (size_t i = 0; i < ntargets; ++i) {
        double x = host_target_coords[i * 3];
        double y = host_target_coords[i * 3 + 1];
        double z = host_target_coords[i * 3 + 2];

        double r = sqrt(x * x + y * y);  // Compute radial distance

        if (x >= min_x && x <= max_x &&
            y >= min_y && y <= max_y &&
            z >= min_z && z <= max_z &&
            r >= 1.4 && r <= 2.6) {  
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

    // Run the neighbor search function
    SupportResults support = findNeighbors(source_data, host_target_data, cutoffDistance);

    auto host_supports_ptr = HostRead<LO>(support.supports_ptr);
    auto host_supports_idx = HostRead<LO>(support.supports_idx);
    auto host_radii2 = HostRead<Real>(support.radii2);

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

            double dx = host_target_data[i][0] - host_source_data[neighbor_index][0];
            double dy = host_target_data[i][1] - host_source_data[neighbor_index][1];
            double dz = host_target_data[i][2] - host_source_data[neighbor_index][2];

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
    Kokkos::View<double**> source_view("source_view", numSources, 3);
    auto host_source_view = Kokkos::create_mirror_view(source_view);

    for (size_t i = 0; i < numSources; ++i) {
        for (int d = 0; d < 3; ++d) {
            host_source_view(i, d) = source_data[i][d];
        }
    }

    // Copy data from host to device
    Kokkos::deep_copy(source_view, host_source_view);

    // Convert source points to flat `Reals` array for MLS
    Write<Real> source_coordinates(dim * numSources, 0.0, "source coordinates");

    Kokkos::parallel_for("convert source_data", numSources, KOKKOS_LAMBDA(const size_t i) {
        source_coordinates[i * dim] = source_view(i, 0);
        source_coordinates[i * dim + 1] = source_view(i, 1);
        source_coordinates[i * dim + 2] = source_view(i, 2);
    });

    // Convert target points to a Kokkos View
    Kokkos::View<double**> target_view("target_view", ntargets_filtered, 3);
    auto host_target_view = Kokkos::create_mirror_view(target_view);

    for (size_t i = 0; i < ntargets_filtered; ++i) {
        for (int d = 0; d < dim; ++d) {
            // host_target_view(i, d) = host_target_coords[i * dim + d];
            host_target_view(i, d) = host_target_data[i][d];
        }
    }

    // Copy from host to device
    Kokkos::deep_copy(target_view, host_target_view);

    SECTION("MLS interpolation with isolated target") {
        // auto host_source_data = generateSourcePoints(gridData);
        auto host_source_data = generateSourcePoints();
        std::vector<std::vector<double>> isolated_target_data = { {100, 100, 100} }; // Far from all sources

        SupportResults support = findNeighbors(host_source_data, isolated_target_data, cutoffDistance);

        CHECK(support.supports_idx.size() == 0); // No neighbors
        CHECK(support.supports_ptr.size() == 2); // Should still be valid
    }

    SECTION("MLS degree 1, poly 0") {
        int degree = 1;
        LO min_num_supports = 10;

        Write<Real> source_values(numSources, 0.0, "source values");
        Write<Real> exact_target_values(ntargets_filtered, 0.0, "exact target values");

        Kokkos::parallel_for("assign source values", numSources, KOKKOS_LAMBDA(const size_t i) {
            Coord p;
            p.x = source_view(i, 0);
            p.y = source_view(i, 1);
            p.z = source_view(i, 2);
            source_values[i] = evaluatePolynomial(p, degree - 1);
        });

        Kokkos::parallel_for("compute exact target values", ntargets_filtered, KOKKOS_LAMBDA(int i) {
            Coord p;
            p.x = target_view(i, 0);
            p.y = target_view(i, 1);
            p.z = target_view(i, 2);
            exact_target_values[i] = evaluatePolynomial(p, degree - 1);
        });

        test_interpolation_point_to_mesh(mesh, cutoffDistance, degree, min_num_supports,
                                         Reals(source_values), Reals(exact_target_values),
                                         Reals(source_coordinates), Reals(target_coordinates_filtered));
    }

    SECTION("MLS degree 1, poly 1") {
        int degree = 1;
        LO min_num_supports = 10;

        Write<Real> source_values(numSources, 0.0, "source values");
        Write<Real> exact_target_values(ntargets_filtered, 0.0, "exact target values");

        Kokkos::parallel_for("assign source values", numSources, KOKKOS_LAMBDA(const size_t i) {
            Coord p;
            p.x = source_view(i, 0);
            p.y = source_view(i, 1);
            p.z = source_view(i, 2);
            source_values[i] = evaluatePolynomial(p, degree);
        });

        Kokkos::parallel_for("compute exact target values", ntargets_filtered, KOKKOS_LAMBDA(int i) {
            Coord p;
            p.x = target_view(i, 0);
            p.y = target_view(i, 1);
            p.z = target_view(i, 2);
            exact_target_values[i] = evaluatePolynomial(p, degree);
        });

        test_interpolation_point_to_mesh(mesh, cutoffDistance, degree, min_num_supports,
                                         Reals(source_values), Reals(exact_target_values),
                                         Reals(source_coordinates), Reals(target_coordinates_filtered));
    }

    SECTION("MLS degree 2, poly 0") {
        int degree = 2;
        LO min_num_supports = 16;

        Write<Real> source_values(numSources, 0.0, "source values");
        Write<Real> exact_target_values(ntargets_filtered, 0.0, "exact target values");

        Kokkos::parallel_for("assign source values", numSources, KOKKOS_LAMBDA(const size_t i) {
            Coord p;
            p.x = source_view(i, 0);
            p.y = source_view(i, 1);
            p.z = source_view(i, 2);
            source_values[i] = evaluatePolynomial(p, degree - 2);  // Polynomial degree 0
        });

        Kokkos::parallel_for("compute exact target values", ntargets_filtered, KOKKOS_LAMBDA(int i) {
            Coord p;
            p.x = target_view(i, 0);
            p.y = target_view(i, 1);
            p.z = target_view(i, 2);
            exact_target_values[i] = evaluatePolynomial(p, degree - 2);  // Polynomial degree 0
        });

        test_interpolation_point_to_mesh(mesh, cutoffDistance, degree, min_num_supports,
                                        Reals(source_values), Reals(exact_target_values),
                                        Reals(source_coordinates), Reals(target_coordinates_filtered));
    }

    SECTION("MLS degree 2, poly 1") {
        int degree = 2;
        LO min_num_supports = 16;

        Write<Real> source_values(numSources, 0.0, "source values");
        Write<Real> exact_target_values(ntargets_filtered, 0.0, "exact target values");

        Kokkos::parallel_for("assign source values", numSources, KOKKOS_LAMBDA(const size_t i) {
            Coord p;
            p.x = source_view(i, 0);
            p.y = source_view(i, 1);
            p.z = source_view(i, 2);
            source_values[i] = evaluatePolynomial(p, degree - 1);  // Polynomial degree 1
        });

        Kokkos::parallel_for("compute exact target values", ntargets_filtered, KOKKOS_LAMBDA(int i) {
            Coord p;
            p.x = target_view(i, 0);
            p.y = target_view(i, 1);
            p.z = target_view(i, 2);
            exact_target_values[i] = evaluatePolynomial(p, degree - 1);  // Polynomial degree 1
        });

        test_interpolation_point_to_mesh(mesh, cutoffDistance, degree, min_num_supports,
                                        Reals(source_values), Reals(exact_target_values),
                                        Reals(source_coordinates), Reals(target_coordinates_filtered));
    }

    SECTION("MLS degree 2, poly 2") {
        int degree = 2;
        LO min_num_supports = 16;

        Write<Real> source_values(numSources, 0.0, "source values");
        Write<Real> exact_target_values(ntargets_filtered, 0.0, "exact target values");

        Kokkos::parallel_for("assign source values", numSources, KOKKOS_LAMBDA(const size_t i) {
            Coord p;
            p.x = source_view(i, 0);
            p.y = source_view(i, 1);
            p.z = source_view(i, 2);
            source_values[i] = evaluatePolynomial(p, degree);  // Polynomial degree 2
        });

        Kokkos::parallel_for("compute exact target values", ntargets_filtered, KOKKOS_LAMBDA(int i) {
            Coord p;
            p.x = target_view(i, 0);
            p.y = target_view(i, 1);
            p.z = target_view(i, 2);
            exact_target_values[i] = evaluatePolynomial(p, degree);  // Polynomial degree 2
        });

        test_interpolation_point_to_mesh(mesh, cutoffDistance, degree, min_num_supports,
                                        Reals(source_values), Reals(exact_target_values),
                                        Reals(source_coordinates), Reals(target_coordinates_filtered));
    }
}