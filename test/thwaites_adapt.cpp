#include <iostream>

#include "Omega_h_adapt.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_build.hpp"
#include "Omega_h_class.hpp"
#include "Omega_h_compare.hpp"
#include "Omega_h_for.hpp"
#include "Omega_h_shape.hpp"
#include "Omega_h_timer.hpp"
#include "Omega_h_recover.hpp" //project_by_fit
#include <Omega_h_file.hpp> //Omega_h::binary
#include <Omega_h_atomics.hpp> //Omega_h::atomic_fetch_add
#include <sstream> //ostringstream
#include <iomanip> //precision
#include <Omega_h_dbg.hpp>
#include <Omega_h_for.hpp>

#include <pcms/interpolator/mls_interpolation.hpp> //mls_interpolation

#include "thwaites_errorEstimator.hpp"

//detect floating point exceptions
#include <fenv.h>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;

using namespace Omega_h;


//DEBUG {
void print_patches(Omega_h::Graph& patches) {
  auto offsets = HostRead(patches.a2ab);
  auto values = HostRead(patches.ab2b);
  std::cout << "num patches " << patches.nnodes() << "\n";
  for(int patch=0; patch<patches.nnodes(); patch++) {
    std::cout << "patch " << patch << " patchElms ";
    for (auto valIdx = offsets[patch]; valIdx < offsets[patch + 1]; ++valIdx) {
      auto patchElm = values[valIdx];
      std::cout << patchElm << " ";
    }
    std::cout << "\n";
  }
}

void print_patch_vals(const Omega_h::Graph& patches_d, Reals vtxCoords_d,
    Reals elmSrcVals_d, Reals elmCentroids, size_t vtx) {
  const auto meshDim = 2;
  const auto offsets = HostRead(patches_d.a2ab);
  const auto values = HostRead(patches_d.ab2b);
  const auto vtxCoords = HostRead(vtxCoords_d);
  const auto elmSrcVals = HostRead(elmSrcVals_d);
  std::cout << std::setprecision (15);
  for(int patch=0; patch<patches_d.nnodes(); patch++) {
    if(patch == vtx) {
      std::cout << "vtxCoords[" << patch << "] " << vtxCoords[patch*meshDim] << " " << vtxCoords[patch*meshDim+1] << "\n";
      std::cout << "<elementIdx> <centroid x> <centroid y> <source_value>\n";
      for (auto valIdx = offsets[patch]; valIdx < offsets[patch + 1]; ++valIdx) {
        auto patchElm = values[valIdx];
        std::cout <<  patchElm << " "
                  << elmCentroids[patchElm*meshDim] << " "
                  << elmCentroids[patchElm*meshDim+1] << " "
                  << elmSrcVals[patchElm] << "\n";
      }
      std::cout << "\n";
    }
  }
}
//END DEBUG }

void setupFieldTransfer(AdaptOpts& opts) {
  opts.xfer_opts.type_map["solution_1"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["solution_2"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["ice_thickness"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["surface_height"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["bed_topography"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["mu_log"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["prescribed_velocity_1"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["prescribed_velocity_2"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["observed_surface_velocity_rms"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["observed_surface_velocity_1"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["observed_surface_velocity_2"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["stiffening_factor_log"] = OMEGA_H_LINEAR_INTERP;
  const int numLayers = 11;
  for(int i=1; i<=numLayers; i++) {
    std::stringstream ss;
    ss << "temperature_" << std::setfill('0') << std::setw(2) << i;
    opts.xfer_opts.type_map[ss.str()] = OMEGA_H_LINEAR_INTERP;
  }
}

Reals getMagVelocity(Mesh& mesh) {
  const auto v_x = mesh.get_array<Real>(VERT, "solution_1");
  const auto v_y = mesh.get_array<Real>(VERT, "solution_2");
  Write<Real> magVelocity(mesh.nverts());
  Kokkos::parallel_for("magVelocity", mesh.nverts(),
      KOKKOS_LAMBDA(int vtx) {
         const auto x = v_x[vtx];
         const auto y = v_y[vtx];
         magVelocity[vtx] = Kokkos::sqrt(x*x+y*y);
      });
  return read(magVelocity);
}

Reals recoverLinearStrain(Mesh& mesh, Reals effectiveStrain) {
  return project_by_fit(&mesh, effectiveStrain);
}

Reals getElementCentroids(Mesh& mesh, Reals vtxCoords) {
  assert(mesh.dim() == 2);
  Write<Real> centroids(
    mesh.dim() * mesh.nelems(), 0, "stores coordinates of cell centroid of each element");

  const auto& faces2nodes = mesh.ask_down(FACE, VERT).ab2b;

  Kokkos::parallel_for(
    "calculate the centroid in each tri element", mesh.nelems(),
    OMEGA_H_LAMBDA(const LO id) {
      const auto current_el_verts = gather_verts<3>(faces2nodes, id);
      const Omega_h::Few<Omega_h::Vector<2>, 3> current_el_vert_coords =
        gather_vectors<3, 2>(vtxCoords, current_el_verts);
      auto centroid = average(current_el_vert_coords);
      int index = 2 * id;
      centroids[index] = centroid[0];
      centroids[index + 1] = centroid[1];
    });

  return read(centroids);
}

Reals min_max_normalization_coordinates(const Reals& coordinates, int dim = 2) {
  int num_points = coordinates.size() / dim;

  int coords_size = coordinates.size();

  Write<Real> x_coordinates(num_points, 0, "x coordinates");

  Write<Real> y_coordinates(num_points, 0, "y coordinates");

  parallel_for(
      "separates x and y coordinates", num_points, KOKKOS_LAMBDA(const int id) {
      int index = id * dim;
      x_coordinates[id] = coordinates[index];
      y_coordinates[id] = coordinates[index + 1];
      });


  const auto min_x = Omega_h::get_min(read(x_coordinates));
  const auto min_y = Omega_h::get_min(read(y_coordinates));

  const auto max_x = Omega_h::get_max(read(x_coordinates));
  const auto max_y = Omega_h::get_max(read(y_coordinates));

  const auto del_x = max_x - min_x;
  const auto del_y = max_y - min_y;

  Write<Real> normalized_coordinates(coords_size, 0,
      "stores scaled coordinates");

  parallel_for(
      "scale coordinates", num_points, KOKKOS_LAMBDA(const int id) {
      int index = id * dim;
      normalized_coordinates[index] = (x_coordinates[id] - min_x) / del_x;
      normalized_coordinates[index + 1] = (y_coordinates[id] - min_y) / del_y;
      });

  return read(normalized_coordinates);
}

Reals recoverMagVelocityPCMS(Mesh& mesh, Reals magVelocity, size_t degree, size_t minPatchSize) {
  assert(degree > 0 && degree < 4);
  assert(minPatchSize > 0 && minPatchSize < 20);
  assert(magVelocity.size() == mesh.nverts());
  assert(mesh.nverts() > 0);
  const auto dim = mesh.dim();
  const Real tolerance = 5e-4;

  const auto target_coordinates = mesh.coords();

  const auto source_coordinates = mesh.coords();

  const auto patches = mesh.get_vtx_patches(minPatchSize, VERT);
  Omega_h::Write<Real> ignored(patches.ab2b.size(), 1);
  SupportResults support{patches.a2ab,patches.ab2b,ignored};

  const auto source_values = magVelocity;
  auto recovered =
    pcms::mls_interpolation (source_values, source_coordinates, target_coordinates,
        support, dim, degree, support.radii2, pcms::RadialBasisFunction::NO_OP);

  return recovered;
}

int main(int argc, char** argv) {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);  // Enable all floating point exceptions but FE_INEXACT
  auto lib = Library(&argc, &argv);
  if( argc != 4 ) {
    fprintf(stderr, "Usage: %s inputMesh.osh outputMeshPrefix minPatchSize\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  auto world = lib.world();
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[1], world, &mesh);
  const auto prefix = std::string(argv[2]);
  const auto minPatchSize = std::atoi(argv[3]);
  std::cout << "input mesh: " << argv[1] << " outputMeshPrefix: " << prefix << " minPatchSize: " << minPatchSize << "\n";
  const auto outname = prefix + "_recoverMagVelocity_deg2";

  const auto magVelocity = getMagVelocity(mesh);
  mesh.add_tag<Real>(VERT, "magVelocity", 1, magVelocity);

  { //write vtk
  const std::string vtkFileName = "beforeMls" + outname + ".vtk";
  Omega_h::vtk::write_parallel(vtkFileName, &mesh, 2);
  }

  const auto recoveredFieldDegree = 2;
  auto recoveredMagVelocityPCMS = recoverMagVelocityPCMS(mesh, magVelocity, recoveredFieldDegree, minPatchSize);
  mesh.add_tag<Real>(VERT, "recoveredMagVelocityPCMS", 1, recoveredMagVelocityPCMS);

  { //write vtk
  const std::string vtkFileName = "beforeAdapt" + outname + ".vtk";
  Omega_h::vtk::write_parallel(vtkFileName, &mesh, 2);
  }

  return 0;
}
