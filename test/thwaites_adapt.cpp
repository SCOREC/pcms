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

#include "thwaites_errorEstimator.hpp"

//detect floating point exceptions
#include <fenv.h>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;

using namespace Omega_h;

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
  const int numLayers = 11;
  for(int i=1; i<=numLayers; i++) {
    std::stringstream ss;
    ss << "temperature_" << std::setfill('0') << std::setw(2) << i;
    opts.xfer_opts.type_map[ss.str()] = OMEGA_H_LINEAR_INTERP;
  }
}

/**
 * retrieve the effective strain rate from the mesh
 *
 * despite the name being 'effective_stress' it is the effective strain:
 * the frobenius norm of  the strain tensor
  */
Reals getEffectiveStrainRate(Mesh& mesh) {
  return mesh.get_array<Real>(2, "effective_stress");
}

Reals recoverLinearStrain(Mesh& mesh, Reals effectiveStrain) {
  return project_by_fit(&mesh, effectiveStrain);
}

Reals computeError(Mesh& mesh, Reals effectiveStrain, Reals recoveredStrain) {
  return Reals(mesh.nelems());
}

template <typename ShapeField>
void setFieldAtVertices(Omega_h::Mesh &mesh, Reals recoveredStrain, ShapeField field) {
  const auto MeshDim = mesh.dim();
  auto setFieldAtVertices = KOKKOS_LAMBDA(const int &vtx) {
    field(0, 0, vtx, MeshField::Vertex) = recoveredStrain[vtx];
  };
  MeshField::parallel_for(ExecutionSpace(), {0}, {mesh.nverts()},
                          setFieldAtVertices, "setFieldAtVertices");
}

void printTriCount(Mesh* mesh) {
  const auto nTri = mesh->nglobal_ents(2);
  if (!mesh->comm()->rank())
    std::cout << "nTri: " << nTri << "\n";
}

int main(int argc, char** argv) {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);  // Enable all floating point exceptions but FE_INEXACT
  auto lib = Library(&argc, &argv);
  if( argc != 4 ) {
    fprintf(stderr, "Usage: %s inputMesh.osh outputMeshPrefix adaptRatio\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  auto world = lib.world();
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[1], world, &mesh);
  const auto prefix = std::string(argv[2]);
  //adaptRatio = 0.1 is used in scorec/core:cws/sprThwaites test/spr_test.cc
  const MeshField::Real adaptRatio = std::stof(argv[3]);
  auto max_size = Real(16); //roughly maintains the boundary shape 
                            //in the downstream parts of the domain
                            //where there is significant coarsening
  auto min_size = Real(1);  //roughly the smallest edge length in the
                            //source mesh
  std::cout << "input mesh: " << argv[1] << " outputMeshPrefix: "
            << prefix << " adaptRatio: " << adaptRatio
            << " max_size: " << max_size
            << " min_size: " << min_size << "\n";
  const auto outname = prefix + "_adaptRatio_" + std::string(argv[3]) +
                       "_maxSz_" + std::to_string(max_size) +
                       "_minSz_" + std::to_string(min_size);

  Omega_h::vtk::write_parallel("beforeClassFix_edges.vtk", &mesh, 1);

  auto effectiveStrain = getEffectiveStrainRate(mesh);
  auto recoveredStrain = recoverLinearStrain(mesh,effectiveStrain);
  mesh.add_tag<Real>(VERT, "recoveredStrain", 1, recoveredStrain);
  auto error = computeError(mesh,effectiveStrain,recoveredStrain);

  MeshField::OmegahMeshField<ExecutionSpace, MeshField::KokkosController> omf(
        mesh);

  const auto MeshDim = 2;
  const auto ShapeOrder = 1;
  auto recoveredStrainField = omf.CreateLagrangeField<Real, ShapeOrder, MeshDim>();
  setFieldAtVertices(mesh, recoveredStrain, recoveredStrainField);

  auto coordField = omf.getCoordField();
  const auto [shp, map] = MeshField::Omegah::getTriangleElement<ShapeOrder>(mesh);
  MeshField::FieldElement coordFe(mesh.nelems(), coordField, shp, map);

  auto estimation = Estimation(mesh, effectiveStrain,
      recoveredStrainField, adaptRatio);

  const auto tgtLength = getSprSizeField(estimation, omf, coordFe);
  Omega_h::Write<MeshField::Real> tgtLength_oh(tgtLength);
  mesh.add_tag<Real>(VERT, "tgtLength", 1, tgtLength_oh);

  { //write vtk
  const std::string vtkFileName = "beforeAdapt" + outname + ".vtk";
  Omega_h::vtk::write_parallel(vtkFileName, &mesh, 2);
  const std::string vtkFileName_edges = "beforeAdapt" + outname + "_edges.vtk";
  Omega_h::vtk::write_parallel(vtkFileName_edges, &mesh, 1);
  }

  //adapt
  auto opts = Omega_h::AdaptOpts(&mesh);
  setupFieldTransfer(opts);
  printTriCount(&mesh);

  auto verbose = true;
  const auto isos = Omega_h::isos_from_lengths(tgtLength_oh);
  auto metric = clamp_metrics(mesh.nverts(), isos, min_size, max_size);
  Omega_h::grade_fix_adapt(&mesh, opts, metric, verbose);

  { //write vtk and osh for adapted mesh
  const std::string outfilename = "afterAdapt" + outname;
  Omega_h::vtk::write_parallel(outfilename + ".vtk", &mesh, 2);
  Omega_h::binary::write(outfilename + ".osh", &mesh);
  std::cout << "wrote adapted mesh: " << outfilename + ".osh" << "\n";
  }


  return 0;
}
