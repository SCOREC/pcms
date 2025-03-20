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

#include <MeshField.hpp>
#include "massMatrixIntegrator.hpp"

#include <petscmat.h>
#include <petscvec_kokkos.hpp>

//detect floating point exceptions
#include <fenv.h>

using ExecutionSpace = Kokkos::DefaultExecutionSpace;
using MemorySpace = Kokkos::DefaultExecutionSpace::memory_space;

using namespace Omega_h;

template <typename ShapeField>
void setFieldAtVertices(Omega_h::Mesh &mesh, ShapeField field, const MeshField::Real val) {
  const auto MeshDim = mesh.dim();
  auto setFieldAtVertices = KOKKOS_LAMBDA(const int &vtx) {
    field(0, 0, vtx, MeshField::Vertex) = val;
  };
  MeshField::parallel_for(ExecutionSpace(), {0}, {mesh.nverts()},
                          setFieldAtVertices, "setFieldAtVertices");
}

//FIXME remove the hard coded 3x3
static PetscErrorCode CreateMatrix(Omega_h::Mesh& mesh, Mat *A) {
  PetscInt *oor, *ooc, cnt = 0;
  PetscFunctionBeginUser;
  PetscCall(MatCreate(PETSC_COMM_WORLD, A));
  PetscCall(MatSetSizes(*A, mesh.nverts(), mesh.nverts(), PETSC_DECIDE, PETSC_DECIDE));
  PetscCall(MatSetFromOptions(*A));
  /* determine for each entry in each element stiffness matrix the global row and column */
  /* since the element is triangular with piecewise linear basis functions there are three degrees of freedom per element, one for each vertex */
  PetscCall(PetscMalloc2(3 * 3 * mesh.nelems(), &oor, 3 * 3 * mesh.nelems(), &ooc));
  for (PetscInt e = 0; e < mesh.nelems(); e++) {
    for (PetscInt vi = 0; vi < 3; vi++) {
      for (PetscInt vj = 0; vj < 3; vj++) {
        oor[cnt]   = fe->vertices[3 * e + vi]; //FIXME - elmToVtx
        ooc[cnt++] = fe->vertices[3 * e + vj]; //FIXME - elmToVtx
      }
    }
  }
  PetscCall(MatSetPreallocationCOO(*A, 3 * 3 * mesh.nelems(), oor, ooc));
  PetscCall(PetscFree2(oor, ooc));
  PetscFunctionReturn(PETSC_SUCCESS);
}


int main(int argc, char** argv) {
  PetscCall(PetscInitialize(&argc,&argv,NULL,NULL));
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);  // Enable all floating point exceptions but FE_INEXACT
  auto lib = Library(&argc, &argv);
  if( argc != 3 ) {
    fprintf(stderr, "Usage: %s inputMesh.osh outputMeshPrefix\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  auto world = lib.world();
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[1], world, &mesh);
  const auto prefix = std::string(argv[2]);
  std::cout << "input mesh: " << argv[1] << " outputMeshPrefix: " << prefix << "\n";

  MeshField::OmegahMeshField<ExecutionSpace, MeshField::KokkosController> omf(
        mesh);

  const auto ShapeOrder = 1;
  auto coordField = omf.getCoordField();
  const auto [shp, map] = MeshField::Omegah::getTriangleElement<ShapeOrder>(mesh);
  MeshField::FieldElement coordFe(mesh.nelems(), coordField, shp, map);

  auto elmMassMatrix = buildMassMatrix(mesh, coordFe);

  mesh.add_tag(2, "elmMassMatrix", 3*3, Omega_h::read(Omega_h::Write<MeshField::Real>(elmMassMatrix)));

  Omega_h::vtk::write_parallel("massMatrix.vtk", &mesh, 2);

  Mat mass;
  PetscCall(CreateMatrix(mesh, mass));
  PetscCall(MatZeroEntries(mass));
  PetscCall(MatSetValuesCOO(mass, elmMassMatrix.data(), INSERT_VALUES));
  if( mesh.nelems() < 10 ) {
    PetscCall(MatView(mass, PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(MatDestroy(&mass));
  PetscCall(PetscFinalize());

  return 0;
}
