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

#include <petscvec_kokkos.hpp>
#include <petscmat.h>

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

static PetscErrorCode CreateMatrix(Omega_h::Mesh& mesh, Mat *A) {
  const auto numNodesPerTri = 3; //FIXME query the mesh
  const auto matSize = numNodesPerTri*numNodesPerTri*mesh.nelems();
  auto elmVerts = Omega_h::HostRead(mesh.ask_elem_verts());
  PetscInt *oor, *ooc, cnt = 0;
  PetscFunctionBeginUser;
  PetscCall(MatCreate(PETSC_COMM_WORLD, A));
  PetscCall(MatSetSizes(*A, mesh.nverts(), mesh.nverts(), PETSC_DECIDE, PETSC_DECIDE));
  PetscCall(MatSetFromOptions(*A));
  /* determine for each entry in each element stiffness matrix the global row and column */
  /* since the element is triangular with piecewise linear basis functions there are three degrees of freedom per element, one for each vertex */
  PetscCall(PetscMalloc2(matSize, &oor, matSize, &ooc));
  for (PetscInt e = 0; e < mesh.nelems(); e++) {
    for (PetscInt vi = 0; vi < numNodesPerTri; vi++) {
      for (PetscInt vj = 0; vj < numNodesPerTri; vj++) {
        oor[cnt]   = elmVerts[numNodesPerTri * e + vi];
        ooc[cnt++] = elmVerts[numNodesPerTri * e + vj];
      }
    }
  }
  PetscCall(MatSetPreallocationCOO(*A, matSize, oor, ooc));
  PetscCall(PetscFree2(oor, ooc));
  PetscFunctionReturn(PETSC_SUCCESS);
}


int main(int argc, char** argv) {
  auto lib = Library(&argc, &argv); //initializes MPI
  PetscCall(PetscInitialize(&argc,&argv,NULL,NULL));
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);  // Enable all floating point exceptions but FE_INEXACT
  if( argc < 3 ) {
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
  PetscCall(CreateMatrix(mesh, &mass));
  PetscBool is_kokkos;
  PetscCall(PetscObjectBaseTypeCompare((PetscObject)mass, MATSEQAIJKOKKOS, &is_kokkos));
  std::cerr << "Matrix type is kokkos: " << is_kokkos << "\n";
  PetscCall(MatZeroEntries(mass));
  PetscCall(MatSetValuesCOO(mass, elmMassMatrix.data(), INSERT_VALUES)); //FIXME fails here on gpu, calls into host implementation... AFAIK, petsc checks the type of the input array of values to decide which backend to use...
  if( mesh.nelems() < 10 ) {
    PetscCall(MatView(mass, PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(MatDestroy(&mass));
  PetscCall(PetscFinalize());

  return 0;
}
