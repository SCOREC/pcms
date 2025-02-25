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

//detect floating point exceptions
#include <fenv.h>

using namespace Omega_h;

/* common base for Scalar Integrator. */
class SInt : public MeshField::Integrator {
  public:
    SInt(int order):
      MeshField::Integrator(order),r(0)
    {}
    void reset() {r=0;}
    MeshField::Real r;
};

struct Estimation {
  Omega_h::Mesh& mesh;
  /* the maximum polynomial order that can be
     integrated exactly using the input integration points.
     not necessarily equal to recovered_order, sometimes
     users give more points than strictly necessary */
  int integration_order;
  /* the polynomial order of the recovered field */
  int recovered_order;
  /* the input field consisting of values of a key
     quantity at integration points (stress or strain, for example).
     this field is related to integration_order */
  Omega_h::Reals eps;
  /* the recovered field, consisting of nodal values
     of the quantity of interest.
     this uses a Lagrange basis of recovered_order */
  Omega_h::Reals eps_star;
  /* the acceptable margin of error, expressed as a factor
     greater than zero.
     setting this equal to zero would request zero element
     size everywhere, i.e. infinite refinement.
     increasing it scales up the desired element size
     throughout.
     basically, this is a (nonlinear) scaling factor on the resulting
     size field */
  Omega_h::Real tolerance;
  /* the uniform linear scaling factor derived from the tolerance
     and integrals of the recovered field over the mesh.
     desired element size = current size * current error * size_factor */
  Omega_h::Real size_factor;
  /* a temporary field storing desired sizes at elements */
  Omega_h::Reals element_size;
  /* the resulting size field, recovered from the element_size field
     (using a local average recovery method much weaker than SPR) */
  Omega_h::Reals size;
};

/* useful for initializing values to quickly
   detect "uninitialized" value bugs */
static double getNaN()
{
  return std::numeric_limits<double>::quiet_NaN();
}

static void setupEstimation(Mesh& mesh, Estimation* e, const Reals eps, Real tolerance)
{
  /* note that getOrder being used to convey this meaning
     is a bit of a hack, but looks decent if you don't
     try to define the FieldShape API too rigorously */
  e->integration_order = apf::getShape(eps)->getOrder(); //FIXME
  e->mesh = mesh;
  /* so far recovery order is directly tied to the
     mesh's coordinate field order, coordinate this
     with field recovery code */
  e->recovered_order = e->mesh->getShape()->getOrder(); //FIXME
  e->eps = eps;
  e->tolerance = tolerance;
  e->size_factor = getNaN();
}

/* computes $\|f\|^2$ */
class SelfProduct : public SInt
{
  public:
    SelfProduct(Estimation* e):
      SInt(e->integration_order), estimation(e), element(0)
    {
      v.setSize(apf::countComponents(e->eps_star)); //FIXME Vector2?
    }
    void atPoint(apf::Vector3 const& p, double w, double dV) //FIXME bulk operation
    {
      apf::getComponents(element, p, &v[0]); //FIXME get components of eps_star
      r += (v * v) * w * dV; //FIXME vector*vector needed
    }
  private:
    Estimation* estimation;
    apf::Element* element;
    apf::DynamicVector v;
};





template <typename T>
void printTagInfo(Omega_h::Mesh& mesh, std::ostringstream& oss, int dim, int tag, std::string type) {
    auto tagbase = mesh.get_tag(dim, tag);
    auto array = Omega_h::as<T>(tagbase)->array();

    Omega_h::Real min = get_min(array);
    Omega_h::Real max = get_max(array);

    oss << std::setw(18) << std::left << tagbase->name().c_str()
        << std::setw(5) << std::left << dim
        << std::setw(7) << std::left << type
        << std::setw(5) << std::left << tagbase->ncomps()
        << std::setw(10) << std::left << min
        << std::setw(10) << std::left << max
        << "\n";
}

void printTags(Mesh& mesh) {
    std::ostringstream oss;
    // always print two places to the right of the decimal
    // for floating point types (i.e., imbalance)
    oss.precision(2);
    oss << std::fixed;

    if (!mesh.comm()->rank()) {
        oss << "\nTag Properties by Dimension: (Name, Dim, Type, Number of Components, Min. Value, Max. Value)\n";
        for (int dim=0; dim <= mesh.dim(); dim++)
        for (int tag=0; tag < mesh.ntags(dim); tag++) {
            auto tagbase = mesh.get_tag(dim, tag);
            if (tagbase->type() == OMEGA_H_I8)
                printTagInfo<Omega_h::I8>(mesh, oss, dim, tag, "I8");
            if (tagbase->type() == OMEGA_H_I32)
                printTagInfo<Omega_h::I32>(mesh, oss, dim, tag, "I32");
            if (tagbase->type() == OMEGA_H_I64)
                printTagInfo<Omega_h::I64>(mesh, oss, dim, tag, "I64");
            if (tagbase->type() == OMEGA_H_F64)
                printTagInfo<Omega_h::Real>(mesh, oss, dim, tag, "F64");
        }

        std::cout << oss.str();
    }

}

void printTriCount(Mesh* mesh) {
  const auto nTri = mesh->nglobal_ents(2);
  if (!mesh->comm()->rank())
    std::cout << "nTri: " << nTri << "\n";
}

void setupFieldTransfer(AdaptOpts& opts) {
  opts.xfer_opts.type_map["velocity"] = OMEGA_H_LINEAR_INTERP;
  opts.xfer_opts.type_map["ice_thickness"] = OMEGA_H_LINEAR_INTERP;
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
  return error;
}

int main(int argc, char** argv) {
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);  // Enable all floating point exceptions but FE_INEXACT
  auto lib = Library(&argc, &argv);
  if( argc != 3 ) {
    fprintf(stderr, "Usage: %s inputMesh.osh outputMeshPrefix\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  auto world = lib.world();
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(argv[1], world, &mesh);

  Omega_h::vtk::write_parallel("beforeClassFix_edges.vtk", &mesh, 1);

  auto effectiveStrain = getEffectiveStrainRate(mesh);
  auto recoveredStrain = recoverLinearStrain(mesh,effectiveStrain);
  mesh.add_tag<Real>(VERT, "recoveredStrain", 1, recoveredStrain);

  const std::string vtkFileName = std::string(argv[2]) + ".vtk";
  Omega_h::vtk::write_parallel(vtkFileName, &mesh, 2);
  const std::string vtkFileName_edges = std::string(argv[2]) + "_edges.vtk";
  Omega_h::vtk::write_parallel(vtkFileName_edges, &mesh, 1);
  return 0;
}
