#ifndef THWAITES_ERROR_ESTIMATOR_H
#define THWAITES_ERROR_ESTIMATOR_H
#include <Omega_h_mesh.hpp>
#include <KokkosController.hpp>
#include <MeshField.hpp>
#include <MeshField_Element.hpp>
#include <MeshField_Integrate.hpp>
#include <MeshField_Fail.hpp>
#include <MeshField_For.hpp>
#include <MeshField_ShapeField.hpp>


/* useful for initializing values to quickly
   detect "uninitialized" value bugs */
static double getNaN()
{
  return std::numeric_limits<double>::quiet_NaN();
}

/* common base for Scalar Integrator. */
class SInt : public MeshField::Integrator {
  public:
    SInt(int order):
      MeshField::Integrator(order),r(0)
    {}
    void reset() {r=0;}
    MeshField::Real r;
};

template <typename ShapeField>
class Estimation {
  private:
  Estimation() {}
  public:
  Estimation(Omega_h::Mesh& mesh_in, Omega_h::Reals eps_in, ShapeField eps_star_in,
      MeshField::Real tolerance_in) 
    : mesh(mesh_in), eps(eps_in), eps_star(eps_star_in), 
      tolerance(tolerance_in), size_factor(getNaN())
  {
    /* note that getOrder being used to convey this meaning
       is a bit of a hack, but looks decent if you don't
       try to define the FieldShape API too rigorously */
    integration_order = -1; //FIXME get from eps_star
    /* so far recovery order is directly tied to the
       mesh's coordinate field order, coordinate this
       with field recovery code */
    recovered_order = -1; //FIXME
  }
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
  ShapeField& eps_star;
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

/* computes $\|f\|^2$ */
template<typename EstimationT, typename OmegahMeshField> 
class SelfProduct : public SInt
{
  public:
    SelfProduct(EstimationT& estimation_in, OmegahMeshField& omf_in):
      SInt(estimation_in.integration_order),
      estimation(estimation_in),
      omf(omf_in)
    {
    }
    void atPoints(Kokkos::View<MeshField::Real**> p,
                  Kokkos::View<MeshField::Real*> w,
                  Kokkos::View<MeshField::Real*> dV)
    {

      std::cerr << "SelfProduct::atPoints(...)\n";
      const size_t numPtsPerElem = p.size()/estimation.mesh.nelems();
      auto eps_star_atPts = omf.triangleLocalPointEval(p,numPtsPerElem,
          estimation.eps_star);

      r = 0;
      Kokkos::parallel_reduce(
        "eval", estimation.mesh.nelems(),
        KOKKOS_LAMBDA(const int &elm, MeshField::Real &r_local) {
          const auto first = elm * numPtsPerElem;
          const auto last = first + numPtsPerElem;
          for (auto pt = first; pt < last; pt++) {
            const auto vPt = eps_star_atPts(pt,0);
            const auto wPt = w(pt);
            const auto dVPt = dV(pt);
            r_local += (vPt * vPt) * wPt * dVPt;
          }
        },
        r);

    }
  private:
    EstimationT& estimation;
    OmegahMeshField& omf;
};



#endif
