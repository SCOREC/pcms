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
    integration_order = 1; //FIXME get from eps_star
    /* so far recovery order is directly tied to the
       mesh's coordinate field order, coordinate this
       with field recovery code */
    recovered_order = 1; //FIXME
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
  Kokkos::View<MeshField::Real*> element_size;
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
      const size_t numPtsPerElem = p.extent(0)/estimation.mesh.nelems();
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

/* computes the integral over the element of the
   sum of the squared differences between the
   original and recovered fields */
template<typename EstimationT, typename OmegahMeshField> 
class Error : public SInt
{
  public:
    Error(EstimationT& estimation_in, OmegahMeshField& omf_in):
      SInt(estimation_in.integration_order),
      estimation(estimation_in),
      omf(omf_in),
      errorNorm("errorNorm", estimation_in.mesh.nelems())
    {
    }
    void atPoints(Kokkos::View<MeshField::Real**> p,
                  Kokkos::View<MeshField::Real*> w,
                  Kokkos::View<MeshField::Real*> dV)
    {
      std::cerr << "SelfProduct::atPoints(...)\n";
      const size_t numPtsPerElem = p.extent(0)/estimation.mesh.nelems();
      //FIXME eps isn't a ShapeField so we can't call 
      //      omf.triangleLocalPointEval.  For now, just get
      //      the value from the element and assert that there is
      //      one integration point per element.
      assert(numPtsPerElem == 1);
      auto eps_star_atPts = omf.triangleLocalPointEval(p,numPtsPerElem,
          estimation.eps_star);
      double meshDim = estimation.mesh.dim();
      double orderP = estimation.recovered_order;

      const auto& eps = estimation.eps;
      r = 0;
      Kokkos::parallel_reduce(
        "eval", estimation.mesh.nelems(),
        KOKKOS_LAMBDA(const int &elm, MeshField::Real &r_local) {
          const auto first = elm * numPtsPerElem;
          const auto last = first + numPtsPerElem;
          const auto eps_elm = eps[elm];
          MeshField::Real sum = 0;
          for (auto pt = first; pt < last; pt++) {
            const auto eps_star_Pt = eps_star_atPts(pt,0);
            const auto diff = eps_elm - eps_star_Pt;
            const auto wPt = w(pt);
            const auto dVPt = dV(pt);
            sum += (diff * diff) * wPt * dVPt;
          }
          r_local += Kokkos::pow(Kokkos::sqrt(sum), ((2 * meshDim) / (2 * orderP + meshDim)));
          errorNorm(elm) = Kokkos::pow(Kokkos::sqrt(sum), -(2 / (2 * orderP + meshDim)));  
        },
        r); // $\sum_{i=1}^n \|e_\epsilon\|^{\frac{2d}{2p+d}}$ 
    }
    EstimationT& estimation;
    OmegahMeshField& omf;
    Kokkos::View<MeshField::Real*> errorNorm; // (||e_eps||_e)^(-2/(2p+d))
};

//TODO move this into Estimation class
template<typename EstimationT, typename OmegahMeshField, typename FieldElement, typename ErrorT>
void computeSizeFactor(EstimationT& e, OmegahMeshField& omf, FieldElement& coordFe, ErrorT& errorIntegrator) {
  SelfProduct sp(e,omf);
  sp.process(coordFe);
  const double epsStarNorm = Kokkos::sqrt(sp.r);
  std::cout << "SelfProduct: " << epsStarNorm << "\n";
  const double a = e.tolerance * e.tolerance * // (n^hat)^2
             epsStarNorm * epsStarNorm; // ||e*||^2
  const double b = a / errorIntegrator.r; // term in parenthesis in section 4 of spr.tex
  const double p = e.recovered_order;
  e.size_factor = Kokkos::pow(b, 1.0 / (2.0 * p));
  std::cout << "size_factor: " << e.size_factor << "\n";

}

/* computes h_e^new from section 4 of spr.tex */
template<typename EstimationT, typename ErrorT>
void getElementSizeField(EstimationT& e, ErrorT& errorIntegrator) {
  Kokkos::View<MeshField::Real*> eSize("eSize", e.mesh.nelems());
  const auto errorNorm = errorIntegrator.errorNorm;
  const auto size_factor = e.size_factor;
  const auto currentElmSize = e.mesh.ask_sizes();
  Kokkos::parallel_for(e.mesh.nelems(), KOKKOS_LAMBDA(const int elm) {
    const double h = currentElmSize[elm]; // h_e^current //FIXME
    eSize(elm) = h * errorNorm(elm) * size_factor;
  });
  e.element_size = eSize;
}

Kokkos::View<MeshField::Real*>
averageToVertex(Omega_h::Mesh& mesh, const Kokkos::View<MeshField::Real*>& elmSize) {
  Kokkos::View<MeshField::Real*> sizeField("sizeField", mesh.nverts());
  const auto v2e = mesh.ask_up(Omega_h::VERT, mesh.dim());
  const auto v2e_offsets = v2e.a2ab;
  const auto v2e_values = v2e.ab2b;
  Kokkos::parallel_for(mesh.nverts(), KOKKOS_LAMBDA(const int vtx) {
    MeshField::Real s = 0;
    for (auto idx = v2e_offsets[vtx]; idx < v2e_offsets[vtx + 1]; ++idx) {
      const auto elm = v2e_values[idx];
      s += elmSize(elm);
    }
    const auto numUpElms = v2e_offsets[vtx+1] - v2e_offsets[vtx];
    sizeField(vtx) =  s / numUpElms;
  });
  return sizeField;
}

template<typename EstimationT, typename OmegahMeshField, typename FieldElement>
Kokkos::View<MeshField::Real*>
getSprSizeField(EstimationT& e, OmegahMeshField& omf, FieldElement& coordFe) {
  Error errorIntegrator(e,omf);
  errorIntegrator.process(coordFe);
  std::cout << "Error: " << errorIntegrator.r << "\n";
  computeSizeFactor(e, omf, coordFe, errorIntegrator);
  getElementSizeField(e, errorIntegrator);
  return averageToVertex(e.mesh, e.element_size);
}

#endif
