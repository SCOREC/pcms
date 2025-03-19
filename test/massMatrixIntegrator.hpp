#ifndef MASS_MATRIX_INTEGRATOR_H
#define MASS_MATRIX_INTEGRATOR_H
#include <Omega_h_mesh.hpp>
#include <KokkosController.hpp>
#include <MeshField.hpp>
#include <MeshField_Element.hpp>
#include <MeshField_Integrate.hpp>
#include <MeshField_Fail.hpp>
#include <MeshField_For.hpp>
#include <MeshField_ShapeField.hpp>

/* computes the mass matrix for each element */
template<typename OmegahMeshField> 
class MassMatrixIntegrator : public MeshField::Integrator {
{
  public:
    MassMatrixIntegrator(Omega_h::Mesh mesh_in, ShapeField& constantField_in, OmegahMeshField& omf_in) :
      mesh(mesh_in),
      constantField(constantField_in),
      omf(omf_in),
      elmMassMatrix("elmMassMatrix", mesh_in.nelems()*3*3) //FIXME - remove hard code size based on numNodes per elm 
    {
      assert(mesh.dim() == 2); //TODO support 1d,2d,3d
      assert(mesh.family() == OMEGA_H_SIMPLEX);
    }
    void atPoints(Kokkos::View<MeshField::Real**> p,
                  Kokkos::View<MeshField::Real*> w,
                  Kokkos::View<MeshField::Real*> dV)
    {
      std::cerr << "MassMatrixIntegrator::atPoints(...)\n";
      const size_t numPtsPerElem = p.extent(0)/estimation.mesh.nelems();
      //FIXME constantField isn't a ShapeField so we can't call 
      //      omf.triangleLocalPointEval.  For now, just get
      //      the value from the element and assert that there is
      //      one integration point per element.
      assert(numPtsPerElem == 1);
      auto shpFn_at_pts = omf.triangleLocalPointEval(p,numPtsPerElem,constantField);
      Kokkos::parallel_for("eval", estimation.mesh.nelems(),
        KOKKOS_CLASS_LAMBDA(const int &elm) {
          const auto first = elm * numPtsPerElem;
          const auto last = first + numPtsPerElem;
          for (auto pt = first; pt < last; pt++) {
            const auto shpFn_pt = shpFn_at_pts(pt,0); //FIXME same as apf::getShapeValues that returns one value per node
            const auto wPt = w(pt);
            const auto dVPt = dV(pt);
            sum += shpFn_pt * wPt * dVPt;
          }
          elmMassMatrix(elm) = mm;//FIXME
        });
    }
    Omega_h::Mesh& mesh;
    ShapeField& constantField;
    OmegahMeshField& omf;
    Kokkos::View<MeshField::Real*> elmMassMatrix; //numNodes^2 entries per element
};

template<typename EstimationT, typename OmegahMeshField, typename FieldElement>
Kokkos::View<MeshField::Real*>
buildMassMatrix(EstimationT& e, OmegahMeshField& omf, FieldElement& coordFe) {
  MassMatrixIntegrator mmi(e,omf);
  mmi.process(coordFe);
  return Kokkos::View<MeshField::Real*>("foo",1); //FIXME
}

#endif
