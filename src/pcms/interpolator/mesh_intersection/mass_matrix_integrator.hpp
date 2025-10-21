#ifndef PCMS_INTERPOLATOR_MASS_MATRIX_INTEGRATOR_H
#define PCMS_INTERPOLATOR_MASS_MATRIX_INTEGRATOR_H

#include <KokkosController.hpp>
#include <MeshField.hpp>
#include <MeshField_Element.hpp>
#include <MeshField_Fail.hpp>
#include <MeshField_For.hpp>
#include <MeshField_Integrate.hpp>
#include <MeshField_ShapeField.hpp>
#include <Omega_h_mesh.hpp>

// computes the mass matrix for each element
template <typename FieldElement>
class MassMatrixIntegrator : public MeshField::Integrator
{
public:
  MassMatrixIntegrator(Omega_h::Mesh mesh_in, FieldElement& fe_in,
                       int order = 2)
    : mesh(mesh_in),
      fe(fe_in),
      subMatrixSize(3 * 3), // FIXME remove hard coded size
      elmMassMatrix("elmMassMatrix", mesh_in.nelems() * 3 * 3),
      Integrator(order)
  {
    Kokkos::deep_copy(elmMassMatrix, 0);
    assert(mesh.dim() == 2); // TODO support 1d,2d,3d
    assert(mesh.family() == OMEGA_H_SIMPLEX);
  }
  void atPoints(Kokkos::View<MeshField::Real**> p,
                Kokkos::View<MeshField::Real*> w,
                Kokkos::View<MeshField::Real*> dV)
  {
    // std::cerr << "MassMatrixIntegrator::atPoints(...)\n";
    const size_t numPtsPerElem = p.extent(0) / mesh.nelems();
    // std::cerr << " Number points per Elem : " << numPtsPerElem << "\n";
    assert(numPtsPerElem >= 1);
    const size_t ptDim = p.extent(1);
    assert(ptDim == fe.MeshEntDim + 1);
    Kokkos::parallel_for(
      "eval", mesh.nelems(), KOKKOS_CLASS_LAMBDA(const int& elm) {
        const auto first = elm * numPtsPerElem;
        const auto last = first + numPtsPerElem;
        for (auto pt = first; pt < last; pt++) {
          // FIXME better way to fill? pass kokkos::subview to getValues?
          Kokkos::Array<MeshField::Real, FieldElement::MeshEntDim + 1>
            localCoord;
          for (auto i = 0; i < localCoord.size(); i++) {
            localCoord[i] = p(pt, i);
          }
          const auto N = fe.shapeFn.getValues(localCoord);
          const auto wPt = w(pt);
          const auto dVPt = dV(pt);
          //       printf("Shape Functions: %f, %f, %f \n", N[0], N[1], N[2]);
          //       printf("wPt, dVPt: %f, %f \n", wPt, dVPt);
          for (auto i = 0; i < N.size(); i++) {
            for (auto j = 0; j < N.size(); j++) {
              elmMassMatrix(elm * subMatrixSize + i * 3 + j) +=
                N[i] * N[j] * wPt * dVPt;
            }
          }
        }
      });
  }
  Omega_h::Mesh& mesh;
  FieldElement& fe;
  const int subMatrixSize;
  Kokkos::View<MeshField::Real*>
    elmMassMatrix; // numNodes^2 entries per element
};

template <typename FieldElement>
Kokkos::View<MeshField::Real*> buildMassMatrix(Omega_h::Mesh& mesh,
                                               FieldElement& coordFe)
{
  MassMatrixIntegrator mmi(mesh, coordFe);
  mmi.process(coordFe);
  return mmi.elmMassMatrix;
}

#endif
