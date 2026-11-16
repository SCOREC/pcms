#ifndef PCMS_TRANSFER_MASS_MATRIX_INTEGRATOR_HPP
#define PCMS_TRANSFER_MASS_MATRIX_INTEGRATOR_HPP

#include <KokkosController.hpp>
#include <MeshField.hpp>
#include <MeshField_Config.hpp>
#include <MeshField_Element.hpp>
#include <MeshField_Fail.hpp>
#include <MeshField_For.hpp>
#include <MeshField_Integrate.hpp>
#include <MeshField_ShapeField.hpp>
#include <Omega_h_mesh.hpp>

namespace pcms
{
// computes the mass matrix for each element
template <typename FieldElement>
class MassMatrixIntegrator : public MeshField::Integrator
{
public:
  // Linear simplex: numNodes = spatial dim + 1 (3 for triangles, 4 for tets).
  static constexpr int numNodes = FieldElement::MeshEntDim + 1;

  MassMatrixIntegrator(Omega_h::Mesh& mesh_in, FieldElement& fe_in,
                       int order = 2)
    : mesh(mesh_in),
      fe(fe_in),
      subMatrixSize(numNodes * numNodes),
      elmMassMatrix("elmMassMatrix", mesh_in.nelems() * numNodes * numNodes),
      Integrator(order)
  {
    Kokkos::deep_copy(elmMassMatrix, 0);
    assert(mesh.dim() == 2 || mesh.dim() == 3);
    assert(mesh.family() == OMEGA_H_SIMPLEX);
  }
  void atPoints(Kokkos::View<MeshField::Real**> p,
                Kokkos::View<MeshField::Real*> w,
                Kokkos::View<MeshField::Real*> dV) override
  {
    // std::cerr << "MassMatrixIntegrator::atPoints(...)\n";
    const size_t numPtsPerElem = p.extent(0) / mesh.nelems();
    // std::cerr << " Number points per Elem : " << numPtsPerElem << "\n";
    assert(numPtsPerElem >= 1);
    const size_t ptDim = p.extent(1);
    assert(ptDim == fe.MeshEntDim);
    // Copy values needed in the kernel to avoid capturing host references
    // (mesh and fe are host objects and cannot be dereferenced on the device)
    const auto numElems = mesh.nelems();
    const auto shapeFn = fe.shapeFn;
    const auto subMat = subMatrixSize;
    auto massMatrix = elmMassMatrix;
    Kokkos::parallel_for(
      "eval", numElems, KOKKOS_LAMBDA(const int& elm) {
        const auto first = elm * numPtsPerElem;
        const auto last = first + numPtsPerElem;
        for (auto pt = first; pt < last; pt++) {
          // FIXME better way to fill? pass kokkos::subview to getValues?
          Kokkos::Array<MeshField::Real, FieldElement::MeshEntDim> localCoord;
          for (auto i = 0; i < localCoord.size(); i++) {
            localCoord[i] = p(pt, i);
          }
          const auto N = shapeFn.getValues(localCoord);
          const auto wPt = w(pt);
          // Use the unsigned volume element: MeshField returns a signed
          // Jacobian determinant, which is negative for tetrahedra whose vertex
          // ordering has negative orientation. A mass matrix integrates against
          // the positive volume measure, so take the magnitude (a no-op in 2D
          // where the differential area is already positive).
          const auto dVPt = Kokkos::fabs(dV(pt));
          //       printf("Shape Functions: %f, %f, %f \n", N[0], N[1], N[2]);
          //       printf("wPt, dVPt: %f, %f \n", wPt, dVPt);
          for (auto i = 0; i < N.size(); i++) {
            for (auto j = 0; j < N.size(); j++) {
              massMatrix(elm * subMat + i * numNodes + j) +=
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
Kokkos::View<MeshField::Real*> buildElementMassMatrix(Omega_h::Mesh& mesh,
                                                      FieldElement& coordFe)
{
  MassMatrixIntegrator mmi(mesh, coordFe);
  mmi.process(coordFe);
  return mmi.elmMassMatrix;
}
} // namespace pcms
#endif // PCMS_TRANSFER_MASS_MATRIX_INTEGRATOR_HPP
