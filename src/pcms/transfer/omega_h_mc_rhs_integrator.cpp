#include "pcms/transfer/omega_h_mc_rhs_integrator.hpp"
#include "pcms/transfer/omega_h_form_integrator_utils.hpp"
#include "pcms/transfer/petsc_utils.hpp"
#include "pcms/utility/assert.h"
#include <Kokkos_Random.hpp>
#include <Omega_h_shape.hpp>
#include <petscksp.h>
#include "pcms/field/coordinate_systems/cartesian.hpp"

namespace pcms
{

namespace
{

// Maps a uniform point on the unit square to uniform barycentric coordinates on
// a triangle. Shape Distributions (ACM Transactions on Graphics, Vol. 21,
// No. 4, October 2002.) page 814 Eq 1.
KOKKOS_INLINE_FUNCTION Omega_h::Vector<3> UniformTriangleBarycentric(Real u,
                                                                     Real v)
{
  const Real s = Kokkos::sqrt(u);
  return {1.0 - s, s * (1.0 - v), s * v};
}

// Maps three uniform draws on the unit cube to uniform barycentric coordinates
// on a tetrahedron via the cut-and-fold method (Rocchini & Cignoni,
// "Generating Random Points in a Tetrahedron", J. Graphics Tools 2000).
KOKKOS_INLINE_FUNCTION Omega_h::Vector<4> UniformTetBarycentric(Real s, Real t,
                                                                Real u)
{
  if (s + t > 1.0) { // fold the cube into a prism
    s = 1.0 - s;
    t = 1.0 - t;
  }
  if (t + u > 1.0) { // fold the prism into a tetrahedron
    const Real tmp = u;
    u = 1.0 - s - t;
    t = 1.0 - tmp;
  } else if (s + t + u > 1.0) {
    const Real tmp = u;
    u = s + t + u - 1.0;
    s = 1.0 - t - tmp;
  }
  const Real a = 1.0 - s - t - u;
  return {a, s, t, u};
}

template <int Dim>
KOKKOS_INLINE_FUNCTION Omega_h::Vector<Dim + 1> UniformSimplexBarycentric(
  const Real r[Dim])
{
  if constexpr (Dim == 3) {
    return UniformTetBarycentric(r[0], r[1], r[2]);
  } else {
    return UniformTriangleBarycentric(r[0], r[1]);
  }
}

// Fills coords, node_gids, and coeffs for all samples of one target element.
// unit_sample_at(s, r) fills the s-th sample's Dim unit-hypercube draws.
template <int Dim, typename UnitSampleAt>
KOKKOS_INLINE_FUNCTION void FillElementSamples(
  const int elm, const int samples_per_element,
  const Omega_h::Reals& mesh_coords, const Omega_h::LOs& elems2nodes,
  const Kokkos::View<const LO*, DeviceMemorySpace>& global_to_local,
  const Kokkos::View<Real**, DeviceMemorySpace>& coords,
  const Kokkos::View<PetscInt*, DeviceMemorySpace>& node_gids,
  const Kokkos::View<Real*, DeviceMemorySpace>& coeffs,
  const UnitSampleAt& unit_sample_at)
{
  constexpr int nv = Dim + 1;
  const auto verts = Omega_h::gather_verts<nv>(elems2nodes, elm);
  const auto vert_coords = Omega_h::gather_vectors<nv, Dim>(mesh_coords, verts);
  Omega_h::Few<Omega_h::Vector<Dim>, Dim> basis;
  for (int d = 0; d < Dim; ++d) {
    basis[d] = vert_coords[d + 1] - vert_coords[0];
  }
  Real measure;
  if constexpr (Dim == 3) {
    measure = Kokkos::fabs(Omega_h::tet_volume_from_basis(basis));
  } else {
    measure = Kokkos::fabs(Omega_h::triangle_area_from_basis(basis));
  }
  const Real weight = measure / samples_per_element;

  for (int s = 0; s < samples_per_element; ++s) {
    Real r[Dim];
    unit_sample_at(s, r);
    const auto bary = UniformSimplexBarycentric<Dim>(r);

    const int i = elm * samples_per_element + s;
    Omega_h::Vector<Dim> x;
    for (int d = 0; d < Dim; ++d) {
      x[d] = 0.0;
    }
    for (int k = 0; k < nv; ++k) {
      for (int d = 0; d < Dim; ++d) {
        x[d] += bary[k] * vert_coords[k][d];
      }
      node_gids(i * nv + k) = static_cast<PetscInt>(global_to_local(verts[k]));
      coeffs(i * nv + k) = bary[k] * weight;
    }
    for (int d = 0; d < Dim; ++d) {
      coords(i, d) = x[d];
    }
  }
}

// Samples samples_per_element from a uniform random distribution over each
// target element, writing their coordinates, node GIDs, and coefficients.
template <int Dim>
void FillElementSamples(
  int nelems, int samples_per_element, const Omega_h::Reals& mesh_coords,
  const Omega_h::LOs& elems2nodes,
  const Kokkos::View<const LO*, DeviceMemorySpace>& global_to_local,
  const Kokkos::View<Real**, DeviceMemorySpace>& coords,
  const Kokkos::View<PetscInt*, DeviceMemorySpace>& node_gids,
  const Kokkos::View<Real*, DeviceMemorySpace>& coeffs, uint64_t seed)
{
  Kokkos::Random_XorShift64_Pool<DefaultExecutionSpace> pool(seed);
  Kokkos::parallel_for(
    "mc_rhs_fill_random", Kokkos::RangePolicy<DefaultExecutionSpace>(0, nelems),
    KOKKOS_LAMBDA(int elm) {
      auto gen = pool.get_state();
      FillElementSamples<Dim>(elm, samples_per_element, mesh_coords,
                              elems2nodes, global_to_local, coords, node_gids,
                              coeffs, [&](int /*s*/, Real r[Dim]) {
                                for (int d = 0; d < Dim; ++d) {
                                  r[d] = gen.drand();
                                }
                              });
      pool.free_state(gen);
    });
  Kokkos::fence();
}

} // namespace

OmegaHMonteCarloRHSIntegrator::OmegaHMonteCarloRHSIntegrator(
  const FunctionSpace& target_space, int samples_per_element,
  MonteCarloSampling sampling, uint64_t seed)
  : OmegaHMonteCarloRHSIntegrator(
      std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(
        target_space.GetLayout()),
      target_space.GetCoordinateSystem(), samples_per_element, sampling, seed)
{
}

OmegaHMonteCarloRHSIntegrator::OmegaHMonteCarloRHSIntegrator(
  std::shared_ptr<const OmegaHLagrangeLayout> target_layout,
  std::shared_ptr<const CoordinateSystem> target_coordinate_system,
  int samples_per_element, MonteCarloSampling /*sampling*/, uint64_t seed)
{
  detail::CheckOmegaHScalarP1Layout(target_coordinate_system, target_layout,
                                    "OmegaHMonteCarloRHSIntegrator", "target");
  if (samples_per_element <= 0) {
    throw pcms_error(
      "OmegaHMonteCarloRHSIntegrator: samples_per_element must be positive");
  }

  Omega_h::Mesh& mesh = target_layout->GetMesh();
  const int dim = mesh.dim();
  nbary_ = dim + 1;

  const int nelems = mesh.nelems();
  const int num_samples = nelems * samples_per_element;
  const auto mesh_coords = mesh.coords();
  const auto elems2nodes = mesh.ask_down(dim, Omega_h::VERT).ab2b;
  const auto global_to_local = target_layout->GetGlobalToLocalPermutation();

  Kokkos::View<Real**, DeviceMemorySpace> coords("mc_rhs_coords", num_samples,
                                                 dim);
  Kokkos::View<PetscInt*, DeviceMemorySpace> node_gids(
    "mc_rhs_node_gids", static_cast<std::size_t>(num_samples) * nbary_);
  Kokkos::View<Real*, DeviceMemorySpace> coeffs(
    "mc_rhs_coeffs", static_cast<std::size_t>(num_samples) * nbary_);

  if (dim == 3) {
    FillElementSamples<3>(nelems, samples_per_element, mesh_coords, elems2nodes,
                          global_to_local, coords, node_gids, coeffs, seed);
  } else {
    FillElementSamples<2>(nelems, samples_per_element, mesh_coords, elems2nodes,
                          global_to_local, coords, node_gids, coeffs, seed);
  }

  coords_ = std::move(coords);
  node_gids_ = std::move(node_gids);
  coeffs_ = std::move(coeffs);

  const PetscInt nnz = static_cast<PetscInt>(node_gids_.extent(0));
  PetscErrorCode ierr =
    createSeqVec(PETSC_COMM_SELF, static_cast<PetscInt>(mesh.nverts()), &vec_);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
  // VecSetPreallocationCOO takes the COO indices on the host
  // TODO: ask Todd/PETSc folks if there is a better way to do this for GPU
  // support
  auto node_gids_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, node_gids_);
  ierr = VecSetPreallocationCOO(vec_, nnz, node_gids_host.data());
  CHKERRABORT(PETSC_COMM_SELF, ierr);
}

OmegaHMonteCarloRHSIntegrator::~OmegaHMonteCarloRHSIntegrator()
{
  if (vec_) {
    VecDestroy(&vec_);
  }
}

// ---------------------------------------------------------------------------
// OmegaHMonteCarloRHSIntegrator public interface
// ---------------------------------------------------------------------------

CoordinateView<DeviceMemorySpace>
OmegaHMonteCarloRHSIntegrator::GetIntegrationPoints() const noexcept
{
  return CoordinateView<DeviceMemorySpace>(
    csys::Cartesian::Create(static_cast<int>(coords_.extent(1))),
    MakeConstRank2View(coords_));
}

Vec OmegaHMonteCarloRHSIntegrator::GetVector() const noexcept
{
  return vec_;
}

void OmegaHMonteCarloRHSIntegrator::Assemble(
  Rank2View<const Real, DeviceMemorySpace> sampled_values)
{
  const int nbary = nbary_;
  const std::size_t num_samples =
    static_cast<std::size_t>(node_gids_.extent(0) / nbary);
  PCMS_ALWAYS_ASSERT(static_cast<std::size_t>(sampled_values.extent(0)) ==
                     num_samples);
  PCMS_ALWAYS_ASSERT(sampled_values.extent(1) >= 1);

  PetscErrorCode ierr = VecZeroEntries(vec_);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  auto sv = sampled_values;
  Kokkos::View<PetscScalar*, DeviceMemorySpace> coo_vals("mc_rhs_coo_vals",
                                                         num_samples * nbary);
  auto coeffs = coeffs_;
  Kokkos::parallel_for(
    "mc_rhs_coo_vals", static_cast<int>(num_samples), KOKKOS_LAMBDA(int i) {
      const PetscScalar f = static_cast<PetscScalar>(sv(i, 0));
      for (int j = 0; j < nbary; ++j) {
        coo_vals(i * nbary + j) =
          static_cast<PetscScalar>(coeffs(i * nbary + j)) * f;
      }
    });

  ierr = VecSetValuesCOO(vec_, coo_vals.data(), ADD_VALUES);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
}

} // namespace pcms
