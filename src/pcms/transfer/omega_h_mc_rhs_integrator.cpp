#include "pcms/transfer/omega_h_mc_rhs_integrator.hpp"
#include "pcms/transfer/omega_h_form_integrator_utils.hpp"
#include "pcms/transfer/petsc_utils.hpp"
#include "pcms/utility/assert.h"
#include <Kokkos_Random.hpp>
#include <Omega_h_shape.hpp>
#include <petscksp.h>

namespace pcms
{

namespace
{

// Maps a uniform point on the unit square to uniform barycentric coordinates
// Shape Distributions (ACM Transactions on Graphics, Vol. 21, No. 4, October
// 2002.) page 814 Eq 1
KOKKOS_INLINE_FUNCTION Omega_h::Vector<3> UniformTriangleBarycentric(Real u,
                                                                     Real v)
{
  const Real s = Kokkos::sqrt(u);
  return {1.0 - s, s * (1.0 - v), s * v};
}

// Fills coords, node_gids, and coeffs for all samples of one target element.
// unit_sample_at(s, u, v) provides the s-th unit-square sample.
template <typename UnitSampleAt>
KOKKOS_INLINE_FUNCTION void FillElementSamples(
  const int elm, const int samples_per_element,
  const Omega_h::Reals& mesh_coords, const Omega_h::LOs& faces2nodes,
  const Kokkos::View<const LO*, DeviceMemorySpace>& global_to_local,
  const Kokkos::View<Real**, DeviceMemorySpace>& coords,
  const Kokkos::View<PetscInt*, DeviceMemorySpace>& node_gids,
  const Kokkos::View<Real*, DeviceMemorySpace>& coeffs,
  const UnitSampleAt& unit_sample_at)
{
  const auto verts = Omega_h::gather_verts<3>(faces2nodes, elm);
  const auto vert_coords = Omega_h::gather_vectors<3, 2>(mesh_coords, verts);
  Omega_h::Few<Omega_h::Vector<2>, 2> basis;
  basis[0] = vert_coords[1] - vert_coords[0];
  basis[1] = vert_coords[2] - vert_coords[0];
  const Real area = Kokkos::fabs(Omega_h::triangle_area_from_basis(basis));
  const Real weight = area / samples_per_element;

  for (int s = 0; s < samples_per_element; ++s) {
    Real u = 0.0;
    Real v = 0.0;
    unit_sample_at(s, u, v);
    const auto bary = UniformTriangleBarycentric(u, v);

    const int i = elm * samples_per_element + s;
    Real x = 0.0;
    Real y = 0.0;
    for (int k = 0; k < 3; ++k) {
      x += bary[k] * vert_coords[k][0];
      y += bary[k] * vert_coords[k][1];
      node_gids(i * 3 + k) = static_cast<PetscInt>(global_to_local(verts[k]));
      coeffs(i * 3 + k) = bary[k] * weight;
    }
    coords(i, 0) = x;
    coords(i, 1) = y;
  }
}

// Samples samples_per_element from a uniform random distribution over the
// target element, writing their coordinates, node GIDs, and coefficients.
void FillElementSamples(
  int nelems, int samples_per_element, const Omega_h::Reals& mesh_coords,
  const Omega_h::LOs& faces2nodes,
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
      FillElementSamples(elm, samples_per_element, mesh_coords, faces2nodes,
                         global_to_local, coords, node_gids, coeffs,
                         [&](int /*s*/, Real& u, Real& v) {
                           u = gen.drand();
                           v = gen.drand();
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
  CoordinateSystem target_coordinate_system, int samples_per_element,
  MonteCarloSampling /*sampling*/, uint64_t seed)
{
  detail::CheckOmegaHScalarP1Layout(target_coordinate_system, target_layout,
                                    "OmegaHMonteCarloRHSIntegrator", "target");
  if (samples_per_element <= 0) {
    throw pcms_error(
      "OmegaHMonteCarloRHSIntegrator: samples_per_element must be positive");
  }

  Omega_h::Mesh& mesh = target_layout->GetMesh();

  const int nelems = mesh.nelems();
  const int num_samples = nelems * samples_per_element;
  const auto mesh_coords = mesh.coords();
  const auto faces2nodes = mesh.ask_down(Omega_h::FACE, Omega_h::VERT).ab2b;
  const auto global_to_local = target_layout->GetGlobalToLocalPermutation();

  Kokkos::View<Real**, DeviceMemorySpace> coords("mc_rhs_coords", num_samples,
                                                 2);
  Kokkos::View<PetscInt*, DeviceMemorySpace> node_gids(
    "mc_rhs_node_gids", static_cast<std::size_t>(num_samples) * 3);
  Kokkos::View<Real*, DeviceMemorySpace> coeffs(
    "mc_rhs_coeffs", static_cast<std::size_t>(num_samples) * 3);

  FillElementSamples(nelems, samples_per_element, mesh_coords, faces2nodes,
                     global_to_local, coords, node_gids, coeffs, seed);

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
  return CoordinateView<DeviceMemorySpace>(CoordinateSystem::Cartesian,
                                           MakeConstRank2View(coords_));
}

Vec OmegaHMonteCarloRHSIntegrator::GetVector() const noexcept
{
  return vec_;
}

void OmegaHMonteCarloRHSIntegrator::Assemble(
  Rank2View<const Real, DeviceMemorySpace> sampled_values)
{
  const std::size_t num_samples =
    static_cast<std::size_t>(node_gids_.extent(0) / 3);
  PCMS_ALWAYS_ASSERT(static_cast<std::size_t>(sampled_values.extent(0)) ==
                     num_samples);
  PCMS_ALWAYS_ASSERT(sampled_values.extent(1) >= 1);

  PetscErrorCode ierr = VecZeroEntries(vec_);
  CHKERRABORT(PETSC_COMM_SELF, ierr);

  auto sv = Kokkos::View<const Real**, Kokkos::LayoutRight, DeviceMemorySpace,
                         Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
    sampled_values.data_handle(), sampled_values.extent(0),
    sampled_values.extent(1));
  Kokkos::View<PetscScalar*, DeviceMemorySpace> coo_vals("mc_rhs_coo_vals",
                                                         num_samples * 3);
  auto coeffs = coeffs_;
  Kokkos::parallel_for(
    "mc_rhs_coo_vals", static_cast<int>(num_samples), KOKKOS_LAMBDA(int i) {
      const PetscScalar f = static_cast<PetscScalar>(sv(i, 0));
      for (int j = 0; j < 3; ++j) {
        coo_vals(i * 3 + j) = static_cast<PetscScalar>(coeffs(i * 3 + j)) * f;
      }
    });

  ierr = VecSetValuesCOO(vec_, coo_vals.data(), ADD_VALUES);
  CHKERRABORT(PETSC_COMM_SELF, ierr);
}

} // namespace pcms
