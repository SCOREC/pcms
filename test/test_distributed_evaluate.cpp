// Distributed evaluation, against the serial evaluator and against the closed
// form.
//
// The distributed evaluator produces values by asking each rank's local
// evaluator about the points that rank owns, so the comparison that matters is
// between a distributed answer and a serial one. The serial evaluator is local:
// it can only answer for points its own partition plus ghosts contain, and
// applies its policy to the rest. So it is asked with a sentinel fill value,
// and a sentinel reading means "not defined here" rather than a disagreement.
// To keep full coverage anyway, every point is additionally checked against the
// closed form, which a linear field on linear elements reproduces exactly.
//
// The last case runs a backend that has no distributed code at all -- the
// uniform-grid evaluator, with a different localization hint and a different
// value rule from the Lagrange one -- through the same distributed evaluator,
// to show that a backend needs nothing but local evaluation to be distributed.
//
// The query points are every owned element centroid of the whole decomposition,
// so each rank is asked about points it does not own and the remote path is
// exercised. Nothing here calls SynchronizeGhosts: a rank that reads only what
// it owns must still produce values for everyone else's points.

#include "distributed_test_utils.h"

#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_mesh.hpp>

#include "pcms/field/evaluator/distributed_point_evaluator.h"
#include "pcms/field/field.h"
#include "pcms/field/field_data.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/mesh_geometry.h"
#include "pcms/utility/mpi_utilities.h"
#include "pcms/utility/mpi_type.h"
#include "pcms/utility/uniform_grid.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <vector>

using distributed_test::require;
using pcms::LO;
using pcms::Real;

namespace
{

/// The fill value the serial reference is asked to use for points it cannot
/// reach. Distinct enough to be unmistakable in a comparison.
constexpr Real kNotDefinedHere = -1.0e300;

/// The affine test field: exactly representable on linear elements.
KOKKOS_INLINE_FUNCTION Real linear_f(Real x, Real y)
{
  return x + 2.0 * y;
}

/// Set the DOF holder coefficients from an analytic function.
///
/// A local copy of the test utility, because this test links neither Catch2 nor
/// the Catch2-based header that utility lives in.
///
/// `scale` multiplies the analytic function, so a caller can change the field's
/// coefficients without touching its geometry. A Lagrange evaluation is linear
/// in the coefficients, so the value at a fixed point scales with them.
void SetLinearField(pcms::FieldData<Real>& data,
                    const pcms::FieldLayout& layout, Real scale = 1.0)
{
  const auto dof_coords = layout.GetDOFHolderCoordinates().GetValues();
  const LO n = static_cast<LO>(dof_coords.extent(0));
  Kokkos::View<Real**, pcms::DeviceMemorySpace> values("field values", n, 1);
  Kokkos::parallel_for(
    "set_linear_field",
    Kokkos::RangePolicy<pcms::DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(LO i) {
      values(i, 0) = scale * linear_f(dof_coords(i, 0), dof_coords(i, 1));
    });
  data.SetDOFHolderData(pcms::MakeRank2View(values));
}

/// Query coordinates on the device with their coordinate-system wrapper.
struct QueryCoords
{
  Kokkos::View<Real**, pcms::DeviceMemorySpace> storage;
  pcms::CoordinateView<pcms::DeviceMemorySpace> view;
};

QueryCoords MakeQueryCoords(const std::vector<Real>& interleaved, int dim)
{
  const LO n = static_cast<LO>(interleaved.size()) / dim;
  Kokkos::View<Real**, pcms::HostMemorySpace> host("query host", n, dim);
  for (LO i = 0; i < n; ++i) {
    for (int d = 0; d < dim; ++d) {
      host(i, d) =
        interleaved[static_cast<size_t>(i) * static_cast<size_t>(dim) +
                    static_cast<size_t>(d)];
    }
  }
  Kokkos::View<Real**, pcms::DeviceMemorySpace> device("query device", n, dim);
  pcms::DeepCopyMismatchLayouts(device, host);
  return QueryCoords{
    device, pcms::CoordinateView<pcms::DeviceMemorySpace>(
              pcms::CoordinateSystem::Cartesian, pcms::MakeRank2View(device))};
}

/// Distributed evaluation must equal the serial evaluator wherever the serial
/// evaluator is defined, and must equal the closed form everywhere.
template <int Dim>
void check_agrees_with_serial(Omega_h::Mesh& mesh, int order)
{
  const auto dim = mesh.dim();
  const auto rank = mesh.comm()->rank();

  auto space = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, order, 1, pcms::CoordinateSystem::Cartesian, "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
  auto field = space->CreateFunction<Real>();
  SetLinearField(field.GetData(), *space->GetLayout());

  const auto gathered = distributed_test::GatherOwnedCentroids(mesh);
  require(gathered.total > 0, "no owned element centroids were gathered");

  auto device_coords = MakeQueryCoords(gathered.coords, dim);
  const LO n_points = static_cast<LO>(gathered.total);

  // Serial reference, through the public evaluator API, with a fill value that
  // marks the points it cannot reach.
  pcms::OutOfBoundsPolicy serial_policy;
  serial_policy.mode = pcms::OutOfBoundsMode::FILL;
  serial_policy.fill_value = kNotDefinedHere;
  auto serial =
    space->CreatePointEvaluator<Real>(pcms::EvaluationRequest::FromCoordinates(
      device_coords.view, serial_policy));
  Kokkos::View<Real**, pcms::DeviceMemorySpace> reference("reference", n_points,
                                                          1);
  serial->Evaluate(field, pcms::MakeRank2View(reference));
  const auto reference_h =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), reference);

  // Distributed, from the same space and points. The space, not a layout, is
  // what is handed over: the distributed evaluator reaches a backend only
  // through the local evaluator the space creates.
  auto distributed = pcms::MakeDistributedPointEvaluator<Real>(
    mesh, space, device_coords.storage, pcms::OutOfBoundsPolicy{});

  Kokkos::View<Real**, pcms::DeviceMemorySpace> got("got", n_points, 1);
  distributed->Evaluate(field, pcms::MakeRank2View(got));
  const std::size_t routed = pcms::LastRoutedPointCount(*distributed);
  const auto got_h =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), got);

  LO n_compared = 0;
  LO n_not_local = 0;
  LO n_mismatch = 0;
  LO n_bad_vs_closed_form = 0;
  Real worst = 0;
  Real worst_vs_closed_form = 0;
  LO shown = 0;
  for (LO i = 0; i < n_points; ++i) {
    const Real value = got_h(i, 0);
    const Real reference_value = reference_h(i, 0);
    const size_t base = static_cast<size_t>(i) * static_cast<size_t>(dim);
    const Real truth =
      linear_f(gathered.coords[base], gathered.coords[base + 1]);

    // Independent of the serial evaluator, so this covers the routed points
    // too, which the serial one may not be able to answer at all.
    const Real closed_form_error = std::fabs(value - truth);
    worst_vs_closed_form = closed_form_error > worst_vs_closed_form
                             ? closed_form_error
                             : worst_vs_closed_form;
    if (closed_form_error > 1e-10) {
      ++n_bad_vs_closed_form;
    }

    if (reference_value == kNotDefinedHere) {
      ++n_not_local;
      continue;
    }
    ++n_compared;
    const Real difference = std::fabs(value - reference_value);
    worst = difference > worst ? difference : worst;
    if (difference > 1e-10) {
      ++n_mismatch;
      if (shown < 4) {
        std::printf(
          "  point (%g,%g): distributed=%.9f serial=%.9f exact=%.9f\n",
          static_cast<double>(gathered.coords[base]),
          static_cast<double>(gathered.coords[base + 1]),
          static_cast<double>(value), static_cast<double>(reference_value),
          static_cast<double>(truth));
        ++shown;
      }
    }
  }

  const LO t_compared =
    mesh.comm()->allreduce(n_compared, Omega_h_Op::OMEGA_H_SUM);
  const LO t_not_local =
    mesh.comm()->allreduce(n_not_local, Omega_h_Op::OMEGA_H_SUM);
  const LO t_mismatch =
    mesh.comm()->allreduce(n_mismatch, Omega_h_Op::OMEGA_H_SUM);
  const LO t_bad =
    mesh.comm()->allreduce(n_bad_vs_closed_form, Omega_h_Op::OMEGA_H_SUM);
  const LO t_points = mesh.comm()->allreduce(n_points, Omega_h_Op::OMEGA_H_SUM);
  const LO t_routed =
    mesh.comm()->allreduce(static_cast<LO>(routed), Omega_h_Op::OMEGA_H_SUM);
  const Real t_worst = mesh.comm()->allreduce(worst, Omega_h_Op::OMEGA_H_MAX);
  const Real t_worst_closed =
    mesh.comm()->allreduce(worst_vs_closed_form, Omega_h_Op::OMEGA_H_MAX);

  if (rank == 0) {
    std::printf(
      "%dD order %d: %lld points, %lld routed to an owner, %lld comparable "
      "locally (%lld beyond it), %lld mismatches (worst %.3g)\n",
      Dim, order, static_cast<long long>(t_points),
      static_cast<long long>(t_routed), static_cast<long long>(t_compared),
      static_cast<long long>(t_not_local), static_cast<long long>(t_mismatch),
      static_cast<double>(t_worst));
    std::printf("  against the closed form: %lld bad, worst %.3g\n",
                static_cast<long long>(t_bad),
                static_cast<double>(t_worst_closed));
  }

  require(t_mismatch == 0,
          "distributed evaluation disagrees with the serial evaluator");
  require(t_bad == 0,
          "distributed evaluation does not reproduce a linear field");
  // A comparison is only meaningful if the two paths actually overlapped, and
  // the remote path is only covered if something had to be routed.
  require(t_compared > 0, "no point was localizable by the serial evaluator");
  require(t_not_local > 0,
          "the serial evaluator reached every point, so the sentinel never "
          "marked anything and the comparison is untested at its boundary");
  require(t_routed > 0,
          "every point was already local, so the remote path was not used");

  // A second Evaluate on the same points has to recompute the values, because
  // the coefficients change even when the location does not. Scaling the field
  // is how the two are told apart: a value left over from the first call, or
  // read out of the payload's zero fill, comes back unscaled. Recomputing is
  // not free for the points another rank owns -- only its owner holds the
  // coefficients -- so the same points are put to their owners again, and the
  // routed count repeats.
  const Real scale = 3.0;
  SetLinearField(field.GetData(), *space->GetLayout(), scale);
  Kokkos::View<Real**, pcms::DeviceMemorySpace> second("second", n_points, 1);
  distributed->Evaluate(field, pcms::MakeRank2View(second));
  const LO routed_again =
    static_cast<LO>(pcms::LastRoutedPointCount(*distributed));
  const auto second_h =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), second);

  LO n_not_scaled = 0;
  for (LO i = 0; i < n_points; ++i) {
    const Real expected = scale * got_h(i, 0);
    if (std::fabs(second_h(i, 0) - expected) >
        1e-10 * (1.0 + std::fabs(expected))) {
      ++n_not_scaled;
    }
  }
  const LO t_not_scaled =
    mesh.comm()->allreduce(n_not_scaled, Omega_h_Op::OMEGA_H_SUM);
  const LO t_routed_again =
    mesh.comm()->allreduce(routed_again, Omega_h_Op::OMEGA_H_SUM);

  if (rank == 0) {
    std::printf(
      "  second pass on the same points: %lld points re-asked of their owners, "
      "%lld stale values after the field was scaled by %g\n",
      static_cast<long long>(t_routed_again),
      static_cast<long long>(t_not_scaled), static_cast<double>(scale));
  }
  // The same points need an owner's answer on both passes: the location did not
  // change, but the value has to come from the coefficients as they are now.
  require(t_routed_again == t_routed,
          "a repeated evaluation did not re-ask the owners of the same points");
  require(t_not_scaled == 0,
          "a repeated evaluation returned a value computed from the old "
          "coefficients");
}

/// Order 0 is piecewise constant, so a point on a shared edge has as many
/// candidate values as there are elements touching it, and they differ. Exactly
/// one answer is correct for a distributed evaluator, and it has to be the same
/// answer on every rank, or two ranks disagree about a shared boundary.
///
/// This checks the arbitration rule on such points: the result is identical on
/// every rank, it equals a real element's coefficient rather than an average or
/// a stale value, and the query set really did contain contested points.
template <int Dim>
void check_boundary_arbitration(Omega_h::Mesh& mesh)
{
  const auto dim = mesh.dim();
  const auto rank = mesh.comm()->rank();
  const int nranks = mesh.comm()->size();
  auto comm = mesh.comm()->get_impl();

  // Order 0: one coefficient per element, set from the element's own centroid,
  // so neighbouring elements hold different values.
  auto space = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, /*order=*/0, 1, pcms::CoordinateSystem::Cartesian, "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
  auto field = space->CreateFunction<Real>();
  SetLinearField(field.GetData(), *space->GetLayout());

  const auto gathered =
    distributed_test::GatherOwnedElementVertices(mesh);
  require(gathered.total > 0, "no owned element vertices were gathered");
  auto device_coords = MakeQueryCoords(gathered.coords, dim);
  const LO n_points = static_cast<LO>(gathered.total);

  // Every rank's own answer, with the sentinel for the points it cannot reach.
  pcms::OutOfBoundsPolicy serial_policy;
  serial_policy.mode = pcms::OutOfBoundsMode::FILL;
  serial_policy.fill_value = kNotDefinedHere;
  auto serial =
    space->CreatePointEvaluator<Real>(pcms::EvaluationRequest::FromCoordinates(
      device_coords.view, serial_policy));
  Kokkos::View<Real**, pcms::DeviceMemorySpace> reference("reference", n_points,
                                                          1);
  serial->Evaluate(field, pcms::MakeRank2View(reference));
  const auto reference_h =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), reference);

  auto distributed = pcms::MakeDistributedPointEvaluator<Real>(
    mesh, space, device_coords.storage, pcms::OutOfBoundsPolicy{});
  Kokkos::View<Real**, pcms::DeviceMemorySpace> got("got", n_points, 1);
  distributed->Evaluate(field, pcms::MakeRank2View(got));
  const auto got_h =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), got);

  // Gather everyone's distributed answer and everyone's serial answer, so both
  // can be compared point by point across ranks.
  const auto gather_all = [&](const Real* local) {
    std::vector<Real> all(static_cast<size_t>(nranks) *
                          static_cast<size_t>(n_points));
    std::vector<int> counts(static_cast<size_t>(nranks), 0);
    std::vector<int> displs(static_cast<size_t>(nranks), 0);
    for (int r = 0; r < nranks; ++r) {
      counts[static_cast<size_t>(r)] = static_cast<int>(n_points);
      displs[static_cast<size_t>(r)] = static_cast<int>(n_points) * r;
    }
    MPI_Allgatherv(local, static_cast<int>(n_points), pcms::GetMPIType(Real{}),
                   all.data(), counts.data(), displs.data(),
                   pcms::GetMPIType(Real{}), comm);
    return all;
  };
  const auto all_distributed = gather_all(got_h.data());
  const auto all_serial = gather_all(reference_h.data());

  LO n_rank_disagreement = 0;
  LO n_not_an_element = 0;
  LO n_contested = 0;
  for (LO i = 0; i < n_points; ++i) {
    const Real value = got_h(i, 0);

    // The same point must be answered identically everywhere.
    for (int r = 1; r < nranks; ++r) {
      const Real other =
        all_distributed[static_cast<size_t>(r) * static_cast<size_t>(n_points) +
                        static_cast<size_t>(i)];
      if (std::fabs(other - value) > 1e-12) {
        ++n_rank_disagreement;
      }
    }

    // The value must be some element's coefficient, so it has to appear among
    // the serial answers of the ranks that could locate it.
    bool found = false;
    LO distinct = 0;
    Real first = 0;
    for (int r = 0; r < nranks; ++r) {
      const Real candidate =
        all_serial[static_cast<size_t>(r) * static_cast<size_t>(n_points) +
                   static_cast<size_t>(i)];
      if (candidate == kNotDefinedHere) {
        continue;
      }
      if (std::fabs(candidate - value) <= 1e-12) {
        found = true;
      }
      if (distinct == 0 || std::fabs(candidate - first) > 1e-12) {
        first = candidate;
        ++distinct;
      }
    }
    if (!found) {
      ++n_not_an_element;
    }
    // Two ranks offering different coefficients for one point is what makes the
    // arbitration rule load-bearing.
    if (distinct > 1) {
      ++n_contested;
    }
  }

  const auto reduce_sum = [&](LO value) {
    return mesh.comm()->allreduce(value, Omega_h_Op::OMEGA_H_SUM);
  };
  const LO t_disagreement = reduce_sum(n_rank_disagreement);
  const LO t_not_an_element = reduce_sum(n_not_an_element);
  const LO t_contested = reduce_sum(n_contested);
  const LO t_points = reduce_sum(n_points);

  if (rank == 0) {
    std::printf(
      "%dD order 0 arbitration: %lld boundary points, %lld contested, %lld "
      "rank disagreements, %lld not an element coefficient\n",
      Dim, static_cast<long long>(t_points),
      static_cast<long long>(t_contested),
      static_cast<long long>(t_disagreement),
      static_cast<long long>(t_not_an_element));
  }

  // A vacuous pass would be a test that never saw a contested point.
  require(t_contested > 0,
          "no query point was contested, so the arbitration rule is untested");
  require(t_disagreement == 0, "ranks disagree about a shared boundary point");
  require(t_not_an_element == 0,
          "a distributed value is not any element's coefficient");
}

/// A backend with no distributed code of its own still gets distributed
/// evaluation, because the distributed evaluator reaches a backend only through
/// the local evaluator its `FunctionSpace` creates.
///
/// The backend here is the uniform-grid one, which is as different from
/// Lagrange as it gets: its localization hint is a cell index plus parametric
/// coordinates, not an element plus barycentric weights, and its value rule is
/// multilinear interpolation over a structured grid, not a Lagrange combination
/// over a simplex. None of that reaches the distributed layer, which hands the
/// coordinates of the points a rank owns to that rank's own evaluator and
/// carries the result back as bytes.
///
/// The grid is global and identical on every rank, so an owner's answer is the
/// answer and the two paths must agree exactly. The remote path is still
/// exercised, and the closed form still holds, because multilinear
/// interpolation reproduces an affine field exactly on a tensor-product grid.
template <int Dim>
void check_backend_without_distributed_code(Omega_h::Mesh& mesh)
{
  const auto dim = mesh.dim();
  const auto rank = mesh.comm()->rank();

  // Both types are spelled out rather than deduced: the grid's type depends on
  // `Dim`, and a dependent `auto` would make `space` dependent too, which stops
  // the compiler from recognizing `CreateFunction<Real>` below as a
  // template-id before instantiation.
  using Grid = pcms::UniformGrid<static_cast<unsigned>(Dim)>;
  const Grid grid =
    pcms::CreateUniformGridFromMesh<static_cast<unsigned>(Dim)>(mesh, 8);
  const std::shared_ptr<pcms::LagrangeFunctionSpace> space =
    pcms::LagrangeFunctionSpace::FromUniformGrid(
      grid, 1, pcms::CoordinateSystem::Cartesian, /*order=*/1);
  auto field = space->CreateFunction<Real>();
  SetLinearField(field.GetData(), *space->GetLayout());

  const auto gathered = distributed_test::GatherOwnedCentroids(mesh);
  require(gathered.total > 0, "no owned element centroids were gathered");
  auto device_coords = MakeQueryCoords(gathered.coords, dim);
  const LO n_points = static_cast<LO>(gathered.total);

  // The serial evaluator for this backend, with a sentinel for the points it
  // cannot reach.
  pcms::OutOfBoundsPolicy serial_policy;
  serial_policy.mode = pcms::OutOfBoundsMode::FILL;
  serial_policy.fill_value = kNotDefinedHere;
  auto serial =
    space->CreatePointEvaluator<Real>(pcms::EvaluationRequest::FromCoordinates(
      device_coords.view, serial_policy));
  Kokkos::View<Real**, pcms::DeviceMemorySpace> reference("reference", n_points,
                                                          1);
  serial->Evaluate(field, pcms::MakeRank2View(reference));
  const auto reference_h =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), reference);

  // The same backend, distributed.
  auto distributed = pcms::MakeDistributedPointEvaluator<Real>(
    mesh, space, device_coords.storage, pcms::OutOfBoundsPolicy{});
  Kokkos::View<Real**, pcms::DeviceMemorySpace> got("got", n_points, 1);
  distributed->Evaluate(field, pcms::MakeRank2View(got));
  const std::size_t routed = pcms::LastRoutedPointCount(*distributed);
  const auto got_h =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), got);

  LO n_mismatch = 0;
  LO n_serial_reached = 0;
  LO n_bad_closed_form = 0;
  Real worst = 0;
  Real worst_closed_form = 0;
  for (LO i = 0; i < n_points; ++i) {
    const Real value = got_h(i, 0);
    const size_t base = static_cast<size_t>(i) * static_cast<size_t>(dim);
    const Real truth =
      linear_f(gathered.coords[base], gathered.coords[base + 1]);
    const Real closed_form_error = std::fabs(value - truth);
    worst_closed_form = closed_form_error > worst_closed_form
                          ? closed_form_error
                          : worst_closed_form;
    if (closed_form_error > 1e-10) {
      ++n_bad_closed_form;
    }

    if (reference_h(i, 0) == kNotDefinedHere) {
      continue;
    }
    ++n_serial_reached;
    const Real difference = std::fabs(value - reference_h(i, 0));
    worst = difference > worst ? difference : worst;
    if (difference > 1e-10) {
      ++n_mismatch;
    }
  }

  const LO t_mismatch =
    mesh.comm()->allreduce(n_mismatch, Omega_h_Op::OMEGA_H_SUM);
  const LO t_serial_reached =
    mesh.comm()->allreduce(n_serial_reached, Omega_h_Op::OMEGA_H_SUM);
  const LO t_bad =
    mesh.comm()->allreduce(n_bad_closed_form, Omega_h_Op::OMEGA_H_SUM);
  const LO t_points = mesh.comm()->allreduce(n_points, Omega_h_Op::OMEGA_H_SUM);
  const LO t_routed =
    mesh.comm()->allreduce(static_cast<LO>(routed), Omega_h_Op::OMEGA_H_SUM);
  const Real t_worst = mesh.comm()->allreduce(worst, Omega_h_Op::OMEGA_H_MAX);
  const Real t_worst_closed =
    mesh.comm()->allreduce(worst_closed_form, Omega_h_Op::OMEGA_H_MAX);

  if (rank == 0) {
    std::printf(
      "%dD uniform-grid backend, no distributed code: %lld points, %lld "
      "routed, %lld comparable, %lld mismatches (worst %.3g)\n",
      Dim, static_cast<long long>(t_points), static_cast<long long>(t_routed),
      static_cast<long long>(t_serial_reached),
      static_cast<long long>(t_mismatch), static_cast<double>(t_worst));
    std::printf("  against the closed form: %lld bad, worst %.3g\n",
                static_cast<long long>(t_bad),
                static_cast<double>(t_worst_closed));
  }

  require(t_mismatch == 0,
          "the distributed path disagrees with the second backend's own "
          "evaluator");
  require(t_bad == 0,
          "the second backend does not reproduce a linear field when "
          "distributed");
  require(t_routed > 0,
          "every point was already local, so the remote path was not used with "
          "the second backend");
  // The sentinel must not be hiding anything: this backend's grid spans the
  // whole mesh, so its own evaluator reaches every query point.
  require(t_serial_reached == t_points,
          "the second backend's serial evaluator could not reach every query "
          "point, so the comparison is incomplete");
}

/// The policy is applied to exactly the points nothing answered for, and it is
/// applied after the exchange, so ERROR cannot leave the other ranks waiting.
///
/// Two things are worth pinning down here. A point no rank can locate is the
/// only one the fill touches, so a located point is never overwritten by it
/// even when the location and the payload disagree (which is what the no-answer
/// marker exists to detect). And ERROR has to be raised once the exchange is
/// finished: the point set is the same on every rank, so every rank raises and
/// none is left blocked in a collective.
template <int Dim>
void check_out_of_bounds_policy(Omega_h::Mesh& mesh)
{
  const auto dim = mesh.dim();

  auto space = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, /*order=*/1, 1, pcms::CoordinateSystem::Cartesian, "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
  auto field = space->CreateFunction<Real>();
  SetLinearField(field.GetData(), *space->GetLayout());

  const auto gathered = distributed_test::GatherOwnedCentroids(mesh);
  require(gathered.total > 0, "no owned element centroids were gathered");

  // The interior points, plus one far outside the mesh. The outside point is
  // outside the routing grid too, so it is never routed.
  std::vector<Real> coords = gathered.coords;
  for (int d = 0; d < dim; ++d) {
    coords.push_back(-10.0);
  }
  auto device_coords = MakeQueryCoords(coords, dim);
  const LO n_points = static_cast<LO>(gathered.total) + 1;
  const LO outside_point = n_points - 1;
  const Real fill_value = -1234.5;

  pcms::OutOfBoundsPolicy fill_policy;
  fill_policy.mode = pcms::OutOfBoundsMode::FILL;
  fill_policy.fill_value = fill_value;
  auto filled = pcms::MakeDistributedPointEvaluator<Real>(
    mesh, space, device_coords.storage, fill_policy);
  Kokkos::View<Real**, pcms::DeviceMemorySpace> got("got", n_points, 1);
  filled->Evaluate(field, pcms::MakeRank2View(got));
  const auto got_h =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), got);

  LO n_interior_bad = 0;
  LO n_outside_filled = 0;
  for (LO i = 0; i < n_points; ++i) {
    if (i == outside_point) {
      if (got_h(i, 0) == fill_value) {
        ++n_outside_filled;
      }
      continue;
    }
    const size_t base = static_cast<size_t>(i) * static_cast<size_t>(dim);
    const Real truth =
      linear_f(gathered.coords[base], gathered.coords[base + 1]);
    if (std::fabs(got_h(i, 0) - truth) > 1e-10) {
      ++n_interior_bad;
    }
  }

  const LO t_interior_bad =
    mesh.comm()->allreduce(n_interior_bad, Omega_h_Op::OMEGA_H_SUM);
  const LO t_outside_filled =
    mesh.comm()->allreduce(n_outside_filled, Omega_h_Op::OMEGA_H_SUM);
  const LO t_points = mesh.comm()->allreduce(n_points, Omega_h_Op::OMEGA_H_SUM);

  if (mesh.comm()->rank() == 0) {
    std::printf(
      "%dD FILL policy: %lld points, %lld filled for the point outside the "
      "mesh, %lld interior answers wrong\n",
      Dim, static_cast<long long>(t_points),
      static_cast<long long>(t_outside_filled),
      static_cast<long long>(t_interior_bad));
  }

  require(t_interior_bad == 0,
          "a point inside the mesh was treated as out of bounds");
  // The query set is replicated, so every rank has the one outside point and
  // every rank must have filled it.
  require(t_outside_filled == static_cast<LO>(mesh.comm()->size()),
          "the point outside the mesh was not filled on every rank");

  // ERROR, on the same point set: every rank has the point, so every rank
  // raises, after the exchange rather than inside it.
  pcms::OutOfBoundsPolicy error_policy;
  error_policy.mode = pcms::OutOfBoundsMode::ERROR;
  bool threw = false;
  try {
    auto strict = pcms::MakeDistributedPointEvaluator<Real>(
      mesh, space, device_coords.storage, error_policy);
    Kokkos::View<Real**, pcms::DeviceMemorySpace> strict_values("strict",
                                                                n_points, 1);
    strict->Evaluate(field, pcms::MakeRank2View(strict_values));
  } catch (const pcms::pcms_error&) {
    threw = true;
  }
  require(threw, "the ERROR policy did not report a point outside the mesh");
}

} // namespace

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  {
    Omega_h::Library library(&argc, &argv);
    auto world = library.world();
    require(world->size() >= 2, "test requires at least 2 MPI ranks");

    {
      auto mesh = distributed_test::BuildGhostedBox<2>(world);
      check_agrees_with_serial<2>(mesh, /*order=*/1);
      check_agrees_with_serial<2>(mesh, /*order=*/0);
      check_boundary_arbitration<2>(mesh);
      check_backend_without_distributed_code<2>(mesh);
      check_out_of_bounds_policy<2>(mesh);
    }

    {
      auto mesh = distributed_test::BuildGhostedBox<3>(world);
      check_agrees_with_serial<3>(mesh, /*order=*/1);
      check_out_of_bounds_policy<3>(mesh);
    }

    if (world->rank() == 0) {
      std::printf("distributed evaluation test passed\n");
    }
  }
  MPI_Finalize();
  return 0;
}
