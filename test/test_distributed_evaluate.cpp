// Distributed evaluation, against the serial evaluator and against the closed
// form.

#include "distributed_test_utils.hpp"

#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_mesh.hpp>

#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/backends/omega_h/omega_h_geometry.hpp"
#include "pcms/field/field.hpp"
#include "pcms/field/distributed_evaluator.hpp"
#include "pcms/geometry/distributed_binding.hpp"
#include "pcms/field/field_data.hpp"
#include "pcms/field/function_spaces.hpp"
#include "pcms/field/query_points.hpp"
#include "pcms/transfer/interpolator.hpp"
#include "pcms/utility/arrays.hpp"
#include "pcms/utility/detail/arrays.hpp"
#include "pcms/utility/coordinate_systems/cartesian.hpp"
#include "pcms/utility/mpi_type.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <vector>

using distributed_test::require;
using pcms::LO;
using pcms::Real;

namespace
{

// The fill value the serial reference uses for points it cannot reach.
constexpr Real kNotDefinedHere = -1.0e300;

// The affine test field: exactly representable on linear elements.
KOKKOS_INLINE_FUNCTION Real linear_f(Real x, Real y)
{
  return x + 2.0 * y;
}

// Set the DOF holder values from linear_f times `scale`, and synchronize.
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
  data.GetDOFAccess()->SetValues(pcms::MakeConstRank2View(values));
  data.Synchronize();
}

// Query coordinates on the device with their coordinate-system wrapper.
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
              pcms::csys::CreateCartesian(dim), pcms::MakeRank2View(device))};
}

// The local geometry of a mesh and the distributed geometry on it.
struct Geometries
{
  std::shared_ptr<const pcms::omega_h::OmegaHGeometry> local;
  std::shared_ptr<const pcms::DistributedGeometry> distributed;
};

Geometries MakeGeometries(Omega_h::Mesh& mesh)
{
  auto local = std::dynamic_pointer_cast<const pcms::omega_h::OmegaHGeometry>(
    pcms::omega_h::CreateGeometry(mesh,
                                  pcms::csys::CreateCartesian(mesh.dim())));
  auto distributed = pcms::omega_h::CreateDistributedGeometry(local);
  return {std::move(local), std::move(distributed)};
}

std::shared_ptr<const pcms::FunctionSpace> MakeSpace(
  std::shared_ptr<const pcms::Geometry> geometry, int order)
{
  pcms::LagrangeOptions options;
  options.backend = pcms::LagrangeBackend::OmegaH;
  return pcms::CreateLagrangeSpace(std::move(geometry), order, 1, options);
}

using HostValues = Kokkos::View<Real**, pcms::DeviceMemorySpace>::HostMirror;

// The values of `evaluator` at its `n` points, on the host.
HostValues EvaluateToHost(const pcms::PointEvaluator<Real>& evaluator,
                          const pcms::Field<Real>& field, LO n, int nc = 1)
{
  Kokkos::View<Real**, pcms::DeviceMemorySpace> got("got", n, nc);
  evaluator.Evaluate(field,
                     pcms::ValueView<Real, pcms::DeviceMemorySpace>(
                       evaluator.OutputBasis(field.GetData().GetValueBasis()),
                       pcms::MakeRank2View(got)));
  return Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), got);
}

// The serial evaluation of `field`, of `space` on the local geometry, with
// kNotDefinedHere at the points this rank cannot reach.
HostValues SerialReference(
  const pcms::FunctionSpace& space, const pcms::Field<Real>& field,
  const pcms::CoordinateView<pcms::DeviceMemorySpace>& coords)
{
  pcms::OutOfBoundsPolicy policy;
  policy.mode = pcms::OutOfBoundsMode::FILL;
  policy.fill_value = kNotDefinedHere;
  const auto bound = space.GetGeometry()->Bind(coords, policy);
  return EvaluateToHost(*space.CreatePointEvaluator<Real>(bound), field,
                        bound->NumPoints());
}

LO RoutedPointCount(const std::shared_ptr<const pcms::BoundPoints>& bound)
{
  const auto distributed =
    std::dynamic_pointer_cast<const pcms::DistributedBinding>(bound);
  require(distributed != nullptr, "not a distributed bound points");
  return static_cast<LO>(distributed->RoutedPointCount());
}

// Distributed evaluation equals the serial evaluator wherever that is defined,
// and the closed form everywhere.
template <int Dim>
void check_agrees_with_serial(Omega_h::Mesh& mesh, const Geometries& g,
                              int order)
{
  const auto comm = mesh.comm();
  auto local_space = MakeSpace(g.local, order);
  auto space = MakeSpace(g.local, order);
  const auto source = pcms::CreateDistributedEvaluator(g.distributed, space);
  auto local_field = local_space->CreateFunction<Real>();
  auto field = space->CreateFunction<Real>();
  SetLinearField(local_field.GetData(), *local_space->GetLayout());
  SetLinearField(field.GetData(), *space->GetLayout());

  const auto gathered = distributed_test::GatherOwnedCentroids(mesh);
  require(gathered.total > 0, "no owned element centroids were gathered");
  auto coords = MakeQueryCoords(gathered.coords, Dim);
  const LO n_points = static_cast<LO>(gathered.total);

  const auto reference =
    SerialReference(*local_space, local_field, coords.view);

  const auto bound = g.distributed->Bind(coords.view);
  auto evaluator = source->CreatePointEvaluator<Real>(bound);
  const auto got = EvaluateToHost(*evaluator, field, n_points);

  LO n_compared = 0, n_not_local = 0, n_mismatch = 0, n_bad = 0;
  for (LO i = 0; i < n_points; ++i) {
    const size_t base = static_cast<size_t>(i) * Dim;
    const Real truth =
      linear_f(gathered.coords[base], gathered.coords[base + 1]);
    n_bad += std::fabs(got(i, 0) - truth) > 1e-10 ? 1 : 0;
    if (reference(i, 0) == kNotDefinedHere) {
      ++n_not_local;
      continue;
    }
    ++n_compared;
    n_mismatch += std::fabs(got(i, 0) - reference(i, 0)) > 1e-10 ? 1 : 0;
  }
  const auto sum = [&](LO v) {
    return comm->allreduce(v, Omega_h_Op::OMEGA_H_SUM);
  };
  const LO t_routed = sum(RoutedPointCount(bound));
  const LO t_points = sum(n_points), t_compared = sum(n_compared),
           t_not_local = sum(n_not_local), t_mismatch = sum(n_mismatch),
           t_bad = sum(n_bad);
  if (comm->rank() == 0) {
    std::printf(
      "%dD order %d: %lld points, %lld routed, %lld comparable "
      "locally (%lld beyond), %lld mismatches, %lld bad\n",
      Dim, order, static_cast<long long>(t_points),
      static_cast<long long>(t_routed), static_cast<long long>(t_compared),
      static_cast<long long>(t_not_local), static_cast<long long>(t_mismatch),
      static_cast<long long>(t_bad));
  }
  require(t_mismatch == 0,
          "distributed evaluation disagrees with the serial evaluator");
  require(t_bad == 0,
          "distributed evaluation does not reproduce a linear field");
  require(t_compared > 0, "no point was reachable by the serial one");
  require(t_not_local > 0, "the serial evaluator reached every point");
  require(t_routed > 0, "no point was routed, so the remote path is untested");

  // A second Evaluate on the same points reads the coefficients as they are.
  const Real scale = 3.0;
  SetLinearField(field.GetData(), *space->GetLayout(), scale);
  const auto second = EvaluateToHost(*evaluator, field, n_points);
  LO n_not_scaled = 0;
  for (LO i = 0; i < n_points; ++i) {
    const Real expected = scale * got(i, 0);
    n_not_scaled +=
      std::fabs(second(i, 0) - expected) > 1e-10 * (1.0 + std::fabs(expected))
        ? 1
        : 0;
  }
  require(sum(n_not_scaled) == 0,
          "a repeated evaluation used the old coefficients");
}

// Order 0 is piecewise constant, so a point on a shared edge has a value per
// element touching it. Every rank must give the same answer, and it must be
// one of those elements' values.
template <int Dim>
void check_boundary_arbitration(Omega_h::Mesh& mesh, const Geometries& g)
{
  const auto comm = mesh.comm();
  const int nranks = comm->size();
  auto local_space = MakeSpace(g.local, 0);
  auto space = MakeSpace(g.local, 0);
  const auto source = pcms::CreateDistributedEvaluator(g.distributed, space);
  auto local_field = local_space->CreateFunction<Real>();
  auto field = space->CreateFunction<Real>();
  SetLinearField(local_field.GetData(), *local_space->GetLayout());
  SetLinearField(field.GetData(), *space->GetLayout());

  const auto gathered = distributed_test::GatherOwnedElementVertices(mesh);
  require(gathered.total > 0, "no owned element vertices were gathered");
  auto coords = MakeQueryCoords(gathered.coords, Dim);
  const LO n_points = static_cast<LO>(gathered.total);

  const auto reference =
    SerialReference(*local_space, local_field, coords.view);
  const auto bound = g.distributed->Bind(coords.view);
  const auto got =
    EvaluateToHost(*source->CreatePointEvaluator<Real>(bound), field, n_points);

  // The query set is the same on every rank.
  const auto gather_all = [&](const Real* local) {
    std::vector<Real> all(static_cast<size_t>(nranks) *
                          static_cast<size_t>(n_points));
    MPI_Allgather(local, static_cast<int>(n_points),
                  pcms::utility::GetMPIType(Real{}), all.data(),
                  static_cast<int>(n_points), pcms::utility::GetMPIType(Real{}),
                  comm->get_impl());
    return all;
  };
  const auto all_got = gather_all(got.data());
  const auto all_serial = gather_all(reference.data());

  LO n_disagree = 0, n_not_an_element = 0, n_contested = 0;
  for (LO i = 0; i < n_points; ++i) {
    const Real value = got(i, 0);
    for (int r = 0; r < nranks; ++r) {
      const Real other =
        all_got[static_cast<size_t>(r) * static_cast<size_t>(n_points) +
                static_cast<size_t>(i)];
      n_disagree += std::fabs(other - value) > 1e-12 ? 1 : 0;
    }
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
      found = found || std::fabs(candidate - value) <= 1e-12;
      if (distinct == 0 || std::fabs(candidate - first) > 1e-12) {
        first = candidate;
        ++distinct;
      }
    }
    n_not_an_element += found ? 0 : 1;
    n_contested += distinct > 1 ? 1 : 0;
  }
  const auto sum = [&](LO v) {
    return comm->allreduce(v, Omega_h_Op::OMEGA_H_SUM);
  };
  const LO t_contested = sum(n_contested), t_disagree = sum(n_disagree),
           t_not_an_element = sum(n_not_an_element);
  if (comm->rank() == 0) {
    std::printf("%dD order 0 arbitration: %lld contested, %lld rank "
                "disagreements, %lld not an element value\n",
                Dim, static_cast<long long>(t_contested),
                static_cast<long long>(t_disagree),
                static_cast<long long>(t_not_an_element));
  }
  require(t_contested > 0, "no point was contested");
  require(t_disagree == 0, "ranks disagree about a shared point");
  require(t_not_an_element == 0,
          "a distributed value is not any element's value");
}

// FILL fills exactly the points no rank found; ERROR throws at Bind on every
// rank.
template <int Dim>
void check_out_of_bounds_policy(Omega_h::Mesh& mesh, const Geometries& g)
{
  const auto comm = mesh.comm();
  auto space = MakeSpace(g.local, 1);
  const auto source = pcms::CreateDistributedEvaluator(g.distributed, space);
  auto field = space->CreateFunction<Real>();
  SetLinearField(field.GetData(), *space->GetLayout());

  const auto gathered = distributed_test::GatherOwnedCentroids(mesh);
  std::vector<Real> points = gathered.coords;
  for (int d = 0; d < Dim; ++d) {
    points.push_back(-10.0);
  }
  auto coords = MakeQueryCoords(points, Dim);
  const LO n_points = static_cast<LO>(gathered.total) + 1;
  const Real fill_value = -1234.5;

  pcms::OutOfBoundsPolicy fill;
  fill.mode = pcms::OutOfBoundsMode::FILL;
  fill.fill_value = fill_value;
  const auto bound = g.distributed->Bind(coords.view, fill);
  require(bound->MissingPoints().extent(0) == 1,
          "the point outside the mesh is not the only missing point");
  const auto got =
    EvaluateToHost(*source->CreatePointEvaluator<Real>(bound), field, n_points);
  LO n_bad = got(n_points - 1, 0) == fill_value ? 0 : 1;
  for (LO i = 0; i + 1 < n_points; ++i) {
    const size_t base = static_cast<size_t>(i) * Dim;
    const Real truth =
      linear_f(gathered.coords[base], gathered.coords[base + 1]);
    n_bad += std::fabs(got(i, 0) - truth) > 1e-10 ? 1 : 0;
  }
  require(comm->allreduce(n_bad, Omega_h_Op::OMEGA_H_SUM) == 0,
          "FILL did not fill exactly the point outside the mesh");

  // Only rank 0 has the outside point; every rank must still throw.
  const auto one_rank =
    MakeQueryCoords(comm->rank() == 0 ? points : gathered.coords, Dim);
  pcms::OutOfBoundsPolicy error;
  error.mode = pcms::OutOfBoundsMode::ERROR;
  bool threw = false;
  try {
    (void)g.distributed->Bind(one_rank.view, error);
  } catch (const pcms::pcms_error&) {
    threw = true;
  }
  require(threw, "ERROR did not throw on every rank");
  if (comm->rank() == 0) {
    std::printf("%dD FILL and ERROR policies: ok\n", Dim);
  }
}

// The evaluator reads ghost values, so it is correct only after Synchronize.
template <int Dim>
void check_ghosts_must_be_synchronized(Omega_h::Mesh& mesh, const Geometries& g)
{
  const auto comm = mesh.comm();
  auto space = MakeSpace(g.local, 1);
  const auto source = pcms::CreateDistributedEvaluator(g.distributed, space);
  const auto& layout = *space->GetLayout();
  require(comm->allreduce(mesh.nverts() - layout.GetNumDofHolder(),
                          Omega_h_Op::OMEGA_H_SUM) > 0,
          "mesh has no ghost DOFs");

  // The owned values only: the ghosts keep zero until Synchronize.
  const LO n_owned = layout.GetNumDofHolder();
  const auto dof_coords = layout.GetDOFHolderCoordinates().GetValues();
  Kokkos::View<Real**, pcms::DeviceMemorySpace> values_d("values", n_owned, 1);
  Kokkos::parallel_for(
    "owned_values",
    Kokkos::RangePolicy<pcms::DeviceMemorySpace::execution_space>(0, n_owned),
    KOKKOS_LAMBDA(LO i) {
      values_d(i, 0) = linear_f(dof_coords(i, 0), dof_coords(i, 1));
    });
  auto field = space->CreateFunction<Real>();
  field.GetData().GetDOFAccess()->SetValues(pcms::MakeConstRank2View(values_d));

  const auto gathered = distributed_test::GatherOwnedCentroids(mesh);
  auto coords = MakeQueryCoords(gathered.coords, Dim);
  const LO n_points = static_cast<LO>(gathered.total);
  auto evaluator =
    source->CreatePointEvaluator<Real>(g.distributed->Bind(coords.view));

  const auto count_wrong = [&]() {
    const auto got = EvaluateToHost(*evaluator, field, n_points);
    LO n_wrong = 0;
    for (LO i = 0; i < n_points; ++i) {
      const size_t base = static_cast<size_t>(i) * Dim;
      const Real truth =
        linear_f(gathered.coords[base], gathered.coords[base + 1]);
      n_wrong += std::fabs(got(i, 0) - truth) > 1e-10 ? 1 : 0;
    }
    return comm->allreduce(n_wrong, Omega_h_Op::OMEGA_H_SUM);
  };
  const LO before = count_wrong();
  field.GetData().Synchronize();
  const LO after = count_wrong();
  if (comm->rank() == 0) {
    std::printf("%dD ghost sync: %lld wrong with stale ghosts, %lld after\n",
                Dim, static_cast<long long>(before),
                static_cast<long long>(after));
  }
  require(after == 0, "distributed evaluation is wrong after Synchronize");
}

// One bound points serves the spaces of the geometry, and gives what a
// one-shot Evaluate gives.
template <int Dim>
void check_shared_bound_source(Omega_h::Mesh& mesh, const Geometries& g)
{
  const auto comm = mesh.comm();
  auto linear_space = MakeSpace(g.local, 1);
  auto constant_space = MakeSpace(g.local, 0);
  const auto linear_source =
    pcms::CreateDistributedEvaluator(g.distributed, linear_space);
  const auto constant_source =
    pcms::CreateDistributedEvaluator(g.distributed, constant_space);
  auto linear = linear_space->CreateFunction<Real>();
  auto constant = constant_space->CreateFunction<Real>();
  SetLinearField(linear.GetData(), *linear_space->GetLayout());
  SetLinearField(constant.GetData(), *constant_space->GetLayout());

  const auto gathered = distributed_test::GatherOwnedCentroids(mesh);
  auto coords = MakeQueryCoords(gathered.coords, Dim);
  const LO n_points = static_cast<LO>(gathered.total);

  const auto bound = g.distributed->Bind(coords.view);
  const auto a = EvaluateToHost(
    *linear_source->CreatePointEvaluator<Real>(bound), linear, n_points);
  const auto c = EvaluateToHost(
    *constant_source->CreatePointEvaluator<Real>(bound), constant, n_points);

  // a bind of its own per evaluation
  const auto one_shot = [&](const pcms::PointEvaluatorFactory& source,
                            const pcms::Function<Real>& f) {
    return EvaluateToHost(
      *source.CreatePointEvaluator<Real>(g.distributed->Bind(coords.view)), f,
      n_points);
  };
  const auto b = one_shot(*linear_source, linear);
  const auto d = one_shot(*constant_source, constant);

  LO n_bad = 0;
  for (LO i = 0; i < n_points; ++i) {
    n_bad += a(i, 0) != b(i, 0) ? 1 : 0;
    n_bad += c(i, 0) != d(i, 0) ? 1 : 0;
  }
  const LO t_bad = comm->allreduce(n_bad, Omega_h_Op::OMEGA_H_SUM);
  if (comm->rank() == 0) {
    std::printf("%dD shared bound points: %lld disagreements\n", Dim,
                static_cast<long long>(t_bad));
  }
  require(t_bad == 0, "a shared bound points disagrees with a one-shot one");
}

// An Interpolator from a distributed evaluator works unchanged.
template <int Dim>
void check_interpolator(Omega_h::Mesh& mesh, const Geometries& g)
{
  const auto comm = mesh.comm();
  auto source_space = MakeSpace(g.local, 1);
  auto target_space = MakeSpace(g.local, 0);
  auto source = source_space->CreateFunction<Real>();
  auto target = target_space->CreateFunction<Real>();
  SetLinearField(source.GetData(), *source_space->GetLayout());

  const auto interpolator = pcms::CreateInterpolator<Real>(
    pcms::CreateDistributedEvaluator(g.distributed, source_space),
    target_space);
  interpolator->Apply(source, target);

  // The owned values and their coordinates, on the host.
  const auto owned = target.GetData().GetDOFAccess()->GetValues();
  const auto src =
    target_space->GetLayout()->GetDOFHolderCoordinates().GetValues();
  const LO n = static_cast<LO>(owned.extent(0));
  Kokkos::View<Real**, pcms::DeviceMemorySpace> values_d("values", n, 1);
  Kokkos::View<Real**, pcms::DeviceMemorySpace> coords_d("coords", n, Dim);
  Kokkos::parallel_for(
    "copy_target",
    Kokkos::RangePolicy<pcms::DeviceMemorySpace::execution_space>(0, n),
    KOKKOS_LAMBDA(LO i) {
      values_d(i, 0) = owned(i, 0);
      for (int d = 0; d < Dim; ++d) {
        coords_d(i, d) = src(i, d);
      }
    });
  const auto values =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), values_d);
  const auto coords =
    Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(), coords_d);
  LO n_bad = 0;
  for (LO i = 0; i < n; ++i) {
    n_bad +=
      std::fabs(values(i, 0) - linear_f(coords(i, 0), coords(i, 1))) > 1e-10
        ? 1
        : 0;
  }
  const LO t_bad = comm->allreduce(n_bad, Omega_h_Op::OMEGA_H_SUM);
  if (comm->rank() == 0) {
    std::printf("%dD interpolator: %lld wrong\n", Dim,
                static_cast<long long>(t_bad));
  }
  require(t_bad == 0, "interpolation from a distributed space is wrong");
}

// The search records of `bound`, an Omega_h local bound points, on the host.
template <int Dim>
auto RecordsHost(const pcms::BoundPoints& bound)
{
  const auto& binding =
    dynamic_cast<const pcms::detail::OmegaHMeshBinding&>(bound);
  return Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace(),
                                             binding.Results<Dim>());
}

// How many of the records of `hinted` differ from those of an unhinted bind of
// the same points on the local geometry.
template <int Dim>
LO CountDifferentRecords(
  const pcms::BoundPoints& hinted,
  Kokkos::View<const Real**, pcms::DeviceMemorySpace> points,
  const pcms::Geometry& local)
{
  const auto unhinted =
    local.Bind(pcms::CoordinateView<pcms::DeviceMemorySpace>(
      local.GetCoordinateSystem(), pcms::MakeConstRank2View(points)));
  const auto a = RecordsHost<Dim>(hinted);
  const auto b = RecordsHost<Dim>(*unhinted);
  require(a.extent(0) == b.extent(0), "the binds hold different point counts");
  LO n_different = 0;
  for (size_t i = 0; i < a.extent(0); ++i) {
    bool same =
      a(i).dimensionality == b(i).dimensionality && a(i).entity == b(i).entity;
    for (int d = 0; d < Dim; ++d) {
      same = same && a(i).parametric[d] == b(i).parametric[d];
    }
    n_different += same ? 0 : 1;
  }
  return n_different;
}

// The local binds of a distributed bind start from the elements the routing
// found: one hint per own and per incoming point, naming the local geometry,
// and the records they give are those an unhinted bind gives.
template <int Dim>
void check_local_binds_use_hints(Omega_h::Mesh& mesh, const Geometries& g)
{
  const auto comm = mesh.comm();
  const auto gathered = distributed_test::GatherOwnedCentroids(mesh);
  auto coords = MakeQueryCoords(gathered.coords, Dim);
  const auto bound = std::dynamic_pointer_cast<const pcms::DistributedBinding>(
    g.distributed->Bind(coords.view));
  require(bound != nullptr, "not a distributed bound points");
  const auto& plan = bound->GetRoutingPlan();
  require(plan.own_hints.geometry == g.local &&
            plan.incoming_hints.geometry == g.local,
          "the hints do not name the local geometry");
  require(plan.own_hints.entities.extent(0) == plan.own_points.extent(0) &&
            plan.incoming_hints.entities.extent(0) == plan.IncomingPointCount(),
          "not one hint per own and per incoming point");

  const auto own_coords =
    pcms::detail::GatherRows(coords.storage, plan.own_points);
  const LO n_different =
    CountDifferentRecords<Dim>(*bound->OwnBound(), own_coords, *g.local) +
    CountDifferentRecords<Dim>(*bound->IncomingBound(), plan.incoming_coords,
                               *g.local);
  const LO n_records =
    static_cast<LO>(plan.own_points.extent(0) + plan.IncomingPointCount());
  const LO t_records = comm->allreduce(n_records, Omega_h_Op::OMEGA_H_SUM);
  const LO t_different = comm->allreduce(n_different, Omega_h_Op::OMEGA_H_SUM);
  if (comm->rank() == 0) {
    std::printf("%dD hinted local binds: %lld records, %lld differ from an "
                "unhinted bind\n",
                Dim, static_cast<long long>(t_records),
                static_cast<long long>(t_different));
  }
  require(t_different == 0, "a hinted local bind differs from an unhinted one");
}

template <int Dim>
void run_all(Omega_h::CommPtr world)
{
  auto mesh = distributed_test::BuildGhostedBox<Dim>(world);
  const auto g = MakeGeometries(mesh);
  check_agrees_with_serial<Dim>(mesh, g, /*order=*/1);
  if constexpr (Dim == 2) {
    check_agrees_with_serial<Dim>(mesh, g, /*order=*/0);
    check_boundary_arbitration<Dim>(mesh, g);
  }
  check_out_of_bounds_policy<Dim>(mesh, g);
  check_ghosts_must_be_synchronized<Dim>(mesh, g);
  check_shared_bound_source<Dim>(mesh, g);
  check_interpolator<Dim>(mesh, g);
  check_local_binds_use_hints<Dim>(mesh, g);
}

} // namespace

int main(int argc, char** argv)
{
  try {
    MPI_Init(&argc, &argv);
    {
      Omega_h::Library library(&argc, &argv);
      auto world = library.world();
      require(world->size() >= 2, "test requires at least 2 MPI ranks");
      run_all<2>(world);
      run_all<3>(world);
      if (world->rank() == 0) {
        std::printf("distributed evaluation test passed\n");
      }
    }
    MPI_Finalize();
    return 0;
  } catch (const std::exception& e) {
    distributed_test::AbortAllRanks(e.what());
  } catch (...) {
    distributed_test::AbortAllRanks("unknown exception");
  }
}
