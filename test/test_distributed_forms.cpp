#include <mpi.h>

#include <Kokkos_Core.hpp>
#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>

#include "distributed_test_utils.hpp"
#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/field/detail/field_ops.hpp"
#include "pcms/field/function_spaces.hpp"
#include "pcms/transfer/forms.hpp"
#include "pcms/transfer/linear_solver.hpp"
#include "pcms/utility/arrays.hpp"

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using pcms::GO;
using pcms::LO;
using pcms::Real;
using pcms::transfer::detail::BilinearForm;
using pcms::transfer::detail::FormOptions;
using Space = std::shared_ptr<const pcms::FunctionSpace>;

namespace
{

using Host2D = Kokkos::View<Real**, pcms::HostMemorySpace>;
using Device2D = Kokkos::View<Real**, pcms::DeviceMemorySpace>;

int failures = 0;

int Rank()
{
  int rank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

Real SumOverRanks(Real local)
{
  Real global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global;
}

// Checks a value that is the same on every rank; rank 0 reports it.
void CheckClose(Real value, Real expected, Real tolerance,
                const std::string& what)
{
  if (!(std::abs(value - expected) <= tolerance)) {
    ++failures;
    if (Rank() == 0) {
      std::fprintf(stderr, "check failed: %s: %.17g != %.17g\n", what.c_str(),
                   value, expected);
    }
  }
}

// Every rank's DOFs depend only on the holders' global ids.
Real DofOfGid(GO gid, int c, Real phase)
{
  return std::sin(1.3 * static_cast<Real>(gid) + 0.7 * c + phase);
}

pcms::Function<Real> FieldOfGids(const Space& space, Real phase)
{
  const auto& layout = *space->GetLayout();
  const auto gids = layout.GetGidsHost();
  Host2D x("x", layout.GetNumDofHolder(), layout.GetNumComponents());
  for (size_t h = 0; h < x.extent(0); ++h) {
    for (size_t c = 0; c < x.extent(1); ++c) {
      x(h, c) = phase < 0 ? 1.0 : DofOfGid(gids(h), c, phase);
    }
  }
  auto field = space->CreateFunction<Real>();
  field.SetDOFHolderDataHost(pcms::ValueView<const Real, pcms::HostMemorySpace>(
    field.GetData().GetValueBasis(), pcms::MakeConstRank2View(x)));
  field.Synchronize();
  return field;
}

pcms::Function<Real> Ones(const Space& space)
{
  return FieldOfGids(space, -1);
}

// sum over the ranks of the owned entries of component c (or of all, c < 0)
Real Total(const pcms::Field<Real>& field, int c = -1)
{
  const auto values = field.GetDOFHolderDataHost().GetValues();
  Real sum = 0;
  for (size_t h = 0; h < values.extent(0); ++h) {
    for (size_t k = 0; k < values.extent(1); ++k) {
      if (c < 0 || static_cast<int>(k) == c) {
        sum += values(h, k);
      }
    }
  }
  return SumOverRanks(sum);
}

// sum over the ranks of the owned x_i y_i
Real GlobalDot(const pcms::Field<Real>& x, const pcms::Field<Real>& y)
{
  const auto a = x.GetDOFHolderDataHost().GetValues();
  const auto b = y.GetDOFHolderDataHost().GetValues();
  Real sum = 0;
  for (size_t h = 0; h < a.extent(0); ++h) {
    for (size_t k = 0; k < a.extent(1); ++k) {
      sum += a(h, k) * b(h, k);
    }
  }
  return SumOverRanks(sum);
}

Real MaxDifference(const pcms::Field<Real>& x, const pcms::Field<Real>& y)
{
  const auto a = x.GetDOFHolderDataHost().GetValues();
  const auto b = y.GetDOFHolderDataHost().GetValues();
  Real local = 0;
  for (size_t h = 0; h < a.extent(0); ++h) {
    for (size_t k = 0; k < a.extent(1); ++k) {
      local = std::max(local, std::abs(a(h, k) - b(h, k)));
    }
  }
  Real global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return global;
}

// DotProd, Apply and Norm2 of the coefficients reduce over the ranks.
void CheckVectorOps(const pcms::Function<Real>& v,
                    const pcms::Function<Real>& w, const std::string& name)
{
  const Real vw = GlobalDot(v, w);
  const Real tolerance = 1e-12 * std::max(1.0, std::abs(vw));
  CheckClose(pcms::detail::DotProd(v, w), vw, tolerance,
             name + ": DotProd over the ranks");
  CheckClose(pcms::detail::Apply(w, v), vw, tolerance,
             name + ": Apply over the ranks");
  const Real vv = GlobalDot(v, v);
  CheckClose(pcms::detail::Norm2(v), std::sqrt(vv), 1e-12 * std::max(1.0, vv),
             name + ": Norm2 over the ranks");
}

Space Lagrange(const std::shared_ptr<const pcms::Geometry>& geometry, int order,
               int nc)
{
  return pcms::CreateLagrangeSpace(
    geometry, order, nc,
    pcms::LagrangeOptions{.backend = pcms::LagrangeBackend::OmegaH});
}

std::string Name(bool ghosted, int dim, const std::string& what)
{
  return std::string(ghosted ? "ghosted" : "element-based") + " dim " +
         std::to_string(dim) + ": " + what;
}

// Assemble of f = 1 over the space's quadrature: the owned entries total the
// measure, per component.
void CheckAssemble(const Space& space, const std::string& name)
{
  const BilinearForm form(space);
  const int nc = space->GetLayout()->GetNumComponents();
  Device2D ones("ones", form.GetQuadrature().weights.extent(0), nc);
  Kokkos::deep_copy(ones, 1.0);
  auto functional = space->CreateFunction<Real>();
  const pcms::ValueView<const Real, pcms::DeviceMemorySpace> samples(
    functional.GetData().GetValueBasis(), pcms::MakeConstRank2View(ones));
  form.Assemble(samples, functional);
  for (int c = 0; c < nc; ++c) {
    CheckClose(Total(functional, c), 1.0, 1e-12,
               name + ": Assemble of 1 totals the measure");
  }
  // overwrites: a second call gives the same result
  form.Assemble(samples, functional);
  CheckClose(Total(functional), nc, 1e-12, name + ": Assemble overwrites");
}

// B u and B^T v of random u and v: the global adjoint identity, Pair summed
// over the ranks, B 1 and B^T 1 totalling the measure, and the overwrite.
void CheckBilinearForm(const BilinearForm& form, const Space& test,
                       const Space& trial, const std::string& name)
{
  const int nc = test->GetLayout()->GetNumComponents();
  const auto u = FieldOfGids(trial, 0.3);
  const auto v = FieldOfGids(test, 1.9);
  auto bu = test->CreateFunction<Real>();
  auto btv = trial->CreateFunction<Real>();
  form.Apply(u, bu);
  form.ApplyTranspose(v, btv);
  const Real vbu = GlobalDot(v, bu);
  const Real ubtv = GlobalDot(u, btv);
  CheckClose(ubtv, vbu, 1e-12 * std::max(1.0, std::abs(vbu)),
             name + ": <v, B u> == <B^T v, u> over the ranks");
  CheckClose(form.Pair(v, u), vbu, 1e-12 * std::max(1.0, std::abs(vbu)),
             name + ": Pair is <v, B u> over the ranks");
  CheckVectorOps(v, bu, name);

  form.Apply(Ones(trial), bu);
  CheckClose(Total(bu), nc, 1e-12, name + ": B 1 totals the measure");
  form.Apply(Ones(trial), bu);
  CheckClose(Total(bu), nc, 1e-12, name + ": Apply overwrites");
  form.ApplyTranspose(Ones(test), btv);
  CheckClose(Total(btv), nc, 1e-12, name + ": B^T 1 totals the measure");
  form.ApplyTranspose(Ones(test), btv);
  CheckClose(Total(btv), nc, 1e-12, name + ": ApplyTranspose overwrites");
}

// The square mass form: the lumped mass totals the measure, and the lumped
// and matrix-free solves recover u.
void CheckMassForm(const Space& space, const std::string& name)
{
  const BilinearForm form(space);
  CheckBilinearForm(form, space, space, name);
  const int nc = space->GetLayout()->GetNumComponents();
  auto lumped = space->CreateFunction<Real>();
  form.Apply(Ones(space), lumped);
  CheckClose(Total(lumped), nc, 1e-12,
             name + ": the lumped mass totals the measure");

  const auto u = FieldOfGids(space, 0.3);
  auto bu = space->CreateFunction<Real>();
  form.Apply(u, bu);
  auto x = space->CreateFunction<Real>();
  const auto solver = pcms::transfer::detail::CreateLinearSolver(
    form, pcms::MassMatrixType::Consistent, {});
  const bool matrix_free =
    dynamic_cast<const pcms::transfer::detail::CgSolver*>(solver.get()) !=
    nullptr;
  CheckClose(matrix_free ? 1 : 0, 1, 0,
             name + ": the default solve of a distributed form is matrix free");
  solver->Solve(bu, x);
  CheckClose(MaxDifference(x, u), 0, 1e-8,
             name + ": the matrix-free solve of B x = B u recovers u");
}

void RunCases(Omega_h::Library& library, int dim, bool ghosted)
{
  auto world = library.world();
  const int n = (dim == 2 ? 4 : 2) * world->size();
  const auto build = [&]() {
    auto mesh =
      Omega_h::build_box(world, OMEGA_H_SIMPLEX, 1.0, 1.0, dim == 3 ? 1.0 : 0.0,
                         n, n, dim == 3 ? n : 0, false);
    if (ghosted) {
      mesh.set_parting(OMEGA_H_GHOSTED, 1, false);
    } else {
      mesh.set_parting(OMEGA_H_ELEM_BASED, 0, false);
    }
    return mesh;
  };
  auto mesh = build();
  auto other_mesh = build();
  const auto geometry =
    pcms::omega_h::CreateGeometry(mesh, pcms::detail::DeferredCartesian());
  const auto other_geometry = pcms::omega_h::CreateGeometry(
    other_mesh, pcms::detail::DeferredCartesian());
  for (int nc : {1, 2}) {
    for (int order : {0, 1}) {
      const auto space = Lagrange(geometry, order, nc);
      const auto what =
        "P" + std::to_string(order) + " components " + std::to_string(nc);
      CheckAssemble(space, Name(ghosted, dim, what));
      CheckMassForm(space, Name(ghosted, dim, what + " mass"));
    }
    // P1 test and P0 trial on one mesh
    const auto p1 = Lagrange(geometry, 1, nc);
    const auto p0 = Lagrange(geometry, 0, nc);
    CheckBilinearForm(
      BilinearForm(p1, p0), p1, p0,
      Name(ghosted, dim, "P1 x P0 components " + std::to_string(nc)));
    // P1 on two meshes with the same partition
    const auto other_p1 = Lagrange(other_geometry, 1, nc);
    CheckBilinearForm(
      BilinearForm(p1, other_p1), p1, other_p1,
      Name(ghosted, dim,
           "P1 x P1 of another mesh components " + std::to_string(nc)));
  }
}

} // namespace

int main(int argc, char** argv)
{
  try {
    MPI_Init(&argc, &argv);
    {
      Omega_h::Library library(&argc, &argv);
      distributed_test::require(library.world()->size() >= 2,
                                "test requires at least 2 MPI ranks");
      for (int dim : {2, 3}) {
        for (bool ghosted : {false, true}) {
          RunCases(library, dim, ghosted);
        }
      }
    }
    int total = 0;
    MPI_Allreduce(&failures, &total, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (Rank() == 0) {
      std::printf("distributed forms test: %d failed checks\n", total);
    }
    MPI_Finalize();
    return total == 0 ? 0 : 1;
  } catch (const std::exception& e) {
    distributed_test::AbortAllRanks(e.what());
  } catch (...) {
    distributed_test::AbortAllRanks("unknown exception");
  }
}
