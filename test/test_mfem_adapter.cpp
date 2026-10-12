#include <catch2/catch_test_macros.hpp>

#include <pcms/coupler/field_serializer.h>
#include <pcms/field/coordinate_systems/cartesian.hpp>
#include <pcms/field/data/mfem.hpp>
#include <pcms/field/function_space/mfem.hpp>
#include <pcms/field/layout/mfem.hpp>

#include <mfem.hpp>

#include <vector>

namespace
{

double LinearField(const mfem::Vector& x)
{
  return 1.0 + 2.0 * x[0] + 3.0 * x[1];
}

} // namespace

TEST_CASE("MFEM vertex-scalar field adapter", "[mfem]")
{
  auto serial = mfem::Mesh::MakeCartesian2D(2, 2, mfem::Element::QUADRILATERAL);
  mfem::ParMesh pmesh(MPI_COMM_WORLD, serial);
  mfem::H1_FECollection fec(1, pmesh.Dimension());
  mfem::ParFiniteElementSpace pfes(&pmesh, &fec);
  mfem::ParGridFunction gf(&pfes);
  gf = 0.0;

  SECTION("layout reports vertex-scalar properties")
  {
    pcms::MFEMLayout layout(pmesh, pfes, pcms::csys::Cartesian::Deferred());

    REQUIRE(layout.GetNumComponents() == 1);
    REQUIRE(layout.GetDimension() == 2);
    REQUIRE(layout.GetNumOwnedDofHolder() == pmesh.GetNV());
    REQUIRE(layout.IsDistributed());

    auto gids = layout.GetGidsHost();
    REQUIRE(static_cast<int>(gids.size()) == pmesh.GetNV());

    auto coords = layout.GetDOFHolderCoordinates().GetValues();
    REQUIRE(static_cast<int>(coords.extent(0)) == pmesh.GetNV());
    REQUIRE(static_cast<int>(coords.extent(1)) == 2);

    // Owned mask covers the global unique vertices across ranks.
    auto owned = layout.GetOwnedHost();
    int local_owned = 0;
    for (size_t i = 0; i < owned.size(); ++i) {
      if (owned[i])
        ++local_owned;
    }
    int total_owned = 0;
    MPI_Allreduce(&local_owned, &total_owned, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    REQUIRE(total_owned == static_cast<int>(layout.GetNumGlobalDofHolder()));
  }

  SECTION("field data round-trips through the grid function")
  {
    mfem::FunctionCoefficient coeff(LinearField);
    gf.ProjectCoefficient(coeff);

    pcms::MFEMLayout layout(pmesh, pfes, pcms::csys::Cartesian::Deferred());
    pcms::MFEMVertexFieldData data(pfes, gf);

    auto host = data.GetDOFHolderDataHost();
    REQUIRE(static_cast<int>(host.extent(0)) == pmesh.GetNV());
    REQUIRE(host.extent(1) == 1);

    std::vector<pcms::Real> captured(host.extent(0));
    for (size_t v = 0; v < captured.size(); ++v) {
      captured[v] = host(v, 0);
    }

    // Write the captured values back and confirm the grid function is restored.
    gf = 0.0;
    Kokkos::View<pcms::Real**, pcms::HostMemorySpace> in("in", captured.size(),
                                                         1);
    for (size_t v = 0; v < captured.size(); ++v) {
      in(v, 0) = captured[v];
    }
    data.SetDOFHolderDataHost(pcms::MakeConstRank2View(in));

    auto host2 = data.GetDOFHolderDataHost();
    auto owned = layout.GetOwnedHost();
    for (size_t v = 0; v < host2.extent(0); ++v) {
      if (owned[v]) {
        REQUIRE(host2(v, 0) == captured[v]);
      }
    }
  }

  SECTION("serializer identity round-trip on a single rank")
  {
    int nproc = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if (nproc != 1) {
      SUCCEED("serializer identity round-trip is only checked on one rank");
      return;
    }

    mfem::FunctionCoefficient coeff(LinearField);
    gf.ProjectCoefficient(coeff);

    pcms::MFEMFieldFactory factory(pmesh, pfes, gf,
                                   pcms::csys::Cartesian::Deferred(), "mfem");
    auto field = factory.CreateField<pcms::Real>("source");
    const auto& layout = field.GetLayout();

    const auto n = static_cast<size_t>(layout.GetNumOwnedDofHolder());

    // Identity permutation: buffer index equals DOF-holder index.
    Kokkos::View<pcms::LO*, pcms::HostMemorySpace> perm("perm", n);
    for (size_t i = 0; i < n; ++i) {
      perm(i) = static_cast<pcms::LO>(i);
    }
    Kokkos::View<pcms::Real*, pcms::HostMemorySpace> buffer("buffer", n);

    pcms::FieldSerializer<pcms::Real> serializer;
    serializer.Serialize(field, pcms::make_array_view(buffer),
                         pcms::make_const_array_view(perm));

    // Deserialize into a second field bound to a fresh grid function.
    mfem::ParGridFunction gf2(&pfes);
    gf2 = 0.0;
    auto field2 = factory.CreateField<pcms::Real>(
      "target", std::make_unique<pcms::MFEMVertexFieldData>(pfes, gf2));
    serializer.Deserialize(field2, pcms::make_const_array_view(buffer),
                           pcms::make_const_array_view(perm));

    auto restored = field2.GetDOFHolderDataHost();
    auto original = field.GetDOFHolderDataHost();
    for (size_t v = 0; v < restored.extent(0); ++v) {
      REQUIRE(restored(v, 0) == original(v, 0));
    }
  }
}
