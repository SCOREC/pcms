#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "pcms/field/function_space/polynomial_reconstruction.hpp"
#include "pcms/field/field_metadata.h"
#include "pcms/discretization/discretization.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/layout/point_cloud.h"
#include "field_test_utils.h"

#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <vector>

using pcms::CoordinateSystem;
using pcms::HostMemorySpace;
using pcms::LO;
using pcms::Rank1View;
using pcms::Rank2View;
using pcms::Real;

namespace
{

std::vector<Real> MakeCoords2D()
{
  return {0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0};
}

} // namespace

TEST_CASE(
  "PolynomialReconstructionFunctionSpace creates point-cloud layout metadata")
{
  auto coords = MakeCoords2D();
  Rank2View<Real, HostMemorySpace> coords_view(coords.data(), 4, 2);

  auto factory = pcms::PolynomialReconstructionFunctionSpace::Create(
    coords_view, CoordinateSystem::Cartesian);
  auto layout = factory.GetLayout();

  REQUIRE(layout->GetNumComponents() == 1);
  REQUIRE(layout->GetNumOwnedDofHolder() == 4);
  REQUIRE(layout->GetNumGlobalDofHolder() == 4);
  REQUIRE_FALSE(layout->IsDistributed());

  auto dof_coords = layout->GetDOFHolderCoordinates().GetValues();
  auto dof_coords_host = pcms::test::CopyCoordinatesToHost(dof_coords, 4, 2);

  REQUIRE(static_cast<int>(dof_coords_host.extent(0)) == 4);
  REQUIRE(static_cast<int>(dof_coords_host.extent(1)) == 2);
  for (int i = 0; i < 4; ++i) {
    REQUIRE(dof_coords_host(i, 0) ==
            Catch::Approx(coords[2 * static_cast<size_t>(i)]));
    REQUIRE(dof_coords_host(i, 1) ==
            Catch::Approx(coords[2 * static_cast<size_t>(i) + 1]));
  }
}

TEST_CASE("PolynomialReconstructionFunctionSpace fields share layout")
{
  auto coords = MakeCoords2D();
  Rank2View<Real, HostMemorySpace> coords_view(coords.data(), 4, 2);

  auto factory = pcms::PolynomialReconstructionFunctionSpace::Create(
    coords_view, CoordinateSystem::Cartesian);
  auto source = factory.CreateField<Real>(pcms::FieldMetadata{});
  auto target = factory.CreateField<Real>(pcms::FieldMetadata{});

  REQUIRE(&source.GetLayout() == &target.GetLayout());
}

TEST_CASE("PolynomialReconstructionFunctionSpace point-cloud field set/get DOF "
          "round-trip")
{
  auto coords = MakeCoords2D();
  Rank2View<Real, HostMemorySpace> coords_view(coords.data(), 4, 2);

  auto field = pcms::PolynomialReconstructionFunctionSpace::Create(
                 coords_view, CoordinateSystem::Cartesian)
                 .CreateField<Real>(pcms::FieldMetadata{});

  std::vector<Real> data{1.0, 2.0, 3.0, 4.0};
  Rank2View<const Real, HostMemorySpace> data_view(
    data.data(), static_cast<LO>(data.size()), 1);
  field.GetData().SetDOFHolderDataHost(data_view);

  auto got = pcms::FlattenToRank1View(field.GetData().GetDOFHolderDataHost());
  REQUIRE(got.size() == data.size());
  for (LO i = 0; i < static_cast<LO>(data.size()); ++i) {
    REQUIRE(got[i] == Catch::Approx(data[i]));
  }
}

TEST_CASE("PolynomialReconstructionFunctionSpace point-cloud field serialize / "
          "deserialize round-trip")
{
  auto coords = MakeCoords2D();
  Rank2View<Real, HostMemorySpace> coords_view(coords.data(), 4, 2);

  auto factory = pcms::PolynomialReconstructionFunctionSpace::Create(
    coords_view, CoordinateSystem::Cartesian);
  auto field = factory.CreateField<Real>(pcms::FieldMetadata{});

  std::vector<Real> data{5.0, 6.0, 7.0, 8.0};
  Rank2View<const Real, HostMemorySpace> data_view(
    data.data(), static_cast<LO>(data.size()), 1);
  field.GetData().SetDOFHolderDataHost(data_view);

  pcms::test::CheckSerializeDeserialize(*factory.GetLayout(), field.GetData());
}

TEST_CASE("PolynomialReconstructionFunctionSpace field keeps layout alive "
          "after temporary factory destruction")
{
  auto coords = MakeCoords2D();
  Rank2View<Real, HostMemorySpace> coords_view(coords.data(), 4, 2);

  auto field = [&]() {
    auto factory = pcms::PolynomialReconstructionFunctionSpace::Create(
      coords_view, CoordinateSystem::Cartesian);
    return factory.CreateField<Real>(pcms::FieldMetadata{});
  }();

  auto point_cloud_layout =
    dynamic_cast<const pcms::PointCloudLayout*>(&field.GetLayout());
  REQUIRE(point_cloud_layout != nullptr);
  REQUIRE(point_cloud_layout->GetNumOwnedDofHolder() == 4);

  std::vector<Real> data{9.0, 10.0, 11.0, 12.0};
  Rank2View<const Real, HostMemorySpace> data_view(
    data.data(), static_cast<LO>(data.size()), 1);
  field.SetDOFHolderDataHost(data_view);

  auto got = pcms::FlattenToRank1View(field.GetDOFHolderDataHost());
  REQUIRE(got[0] == Catch::Approx(9.0));
  REQUIRE(got[3] == Catch::Approx(12.0));
}

TEST_CASE("Different layouts on the same mesh report SameEntities")
{
  Omega_h::Library lib;
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 2,
                                 2, 0, false);

  auto nodal = pcms::PolynomialReconstructionFunctionSpace::FromMesh(
    mesh, pcms::Face, CoordinateSystem::Cartesian);
  auto lagrange = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, CoordinateSystem::Cartesian);

  auto nodal_disc = nodal.GetLayout()->GetDiscretization();
  auto lagrange_disc = lagrange.GetLayout()->GetDiscretization();

  REQUIRE(nodal_disc != nullptr);
  REQUIRE(lagrange_disc != nullptr);
  REQUIRE(nodal_disc->SameEntities(*lagrange_disc));

  REQUIRE(nodal.GetLayout()->GetNumOwnedDofHolder() == mesh.nfaces());
  REQUIRE(lagrange.GetLayout()->GetNumOwnedDofHolder() == mesh.nverts());

  REQUIRE(nodal_disc->GetNumEntities(pcms::Face) == mesh.nfaces());
  REQUIRE(lagrange_disc->GetNumEntities(pcms::Vertex) == mesh.nverts());
}

TEST_CASE("Layouts on different meshes do not report SameEntities")
{
  Omega_h::Library lib;
  auto mesh_a = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0,
                                   2, 2, 0, false);
  auto mesh_b = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0,
                                   3, 3, 0, false);

  auto nodal_a = pcms::PolynomialReconstructionFunctionSpace::FromMesh(
    mesh_a, pcms::Vertex, CoordinateSystem::Cartesian);
  auto nodal_b = pcms::PolynomialReconstructionFunctionSpace::FromMesh(
    mesh_b, pcms::Vertex, CoordinateSystem::Cartesian);

  auto disc_a = nodal_a.GetLayout()->GetDiscretization();
  auto disc_b = nodal_b.GetLayout()->GetDiscretization();

  REQUIRE_FALSE(disc_a->SameEntities(*disc_b));
}

TEST_CASE(
  "Standalone point-cloud layout does not report SameEntities with mesh layout")
{
  Omega_h::Library lib;
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 2,
                                 2, 0, false);

  auto lagrange = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, 1, CoordinateSystem::Cartesian);

  auto coords = MakeCoords2D();
  Rank2View<Real, HostMemorySpace> coords_view(coords.data(), 4, 2);
  auto standalone = pcms::PolynomialReconstructionFunctionSpace::Create(
    coords_view, CoordinateSystem::Cartesian);

  auto mesh_disc = lagrange.GetLayout()->GetDiscretization();
  auto point_cloud_disc = standalone.GetLayout()->GetDiscretization();

  REQUIRE_FALSE(mesh_disc->SameEntities(*point_cloud_disc));
  REQUIRE_FALSE(point_cloud_disc->SameEntities(*mesh_disc));
}
