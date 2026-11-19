#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "pcms/field/function_spaces.hpp"
#include <pcms/backends/point_cloud/geometry.hpp>
#include "pcms/field/value_view.hpp"
#include "pcms/backends/omega_h/entity_layout.hpp"
#include "pcms/backends/point_cloud/layout.hpp"
#include "field_test_utils.hpp"

#include <Omega_h_build.hpp>
#include <Omega_h_library.hpp>
#include <vector>
#include "pcms/utility/coordinate_systems/cartesian.hpp"

namespace pcms
{
using detail::PointCloudLayout;
using omega_h::OmegaHEntityLayout;
} // namespace pcms

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

  auto factory = pcms::CreatePolynomialReconstructionSpace(
    pcms::point_cloud::CreateGeometry(pcms::CoordinateView<HostMemorySpace>(
      pcms::csys::Cartesian::Deferred(), coords_view)));
  auto layout = factory->GetLayout();

  REQUIRE(layout->GetNumComponents() == 1);
  REQUIRE(layout->GetNumDofHolder() == 4);
  REQUIRE(layout->GetNumGlobalDofHolder() == 4);
  REQUIRE_FALSE(layout->IsDistributed());

  auto dof_coords = layout->GetDOFHolderCoordinates().GetValues();
  auto dof_coords_host = pcms::test::CopyCoordinatesToHost(dof_coords);

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

  auto factory = pcms::CreatePolynomialReconstructionSpace(
    pcms::point_cloud::CreateGeometry(pcms::CoordinateView<HostMemorySpace>(
      pcms::csys::Cartesian::Deferred(), coords_view)));
  auto source = factory->CreateFunction<Real>();
  auto target = factory->CreateFunction<Real>();

  REQUIRE(&source.GetLayout() == &target.GetLayout());
}

TEST_CASE("PolynomialReconstructionFunctionSpace point-cloud field set/get DOF "
          "round-trip")
{
  auto coords = MakeCoords2D();
  Rank2View<Real, HostMemorySpace> coords_view(coords.data(), 4, 2);

  auto field =
    pcms::CreatePolynomialReconstructionSpace(
      pcms::point_cloud::CreateGeometry(pcms::CoordinateView<HostMemorySpace>(
        pcms::csys::Cartesian::Deferred(), coords_view)))
      ->CreateFunction<Real>();

  std::vector<Real> data{1.0, 2.0, 3.0, 4.0};
  Rank2View<const Real, HostMemorySpace> data_view(
    data.data(), static_cast<LO>(data.size()), 1);
  field.SetDOFHolderDataHost(
    pcms::ValueView<const pcms::Real, pcms::HostMemorySpace>(pcms::ValueBasis{},
                                                             data_view));

  auto got =
    pcms::detail::FlattenToRank1View(field.GetDOFHolderDataHost().GetValues());
  REQUIRE(got.size() == data.size());
  for (LO i = 0; i < static_cast<LO>(data.size()); ++i) {
    REQUIRE(got[i] == Catch::Approx(data[i]));
  }
}

TEST_CASE("PolynomialReconstructionFunctionSpace field keeps layout alive "
          "after temporary factory destruction")
{
  auto coords = MakeCoords2D();
  Rank2View<Real, HostMemorySpace> coords_view(coords.data(), 4, 2);

  auto field = [&]() {
    auto factory = pcms::CreatePolynomialReconstructionSpace(
      pcms::point_cloud::CreateGeometry(pcms::CoordinateView<HostMemorySpace>(
        pcms::csys::Cartesian::Deferred(), coords_view)));
    return factory->CreateFunction<Real>();
  }();

  auto point_cloud_layout =
    dynamic_cast<const pcms::PointCloudLayout*>(&field.GetLayout());
  REQUIRE(point_cloud_layout != nullptr);
  REQUIRE(point_cloud_layout->GetNumDofHolder() == 4);

  std::vector<Real> data{9.0, 10.0, 11.0, 12.0};
  Rank2View<const Real, HostMemorySpace> data_view(
    data.data(), static_cast<LO>(data.size()), 1);
  field.SetDOFHolderDataUncheckedHost(data_view);

  auto got =
    pcms::detail::FlattenToRank1View(field.GetDOFHolderDataHost().GetValues());
  REQUIRE(got[0] == Catch::Approx(9.0));
  REQUIRE(got[3] == Catch::Approx(12.0));
}

TEST_CASE("Layouts on one geometry share it")
{
  Omega_h::Library lib;
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 2,
                                 2, 0, false);

  auto geometry = pcms::test::MakeOmegaHGeometry(mesh);
  auto nodal =
    std::make_shared<pcms::OmegaHEntityLayout>(geometry, pcms::Face, 1);
  auto lagrange = pcms::CreateLagrangeSpace(geometry, 1, 1);

  REQUIRE(nodal->GetGeometry() == geometry);
  REQUIRE(lagrange->GetLayout()->GetGeometry() == geometry);
  REQUIRE(lagrange->GetGeometry() == geometry);

  REQUIRE(nodal->GetNumDofHolder() == mesh.nfaces());
  REQUIRE(lagrange->GetLayout()->GetNumDofHolder() == mesh.nverts());
}

TEST_CASE("Layouts on different geometries do not share one")
{
  Omega_h::Library lib;
  auto mesh_a = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0,
                                   2, 2, 0, false);
  auto mesh_b = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0,
                                   3, 3, 0, false);

  const auto geometry_a = pcms::test::MakeOmegaHGeometry(mesh_a);
  const auto a =
    pcms::OmegaHEntityLayout(geometry_a, pcms::Vertex, 1).GetGeometry();
  const auto b = pcms::OmegaHEntityLayout(
                   pcms::test::MakeOmegaHGeometry(mesh_b), pcms::Vertex, 1)
                   .GetGeometry();
  // another geometry of the same mesh has the same entities but is another
  // geometry: owner hints name geometries by identity
  const auto same_mesh =
    pcms::OmegaHEntityLayout(pcms::test::MakeOmegaHGeometry(mesh_a),
                             pcms::Vertex, 1)
      .GetGeometry();

  REQUIRE(a != b);
  REQUIRE(a != same_mesh);
  REQUIRE(a == geometry_a);
}

TEST_CASE("A standalone point-cloud layout has no geometry")
{
  Omega_h::Library lib;
  auto mesh = Omega_h::build_box(lib.world(), OMEGA_H_SIMPLEX, 1.0, 1.0, 0.0, 2,
                                 2, 0, false);

  auto geometry = pcms::test::MakeOmegaHGeometry(mesh);
  auto lagrange = pcms::CreateLagrangeSpace(geometry, 1, 1);

  auto coords = MakeCoords2D();
  Rank2View<Real, HostMemorySpace> coords_view(coords.data(), 4, 2);
  auto standalone = pcms::CreatePolynomialReconstructionSpace(
    pcms::point_cloud::CreateGeometry(pcms::CoordinateView<HostMemorySpace>(
      pcms::csys::Cartesian::Deferred(), coords_view)));

  REQUIRE(lagrange->GetLayout()->GetGeometry() != nullptr);
  REQUIRE(standalone->GetLayout()->GetGeometry() == nullptr);
}
