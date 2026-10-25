#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "pcms/configuration.h"
#include "pcms/coupler/field_serializer.h"
#include "pcms/field/coordinate_systems/cartesian.hpp"
#include "pcms/utility/arrays.h"

#include <memory>
#include <numeric>
#include <vector>

#ifdef PCMS_ENABLE_OMEGA_H
#include "field_test_utils.h"
#include "pcms/field/function_space/lagrange.h"
#include "pcms/field/function_space/polynomial_reconstruction.hpp"
#include "pcms/utility/uniform_grid.h"
#include <Omega_h_library.hpp>
#endif

#ifdef PCMS_ENABLE_XGC
#include "pcms/coupler/overlap_mask.h"
#include "pcms/field/data/xgc.h"
#include "pcms/field/function_space/xgc.h"
#endif

namespace
{

#ifdef PCMS_ENABLE_OMEGA_H
void CheckSerializeDeserialize(pcms::Field<pcms::Real>& field)
{
  auto data_before = field.GetDOFHolderDataHost().GetValues();
  const auto num_dof = static_cast<pcms::LO>(data_before.extent(0));
  const auto num_components = static_cast<pcms::LO>(data_before.extent(1));
  const auto num_values = static_cast<size_t>(data_before.size());

  std::vector<pcms::Real> expected(num_values);
  for (pcms::LO i = 0; i < num_dof; ++i) {
    for (pcms::LO c = 0; c < num_components; ++c) {
      expected[static_cast<size_t>(i) * num_components + c] = data_before(i, c);
    }
  }

  std::vector<pcms::Real> buffer(num_values);
  std::vector<pcms::LO> permutation(num_dof);
  std::iota(permutation.begin(), permutation.end(), pcms::LO{0});

  pcms::FieldSerializer<pcms::Real> serializer;
  REQUIRE(serializer.Serialize(field, pcms::make_array_view(buffer),
                               pcms::make_const_array_view(permutation)) ==
          static_cast<int>(num_values));

  std::vector<pcms::Real> cleared(num_values, -1.0);
  field.SetDOFHolderDataUncheckedHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
      cleared.data(), num_dof, num_components));
  serializer.Deserialize(field, pcms::make_const_array_view(buffer),
                         pcms::make_const_array_view(permutation));

  auto data_after = field.GetDOFHolderDataHost().GetValues();
  REQUIRE(data_after.size() == expected.size());
  for (pcms::LO i = 0; i < num_dof; ++i) {
    for (pcms::LO c = 0; c < num_components; ++c) {
      REQUIRE(
        data_after(i, c) ==
        Catch::Approx(expected[static_cast<size_t>(i) * num_components + c]));
    }
  }
}
#endif

#ifdef PCMS_ENABLE_XGC
pcms::ReverseClassificationVertex CreateDummyReverseClassification(int size)
{
  pcms::ReverseClassificationVertex rc;
  for (int i = 0; i < size; ++i) {
    rc.Insert(i % 4 == 0 ? pcms::DimID{0, 0} : pcms::DimID{0, 1}, i);
  }
  return rc;
}

bool InXGCOverlap(int, int id)
{
  return id == 0;
}
#endif

} // namespace

#ifdef PCMS_ENABLE_XGC
TEST_CASE("Serializing an XGC field preserves inactive field entries")
{
  static constexpr int data_size = 16;
  auto rc = CreateDummyReverseClassification(data_size);
  pcms::XGCFieldFactory factory(rc, data_size);

  std::vector<pcms::Real> data(data_size);
  std::iota(data.begin(), data.end(), 0.0);
  const auto original = data;
  auto field = factory.CreateField<pcms::Real>(
    "", std::make_unique<pcms::XGCFieldData<pcms::Real>>(
          factory.GetXGCLayout(), pcms::ValueBasis{},
          pcms::make_array_view(data), MPI_COMM_SELF));
  pcms::FieldSerializer<pcms::Real> serializer;

  // Entries outside the overlap region are not exchanged.
  const pcms::OverlapMask overlap_mask(InXGCOverlap);
  auto active = overlap_mask.GetMask(field.GetLayout());
  std::vector<pcms::LO> permutation(data_size, -1);
  int entry = 0;
  for (int i = 0; i < data_size; ++i) {
    if (active[i]) {
      permutation[i] = entry++;
    }
  }
  const int num_active = entry;
  REQUIRE(num_active == 4);

  std::vector<pcms::Real> buffer(num_active, -1.0);
  REQUIRE(serializer.Serialize(field, pcms::make_array_view(buffer),
                               pcms::make_const_array_view(permutation)) ==
          num_active);

  for (int i = 0; i < data_size; ++i) {
    if (active[i]) {
      REQUIRE(buffer[permutation[i]] == data[i]);
      buffer[permutation[i]] += 100.0;
    } else {
      REQUIRE(permutation[i] == -1);
    }
  }

  serializer.Deserialize(field, pcms::make_const_array_view(buffer),
                         pcms::make_const_array_view(permutation));

  auto after =
    pcms::FlattenToRank1View(field.GetDOFHolderDataHost().GetValues());
  for (int i = 0; i < data_size; ++i) {
    if (active[i]) {
      REQUIRE(after[i] == Catch::Approx(original[i] + 100.0));
    } else {
      REQUIRE(after[i] == Catch::Approx(original[i]));
    }
  }
}
#endif

#ifdef PCMS_ENABLE_OMEGA_H
TEST_CASE("FieldSerializer round-trips a polynomial-reconstruction field")
{
  std::vector<pcms::Real> coords{0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  pcms::Rank2View<pcms::Real, pcms::HostMemorySpace> coords_view(coords.data(),
                                                                 4, 2);
  auto space = pcms::PolynomialReconstructionFunctionSpace::Create(
    pcms::CoordinateView<pcms::HostMemorySpace>(
      pcms::csys::Cartesian::Deferred(), coords_view));
  auto field = space->CreateFunction<pcms::Real>();

  std::vector<pcms::Real> data{5.0, 6.0, 7.0, 8.0};
  field.SetDOFHolderDataUncheckedHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(data.data(), 4,
                                                             1));

  CheckSerializeDeserialize(field);
}

TEST_CASE("FieldSerializer round-trips a uniform-grid field")
{
  pcms::UniformGrid<2> grid;
  grid.bot_left = {0.0, 0.0};
  grid.edge_length = {10.0, 10.0};
  grid.divisions = {3, 3};
  auto space = pcms::LagrangeFunctionSpace::FromUniformGrid(
    grid, 1, pcms::csys::Cartesian::Deferred());
  auto field = space->CreateFunction<pcms::Real>();

  std::vector<pcms::Real> data(16);
  for (size_t i = 0; i < data.size(); ++i) {
    data[i] = static_cast<pcms::Real>(i * 10);
  }
  field.SetDOFHolderDataUncheckedHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(data.data(), 16,
                                                             1));

  CheckSerializeDeserialize(field);
}

TEST_CASE("FieldSerializer round-trips an order-1 Omega_h field")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto space = pcms::test::MakeP1Space(mesh);
  auto field = space->CreateFunction<pcms::Real>();

  pcms::test::SetField(
    field, OMEGA_H_LAMBDA(pcms::Real x, pcms::Real y) {
      return pcms::test::linear_f(x, y);
    });

  CheckSerializeDeserialize(field);
}

TEST_CASE("FieldSerializer round-trips a multi-component Omega_h field")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);
  constexpr int num_components = 3;
  auto space = pcms::LagrangeFunctionSpace::FromMesh(
    mesh, 1, num_components, pcms::csys::Cartesian::Deferred(), "global",
    pcms::LagrangeFunctionSpace::Backend::OmegaH);
  auto field = space->CreateFunction<pcms::Real>();

  const int num_dof = field.GetLayout().GetNumOwnedDofHolder();
  std::vector<pcms::Real> data(static_cast<size_t>(num_dof) * num_components);
  for (int i = 0; i < num_dof; ++i) {
    for (int c = 0; c < num_components; ++c) {
      data[static_cast<size_t>(i) * num_components + c] = i + 0.25 * c;
    }
  }
  field.SetDOFHolderDataUncheckedHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(
      data.data(), num_dof, num_components));

  CheckSerializeDeserialize(field);
}

TEST_CASE("FieldSerializer round-trips an order-0 Omega_h field")
{
  auto lib = Omega_h::Library{};
  auto mesh = pcms::test::BuildUnitSquare(lib, 0);
  auto space = pcms::test::MakeP0Space(mesh);
  auto field = space->CreateFunction<pcms::Real>();

  const int num_dof = mesh.nelems();
  std::vector<pcms::Real> data(num_dof);
  std::iota(data.begin(), data.end(), pcms::Real{0});
  field.SetDOFHolderDataUncheckedHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(data.data(),
                                                             num_dof, 1));

  CheckSerializeDeserialize(field);
}
TEST_CASE("FieldSerializer preserves holders outside the exchange on receive")
{
  std::vector<pcms::Real> coords{0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  pcms::Rank2View<pcms::Real, pcms::HostMemorySpace> coords_view(coords.data(),
                                                                 4, 2);
  auto space = pcms::PolynomialReconstructionFunctionSpace::Create(
    pcms::CoordinateView<pcms::HostMemorySpace>(
      pcms::csys::Cartesian::Deferred(), coords_view));
  auto field = space->CreateFunction<pcms::Real>();

  std::vector<pcms::Real> data{5.0, 6.0, 7.0, 8.0};
  field.SetDOFHolderDataUncheckedHost(
    pcms::Rank2View<const pcms::Real, pcms::HostMemorySpace>(data.data(), 4,
                                                             1));

  // Holders 1 and 3 have no slot in the received buffer (e.g. they lie outside
  // the sender's overlap mask); they must keep their current values.
  std::vector<pcms::LO> permutation{0, -1, 1, -1};
  std::vector<pcms::Real> buffer{50.0, 70.0};
  pcms::FieldSerializer<pcms::Real> serializer;
  serializer.Deserialize(field, pcms::make_const_array_view(buffer),
                         pcms::make_const_array_view(permutation));

  auto after = field.GetDOFHolderDataHost().GetValues();
  REQUIRE(after(0, 0) == 50.0);
  REQUIRE(after(1, 0) == 6.0);
  REQUIRE(after(2, 0) == 70.0);
  REQUIRE(after(3, 0) == 8.0);
}
#endif
