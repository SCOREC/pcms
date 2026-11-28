#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <Kokkos_Core.hpp>
#include "pcms/field/point_rows.hpp"
#include "field_test_utils.hpp"
#include <pcms/utility/coordinate_systems/cartesian.hpp>
#include <vector>

namespace
{

using pcms::LO;
using pcms::Real;

template <typename T>
Kokkos::View<T*, pcms::DeviceMemorySpace> Device(const std::vector<T>& v)
{
  Kokkos::View<T*, pcms::HostMemorySpace> host("host", v.size());
  for (std::size_t i = 0; i < v.size(); ++i) {
    host(i) = v[i];
  }
  return Kokkos::create_mirror_view_and_copy(pcms::DeviceMemorySpace{}, host);
}

auto Host(const Kokkos::View<Real**, pcms::DeviceMemorySpace>& view)
{
  return Kokkos::create_mirror_view_and_copy(pcms::HostMemorySpace{}, view);
}

} // namespace

TEST_CASE("PointRows: apply, transpose and matrix of several components",
          "[point_rows]")
{
  // rows of 2, 2 and 1 entries over 4 holders
  const pcms::detail::PointRows rows{Device<LO>({0, 2, 4, 5}),
                                     Device<LO>({0, 1, 1, 2, 3}),
                                     Device<Real>({0.5, 0.5, 0.25, 0.75, 2})};
  Kokkos::View<Real**, pcms::HostMemorySpace> host_dofs("dofs", 4, 2);
  Kokkos::View<Real**, pcms::HostMemorySpace> host_in("in", 3, 2);
  for (int c = 0; c < 2; ++c) {
    for (int h = 0; h < 4; ++h) {
      host_dofs(h, c) = h + 10 * c;
    }
    for (int p = 0; p < 3; ++p) {
      host_in(p, c) = 1 + p - 2 * c;
    }
  }
  const auto dofs = pcms::test::ToDevice(host_dofs);
  const auto in = pcms::test::ToDevice(host_in);

  Kokkos::View<Real**, pcms::DeviceMemorySpace> out("out", 3, 2);
  pcms::detail::Apply(rows, dofs, out);
  const auto host_out = Host(out);
  for (int c = 0; c < 2; ++c) {
    CHECK(host_out(0, c) == Catch::Approx(0.5 + 10 * c));
    CHECK(host_out(1, c) == Catch::Approx(1.75 + 10 * c));
    CHECK(host_out(2, c) == Catch::Approx(6 + 20 * c));
  }

  Kokkos::View<Real**, pcms::DeviceMemorySpace> total("total", 4, 2);
  pcms::detail::ApplyTranspose(rows, in, total);
  const auto host_total = Host(total);
  Real lhs = 0;
  Real rhs = 0;
  for (int c = 0; c < 2; ++c) {
    for (int p = 0; p < 3; ++p) {
      lhs += host_out(p, c) * host_in(p, c);
    }
    for (int h = 0; h < 4; ++h) {
      rhs += host_dofs(h, c) * host_total(h, c);
    }
    CHECK(host_total(1, c) ==
          Catch::Approx(0.5 * host_in(0, c) + 0.25 * host_in(1, c)));
  }
  CHECK(lhs == Catch::Approx(rhs));
}

TEST_CASE("PointRows: the outputs of a point are component-major",
          "[point_rows]")
{
  // one point, two outputs: holder 0 and holder 1
  const pcms::detail::PointRows rows{Device<LO>({0, 1, 2}), Device<LO>({0, 1}),
                                     Device<Real>({1, 1}), 2};
  Kokkos::View<Real**, pcms::HostMemorySpace> host_dofs("dofs", 2, 2);
  host_dofs(0, 0) = 1;
  host_dofs(1, 0) = 2;
  host_dofs(0, 1) = 3;
  host_dofs(1, 1) = 4;
  Kokkos::View<Real**, pcms::DeviceMemorySpace> out("out", 1, 4);
  pcms::detail::Apply(rows, pcms::test::ToDevice(host_dofs), out);
  const auto host_out = Host(out);
  CHECK(host_out(0, 0) == 1);
  CHECK(host_out(0, 1) == 2);
  CHECK(host_out(0, 2) == 3);
  CHECK(host_out(0, 3) == 4);
}
