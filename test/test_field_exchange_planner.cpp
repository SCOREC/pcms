#include <catch2/catch_test_macros.hpp>

#include <pcms/coupler/field_exchange_planner.h>
#include <pcms/field/layout/point_cloud.h>
#include "pcms/field/coordinate_systems/cartesian.hpp"

TEST_CASE("Field exchange plans exclude GID message headers",
          "[field_exchange]")
{
  Kokkos::View<pcms::Real**> coords("coords", 4, 2);
  pcms::PointCloudLayout layout(2, coords, pcms::csys::Cartesian::Deferred());

  // Two incoming messages. Each contains a five-entry entity-offset header
  // followed by two vertex GIDs.
  Kokkos::View<pcms::GO*, pcms::HostMemorySpace> received(
    "received_gids", 2 * static_cast<std::size_t>(pcms::ent_offsets_len + 2));
  const pcms::GO message[] = {
    0, 2, 2, 2, 2, 2, 0, 0, 2, 2, 2, 2, 3, 1,
  };
  for (size_t i = 0; i < received.size(); ++i)
    received(i) = message[i];

  redev::InMessageLayout in;
  in.srcRanks = {0, 7};
  in.offset = {0, 14};

  pcms::GenericFieldExchangePlanner planner;
  const auto plan = planner.BuildReceivePlan(
    layout,
    pcms::GlobalIDView<pcms::HostMemorySpace>(received.data(), received.size()),
    0, 1, in);

  REQUIRE(plan.dest_ranks == redev::LOs{0, 1});
  REQUIRE(plan.offsets == redev::LOs{0, 2, 4});
  REQUIRE(plan.permutation == std::vector<pcms::LO>{1, 3, 0, 2});
  REQUIRE(plan.msg_size == 4);
}

TEST_CASE("GID messages insert headers around compact field payloads",
          "[field_exchange]")
{
  Kokkos::View<pcms::Real**> coords("coords", 4, 2);
  pcms::PointCloudLayout layout(2, coords, pcms::csys::Cartesian::Deferred());

  pcms::ExchangePlan plan;
  plan.dest_ranks = {0, 1};
  plan.offsets = {0, 2, 4};
  plan.permutation = {0, 2, 1, 3};
  plan.msg_size = 4;

  Kokkos::View<pcms::GO*, pcms::HostMemorySpace> message(
    "gid_message",
    plan.msg_size + 2 * static_cast<std::size_t>(pcms::ent_offsets_len));
  pcms::GenericFieldExchangePlanner planner;
  planner.FillGidMessage(layout, plan, pcms::MakeRank1View(message));

  const pcms::GO expected[] = {
    0, 2, 2, 2, 2, 0, 2, 0, 2, 2, 2, 2, 1, 3,
  };
  for (size_t i = 0; i < message.size(); ++i) {
    CAPTURE(i);
    REQUIRE(message(i) == expected[i]);
  }
}

// Verify the permutation documented on ExchangePlan::permutation:
//   - permutation[i] < 0 for holders excluded from this exchange.
//   - If all holders are owned, a partial GID message leaves the unmatched
//     (owned but not participating) holders with permutation[i] < 0.
TEST_CASE("ExchangePlan permutation is negative for non-participating holders",
          "[field_exchange][permutation]")
{
  // Create 6 points (all owned), GIDs 0-5. Only GIDs 0-3 will appear in
  // the received message, so holders 4 and 5 should have negative
  // permutation entries.
  Kokkos::View<pcms::Real**> coords("coords", 6, 2);
  pcms::PointCloudLayout layout(2, coords, pcms::csys::Cartesian::Deferred());

  const std::size_t header_len =
    static_cast<std::size_t>(pcms::ent_offsets_len);
  const std::size_t msg_total = header_len + 4;
  Kokkos::View<pcms::GO*, pcms::HostMemorySpace> received("received_gids",
                                                          msg_total);
  const pcms::GO message[] = {0, 4, 4, 4, 4, 0, 1, 2, 3};
  for (size_t i = 0; i < received.size(); ++i)
    received(i) = message[i];

  redev::InMessageLayout in;
  in.srcRanks = {0};
  in.offset = {0, static_cast<redev::LO>(msg_total)};

  pcms::GenericFieldExchangePlanner planner;
  const auto plan = planner.BuildReceivePlan(
    layout,
    pcms::GlobalIDView<pcms::HostMemorySpace>(received.data(), received.size()),
    0, 1, in);

  REQUIRE(plan.permutation.size() == 6);
  // GIDs 0-3 are in the received message → non-negative permutation entries.
  REQUIRE(plan.permutation[0] == 0);
  REQUIRE(plan.permutation[1] == 1);
  REQUIRE(plan.permutation[2] == 2);
  REQUIRE(plan.permutation[3] == 3);
  // GIDs 4-5 are not in the message → negative sentinel.
  REQUIRE(plan.permutation[4] < 0);
  REQUIRE(plan.permutation[5] < 0);
  // Payload length excludes headers.
  REQUIRE(plan.msg_size == 4);
}
