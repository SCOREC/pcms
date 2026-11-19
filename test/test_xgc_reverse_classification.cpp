#include "pcms/backends/xgc/reverse_classification.hpp"
#include <catch2/catch_test_macros.hpp>
#include <sstream>

static constexpr auto test_data = R"(
17
0 1
13
0 2
37
1 11
734 735 736 737
1 12
871 872 873 874
1 13
1020 1021 1022 1023
1 22
-1
2 20
-1
2 1
5 6 7
)";

TEST_CASE("reverse classification")
{
  std::stringstream ss{test_data};
  auto rc = pcms::xgc::ReadReverseClassificationVertex(ss);
  {
    const auto* q =
      pcms::xgc::detail::ReverseClassificationAccess::Query(rc, {2, 20});
    REQUIRE(q == nullptr);
  }
  {
    const auto* q =
      pcms::xgc::detail::ReverseClassificationAccess::Query(rc, {2, 1});
    REQUIRE(q != nullptr);
    REQUIRE(q->size() == 3);
    std::set<pcms::LO> s{4, 5, 6};
    REQUIRE(s == *q);
  }
  auto vec = pcms::xgc::detail::ReverseClassificationAccess::Serialize(rc);
  pcms::xgc::ReverseClassificationVertex rc_deserialized;
  auto av = pcms::detail::make_array_view(vec);
  pcms::xgc::detail::ReverseClassificationAccess::Deserialize(rc_deserialized,
                                                              av);
  REQUIRE(
    pcms::xgc::detail::ReverseClassificationAccess::Equal(rc_deserialized, rc));
}
