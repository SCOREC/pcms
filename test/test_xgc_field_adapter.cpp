#include <catch2/catch.hpp>
#include <wdmcpl/xgc_field_adapter.h>
#include <algorithm>

using wdmcpl::DimID;
using wdmcpl::make_array_view;
using wdmcpl::make_const_array_view;
using wdmcpl::ReverseClassificationVertex;
using wdmcpl::XGCFieldAdapter;

ReverseClassificationVertex create_dummy_rc(int size)
{
  ReverseClassificationVertex rc;
  for (int i = 0; i < size; ++i) {
    if (i % 4 == 0) {
      rc.Insert({0, 0}, i);
    } else {
      rc.Insert({0, 1}, i);
    }
  }
  return rc;
}
template <typename T, typename T2>
bool is_close(T val1, T2 val2)
{
  static_assert(std::is_integral_v<T2>, "T2 should be counter/integral");
  if constexpr (std::is_integral_v<T>) {
    return val1 == val2;
  }
  return fabs(val1 - val2) < 1E-16;
}
template <typename T, typename Func>
int check_data(const std::vector<T>& data,
               const ReverseClassificationVertex& rc, const Func& in_overlap,
               int offset = 0)
{

  for (const auto& [geom, verts] : rc) {
    auto loffset = in_overlap(geom.dim, geom.id) ? offset : 0;
    for (auto v : verts) {
      if (!is_close(data[v], v + loffset))
        return 1;
    }
  }
  return 0;
}
template <typename T, typename Func>
int check_gids(const std::vector<T>& gids, int size,
               const ReverseClassificationVertex& rc, const Func& in_overlap,
               int offset = 0)
{
  size_t cnt = 0;
  for (const auto& [geom, verts] : rc) {
    if (in_overlap(geom.dim, geom.id)) {
      if (cnt > verts.size()) {
        return 2;
      }
      std::vector<wdmcpl::GO> sorted_verts{verts.begin(), verts.end()};
      std::sort(sorted_verts.begin(), sorted_verts.end());
      for (auto v : sorted_verts) {
        if (!is_close(gids[cnt++], v + offset)) {
          return 1;
        }
      }
    }
  }
  return 0;
}

bool in_overlap(int, int id)
{
  return (id == 0);
}

TEMPLATE_TEST_CASE("XGC Field Adapter", "[adapter]", wdmcpl::LO, wdmcpl::Real)
{
  static constexpr auto data_size = 100;
  // serialize, deserialize, getgids, reverse_partition_map
  std::vector<TestType> dummy_data(data_size);
  std::iota(dummy_data.begin(), dummy_data.end(), 0);

  const auto reverse_classification = create_dummy_rc(data_size);
  const auto num_in_overlap = std::accumulate(
    reverse_classification.begin(), reverse_classification.end(), 0UL,
    [](auto cur, const auto& it) {
      return in_overlap(it.first.dim, it.first.id) ? cur + it.second.size()
                                                   : cur;
    });
  std::cout << "Num in overlap: " << num_in_overlap << "\n";

  XGCFieldAdapter<TestType> field_adapter("fa", make_array_view(dummy_data),
                                          reverse_classification, in_overlap);
  auto gids = field_adapter.GetGids();
  REQUIRE(gids.size() == static_cast<size_t>(num_in_overlap));
  REQUIRE(check_gids(gids, data_size, reverse_classification, in_overlap) == 0);
  std::vector<TestType> buffer;
  std::vector<wdmcpl::LO> permutation;
  auto serialize_size = field_adapter.Serialize(
    make_array_view(buffer), make_const_array_view(permutation));
  REQUIRE(serialize_size == num_in_overlap);
  buffer.resize(serialize_size);
  field_adapter.Serialize(make_array_view(buffer),
                          make_const_array_view(permutation));
  // TODO test with not empty permutation array
  REQUIRE(check_gids(buffer, data_size, reverse_classification, in_overlap) ==
          0);

  // verify that deserializing data writes same values back into dummy data
  field_adapter.Deserialize(make_const_array_view(buffer),
                            make_const_array_view(permutation));
  REQUIRE(check_data(dummy_data, reverse_classification, in_overlap) == 0);
  // modify the buffer and verify that the correct data is written into the
  // dummy_data
  for (auto& val : buffer) {
    val += 5;
  }
  field_adapter.Deserialize(make_const_array_view(buffer),
                            make_const_array_view(permutation));
  REQUIRE(check_data(dummy_data, reverse_classification, in_overlap, 5) == 0);
}
