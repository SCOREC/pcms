#include <catch2/catch_test_macros.hpp>
#include "pcms/field/data/xgc.h"
#include "pcms/field/function_space/xgc.h"
#include <numeric>

namespace
{

pcms::ReverseClassificationVertex create_dummy_rc(int size)
{
  pcms::ReverseClassificationVertex rc;
  for (int i = 0; i < size; ++i) {
    if (i % 4 == 0) {
      rc.Insert({0, 0}, i);
    } else {
      rc.Insert({0, 1}, i);
    }
  }
  return rc;
}

bool in_overlap(int, int id)
{
  return id == 0;
}

} // namespace

TEST_CASE("XGC FieldLayout marks overlap entries and gids")
{
  static constexpr int data_size = 16;
  auto rc = create_dummy_rc(data_size);
  pcms::XGCFieldLayout layout(rc, in_overlap, data_size);

  auto owned = layout.GetOwnedHost();
  auto gids = layout.GetGidsHost();
  auto class_dims = layout.GetDOFHolderClassificationDimensionsHost();
  auto class_ids = layout.GetDOFHolderClassificationIdsHost();

  for (int i = 0; i < data_size; ++i) {
    REQUIRE(gids[i] == i + 1);
    if (i % 4 == 0) {
      REQUIRE(owned[i]);
      REQUIRE(class_dims[i] == 0);
      REQUIRE(class_ids[i] == 0);
    } else {
      REQUIRE(!owned[i]);
      REQUIRE(class_dims[i] == -1);
      REQUIRE(class_ids[i] == -1);
    }
  }
}

TEST_CASE("XGCFieldFactory creates fields and rejects evaluator access")
{
  static constexpr int data_size = 16;
  auto rc = create_dummy_rc(data_size);
  pcms::XGCFieldFactory function_space(rc, in_overlap, data_size);

  std::vector<pcms::Real> data(data_size);
  std::iota(data.begin(), data.end(), 0.0);
  auto field = function_space.CreateField<pcms::Real>(
    "", std::make_unique<pcms::XGCFieldData<pcms::Real>>(
          function_space.GetXGCLayout(), pcms::FieldMetadata{},
          pcms::make_array_view(data)));

  REQUIRE(&field.GetLayout() == function_space.GetLayout().get());
}
