#ifndef PCMS_SRC_PCMS_DUMMY_FIELD_ADAPTER_H
#define PCMS_SRC_PCMS_DUMMY_FIELD_ADAPTER_H
#include "pcms/field.h"
#include "pcms/partition.h"
#include <vector>
namespace pcms
{
class DummyFieldAdapter
{
public:
  using value_type = int;
  [[nodiscard]] std::vector<GO> GetGids() const { return {}; }
  [[nodiscard]] ReversePartitionMap GetReversePartitionMap(
    const Partition& partition) const
  {
    return {};
  }
  template <typename T1, typename T2>
  int Serialize(T1, T2) const noexcept
  {
    return 0;
  }
  // function so that call to Serialize({},{}) works.
  int Serialize(int, int) const noexcept { return 0; }
  template <typename T1, typename T2>
  void Deserialize(T1, T2) const noexcept
  {
  }
};

} // namespace pcms

#endif // PCMS_SRC_PCMS_DUMMY_FIELD_ADAPTER_H
