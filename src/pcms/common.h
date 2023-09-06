#ifndef WDM_COUPLING_COMMON_H
#define WDM_COUPLING_COMMON_H
#include <redev.h>
#include "pcms/transfer_field.h"
#include "pcms/assert.h"
#include <map>
namespace pcms
{
using ProcessType = redev::ProcessType;
class GatherOperation;
class ScatterOperation;

namespace detail
{
// helper function for dealing with field maps
template <typename T, typename U>
auto& find_or_error(const std::string& name,
                    const std::map<T, U>& map)
{
  auto it = map.find(name);
  WDMCPL_ALWAYS_ASSERT(it != map.end());
  return it->second;
}

template <typename T, typename U>
auto& find_or_error(const std::string& name, std::map<T, U>& map)
{
  auto it = map.find(name);
  WDMCPL_ALWAYS_ASSERT(it != map.end());
  return it->second;
}
template <typename T, typename U>
auto find_many_or_error(const std::vector<T>& keys,
                        const std::map<T, U>& map)
{

  std::vector<std::reference_wrapper<U>> results;
  results.reserve(keys.size());
  std::transform(keys.begin(), keys.end(), std::back_inserter(results),
                 [&map](const std::string& key) {
                   return std::ref(detail::find_or_error(key, map));
                 });
  return results;
}
template <typename T, typename U>
auto find_many_or_error(const std::vector<T>& keys,
                        std::map<T, U>& map)
{

  std::vector<std::reference_wrapper<U>> results;
  results.reserve(keys.size());
  std::transform(keys.begin(), keys.end(), std::back_inserter(results),
                 [&map](const std::string& key) {
                   return std::ref(detail::find_or_error(key, map));
                 });
  return results;
}
} // namespace detail

struct TransferOptions
{
  FieldTransferMethod transfer_method;
  FieldEvaluationMethod evaluation_method;
};
} // namespace pcms

#endif // WDM_COUPLING_COMMON_H
