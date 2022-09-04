#ifndef WDM_COUPLING_FIELD_H
#define WDM_COUPLING_FIELD_H
#include "wdmcpl/external/span.h"
#include "wdmcpl/types.h"
#include "wdmcpl/arrays.h"
#include "wdmcpl/memory_spaces.h"
#include <map>
#include <redev.h>
#include "wdmcpl/field_evaluation_methods.h"
#include <any>
namespace wdmcpl
{

namespace detail
{
template <typename...>
using VoidT = void;

template <typename, typename = void>
struct HasCoordinateSystem : public std::false_type
{};

template <typename T>
struct HasCoordinateSystem<T, VoidT<typename T::coordinate_system>>
  : public std::true_type
{};

} // namespace detail

/**
 * Key: result of partition object i.e. rank that the data is sent to on host
 * Value: Vector of local index (ordered)
 */
using ReversePartitionMap = std::map<wdmcpl::LO, std::vector<wdmcpl::LO>>;

// template <typename T, typename MemorySpace = HostMemorySpace, typename
// CoordinateElementType = Real>
template <typename T, typename MemorySpace = HostMemorySpace>
class FieldAdapter
{
  using InternalCoordinateType = Real;

public:
  using memory_space = MemorySpace;
  using value_type = T;
  virtual const std::string& GetName() const noexcept = 0;
  virtual int Serialize(
    ScalarArrayView<T, MemorySpace> buffer,
    ScalarArrayView<const wdmcpl::LO, MemorySpace> permutation) const = 0;
  virtual void Deserialize(
    ScalarArrayView<T, MemorySpace> buffer,
    ScalarArrayView<const wdmcpl::LO, MemorySpace> permutation) const = 0;
  virtual std::vector<GO> GetGids() const = 0;
  virtual ReversePartitionMap GetReversePartitionMap(
    const redev::Partition& partition) const = 0;
};
} // namespace wdmcpl

#endif // WDM_COUPLING_FIELD_H
