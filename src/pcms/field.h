#ifndef PCMS_COUPLING_FIELD_H
#define PCMS_COUPLING_FIELD_H
#include "pcms/types.h"
#include "pcms/arrays.h"
#include "pcms/memory_spaces.h"
#include <map>
#include <redev.h>
#include "pcms/field_evaluation_methods.h"
#include <any>
namespace pcms
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
 * Key: result of partition object i.e. rank that the data is sent to on coupling server
 * Value: Vector of local index (ordered)
 */
using ReversePartitionMap = std::map<pcms::LO, std::vector<pcms::LO>>;

// This is a model interface. The current pcms design does not require that
// the FieldAdapter must be inherited from this class
/*
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
    ScalarArrayView<const pcms::LO, MemorySpace> permutation) const = 0;
  virtual void Deserialize(
    ScalarArrayView<T, MemorySpace> buffer,
    ScalarArrayView<const pcms::LO, MemorySpace> permutation) const = 0;
  virtual std::vector<GO> GetGids() const = 0;
  virtual ReversePartitionMap GetReversePartitionMap(
    const redev::Partition& partition) const = 0;
};
*/
} // namespace pcms

#endif // PCMS_COUPLING_FIELD_H
