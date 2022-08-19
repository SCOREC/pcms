#ifndef WDM_COUPLING_FIELD_H
#define WDM_COUPLING_FIELD_H
#include "wdmcpl/external/span.h"
#include "wdmcpl/types.h"
#include "wdmcpl/arrays.h"
#include "wdmcpl/memory_spaces.h"
#include <map>
#include <redev.h>
#include "wdmcpl/field_evaluation_methods.h"
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
template <typename T, typename CoordinateElementType>
class OmegaHField;

//template <typename T, typename MemorySpace = HostMemorySpace, typename CoordinateElementType = Real>
template <typename T, typename MemorySpace = HostMemorySpace>
class FieldShim
{
public:
  using memory_space = MemorySpace;
  using value_type = T;
  virtual const std::string& GetName() const noexcept= 0;
  virtual int Serialize(
    ScalarArrayView<T, MemorySpace> buffer,
    ScalarArrayView<const wdmcpl::LO, MemorySpace> permutation) const = 0;
  virtual void Deserialize(
    ScalarArrayView<T, MemorySpace> buffer,
    ScalarArrayView<const wdmcpl::LO, MemorySpace> permutation) const = 0;
  virtual std::vector<GO> GetGids() const = 0;
  virtual ReversePartitionMap GetReversePartitionMap(
    const redev::Partition& partition) const = 0;
  // FIXME these functions are only needed on the server and we need to make sure we
  // can compile the applications without linking against OmegaH.
  virtual void ToOmegaH(OmegaHField<T, Real>& internal_field, FieldTransferMethod, FieldEvaluationMethod) = 0;
  virtual void FromOmegaH(const OmegaHField<T, Real>& internal_field, FieldTransferMethod, FieldEvaluationMethod) = 0;
};

} // namespace wdmcpl

#endif // WDM_COUPLING_FIELD_H
