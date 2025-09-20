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
 * Key: result of partition object i.e. rank that the data is sent to on
 * coupling server Value: Vector of local index (ordered)
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
    Rank1View<T, MemorySpace> buffer,
    Rank1View<const pcms::LO, MemorySpace> permutation) const = 0;
  virtual void Deserialize(
    Rank1View<T, MemorySpace> buffer,
    Rank1View<const pcms::LO, MemorySpace> permutation) const = 0;
  virtual std::vector<GO> GetGids() const = 0;
  virtual ReversePartitionMap GetReversePartitionMap(
    const redev::Partition& partition) const = 0;
};
*/

struct LocalizationHint
{
  void* data = nullptr;
};

/*
 * A field expresses the highest level view of operations
 * that a user can perform on a field
 * Note results must be tagged by the coordinate system type
 * Shape functions can be thought of as a particular field type.
 */
template <typename T>
class Field
{
  /*
   * specify the coordinates to evaluate the field
   */
  /*
  virtual void SetEvaluationCoordinates(CoordinateView coordinates,
  LocalizationHint hint = {}) = 0;

  // returns a hint that can be given to the Evaluate method
  // this can be useful to cache data if you have multiple sets of coordinates
  you may evaluate virtual LocalizationHint GetLocalizationHint(CoordinateView
  coordinates) = 0;

  // always takes 3D view, dof holder #, dimension, component
  // underlying allocated buffer needs to be #dof holder * # components
  // We return a FieldDataView to make sure we get both the data, and the
  coordinate system that the data is in virtual void Evaluate(FieldDataView<T>
  results) = 0;

  // should offer component wise version?
  // if data is scalar results are vector, if data is
  // Results should use same coordinate frame as Coordinates passed in
  virtual void EvaluateGradient(FieldDataView<T> results) = 0;

  virtual void SetDOFHolderData(FieldDataView<T> data) = 0;
  virtual CoordinateView GetDOFHolderCoordinates() = 0;

  virtual const FieldLayout &GetLayout() = 0;
  // number of physical dimensions (typically 1-6)
  //int GetDimension();
  virtual bool CanEvaluateGradient() = 0;
  */
  int Serialize(
    Rank1View<T, pcms::HostMemorySpace> buffer,
    Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const = 0;

  void Deserialize(
    Rank1View<const T, pcms::HostMemorySpace> buffer,
    Rank1View<const pcms::LO, pcms::HostMemorySpace> permutation) const = 0;

  virtual ~Field() noexcept = default;
};

} // namespace pcms

#endif // PCMS_COUPLING_FIELD_H
