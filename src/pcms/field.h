#ifndef PCMS_COUPLING_FIELD_H
#define PCMS_COUPLING_FIELD_H
#include "adapter/omega_h/omega_h_field.h"
#include "pcms/types.h"
#include "pcms/arrays.h"
#include "pcms/memory_spaces.h"
#include <map>
#include <string>
#include <Kokkos_Core.hpp>
#include <redev.h> // TODO remove this include
#include "pcms/field_evaluation_methods.h" // TODO remove this include
#include "pcms/coordinate_system.h"
#include "pcms/field_layout.h"

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

template <typename T, typename MemorySpace>
class FieldDataView {
public:
  LO Size() {return values.size(); }
  CoordinateSystem GetCoordinateSystem() { return coordinate_system; }

  [[nodiscard]] Rank2View<const T, MemorySpace> GetValues() const noexcept { return values; }
  [[nodiscard]] Rank2View<T, MemorySpace> GetValues() noexcept { return values; }

  // Note: currently don't beleive we should allow changing the coordinate system

private:
  Rank2View<T, MemorySpace> values;
  CoordinateSystem coordinate_system;
};


/*
* The LocalizationHint can hold any data that the underlying field finds useful.
* This is essentially a method to store an external cache of localization information.
* The API of the Field does allow for internal cacheing as well, however some wrapped
* Fields may use a C interface. This avoids the need for additional wrapping of the C
* interface. May re-evaluate the need
*/
struct LocalizationHint {
  void* data = nullptr;
};

/*
* A field expresses the highest level view of operations
* that a user can perform on a field
* Note results must be tagged by the coordinate system type
* Shape functions can be thought of as a particular field type.
 */
template <typename T>
class FieldT {
public:
  virtual const std::string& GetName() const = 0;

  virtual mesh_entity_type GetEntityType() const = 0;

  virtual CoordinateSystem GetCoordinateSystem() const = 0;

  virtual Kokkos::View<const T*> GetNodalData() const = 0;

  virtual Kokkos::View<const T*> GetNodalCoordinates() const = 0;

  virtual void SetNodalData(Kokkos::View<const T*> data) = 0;

  /*
  * specify the coordinates to evaluate the field
   */
  virtual void SetEvaluationCoordinates(CoordinateView<HostMemorySpace> coordinates, LocalizationHint hint = {}) = 0;

  // returns a hint that can be given to the Evaluate method
  // this can be useful to cache data if you have multiple sets of coordinates you may evaluate
  virtual LocalizationHint GetLocalizationHint(CoordinateView<HostMemorySpace> coordinates) = 0;

  // always takes 3D view, dof holder #, dimension, component
  // underlying allocated buffer needs to be #dof holder * # components
  // We return a FieldDataView to make sure we get both the data, and the coordinate system that the data is in
  virtual void Evaluate(FieldDataView<T, HostMemorySpace> results) = 0;

  // should offer component wise version?
  // if data is scalar results are vector, if data is
  // Results should use same coordinate frame as Coordinates passed in
  virtual void EvaluateGradient(FieldDataView<T, HostMemorySpace> results) = 0;

  virtual void SetDOFHolderData(FieldDataView<T, HostMemorySpace> data) = 0;
  virtual CoordinateView<HostMemorySpace> GetDOFHolderCoordinates() = 0;

  virtual const FieldLayout &GetLayout() = 0;
  // number of physical dimensions (typically 1-6)
  //int GetDimension();
  virtual bool CanEvaluateGradient() = 0;

  virtual int Serialize(
    Rank1View<T, pcms::HostMemorySpace> buffer,
    Rank1View<const pcms::LO, pcms::HostMemorySpace>
      permutation) const = 0;

  virtual void Deserialize(
    Rank1View<const T, pcms::HostMemorySpace> buffer,
    Rank1View<const pcms::LO, pcms::HostMemorySpace>
      permutation) const = 0;

  virtual ~FieldT() noexcept = default;
};
// Should statically instantiate types
using FieldPtr = std::variant< FieldT<int8_t>*,
                               FieldT<int32_t>*,
                               FieldT<int64_t>*,
                               FieldT<float>*,
                               FieldT<double>*>;



} // namespace pcms

#endif // PCMS_COUPLING_FIELD_H
