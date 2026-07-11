#ifndef PCMS_FIELD_DATA_H
#define PCMS_FIELD_DATA_H

#include "field_metadata.h"
#include "pcms/utility/arrays.h"
#include "pcms/utility/memory_spaces.h"
#include <memory>
#include <variant>

namespace pcms
{

// FieldData<T> stores the data that parameterizes a field function, plus
// metadata. The interpretation of the stored data is defined entirely by the
// FieldEvaluatorFactory that consumes it:
//
//   - Mesh-backed fields:    DOF holder coefficients (nodal values, modal
//                            coefficients, etc.) at locations defined by the
//                            layout.
//   - Analytic functions:    function parameters (e.g. wavenumber, amplitude).
//   - ML surrogates:         model weights or a handle to model state.
//
// This generality means FieldData is not restricted to fields with spatial DOF
// locations. HasDOFHolderCoordinates() on the corresponding
// FieldEvaluatorFactory distinguishes spatially-located fields from those
// whose parameter data has no physical-space interpretation.
//
// Concrete backends manage the ownership and storage policy for the underlying
// data.
//
// Contract:
//   - GetMetadata() describes the stored coefficient data.
//   - DOF-holder data is exposed as a 2D [num_dof_holder][num_components] view
//     (the canonical shape); the DOF-holder ordering must match the layout.
//     The default layout_right means it is contiguous node-major (component
//     fastest, i.e. flat index i*num_components + c) on every backend, so a
//     flat/wire representation is just the contiguous data_handle().
//   - Views returned by GetDOFHolderDataHost()/GetDOFHolderData() remain
//     valid until this FieldData object is mutated or destroyed.
//   - SetDOFHolderDataHost()/SetDOFHolderData() replace the entire stored
//     coefficient array; partial updates are not part of this interface.
//   - Basic coefficient operations should be expressible from FieldData and
//     its associated FieldLayout; callers should not need the originating
//     factory for routine get/set workflows.
//
// NOTE: FieldData<T> is being introduced incrementally alongside the existing
// FieldT<T>. Concrete backends will migrate from FieldT<T> to FieldData<T>
// over time.
template <typename T>
class FieldData
{
public:
  using value_type = T;

  virtual const FieldMetadata& GetMetadata() const = 0;

  // Direct access to DOF-holder coefficient data as a 2D
  // [num_dof_holder][num_components] view. The returned view remains valid
  // until the FieldData is mutated or destroyed.
  virtual Rank2View<const T, HostMemorySpace> GetDOFHolderDataHost() const = 0;
  virtual void SetDOFHolderDataHost(
    Rank2View<const T, HostMemorySpace> values) = 0;

  // The returned view remains valid until the FieldData is mutated or
  // destroyed.
  virtual Rank2View<const T, DeviceMemorySpace> GetDOFHolderData() const = 0;
  virtual void SetDOFHolderData(
    Rank2View<const T, DeviceMemorySpace> values) = 0;

  virtual ~FieldData() noexcept = default;
};

using FieldDataVariant = std::variant<
  std::unique_ptr<FieldData<int8_t>>, std::unique_ptr<FieldData<int32_t>>,
  std::unique_ptr<FieldData<int64_t>>, std::unique_ptr<FieldData<float>>,
  std::unique_ptr<FieldData<double>>>;

} // namespace pcms

#endif // PCMS_FIELD_DATA_H
