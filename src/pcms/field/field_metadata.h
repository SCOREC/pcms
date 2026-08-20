#ifndef PCMS_FIELD_METADATA_H
#define PCMS_FIELD_METADATA_H

#include "pcms/field/coordinate_system.hpp"
#include <memory>

namespace pcms
{

enum class FieldValueType
{
  Scalar,
  Vector,
  Tensor
  // Tensor variance and transformation rules are intentionally deferred to a
  // future richer metadata model.
};

struct FieldMetadata
{
  FieldValueType value_type = FieldValueType::Scalar;
  std::shared_ptr<const CoordinateSystem> value_coordinate_system;
};

} // namespace pcms

#endif // PCMS_FIELD_METADATA_H
