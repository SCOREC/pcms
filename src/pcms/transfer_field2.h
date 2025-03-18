#ifndef PCMS_TRANSFER_FIELD2_H_
#define PCMS_TRANSFER_FIELD2_H_
#include <utility>
#include <typeinfo>
#include "pcms/arrays.h"
#include "pcms/field_evaluation_methods.h"
#include "pcms/profile.h"
#include "pcms/field.h"

namespace pcms
{

void copy_field2(const FieldT& source, FieldT& target)
{
  if (typeid(source) != typeid(target)) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Mismatched types");
  }

  target.SetNodalData(source.GetNodalData());
}

void interpolate_field2(const FieldT& source, FieldT& target)
{
  if (source.GetCoordinateSystem() != target.GetCoordinateSystem()) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Coordinate system mismatch");
  }
}

} // namespace pcms

#endif // PCMS_TRANSFER_FIELD2_H_
