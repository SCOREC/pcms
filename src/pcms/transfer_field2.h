#ifndef PCMS_TRANSFER_FIELD2_H_
#define PCMS_TRANSFER_FIELD2_H_
#include <utility>
#include <typeinfo>
#include "arrays.h"
#include "field.h"
#include "pcms/arrays.h"
#include "pcms/field_evaluation_methods.h"
#include "pcms/profile.h"
#include "pcms/field.h"

namespace pcms
{

template <typename T>
void copy_field2(const FieldT<T>& source, FieldT<T>& target)
{
  PCMS_FUNCTION_TIMER;
  if (typeid(source) != typeid(target)) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Mismatched types");
  }

  target.SetDOFHolderData(source.GetDOFHolderData());
}

template <typename T>
void interpolate_field2(const FieldT<T>& source, FieldT<T>& target)
{
  PCMS_FUNCTION_TIMER;
  if (source.GetCoordinateSystem() != target.GetCoordinateSystem()) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Coordinate system mismatch");
  }

  auto coords = source.GetDOFHolderCoordinates();
  std::vector<double> evaluation(coords.GetCoordinates().size() / 2);
  FieldDataView<Real, HostMemorySpace> data_view{make_array_view(evaluation), source.GetCoordinateSystem()};
  auto locale = source.GetLocalizationHint(coords);
  source.Evaluate(locale, data_view);
  Rank1View<const double, HostMemorySpace> array_view{evaluation.data(), evaluation.size()};
  FieldDataView<const Real, HostMemorySpace> data_const_view{array_view, source.GetCoordinateSystem()};
  target.SetDOFHolderData(data_const_view);
}

} // namespace pcms

#endif // PCMS_TRANSFER_FIELD2_H_
