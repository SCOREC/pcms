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
  if (typeid(source) != typeid(target)) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Mismatched types");
  }

  target.SetNodalData(source.GetNodalData());
}

template <typename T>
void interpolate_field2(const FieldT<T>& source, FieldT<T>& target)
{
  if (source.GetCoordinateSystem() != target.GetCoordinateSystem()) {
    // TODO when moved to PCMS throw PCMS exception
    throw std::runtime_error("Coordinate system mismatch");
  }

  auto coordinates = target.GetNodalCoordinates();
  CoordinateView view(source.GetCoordinateSystem(), make_const_array_view(coordinates));
  std::vector<double> evaluation(coordinates.size() / 2);
  FieldDataView<Real, HostMemorySpace> data_view{make_array_view(evaluation), source.GetCoordinateSystem()};
  auto locale = source.GetLocalizationHint(view);
  source.Evaluate(locale, data_view);
  target.SetNodalData(data_view);
}

} // namespace pcms

#endif // PCMS_TRANSFER_FIELD2_H_
