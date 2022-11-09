#include "xgc_field_adapter.h"

wdmcpl::ReadXGCNodeClassificationResult wdmcpl::ReadXGCNodeClassification(
  std::istream& in)
{
  std::int8_t dim;
  wdmcpl::LO id;
 wdmcpl::ReadXGCNodeClassificationResult result;
  while(in >> dim >> id) {
    result.dimension.push_back(dim);
    result.geometric_id.push_back(id);
  }
  return result;
}