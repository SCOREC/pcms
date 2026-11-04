#include "pcms/backends/xgc/xgc.hpp"
#include "pcms/backends/xgc/field_data.hpp"
#include "pcms/backends/xgc/field_factory.hpp"
#include "pcms/utility/assert.hpp"

#include <cstdint>

namespace pcms::backend::xgc
{

std::shared_ptr<FieldFactory> CreateFieldFactory(
  const ReverseClassificationVertex& reverse_classification,
  LO num_plane_nodes, std::string layout_name)
{
  return std::make_shared<XGCFieldFactory>(
    reverse_classification, num_plane_nodes, std::move(layout_name));
}

template <typename T>
Field<T> CreateField(const FieldFactory& factory, std::string name,
                     Rank1View<T, HostMemorySpace> data, MPI_Comm plane_comm,
                     ValueBasis basis)
{
  const auto* xgc_factory = dynamic_cast<const XGCFieldFactory*>(&factory);
  if (xgc_factory == nullptr) {
    throw pcms_error("pcms::backend::xgc::CreateField: the factory was not "
                     "created by pcms::backend::xgc::CreateFieldFactory");
  }
  return xgc_factory->CreateField<T>(
    std::move(name),
    std::make_unique<XGCFieldData<T>>(xgc_factory->GetXGCLayout(),
                                      std::move(basis), data, plane_comm));
}

#define PCMS_XGC_INSTANTIATE_CREATE_FIELD(T)                                   \
  template Field<T> CreateField<T>(const FieldFactory&, std::string,           \
                                   Rank1View<T, HostMemorySpace>, MPI_Comm,    \
                                   ValueBasis);
PCMS_XGC_INSTANTIATE_CREATE_FIELD(int8_t)
PCMS_XGC_INSTANTIATE_CREATE_FIELD(int32_t)
PCMS_XGC_INSTANTIATE_CREATE_FIELD(int64_t)
PCMS_XGC_INSTANTIATE_CREATE_FIELD(float)
PCMS_XGC_INSTANTIATE_CREATE_FIELD(double)
#undef PCMS_XGC_INSTANTIATE_CREATE_FIELD

} // namespace pcms::backend::xgc
