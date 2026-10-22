#ifndef PCMS_XGC_FIELD_FACTORY_H
#define PCMS_XGC_FIELD_FACTORY_H

#include "pcms/field/data/xgc.h"
#include "pcms/field/field.h"
#include "pcms/field/field_factory.h"
#include "pcms/field/value_view.hpp"
#include "pcms/utility/common.h"
#include <functional>
#include <memory>

namespace pcms
{

// XGC is a comm-only backend
class XGCFieldFactory : public FieldFactory
{
public:
  XGCFieldFactory(const ReverseClassificationVertex& reverse_classification,
                  std::function<int8_t(int, int)> in_overlap,
                  LO num_plane_nodes, std::string layout_name = "")
    : layout_(std::make_shared<XGCFieldLayout>(
        reverse_classification, std::move(in_overlap), num_plane_nodes))
  {
    layout_->SetName(std::move(layout_name));
  }

  void SetLayoutName(std::string name) const
  {
    layout_->SetName(std::move(name));
  }

  [[nodiscard]] std::shared_ptr<const FieldLayout> GetLayout()
    const noexcept override
  {
    return layout_;
  }

  [[nodiscard]] std::shared_ptr<const XGCFieldLayout> GetXGCLayout()
    const noexcept
  {
    return layout_;
  }

protected:
  [[nodiscard]] FieldVariant CreateFieldImpl(Type storage_type,
                                             ValueBasis basis) const override
  {
    return apply_to_type(storage_type, [&](auto tag) -> FieldVariant {
      using T = typename decltype(tag)::type;
      return WrapField<T>(layout_,
                          std::make_unique<XGCFieldData<T>>(layout_, basis));
    });
  }

  [[nodiscard]] FieldVariant CreateFieldImpl(
    FieldDataVariant data) const override
  {
    return std::visit(
      [this](auto&& fd) -> FieldVariant {
        using FD = std::decay_t<decltype(fd)>;
        using T = typename FD::element_type::value_type;
        PCMS_ALWAYS_ASSERT(fd != nullptr);
        if (dynamic_cast<const XGCFieldData<T>*>(fd.get()) == nullptr) {
          throw pcms_error(
            "XGCFieldFactory::CreateField: requires XGCFieldData");
        }
        if (fd->GetDOFHolderDataHost().size() !=
            static_cast<size_t>(layout_->GetFullDataSize())) {
          throw pcms_error(
            "XGCFieldFactory::CreateField: field data size does not match "
            "layout");
        }
        return WrapField<T>(layout_, std::move(fd));
      },
      std::move(data));
  }

private:
  std::shared_ptr<XGCFieldLayout> layout_;
};

} // namespace pcms

#endif // PCMS_XGC_FIELD_FACTORY_H
