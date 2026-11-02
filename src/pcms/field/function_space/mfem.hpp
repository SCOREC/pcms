#ifndef PCMS_MFEM_FIELD_FACTORY_HPP
#define PCMS_MFEM_FIELD_FACTORY_HPP

#include "pcms/field/coordinate_system.h"
#include "pcms/field/data/mfem.hpp"
#include "pcms/field/field_factory.h"
#include "pcms/field/layout/mfem.hpp"
#include "pcms/utility/assert.h"
#include "pcms/utility/common.h"

#include <mfem.hpp>

#include <memory>
#include <string>
#include <variant>

namespace pcms
{

// Communication-only field factory for an MFEM order-1 H1 (vertex) scalar
// field. The produced field's data is bound to the live mfem::ParGridFunction,
// so coupled get/set operate on the solver's state. Point evaluation is not
// supported (mirrors XGCFieldFactory).
// Lifetime: the factory and any Field it produces must not outlive
// pmesh / pfes / gf.
class MFEMFieldFactory : public FieldFactory
{
public:
  MFEMFieldFactory(mfem::ParMesh& pmesh, mfem::ParFiniteElementSpace& pfes,
                   mfem::ParGridFunction& gf,
                   CoordinateSystem coordinate_system,
                   std::string layout_name = "")
    : layout_(std::make_shared<MFEMLayout>(pmesh, pfes, coordinate_system)),
      pfes_(pfes),
      gf_(gf)
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

  [[nodiscard]] std::shared_ptr<const MFEMLayout> GetMFEMLayout() const noexcept
  {
    return layout_;
  }

protected:
  [[nodiscard]] FieldVariant CreateFieldImpl(
    Type value_type, FieldMetadata metadata) const override
  {
    if (value_type != Type::Real) {
      throw pcms_error("MFEM adapter only supports Real (double) fields");
    }
    return WrapField<Real>(
      layout_, std::make_unique<MFEMVertexFieldData>(pfes_, gf_, metadata));
  }

  [[nodiscard]] FieldVariant CreateFieldImpl(
    FieldDataVariant data) const override
  {
    auto* real = std::get_if<std::unique_ptr<FieldData<Real>>>(&data);
    if (real == nullptr || *real == nullptr ||
        dynamic_cast<const MFEMVertexFieldData*>(real->get()) == nullptr) {
      throw pcms_error(
        "MFEMFieldFactory::CreateField: requires MFEMVertexFieldData");
    }
    if ((*real)->GetDOFHolderDataHost().size() !=
        detail::ExpectedFlatFieldDataSize(*layout_)) {
      throw pcms_error("MFEMFieldFactory::CreateField: field data size does "
                       "not match layout");
    }
    return WrapField<Real>(layout_, std::move(*real));
  }

private:
  std::shared_ptr<MFEMLayout> layout_;
  mfem::ParFiniteElementSpace& pfes_;
  mfem::ParGridFunction& gf_;
};

} // namespace pcms

#endif // PCMS_MFEM_FIELD_FACTORY_HPP
