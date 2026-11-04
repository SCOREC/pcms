#include "pcms/field/function_space.hpp"
#include "pcms/configuration.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/common.hpp"
#ifdef PCMS_ENABLE_MESHFIELDS
#include "pcms/field/layout/mesh_fields.hpp"
#include "pcms/field/data/mesh_fields.hpp"
#endif
#ifdef PCMS_ENABLE_OMEGA_H
#include "pcms/field/layout/omega_h_lagrange.hpp"
#include "pcms/field/data/omega_h.hpp"
#endif

namespace pcms
{

namespace
{

// Checks that `data` is what the backend of `layout` stores and that it has the
// layout's size. Layouts of other backends accept any data of the right size;
// data without DOFAccess (values that live outside PCMS) has no size to check.
template <typename T>
void ValidateWrappedFieldData(const FieldLayout& layout,
                              const FieldData<T>& data)
{
#ifdef PCMS_ENABLE_MESHFIELDS
  if (dynamic_cast<const MeshFieldsAdapterLayout*>(&layout) != nullptr &&
      dynamic_cast<const MeshFieldsFieldData<T>*>(&data) == nullptr) {
    throw pcms_error("FunctionSpace::CreateField: MeshFields layout requires "
                     "MeshFieldsFieldData");
  }
#endif
#ifdef PCMS_ENABLE_OMEGA_H
  if (dynamic_cast<const OmegaHLagrangeLayout*>(&layout) != nullptr &&
      dynamic_cast<const OmegaHFieldData<T>*>(&data) == nullptr) {
    throw pcms_error("FunctionSpace::CreateField: OmegaH layout requires "
                     "OmegaHFieldData");
  }
#endif
  if (data.GetDOFAccess() != nullptr &&
      data.GetDOFAccess()->GetValues().size() !=
        detail::ExpectedFlatFieldDataSize(layout)) {
    throw pcms_error(
      "FunctionSpace::CreateField: field data size does not match layout");
  }
}

} // namespace

bool FunctionSpace::Supports(DifferentialOperator op) const noexcept
{
  return evaluator_ != nullptr
           ? evaluator_->Supports(op)
           : (geometry_ != nullptr && op == DifferentialOperator::Value);
}

std::shared_ptr<const FunctionSpace> FunctionSpace::Owner() const
{
  if (weak_from_this().expired()) {
    throw pcms_error("FunctionSpace::CreateFunction: function spaces must be "
                     "owned by std::shared_ptr");
  }
  return shared_from_this();
}

PointEvaluatorVariant FunctionSpace::CreatePointEvaluatorImpl(
  Type value_type, std::shared_ptr<const BoundSource> bound,
  DifferentialOperator op) const
{
  if (evaluator_ == nullptr) {
    throw pcms_error("CreatePointEvaluator: this function space does not "
                     "evaluate fields at points");
  }
  return apply_to_type(value_type, [&](auto type) -> PointEvaluatorVariant {
    using T = typename decltype(type)::type;
    return CreateUntransformedPointEvaluator<T>(*evaluator_, std::move(bound),
                                                op);
  });
}

FieldVariant FunctionSpace::CreateFieldImpl(Type storage_type,
                                            ValueBasis basis) const
{
  auto field_data = create_field_data_(storage_type, std::move(basis));
  return std::visit(
    [this](auto&& fd) -> FieldVariant {
      using FD = std::decay_t<decltype(fd)>;
      using T = typename FD::element_type::value_type;
      return WrapField<T>(layout_, std::forward<decltype(fd)>(fd));
    },
    std::move(field_data));
}

FieldVariant FunctionSpace::CreateFieldImpl(FieldDataVariant data) const
{
  return std::visit(
    [this](auto&& fd) -> FieldVariant {
      using FD = std::decay_t<decltype(fd)>;
      using T = typename FD::element_type::value_type;
      PCMS_ALWAYS_ASSERT(fd != nullptr);
      if constexpr (std::is_same_v<T, int8_t>) {
        throw pcms_error(
          "FunctionSpace::CreateField: int8_t is not a supported field type");
      } else {
        ValidateWrappedFieldData(*layout_, *fd);
        return WrapField<T>(layout_, std::forward<decltype(fd)>(fd));
      }
    },
    std::move(data));
}

} // namespace pcms
