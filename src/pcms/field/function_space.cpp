#include "pcms/field/function_space.hpp"
#include "pcms/configuration.hpp"
#include "pcms/utility/assert.hpp"
#include "pcms/utility/detail/common.hpp"
#include "pcms/field/detail/simple_field_data.hpp"

namespace pcms
{

FieldDataVariant FunctionSpace::CreateFieldData(Type storage_type,
                                                ValueBasis basis) const
{
  return detail::apply_to_type(storage_type, [&](auto tag) -> FieldDataVariant {
    using T = typename decltype(tag)::type;
    if constexpr (std::is_same_v<T, int8_t>) {
      throw pcms_error("FunctionSpace::CreateField: int8_t is not a supported "
                       "field type");
    } else {
      return std::make_unique<detail::SimpleFieldData<T>>(layout_,
                                                          std::move(basis));
    }
  });
}

void FunctionSpace::ValidateFieldData(const FieldDataVariant& data) const
{
  std::visit(
    [this](const auto& fd) {
      PCMS_ALWAYS_ASSERT(fd != nullptr);
      if (fd->GetDOFAccess() != nullptr &&
          fd->GetDOFAccess()->GetValues().size() !=
            detail::ExpectedFlatFieldDataSize(*layout_)) {
        throw pcms_error("FunctionSpace::CreateField: field data size does "
                         "not match layout");
      }
    },
    data);
}

Quadrature FunctionSpace::CreateQuadrature() const
{
  if (!natural_) {
    throw pcms_error("FunctionSpace::CreateQuadrature: this space has no "
                     "quadrature");
  }
  return natural_(geometry_);
}

bool FunctionSpace::Supports(DifferentialOperator op) const noexcept
{
  return evaluator_ != nullptr ? evaluator_->Supports(op)
                               : op == DifferentialOperator::Value;
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
  Type value_type, std::shared_ptr<const BoundPoints> bound,
  DifferentialOperator op) const
{
  if (evaluator_ == nullptr) {
    throw pcms_error("CreatePointEvaluator: this function space does not "
                     "evaluate fields at points");
  }
  return detail::apply_to_type(value_type,
                               [&](auto type) -> PointEvaluatorVariant {
                                 using T = typename decltype(type)::type;
                                 return CreateUntransformedPointEvaluator<T>(
                                   *evaluator_, std::move(bound), op);
                               });
}

FieldVariant FunctionSpace::CreateFieldImpl(Type storage_type,
                                            ValueBasis basis) const
{
  auto field_data = CreateFieldData(storage_type, std::move(basis));
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
  ValidateFieldData(data);
  return std::visit(
    [this](auto&& fd) -> FieldVariant {
      using FD = std::decay_t<decltype(fd)>;
      using T = typename FD::element_type::value_type;
      PCMS_ALWAYS_ASSERT(fd != nullptr);
      if constexpr (std::is_same_v<T, int8_t>) {
        throw pcms_error(
          "FunctionSpace::CreateField: int8_t is not a supported field type");
      } else {
        return WrapField<T>(layout_, std::forward<decltype(fd)>(fd));
      }
    },
    std::move(data));
}

} // namespace pcms
