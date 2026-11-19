// Part of pcms/field/value_view.hpp, which includes it after the declarations
// it needs.
#include "pcms/field/value_view.hpp"
#ifndef PCMS_FIELD_DETAIL_VALUE_VIEW_HPP
#define PCMS_FIELD_DETAIL_VALUE_VIEW_HPP

#include <memory>
#include <optional>
#include <string>

namespace pcms::detail
{

[[nodiscard]] inline FieldValueType ValueTypeOfRank(int rank)
{
  switch (rank) {
    case 0: return FieldValueType::Scalar;
    case 1: return FieldValueType::Vector;
    case 2: return FieldValueType::Tensor;
    default: break;
  }
  throw pcms_error("ValueTypeOfRank: no FieldValueType label for rank " +
                   std::to_string(rank));
}

[[nodiscard]] inline std::optional<ComponentScaling> DerivedComponentScaling(
  const CoordinateSystem& system) noexcept
{
  if (!system.HasOrthogonalBasis()) {
    return ComponentScaling::Natural;
  }
  if (system.HasUnitScaleFactors()) {
    return ComponentScaling::Physical;
  }
  return std::nullopt;
}

[[nodiscard]] inline bool SameValueBasis(const ValueBasis& a,
                                         const ValueBasis& b) noexcept
{
  if (a.Rank() == 0 && b.Rank() == 0) {
    return true;
  }
  if (a.Rank() != b.Rank() ||
      !detail::SameCoordinateSystem(a.system, b.system) ||
      a.scaling != b.scaling) {
    return false;
  }
  return a.scaling == ComponentScaling::Physical || a.variance == b.variance;
}

/// Gate for a tagged write: throws pcms_error unless the basis the caller
/// claimed on the output view is the basis the callee is about to write.
/// @param context callee name used to prefix the error message
/// @param claimed basis tagged onto the caller-supplied output view
/// @param written basis the callee actually writes
inline void CheckWrittenValueBasis(const char* context,
                                   const ValueBasis& claimed,
                                   const ValueBasis& written)
{
  if (!detail::SameValueBasis(claimed, written)) {
    throw pcms_error(std::string(context) +
                     ": the output view's declared value type/basis does not "
                     "match the basis this call writes");
  }
}

inline void ValidateValueBasis(const ValueBasis& basis)
{
  if (basis.Rank() == 0) {
    return;
  }
  if (basis.Rank() > 2) {
    throw pcms_error("ValidateValueBasis: rank-" +
                     std::to_string(basis.Rank()) +
                     " values are not yet supported (maximum supported rank "
                     "is 2)");
  }
  if (basis.system == nullptr) {
    throw pcms_error(
      "ValidateValueBasis: Vector/Tensor values require a basis coordinate "
      "system");
  }
  if (basis.scaling == ComponentScaling::Physical &&
      !basis.system->HasOrthogonalBasis()) {
    throw pcms_error(
      "ValidateValueBasis: coordinate system '" +
      std::string(basis.system->Kind()) +
      "' has a non-orthogonal basis, so physical components are undefined; "
      "use ComponentScaling::Natural, or borrow an orthogonal system's basis");
  }
}

inline void ValidateValueSemantics(const ValueBasis& basis, int num_components)
{
  ValidateValueBasis(basis);
  const int rank = basis.Rank();
  if (rank == 0) {
    return;
  }
  int expected = 1;
  for (int i = 0; i < rank; ++i) {
    expected *= basis.system->Dimension();
  }
  if (num_components != expected) {
    throw pcms_error("ValidateValueSemantics: rank-" + std::to_string(rank) +
                     " values on coordinate system '" +
                     std::string(basis.system->Kind()) + "' (dimension " +
                     std::to_string(basis.system->Dimension()) + ") require " +
                     std::to_string(expected) + " components, got " +
                     std::to_string(num_components));
  }
}

inline ValueBasis MakeValueBasis(
  const std::shared_ptr<const CoordinateSystem>& space_system,
  VarianceSignature variance,
  std::shared_ptr<const CoordinateSystem> system, // null == the space's own
  std::optional<ComponentScaling> scaling)
{
  if (variance.Rank() == 0) {
    return ValueBasis{};
  }
  if (system == nullptr) {
    system = space_system;
  }
  ComponentScaling resolved;
  if (scaling.has_value()) {
    resolved = *scaling;
  } else if (auto derived = detail::DerivedComponentScaling(*system)) {
    resolved = *derived;
  } else {
    throw pcms_error(
      "CreateFunction: coordinate system '" + std::string(system->Kind()) +
      "' admits both physical and natural components and they differ, so "
      "the scaling must be stated (pass ComponentScaling::Physical or "
      "ComponentScaling::Natural)");
  }
  return ValueBasis{std::move(system), resolved, variance};
}

/// @brief Basis of the gradient of values stored in `stored`.
///
/// Component (j..., i) is the derivative along coordinate i of stored
/// component (j...). Components are physical on orthogonal systems and
/// natural otherwise.
/// @param stored basis of the differentiated values
/// @param system coordinate system of the derivatives
/// @return the stored basis with one covariant index appended; the parity is
///         the stored basis's
/// @throws pcms_error if `system` is null, or if `stored` has rank above 0
///         and `stored` and `system` are not one Cartesian system (the
///         gradient of a rank-N field needs the coordinate system's
///         connection)
[[nodiscard]] inline ValueBasis GradientBasis(
  const ValueBasis& stored, std::shared_ptr<const CoordinateSystem> system)
{
  if (system == nullptr) {
    throw pcms_error("GradientBasis: derivatives require a coordinate system");
  }
  const ComponentScaling scaling = system->HasOrthogonalBasis()
                                     ? ComponentScaling::Physical
                                     : ComponentScaling::Natural;
  if (stored.Rank() == 0) {
    return ValueBasis{std::move(system), scaling, values::Covector,
                      stored.parity};
  }
  const bool cartesian = detail::IsCartesian(*system) &&
                         detail::SameCoordinateSystem(stored.system, system);
  if (!cartesian) {
    throw pcms_error("GradientBasis: the gradient of a rank-" +
                     std::to_string(stored.Rank()) +
                     " field needs the coordinate system's connection, which "
                     "is not implemented");
  }
  return ValueBasis{std::move(system), scaling,
                    stored.variance.WithCovariantIndexAppended(),
                    stored.parity};
}

} // namespace pcms::detail

#endif // PCMS_FIELD_DETAIL_VALUE_VIEW_HPP
