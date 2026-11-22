// Part of pcms/field/value_view.hpp, which includes it after the declarations
// it needs.
#include "pcms/field/value_view.hpp"
#ifndef PCMS_FIELD_DETAIL_VALUE_VIEW_HPP
#define PCMS_FIELD_DETAIL_VALUE_VIEW_HPP

#include <memory>
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

/// @return "" unless both bases have units that differ, else why they differ
[[nodiscard]] inline std::string UnitsMismatch(const ValueBasis& a,
                                               const ValueBasis& b)
{
  if (!a.units || !b.units || SameUnits(*a.units, *b.units)) {
    return "";
  }
  return std::string(": units [") + a.units->ToString() + "] and [" +
         b.units->ToString() +
         (SameDimensions(*a.units, *b.units)
            ? "] differ in scale or offset; PCMS does not convert units, so "
              "convert one of the fields"
            : "] have different dimensions");
}

[[nodiscard]] inline bool SameValueBasis(const ValueBasis& a,
                                         const ValueBasis& b)
{
  if (a.units && b.units && !SameUnits(*a.units, *b.units)) {
    return false;
  }
  if (a.Rank() == 0 && b.Rank() == 0) {
    return a.kind == b.kind &&
           (a.kind != ValueKind::Coordinates ||
            detail::SameCoordinateSystem(a.system, b.system));
  }
  if (a.Rank() != b.Rank() ||
      !detail::SameCoordinateSystem(a.system, b.system) ||
      a.scaling != b.scaling) {
    return false;
  }
  return a.scaling == ComponentScaling::Physical || a.variance == b.variance;
}

/// @brief Whether `a` and `b` are the same basis but for their scaling
[[nodiscard]] inline bool SameValueBasisExceptScaling(const ValueBasis& a,
                                                      ValueBasis b) noexcept
{
  b.scaling = a.scaling;
  return SameValueBasis(a, b);
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
                     "match the basis this call writes" +
                     detail::UnitsMismatch(claimed, written));
  }
}

inline void ValidateValueBasis(const ValueBasis& basis)
{
  if (basis.kind != ValueKind::Tensor && basis.Rank() != 0) {
    throw pcms_error("ValidateValueBasis: coordinates and arrays have no "
                     "variance");
  }
  if (basis.kind == ValueKind::Coordinates && basis.system == nullptr) {
    throw pcms_error("ValidateValueBasis: coordinates are points of a "
                     "coordinate system");
  }
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
}

/// @param blocks whether the values may hold several blocks of the basis's
///        components per point (values at points, e.g. the gradient of a
///        field with several scalar components) rather than exactly one (a
///        field's data)
inline void ValidateValueSemantics(const ValueBasis& basis, int num_components,
                                   bool blocks = false)
{
  ValidateValueBasis(basis);
  const int rank = basis.Rank();
  if (basis.kind == ValueKind::Coordinates &&
      num_components != basis.system->Dimension()) {
    throw pcms_error("ValidateValueSemantics: coordinates of '" +
                     std::string(basis.system->Kind()) + "' require " +
                     std::to_string(basis.system->Dimension()) +
                     " components, got " + std::to_string(num_components));
  }
  if (basis.kind != ValueKind::Tensor) {
    return;
  }
  if (rank == 0) {
    if (!blocks && num_components != 1) {
      throw pcms_error("ValidateValueSemantics: a scalar has one component, "
                       "got " +
                       std::to_string(num_components) +
                       " (use ValueKind::Array for independent components)");
    }
    return;
  }
  int expected = 1;
  for (int i = 0; i < rank; ++i) {
    expected *= basis.system->Dimension();
  }
  const bool valid = blocks
                       ? num_components > 0 && num_components % expected == 0
                       : num_components == expected;
  if (!valid) {
    throw pcms_error("ValidateValueSemantics: rank-" + std::to_string(rank) +
                     " values on coordinate system '" +
                     std::string(basis.system->Kind()) + "' (dimension " +
                     std::to_string(basis.system->Dimension()) + ") require " +
                     (blocks ? "a multiple of " : "") +
                     std::to_string(expected) + " components, got " +
                     std::to_string(num_components));
  }
}

/// @return the basis of a new field given none: scalar for one component,
///         independent components (ValueKind::Array) for several
inline ValueBasis DefaultValueBasis(int num_components)
{
  return num_components > 1 ? ValueBasis{.kind = ValueKind::Array}
                            : ValueBasis{};
}

/// @brief The basis of a new field: an unset system is the space's, an unset
///        scaling is Physical for scalars and where the metric is the
///        identity.
/// @throws pcms_error if a vector or tensor on another system does not state
///         its scaling
inline ValueBasis ResolveValueBasis(
  const std::shared_ptr<const CoordinateSystem>& space_system, ValueBasis basis)
{
  if (basis.Rank() > 0 && basis.system == nullptr) {
    basis.system = space_system;
  }
  if (basis.scaling) {
    return basis;
  }
  if (basis.Rank() > 0 && !basis.system->IsOrthonormal()) {
    throw pcms_error(
      "CreateFunction: physical and natural components differ on coordinate "
      "system '" +
      std::string(basis.system->Kind()) +
      "', so the scaling must be stated (set ValueBasis::scaling to "
      "ComponentScaling::Physical or ComponentScaling::Natural)");
  }
  basis.scaling = ComponentScaling::Physical;
  return basis;
}

/// @brief Basis of the gradient of values stored in `stored`.
///
/// Component (j..., i) is the derivative along coordinate i of stored
/// component (j...). Components are natural unless the system is
/// orthonormal.
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
  const ComponentScaling scaling = system->IsOrthonormal()
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
