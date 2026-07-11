#include "pcms/transfer/omega_h_conservative_projection.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include <Omega_h_array.hpp>
#include <memory>

namespace pcms
{

namespace
{

void CheckSupportedLayout(
  const FunctionSpace& space,
  const std::shared_ptr<const OmegaHLagrangeLayout>& layout, const char* role)
{
  if (layout == nullptr) {
    throw pcms_error(std::string("OmegaHConservativeProjection: ") + role +
                     " space must use OmegaHLagrangeLayout");
  }
  if (layout->GetOrder() != 1) {
    throw pcms_error(std::string("OmegaHConservativeProjection: ") + role +
                     " space must be order-1");
  }
  if (layout->GetNumComponents() != 1) {
    throw pcms_error(std::string("OmegaHConservativeProjection: ") + role +
                     " space must have exactly one component");
  }
  if (space.GetCoordinateSystem() != CoordinateSystem::Cartesian) {
    throw pcms_error(std::string("OmegaHConservativeProjection: ") + role +
                     " space must use Cartesian coordinates");
  }
}

void CheckApplyCompatible(const Field<Real>& source, const Field<Real>& target,
                          const OmegaHLagrangeLayout& source_layout,
                          const OmegaHLagrangeLayout& target_layout)
{
  if (&source.GetLayout() != &source_layout) {
    throw pcms_error(
      "OmegaHConservativeProjection::Apply: source field layout mismatch");
  }
  if (&target.GetLayout() != &target_layout) {
    throw pcms_error(
      "OmegaHConservativeProjection::Apply: target field layout mismatch");
  }

  const auto& source_md = source.GetData().GetMetadata();
  const auto& target_md = target.GetData().GetMetadata();
  if (source_md.value_type != FieldValueType::Scalar ||
      target_md.value_type != FieldValueType::Scalar) {
    throw pcms_error(
      "OmegaHConservativeProjection::Apply: only scalar fields are supported");
  }
  if (source_md.value_coordinate_system != target_md.value_coordinate_system) {
    throw pcms_error("OmegaHConservativeProjection::Apply: source and target "
                     "value coordinate systems differ");
  }
}

Omega_h::Reals MakeOmegaHReals(Rank1View<const Real, HostMemorySpace> values)
{
  Omega_h::HostWrite<Omega_h::Real> values_host(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    values_host[i] = values[i];
  }
  return Omega_h::Reals(values_host);
}

} // namespace

OmegaHConservativeProjection::OmegaHConservativeProjection(
  const FunctionSpace& source_space, const FunctionSpace& target_space)
  : source_layout_(std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(
      source_space.GetLayout())),
    target_layout_(std::dynamic_pointer_cast<const OmegaHLagrangeLayout>(
      target_space.GetLayout()))
{
  CheckSupportedLayout(source_space, source_layout_, "source");
  CheckSupportedLayout(target_space, target_layout_, "target");

  intersections_ =
    intersectTargets(source_layout_->GetMesh(), target_layout_->GetMesh());
}

void OmegaHConservativeProjection::Apply(const Field<Real>& source,
                                         Field<Real>& target) const
{
  CheckApplyCompatible(source, target, *source_layout_, *target_layout_);

  const auto source_values =
    MakeOmegaHReals(FlattenToRank1View(source.GetDOFHolderDataHost()));
  const auto target_values = solveGalerkinProjection(
    target_layout_->GetMesh(), source_layout_->GetMesh(), intersections_,
    source_values);
  auto target_values_h = Omega_h::HostRead<Omega_h::Real>(target_values);
  target.SetDOFHolderDataHost(Rank2View<const Real, HostMemorySpace>(
    target_values_h.data(), target_layout_->GetNumOwnedDofHolder(),
    target_layout_->GetNumComponents()));
}

} // namespace pcms
