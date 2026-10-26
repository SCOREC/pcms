#include "pcms/field/function_space/fusion_io.hpp"
#include "pcms/field/function_space/fusion_io_source.hpp"
#include "pcms/field/coordinate_systems/cylindrical.hpp"
#include "pcms/utility/arrays.h"
#include "pcms/utility/assert.h"
#include "pcms/utility/profile.h"
#include <ranges>

namespace pcms
{

namespace
{

// Names one quantity of a source. PCMS holds none of its values, so it
// provides no DOFAccess.
class FusionIOFieldData final : public FieldData<Real>
{
public:
  FusionIOFieldData(std::shared_ptr<FusionIOSource> source,
                    FusionIOFieldPtr field, ValueBasis basis)
    : source(std::move(source)),
      field(std::move(field)),
      basis_(std::move(basis))
  {
  }

  [[nodiscard]]
  FieldValueType GetValueType() const override
  {
    return ValueTypeOfRank(basis_.Rank());
  }

  [[nodiscard]]
  const ValueBasis& GetValueBasis() const override
  {
    return basis_;
  }

  std::shared_ptr<FusionIOSource> source;
  FusionIOFieldPtr field;

private:
  ValueBasis basis_;
};

// Fusion-IO evaluation is not thread safe, don't call Evaluate concurrently
class FusionIOEvaluator final : public PointEvaluator<Real>
{
public:
  FusionIOEvaluator(std::shared_ptr<FusionIOSource> source,
                    const EvaluationRequest& request)
    : source_(std::move(source))
  {
    if (request.policy.mode != OutOfBoundsMode::ERROR) {
      throw pcms_error("FusionIOFunctionSpace: only OutOfBoundsMode::ERROR is "
                       "supported");
    }
    // Copy to a LayoutRight host view so each point's coordinates are
    // contiguous.
    const auto coords = request.coords.GetValues();
    Kokkos::View<Real**, DeviceMemorySpace> staged(
      "fusion_io_points", coords.extent(0), coords.extent(1));
    CopyDeviceRank2ViewToRank2View(MakeRank2View(staged), coords);
    points_ =
      decltype(points_)("fusion_io_points", coords.extent(0), coords.extent(1));
    DeepCopyMismatchLayouts(points_, staged);
    // Search hints let Fusion-IO start each point's search where it last
    // found it; they belong to the source, so every field shares them. As in
    // Fusion-IO's korc_example, they are one zero-initialized block of
    // sizeof_search_hint() bytes per point. For M3D-C1 a hint is an element to
    // try first, so zero is a valid starting guess.
    hint_size_ = source_->Get().sizeof_search_hint();
    hints_.assign(hint_size_ * points_.extent(0), 0);
  }

  void Evaluate(const Field<Real>& field,
                ValueView<Real, DeviceMemorySpace> values) const override
  {
    PCMS_FUNCTION_TIMER;
    CheckEvaluateWriteTag("FusionIOFunctionSpace", field, values);
    const auto* data = dynamic_cast<const FusionIOFieldData*>(&field.GetData());
    if (data == nullptr || data->source != source_) {
      throw pcms_error(
        "FusionIOFunctionSpace: field is not a quantity of this source");
    }
    const auto out = values.GetValues();
    // Fusion-IO writes the field's components into each row.
    const auto num_components =
      static_cast<size_t>(field.GetLayout().GetNumComponents());
    if (out.extent(0) != points_.extent(0) || out.extent(1) != num_components) {
      throw pcms_error("FusionIOFunctionSpace: output is shaped [" +
                       std::to_string(out.extent(0)) + "][" +
                       std::to_string(out.extent(1)) + "] but expected [" +
                       std::to_string(points_.extent(0)) + "][" +
                       std::to_string(num_components) + "]");
    }
    Kokkos::View<Real**, Kokkos::LayoutRight, HostMemorySpace> host_out(
      "fusion_io_values", out.extent(0), out.extent(1));
    for (size_t i = 0; i < points_.extent(0); ++i) {
      const fio_hint hint = hint_size_ > 0 ? &hints_[i * hint_size_] : nullptr;
      const int status =
        data->field->eval(&points_(i, 0), &host_out(i, 0), hint);
      if (status == FIO_OUT_OF_BOUNDS) {
        throw pcms_error("FusionIOFunctionSpace: query point " +
                         std::to_string(i) + " is outside the source's domain");
      }
      if (status != FIO_SUCCESS) {
        throw pcms_error("FusionIOFunctionSpace: evaluation failed at query "
                         "point " +
                         std::to_string(i) + " (Fusion-IO error " +
                         std::to_string(status) + ")");
      }
    }
    Kokkos::View<Real**, DeviceMemorySpace> device_out(
      "fusion_io_values", out.extent(0), out.extent(1));
    DeepCopyMismatchLayouts(device_out, host_out);
    CopyDeviceRank2ViewToRank2View(out, MakeConstRank2View(device_out));
  }

  Kokkos::View<const LO*, DeviceMemorySpace> FilledPoints() const override
  {
    return {};
  }

private:
  std::shared_ptr<FusionIOSource> source_;
  Kokkos::View<Real**, Kokkos::LayoutRight, HostMemorySpace> points_;
  size_t hint_size_ = 0;
  mutable std::vector<char> hints_;
};

} // namespace

FusionIOSource::FusionIOSource(FusionIOSourcePtr source)
  : source_(std::move(source))
{
  // Points are passed to Fusion-IO as (R, phi, Z), which only holds for
  // cylindrical sources.
  int geometry = -1;
  if (source_->get_int_parameter(FIO_GEOMETRY, &geometry) != FIO_SUCCESS ||
      geometry != FIO_CYLINDRICAL) {
    throw pcms_error("FusionIOSource: only sources with cylindrical geometry "
                     "are supported");
  }
  fio_source::fio_field_list fields;
  source_->get_available_fields(&fields);
  for (const field_type type : fields) {
    std::string name;
    if (fio_get_field_name(type, &name) == FIO_SUCCESS) {
      quantities_.emplace(std::move(name), type);
    }
  }
}

FusionIOFieldPtr FusionIOSource::GetField(field_type type,
                                          const fio_option_list& options) const
{
  fio_field* field = nullptr;
  const int status = source_->get_field(type, &field, &options);
  FusionIOFieldPtr owned(field);
  if (status != FIO_SUCCESS) {
    std::string name;
    fio_get_field_name(type, &name);
    throw pcms_error("FusionIOSource: could not read \"" + name +
                     "\" (Fusion-IO error " + std::to_string(status) + ")");
  }
  return owned;
}

std::shared_ptr<FusionIOSource> OpenFusionIOSource(FusionIOFormat format,
                                                   const std::string& filename)
{
  int type = 0;
  switch (format) {
    case FusionIOFormat::M3DC1: type = FIO_M3DC1_SOURCE; break;
    case FusionIOFormat::GEQDSK: type = FIO_GEQDSK_SOURCE; break;
    case FusionIOFormat::GATO: type = FIO_GATO_SOURCE; break;
    case FusionIOFormat::MARS: type = FIO_MARS_SOURCE; break;
    case FusionIOFormat::GPEC: type = FIO_GPEC_SOURCE; break;
  }
  fio_source* opened = nullptr;
  const int status = fio_open_source(&opened, type, filename.c_str());
  // On failure fio_open_source has already deleted the source.
  if (status != FIO_SUCCESS) {
    throw pcms_error("OpenFusionIOSource: could not open '" + filename +
                     "' (Fusion-IO error " + std::to_string(status) + ")");
  }
  return std::make_shared<FusionIOSource>(FusionIOSourcePtr(opened));
}

std::vector<std::string> AvailableQuantities(const FusionIOSource& source)
{
  const auto names = std::views::keys(source.Quantities());
  return {names.begin(), names.end()};
}

std::shared_ptr<FusionIOFunctionSpace> FusionIOFunctionSpace::Create(
  std::shared_ptr<FusionIOSource> source, int num_components)
{
  return std::shared_ptr<FusionIOFunctionSpace>(
    new FusionIOFunctionSpace(std::move(source), num_components));
}

FusionIOFunctionSpace::FusionIOFunctionSpace(
  std::shared_ptr<FusionIOSource> source, int num_components)
  : source_(std::move(source)),
    layout_(std::make_shared<EmptyFieldLayout>(
      csys::CylindricalRThetaZ::Create(), num_components))
{
  if (source_ == nullptr) {
    throw pcms_error("FusionIOFunctionSpace: source must not be null");
  }
}

Function<Real> FusionIOFunctionSpace::CreateFunction(
  std::string name, const std::string& quantity,
  const FusionIOOptions& options) const
{
  const auto& quantities = source_->Quantities();
  const auto found = quantities.find(quantity);
  if (found == quantities.end()) {
    std::string available;
    for (const auto& known : AvailableQuantities(*source_)) {
      available += (available.empty() ? "" : ", ") + known;
    }
    throw pcms_error(
      "FusionIOFunctionSpace: \"" + quantity +
      "\" is not a quantity of this source; it has: " + available);
  }
  fio_source& source = source_->Get();
  fio_option_list fio_options;
  source.get_field_options(&fio_options);
  if (options.timeslice) {
    fio_options.set_option(FIO_TIMESLICE, *options.timeslice);
  }
  if (options.part) {
    const int parts[] = {FIO_TOTAL, FIO_PERTURBED_ONLY, FIO_EQUILIBRIUM_ONLY};
    fio_options.set_option(FIO_PART, parts[static_cast<int>(*options.part)]);
  }
  if (options.species) {
    const int species[] = {FIO_MAIN_ION, FIO_ELECTRON};
    fio_options.set_option(FIO_SPECIES,
                           species[static_cast<int>(*options.species)]);
  }
  if (options.linear_scale) {
    fio_options.set_option(FIO_LINEAR_SCALE, *options.linear_scale);
  }
  if (options.phase) {
    fio_options.set_option(FIO_PHASE, *options.phase);
  }
  auto field = source_->GetField(found->second, fio_options);
  const int dimension = field->dimension();
  if (dimension != layout_->GetNumComponents()) {
    throw pcms_error("FusionIOFunctionSpace: \"" + quantity + "\" has " +
                     std::to_string(dimension) +
                     " components but the space holds " +
                     std::to_string(layout_->GetNumComponents()));
  }
  // Fusion-IO quantities are scalars or (R, phi, Z) vectors.
  ValueBasis basis = dimension == 3
                       ? ValueBasis{csys::CylindricalRThetaZ::Create(),
                                    ComponentScaling::Physical, values::Vector}
                       : ValueBasis{};
  return FunctionSpace::CreateFunction<Real>(
    std::move(name), std::make_unique<FusionIOFieldData>(
                       source_, std::move(field), std::move(basis)));
}

FieldVariant FusionIOFunctionSpace::CreateFieldImpl(Type /*storage_type*/,
                                                    ValueBasis /*basis*/) const
{
  throw pcms_error("FusionIOFunctionSpace: a value type does not name a "
                   "quantity; use CreateFunction(name, quantity)");
}

FieldVariant FusionIOFunctionSpace::CreateFieldImpl(FieldDataVariant data) const
{
  return WrapField<Real>(
    layout_, std::get<std::unique_ptr<FieldData<Real>>>(std::move(data)));
}

PointEvaluatorVariant FusionIOFunctionSpace::CreatePointEvaluatorImpl(
  Type value_type, const EvaluationRequest& request) const
{
  if (value_type != Type::Real) {
    throw pcms_error("FusionIOFunctionSpace: only Real evaluation is "
                     "supported");
  }
  return std::unique_ptr<PointEvaluator<Real>>(
    std::make_unique<FusionIOEvaluator>(source_, request));
}

} // namespace pcms
