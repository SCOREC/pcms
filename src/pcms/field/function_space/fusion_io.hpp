#ifndef PCMS_FIELD_FUNCTION_SPACE_FUSION_IO_HPP
#define PCMS_FIELD_FUNCTION_SPACE_FUSION_IO_HPP

#include "pcms/field/function_space.h"
#include "pcms/field/layout/empty.h"
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace pcms
{

/// An open Fusion-IO file. Opaque: Fusion-IO is a private dependency.
class FusionIOSource;

enum class FusionIOFormat
{
  M3DC1,
  GEQDSK,
  GATO,
  MARS,
  GPEC
};

/// Fusion-IO's field options. Unset members keep the source's default.
struct FusionIOOptions
{
  enum class Part
  {
    Total,
    Perturbed,
    Equilibrium
  };
  enum class Species
  {
    MainIon,
    Electron
  };

  /// -1 selects the equilibrium.
  std::optional<int> timeslice;
  std::optional<Part> part;
  /// For species-dependent quantities such as density and temperature.
  std::optional<Species> species;
  std::optional<Real> linear_scale;
  /// In degrees.
  std::optional<Real> phase;
};

/// Opens a Fusion-IO file, e.g. an M3D-C1 C1.h5 or a geqdsk file. Only
/// sources with cylindrical geometry are supported. The spaces and fields
/// created from it share and keep it open.
std::shared_ptr<FusionIOSource> OpenFusionIOSource(FusionIOFormat format,
                                                   const std::string& filename);

/// Fusion-IO's names for the quantities `source` provides, e.g.
/// "poloidal flux" or "magnetic field".
std::vector<std::string> AvailableQuantities(const FusionIOSource& source);

/// Fields of a Fusion-IO source, evaluated by Fusion-IO at the query points.
/// Points are (R, phi, Z), which is CylindricalRThetaZ, and vector values are
/// physical cylindrical components. PCMS holds none of the values: functions
/// can be evaluated and used as transfer sources, but not read or written.
///
/// Fusion-IO's (R, phi, Z) is right-handed, like CylindricalRThetaZ, so phi is
/// theta.
///
/// A space holds quantities with one component count (1 for scalars, 3 for
/// vectors); create one space per count over the same source.
class FusionIOFunctionSpace final : public FunctionSpace
{
public:
  static std::shared_ptr<FusionIOFunctionSpace> Create(
    std::shared_ptr<FusionIOSource> source, int num_components);

  std::shared_ptr<const FieldLayout> GetLayout() const noexcept override
  {
    return layout_;
  }

  /// Function for a Fusion-IO quantity, whose component count must be the
  /// space's.
  /// @param quantity Fusion-IO's name for it (see AvailableQuantities)
  Function<Real> CreateFunction(std::string name, const std::string& quantity,
                                const FusionIOOptions& options = {}) const;

protected:
  FieldVariant CreateFieldImpl(Type storage_type,
                               ValueBasis basis) const override;
  FieldVariant CreateFieldImpl(FieldDataVariant data) const override;
  PointEvaluatorVariant CreatePointEvaluatorImpl(
    Type value_type, const EvaluationRequest& request) const override;

private:
  FusionIOFunctionSpace(std::shared_ptr<FusionIOSource> source,
                        int num_components);

  std::shared_ptr<FusionIOSource> source_;
  std::shared_ptr<const EmptyFieldLayout> layout_;
};

} // namespace pcms

#endif // PCMS_FIELD_FUNCTION_SPACE_FUSION_IO_HPP
