#ifndef PCMS_FIELD_FUNCTION_SPACE_FUSION_IO_SOURCE_HPP
#define PCMS_FIELD_FUNCTION_SPACE_FUSION_IO_SOURCE_HPP

// Private (not installed): the definition of the opaque FusionIOSource, for
// the Fusion-IO function space and its tests.

#include <fusion_io.h>
#include <map>
#include <memory>
#include <string>

namespace pcms
{

/// Closes a Fusion-IO source. Unlike deleting it, fio_close_source also closes
/// the file, which some sources (M3D-C1) do not do in their destructor.
struct FusionIOClose
{
  void operator()(fio_source* source) const { fio_close_source(&source); }
};

using FusionIOSourcePtr = std::unique_ptr<fio_source, FusionIOClose>;
using FusionIOFieldPtr = std::unique_ptr<fio_field>;

/// An open Fusion-IO source and Fusion-IO's names for its quantities.
class FusionIOSource
{
public:
  /// @param source must have cylindrical geometry
  explicit FusionIOSource(FusionIOSourcePtr source);

  /// Fusion-IO's interface is not const-correct, so this is shallow const.
  [[nodiscard]] fio_source& Get() const noexcept { return *source_; }

  /// Reads a field of the source; throws if Fusion-IO cannot.
  [[nodiscard]] FusionIOFieldPtr GetField(field_type type,
                                          const fio_option_list& options) const;

  /// Fusion-IO's field type for each quantity name the source provides.
  [[nodiscard]] const std::map<std::string, field_type>& Quantities()
    const noexcept
  {
    return quantities_;
  }

private:
  FusionIOSourcePtr source_;
  std::map<std::string, field_type> quantities_;
};

} // namespace pcms

#endif // PCMS_FIELD_FUNCTION_SPACE_FUSION_IO_SOURCE_HPP
