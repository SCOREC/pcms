#ifndef PCMS_EQDSK_FIELD_H
#define PCMS_EQDSK_FIELD_H

#include "pcms/utility/eqdsk.h"
#include "pcms/field/layout/uniform_grid.h"
#include "pcms/field/function_space/spline.h"

#include <string>

namespace pcms
{

struct EQDSKField
{
  SplineFunctionSpace space;
  Field<Real> field;
};

// TODD: @jacobmerson replace with Function.
struct EQDSKFieldWithData
{
  EQDSKData data;
  SplineFunctionSpace space;
  Field<Real> field;
};

/**
 * @brief Create a SplineFunctionSpace and Field from EQDSKData, combining
 * function space creation, field allocation, and DOF holder data population
 * into a single call.
 *
 * Performs three steps:
 * 1. Creates a SplineFunctionSpace from the EQDSK data
 * 2. Creates a Field<Real> on that function space
 * 3. Populates the field's DOF holder with the PSIZR flux data
 *
 * @param eqdsk_data The EQDSK equilibrium data containing the grid and flux
 *                   values.
 * @return EQDSKField A struct containing the function space and populated
 *                    field, ready for evaluation.
 */
inline EQDSKField MakeEQDSKField(const EQDSKData& eqdsk_data)
{
  auto space = SplineFunctionSpace::FromUniformGrid(
    eqdsk_data.grid, CoordinateSystem::Cartesian);
  auto field = space.CreateField<Real>();
  field.GetData().SetDOFHolderData(Rank2View<const Real, DeviceMemorySpace>(
    eqdsk_data.PSIZR.data(), eqdsk_data.PSIZR.extent(0), 1));
  return {std::move(space), std::move(field)};
}

/**
 * @brief Read an EQDSK file and construct a SplineFunctionSpace and Field from
 * it, combining file I/O, function space creation, field allocation, and DOF
 * holder data population into a single call.
 *
 * Performs four steps:
 * 1. Reads the EQDSK file into an EQDSKData
 * 2. Creates a SplineFunctionSpace from the data
 * 3. Creates a Field<Real> on that function space
 * 4. Populates the field's DOF holder with the PSIZR flux data
 *
 * @param filename Path to the EQDSK equilibrium file.
 * @return EQDSKFieldWithData A struct containing the data, function space,
 *                            and populated field, ready for evaluation.
 */
inline EQDSKFieldWithData MakeEQDSKField(const std::string& filename)
{
  auto eqdsk_data = EQDSKData(filename);
  auto [space, field] = MakeEQDSKField(eqdsk_data);
  return {std::move(eqdsk_data), std::move(space), std::move(field)};
}

} // namespace pcms

#endif // PCMS_EQDSK_FIELD_H
