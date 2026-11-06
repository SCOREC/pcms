#include "pcms/backends/omega_h/omega_h.hpp"
#include "pcms/backends/omega_h/uniform_grid_binary_field.hpp"

namespace pcms::backend::omega_h
{

Field<Real> CreateUniformGridBinaryField(Omega_h::Mesh& mesh,
                                         const UniformGrid<2>& grid)
{
  return CreateUniformGridBinaryField<2>(mesh, grid).second;
}

Field<Real> CreateUniformGridBinaryField(Omega_h::Mesh& mesh,
                                         std::array<LO, 2> divisions)
{
  return CreateUniformGridBinaryField<2>(mesh, divisions).second;
}

} // namespace pcms::backend::omega_h
