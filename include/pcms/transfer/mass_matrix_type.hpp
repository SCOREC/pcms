#ifndef PCMS_TRANSFER_MASS_MATRIX_TYPE_HPP
#define PCMS_TRANSFER_MASS_MATRIX_TYPE_HPP

namespace pcms
{

// Mass-matrix formulation used by Galerkin projections.
enum class MassMatrixType
{
  Consistent, // full Galerkin mass matrix
  Lumped,     // row-sum lumped diagonal mass matrix
};

} // namespace pcms

#endif // PCMS_TRANSFER_MASS_MATRIX_TYPE_HPP
