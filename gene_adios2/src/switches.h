!/* Temporary file for some precompiler switches */
!/* which should only be important for development */
!/* Please don't touch them for production runs! */

!/* Switch to Gyromatrix^dagger * profiles * Gyromatrix */
!/* instead of profiles * Gyromatrix^2 */
#undef GDAGGER_G

!/* Switch to Gyromatrix^dagger in <f(x-r)> gyro averages */
!/* this currently causes problems with Neumann b.c. */
#undef G_DAGGER

!/* Test with artificially 'hermitized' Laplacian */
#define APPROX_METRIC

!preconditioner
#define prec_x
#define prec_ara
#undef prec_coll

#undef hyp_over_j

#define with_extended_diags

#define COLLSIX

#define GY_WITH_MAT
