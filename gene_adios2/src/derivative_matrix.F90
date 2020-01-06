#include "redef.h"
#include "intrinsic_sizes.h"
MODULE DerivativeMatrixModule
  use Grid1DModule
  USE BandedMatrixModule
  implicit none

  private

  public :: static_size_of_DerivativeMatrix

  type,extends(BandedMatrix),public ::  DerivativeMatrix
     integer :: derivative_order !the order of the finite differences formulas for the derivatives
     logical :: periodic_boundary
     TYPE(Grid1D), pointer :: grid
   contains
     procedure,private ::  dm_initialize_matrix
     generic,public :: initialize => dm_initialize_matrix
     final :: dm_finalize_matrix

     procedure,private :: dm_assign_matrix
     generic,public :: assignment(=) => dm_assign_matrix

     procedure :: calculate
  END TYPE DerivativeMatrix

  !interface DerivativeMatrix
  !   module procedure dm_initialize_matrix
  !end interface DerivativeMatrix
CONTAINS

  FUNCTION static_size_of_DerivativeMatrix() RESULT(memory_need)
    integer :: memory_need

    memory_need = static_size_of_BandedMatrix()&
         & + SIZE_OF_INTEGER&
         & + SIZE_OF_LOGICAL

  END FUNCTION static_size_of_DerivativeMatrix

  SUBROUTINE mp_initialize_module
  END SUBROUTINE mp_initialize_module

  !> Initialize Derivative matrices
  !function dm_initialize_matrix(grid,p_derivative_order,transposed) result(dmat)
  subroutine dm_initialize_matrix(this,grid,p_derivative_order,transposed)
    class(DerivativeMatrix),intent(inout) :: this !> derivative matrix
    TYPE(Grid1D),TARGET,intent(in) :: grid     !> grid definition
    INTEGER,intent(in) :: p_derivative_order   !> derivative order
    LOGICAL, OPTIONAL :: transposed
    !type(DerivativeMatrix) :: dmat

    ! Local variables
    integer :: n_points

    n_points = get_number_of_nodes(grid)

    IF (PRESENT(transposed)) THEN
      call this%BandedMatrix%initialize(n_points,n_points,transposed)
      !this%BandedMatrix = BandedMatrix(n_points,n_points,transposed)
    ELSE
      call this%BandedMatrix%initialize(n_points,n_points)
      !this%BandedMatrix = BandedMatrix(n_points,n_points)
    END IF

    call this%BandedMatrix%allocate()
    this%derivative_order = p_derivative_order
    this%periodic_boundary = is_periodic_boundary(grid)

    this%grid => grid
  end subroutine dm_initialize_matrix

  subroutine dm_finalize_matrix(this)
    type(DerivativeMatrix) :: this
  end subroutine dm_finalize_matrix

  SUBROUTINE dm_assign_matrix(lmat,rmat)
    Class(DerivativeMatrix), intent(INOUT) :: lmat
    type(DerivativeMatrix), intent(IN) :: rmat

    lmat%BandedMatrix = rmat%BandedMatrix
    lmat%derivative_order = rmat%derivative_order
    lmat%periodic_boundary = rmat%periodic_boundary
    if (associated(rmat%grid)) then
       lmat%grid => rmat%grid
    else
       nullify(lmat%grid)
    end if
   end SUBROUTINE dm_assign_matrix

  !> Calculate derivative matrices
  SUBROUTINE calculate(this, which_derivative, rad_bc_type)
    Class(DerivativeMatrix) :: this             !> matrix to contain the derivative matrix
    INTEGER, intent(IN) :: which_derivative    !> number of the derivative to calculate the matrix for
    !                   can be 1 or 2 for the first or second derivative
    INTEGER,intent(in) :: rad_bc_type          !> radial boundary condition

    !local variables
    REAL :: dx
    REAL, DIMENSION(:),ALLOCATABLE :: sten
    INTEGER :: i, s, n_points

    PERFON('deriv_mat')

    call this%set_zero()

    IF (this%grid%isEquidistant) THEN

       ALLOCATE(sten(-this%derivative_order/2:this%derivative_order/2))
       dx = get_gridspacing(this%grid)

       SELECT CASE (this%derivative_order)
       CASE(2)
          !PRINT*,"in calculate, derivative_order=2"
          SELECT CASE (which_derivative)
          CASE(1) ! first derivative
             !PRINT*,"in calculate, derivative_order=2, first derivative"
             sten = (/-0.5,0.0,0.5/) / dx
          CASE(2) ! second derivative
             sten = (/1.0,-2.0,1.0/) / dx**2
          END SELECT
       CASE(4)
          SELECT CASE (which_derivative)
          CASE(1)
             sten = (/1.0,-8.0,0.0,8.0,-1.0/) / dx
          CASE(2)
             sten = (/-1.0,16.0,-30.0,16.0,-1.0/) / dx**2
          END SELECT
          sten=sten/12.0
       CASE(6)
          SELECT CASE (which_derivative)
          CASE(1)
             sten = (/-1.0,9.0,-45.0,0.0,45.0,-9.0,1.0/) / dx
          case default
             PRINT*,"for sixth order derivative, the second derivative is not yet implemented."
             STOP
          END SELECT
          sten = sten/60.0
       CASE default
          PRINT*,"only derivation order 2, 4 and 6 implemented"
          STOP
       END SELECT

       !PRINT*, sten
       n_points = get_number_of_nodes(this%grid)
       DO s=-this%derivative_order/2,this%derivative_order/2
          DO i=1,n_points
             IF (((i+s).LT.1).OR.((i+s).GT.n_points)) THEN !outside box
                IF (rad_bc_type.EQ.0) THEN !periodic b.c.
                   IF (is_periodic_boundary(this%grid)) THEN
                      ! The last and the first point are identical and lie on the
                      ! left or right boundary, respectively. So the derivatives
                      ! only run over the nintervals columns and rows of the matrix.
                      CALL this%add_value(i,MOD(i+s-1+n_points,n_points)+1, sten(s))
                   ELSE
                      STOP 'rad_bc_type is 0, but the grid is nonperiodic'
                   END IF
                ELSEIF (rad_bc_type.EQ.1) THEN !Dirichlet b.c. at both boundaries
                   !do nothing (zero's outside box)
                ELSEIF (rad_bc_type.EQ.2) THEN !Neumann at inner and Dirichlet b.c.
                   !at outer boundary, i.e. mirror grid points at inner boundary
                   IF ((i+s).LT.1) THEN
                      Call this%add_value(i,1-(i+s), sten(s))
                   ENDIF
                ELSE
                   STOP 'rad_bc_type>2 not implemented in derivative matrix construction '
                ENDIF
             ELSE !inside box
                CALL this%add_value(i,i+s,sten(s))
             ENDIF
          END DO
       END DO
       call this%commit_values()
       DEALLOCATE(sten)

    ELSE
       ! not equidistant underlying grid
       ! this is for the future and has to be implemented
    END IF

    PERFOFF
  END SUBROUTINE calculate

END MODULE DerivativeMatrixModule
