! This routine written by P. Worley (March 18, 2015)
!

module lbal_module
  implicit none

  private

  type performance_history
     real (8) :: constraint
     real (8) :: cost
     integer  :: age
  end type performance_history

  type(performance_history) :: perf_history(3)

  real (8), parameter :: golden_ratio    = 1.618033988D0
  real (8), parameter :: phi             = 0.618033988D0
  real (8), parameter :: min_lower_bound = 1.0D0   ! collision load-balancing-specific
  real (8), parameter :: def_upper_bound = 2.0D0   ! arbitrary (> min_lower_bound), but need something
  real (8), parameter :: rel_tol         = 0.001D0 ! arbitrary, but need something
  integer, parameter  :: age_limit       = 3       ! arbitrary, but need something 
                                                   ! (= 3 means around half are recomputes)

  real (8) :: upper_bound      = -1.0D0
  real (8) :: lower_bound      =  1.0D0
  real (8) :: tuning_parameter =  1.0D0
  real (8) :: abs_tol          =  1.0E-10          ! arbitrary, but need something

  public :: update_ptl_constraint_init
  public :: update_ptl_constraint 

contains

subroutine update_ptl_constraint_init(min_tuning_parameter, max_tuning_parameter, init_tuning_parameter)
  use sml_module
  implicit none
  real (8), optional, intent(in) :: min_tuning_parameter
  real (8), optional, intent(in) :: max_tuning_parameter
  real (8), optional, intent(in) :: init_tuning_parameter

  ! initialize history
  perf_history(:)%constraint = -1.0D0
  perf_history(:)%cost  = -1.0D0
  perf_history(:)%age   = -1

  ! Currently only makes sense to call this from "master"
  if (sml_mype .ne. 0) return

  ! set bounds on interval over which to optimize tuning parameter
  if (present(min_tuning_parameter)) then
    lower_bound = max(min_tuning_parameter, min_lower_bound)
  else
    lower_bound = min_lower_bound
  endif

  if (present(max_tuning_parameter)) then
    upper_bound = max(max_tuning_parameter, lower_bound)
  else
    upper_bound = def_upper_bound
  endif

  ! set initial tuning parameter
  if (present(init_tuning_parameter)) then
    tuning_parameter = init_tuning_parameter
  else
    tuning_parameter = (upper_bound + lower_bound)/2.0D0
  endif
  tuning_parameter = min(tuning_parameter,upper_bound)
  tuning_parameter = max(tuning_parameter,lower_bound)

  abs_tol = max(rel_tol*(upper_bound - lower_bound), abs_tol)

end subroutine update_ptl_constraint_init

subroutine update_ptl_constraint(cost, new_tuning_parameter)
  use sml_module
  implicit none
  real (8), intent(in)  :: cost
  real (8), intent(out) :: new_tuning_parameter

  real (8) :: a, b, c
  integer  :: min_cost, max_age

  ! Currently only makes sense to call this from "master"
  if (sml_mype .ne. 0) return

  ! Here recent performance history is saved and used to update the 
  ! tuning parameter (particle load imbalance constraint)
  ! via a Golden-Section-type search algorithm.

  ! This formulation is pretty generic, but the tuning_parameter
  ! is constrained to the interval
  ! [lower_bound, upper_bound] 
  ! where lower_bound >= 1.0

  ! First, call the init routine if it has not been called yet
  if (upper_bound < lower_bound) then
     call update_ptl_constraint_init()
  endif

  min_cost = min(perf_history(1)%cost, &
                 perf_history(2)%cost, &
                 perf_history(3)%cost)
  if (min_cost < 0.0D0) then

    ! Second, if this is one of the first three calls, implement
    ! predetermined strategy.

    if (perf_history(1)%cost < 0.0D0) then

      ! For the first time there has not yet been any collision
      ! load balancing, so initial constraint is 1.0D0 (particle load 
      ! balancing only)
      perf_history(1)%constraint = 1.0D0
      perf_history(1)%cost = cost
      perf_history(1)%age  = 0

      new_tuning_parameter = tuning_parameter

      perf_history(1)%age = perf_history(1)%age + 1

      a = perf_history(1)%constraint
      b = perf_history(2)%constraint
      c = perf_history(3)%constraint

      write(6,21) a, perf_history(1)%age, perf_history(1)%cost, &
                  b, perf_history(2)%age, perf_history(2)%cost, &
                  c, perf_history(3)%age, perf_history(3)%cost
      write(6,22) tuning_parameter

      return

    else if (perf_history(3)%cost < 0.0D0) then

      ! For the second time through initialize (3)
      ! and define the final initial sample location
      perf_history(3)%constraint = tuning_parameter
      perf_history(3)%cost = cost
      perf_history(3)%age  = 0

      a = perf_history(1)%constraint
      c = perf_history(3)%constraint

      ! Sample between lower_bound and the initial tuning parameter, 
      ! at the Golden Section 'x1' point:
      ! tuning_parameter = max((c + b*golden_ratio)/(1.0d0 + golden_ratio), lower_bound)
      tuning_parameter = max((c - (c-a)*phi), lower_bound)
      new_tuning_parameter = tuning_parameter

      perf_history(1)%age = perf_history(1)%age + 1
      perf_history(3)%age = perf_history(3)%age + 1

      b = perf_history(2)%constraint

      write(6,21) a, perf_history(1)%age, perf_history(1)%cost, &
                  b, perf_history(2)%age, perf_history(2)%cost, &
                  c, perf_history(3)%age, perf_history(3)%cost
      write(6,22) tuning_parameter

      return

    else

      ! For the third time through initialize (2)
      ! but use standard logic for choosing next sample
      ! location
      perf_history(2)%constraint = tuning_parameter
      perf_history(2)%cost = cost
      perf_history(2)%age  = 0

    endif

  else

    ! If not one of the first three calls, then insert new data, overwriting old
    if (tuning_parameter > (perf_history(3)%constraint - 0.1*abs_tol)) then
      ! larger than current sample points
      if (abs(tuning_parameter - perf_history(3)%constraint) <= 0.1*abs_tol) then
        perf_history(3)%constraint = tuning_parameter
        perf_history(3)%cost       = cost
        perf_history(3)%age        = 0
      else
        perf_history(1)%constraint = perf_history(2)%constraint
        perf_history(1)%cost       = perf_history(2)%cost
        perf_history(1)%age        = perf_history(2)%age

        perf_history(2)%constraint = perf_history(3)%constraint
        perf_history(2)%cost       = perf_history(3)%cost
        perf_history(2)%age        = perf_history(3)%age

        perf_history(3)%constraint = tuning_parameter
        perf_history(3)%cost       = cost
        perf_history(3)%age        = 0
      endif
    else if (tuning_parameter < (perf_history(1)%constraint + 0.1*abs_tol)) then
      ! smaller than current sample points
      if (abs(tuning_parameter - perf_history(1)%constraint) <= 0.1*abs_tol) then
        perf_history(1)%constraint = tuning_parameter
        perf_history(1)%cost       = cost
        perf_history(1)%age        = 0
      else
        perf_history(3)%constraint = perf_history(2)%constraint
        perf_history(3)%cost       = perf_history(2)%cost
        perf_history(3)%age        = perf_history(2)%age

        perf_history(2)%constraint = perf_history(1)%constraint
        perf_history(2)%cost       = perf_history(1)%cost
        perf_history(2)%age        = perf_history(1)%age

        perf_history(1)%constraint = tuning_parameter
        perf_history(1)%cost       = cost
        perf_history(1)%age        = 0
      endif
    else if (tuning_parameter > (perf_history(2)%constraint - 0.1*abs_tol)) then
      ! between samples 2 and 3
      if (abs(tuning_parameter - perf_history(2)%constraint) <= 0.1*abs_tol) then
        perf_history(2)%constraint = tuning_parameter
        perf_history(2)%cost       = cost
        perf_history(2)%age        = 0
      else
        if (cost < perf_history(2)%cost) then
          ! keep samples 2 and 3
          perf_history(1)%constraint = perf_history(2)%constraint
          perf_history(1)%cost       = perf_history(2)%cost
          perf_history(1)%age        = perf_history(2)%age

          perf_history(2)%constraint = tuning_parameter
          perf_history(2)%cost       = cost
          perf_history(2)%age        = 0
        else
          ! keep samples 1 and 2
          perf_history(3)%constraint = tuning_parameter
          perf_history(3)%cost       = cost
          perf_history(3)%age        = 0
        endif
      endif
    else
      ! between samples 1 and 2
      if (cost < perf_history(2)%cost) then
        ! keep samples 1 and 2
        perf_history(3)%constraint = perf_history(2)%constraint
        perf_history(3)%cost       = perf_history(2)%cost
        perf_history(3)%age        = perf_history(2)%age

        perf_history(2)%constraint = tuning_parameter
        perf_history(2)%cost       = cost
        perf_history(2)%age        = 0
      else
        ! keep samples 2 and 3
        perf_history(1)%constraint = tuning_parameter
        perf_history(1)%cost       = cost
        perf_history(1)%age        = 0
      endif
    endif

    ! make sure that duplicate sample points are also updated 
    if ((abs(perf_history(2)%constraint - perf_history(3)%constraint) <= 0.1*abs_tol) .or. &
        (abs(perf_history(1)%constraint - perf_history(2)%constraint) <= 0.1*abs_tol)) then

      if (abs(perf_history(1)%constraint - perf_history(3)%constraint) <= 0.1*abs_tol) then

        if ((perf_history(3)%age >= perf_history(1)%constraint) .and. &
            (perf_history(2)%age >= perf_history(1)%constraint)) then
          perf_history(2)%constraint = perf_history(1)%constraint
          perf_history(2)%cost       = perf_history(1)%cost
          perf_history(2)%age        = perf_history(1)%age

          perf_history(3)%constraint = perf_history(1)%constraint
          perf_history(3)%cost       = perf_history(1)%cost
          perf_history(3)%age        = perf_history(1)%age
        else if ((perf_history(3)%age >= perf_history(2)%constraint) .and. &
                 (perf_history(1)%age >= perf_history(2)%constraint)) then
          perf_history(1)%constraint = perf_history(2)%constraint
          perf_history(1)%cost       = perf_history(2)%cost
          perf_history(1)%age        = perf_history(2)%age

          perf_history(3)%constraint = perf_history(2)%constraint
          perf_history(3)%cost       = perf_history(2)%cost
          perf_history(3)%age        = perf_history(2)%age
        else
          perf_history(1)%constraint = perf_history(3)%constraint
          perf_history(1)%cost       = perf_history(3)%cost
          perf_history(1)%age        = perf_history(3)%age

          perf_history(2)%constraint = perf_history(3)%constraint
          perf_history(2)%cost       = perf_history(3)%cost
          perf_history(2)%age        = perf_history(3)%age
        endif

      else

        if (abs(perf_history(2)%constraint - perf_history(3)%constraint) <= 0.1*abs_tol) then
          perf_history(2)%constraint = perf_history(3)%constraint
          if (perf_history(3)%age >= perf_history(2)%constraint) then
            perf_history(3)%cost       = perf_history(2)%cost
            perf_history(3)%age        = perf_history(2)%age
          else
            perf_history(2)%cost       = perf_history(3)%cost
            perf_history(2)%age        = perf_history(3)%age
          endif
        endif

        if (abs(perf_history(1)%constraint - perf_history(2)%constraint) <= 0.1*abs_tol) then
          perf_history(2)%constraint = perf_history(1)%constraint
          if (perf_history(2)%age > perf_history(1)%constraint) then
            perf_history(2)%cost       = perf_history(1)%cost
            perf_history(2)%age        = perf_history(1)%age
          else
            perf_history(1)%cost       = perf_history(2)%cost
            perf_history(1)%age        = perf_history(2)%age
          endif
        endif

      endif
  
    endif

  endif

  ! update age of data
  perf_history(1)%age = perf_history(1)%age + 1
  perf_history(2)%age = perf_history(2)%age + 1
  perf_history(3)%age = perf_history(3)%age + 1

  ! adjust tuning_parameter
  a = perf_history(1)%constraint
  b = perf_history(2)%constraint
  c = perf_history(3)%constraint

  max_age = max(perf_history(1)%age, perf_history(2)%age, perf_history(3)%age)

  if (abs(upper_bound-lower_bound) <= abs_tol) then
    ! If lower_bound == upper_bound then tuning_parameter == lower_bound.
    tuning_parameter = lower_bound
    !
  else if (max_age > age_limit) then
    ! if data being used to determine new sample location is 'old', sample again
    ! at that location to make sure that an anomalous value is not interfering
    ! with the optimization
    if (perf_history(1)%age >= max(perf_history(2)%age, perf_history(3)%age)) then
       tuning_parameter = a
    else if (perf_history(3)%age >= max(perf_history(1)%age, perf_history(2)%age)) then
       tuning_parameter = c
    else
       tuning_parameter = b
    endif
    !
  else if ((c-b) < abs_tol) then
    ! if 'converged', re-expand the interval (to avoid getting trapped at
    ! nonoptimal setting if simulation costs change). This may be addressed
    ! naturally by the previous age limit requirement, but including this to be
    ! conservative
    tuning_parameter = (c - b*phi)/(1.0D0-phi)
    if (((tuning_parameter-b) < abs_tol) .or. (tuning_parameter >= upper_bound)) then
      tuning_parameter = max((a - b*phi)/(1.0D0-phi), lower_bound)
    endif
    !
  else
    if (perf_history(2)%cost < min(perf_history(1)%cost,perf_history(3)%cost)) then
      ! always sample in the larger region
      if ((b-a) < (c-b)) then
        ! sample between b and c, at 'x2' Golden Section sample point
        ! tuning_parameter = max((c + b*golden_ratio)/(1.0d0 + golden_ratio), lower_bound)
        tuning_parameter = min(a + (c-a)*phi, upper_bound)
      else
        ! sample between a and b, , at 'x1' Golden Section sample point
        ! tuning_parameter = max((a + b*golden_ratio)/(1.0d0 + golden_ratio), lower_bound)
        tuning_parameter = max(c - (c-a)*phi, lower_bound)
      endif
      !
    else if (perf_history(3)%cost < min(perf_history(1)%cost,perf_history(2)%cost)) then
      ! try a larger sample point, subject to the upper bound constraint
      if ((b-a) < (c-b)) then
        ! take a small step
        tuning_parameter = min((c - b*(1.0D0-phi))/phi, upper_bound)
      else
        ! take a large step
        tuning_parameter = min((c - b*phi)/(1.0D0-phi), upper_bound)
      endif
      !
    else
      ! try a smaller sample point, subject to the lower bound constraint
      if ((b-a) < (c-b)) then
        ! take a large step
        tuning_parameter = max((a - b*phi)/(1.0D0-phi), lower_bound)
      else
        ! take a small step
        tuning_parameter = max((a - b*(1.0D0-phi))/phi, lower_bound)
      endif
    endif
  endif

  new_tuning_parameter = tuning_parameter

  write(6,21) a, perf_history(1)%age, perf_history(1)%cost, &
              b, perf_history(2)%age, perf_history(2)%cost, &
              c, perf_history(3)%age, perf_history(3)%cost
21 format("Performance History: (",f6.3,")[",i1,"] ",e12.3, &
                             ", (",f6.3,")[",i1,"] ",e12.3, &
                             ", (",f6.3,")[",i1,"] ",e12.3  )
  write(6,22) tuning_parameter
22 format("New collision load balance tuning parameter: ",f12.6)

end subroutine update_ptl_constraint

end module lbal_module
