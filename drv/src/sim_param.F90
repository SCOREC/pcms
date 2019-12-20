module sim_param

   IMPLICIT NONE

   integer :: sml_mype, sml_intpl_mype, sml_intpl_comm, & 
               sml_time, sml_gstep, sml_nphi_total, mpi_comm_x

       ! process coordinates as seen in discretization.F90
   integer :: my_pespec, my_pew, my_pez, my_pev, my_mpi_comm_world, mype
   real :: time =0.     !current time
   integer :: itime = 0
   real :: cref
   real :: Lref = 0
contains
  subroutine update_sim_param(mype,intpl_mype,intpl_comm,time,gstep,nphi_total)

  integer :: mype, intpl_mype, intpl_comm, time, gstep, nphi_total
                



  end subroutine update_sim_param
end module sim_param
