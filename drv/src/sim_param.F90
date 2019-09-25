module sim_param

   IMPLICIT NONE

   integer :: sml_mype, sml_intpl_mype, sml_intpl_comm, & 
               sml_time, sml_gstep, sml_nphi_total
contains
  subroutine update_sim_param(mype,intpl_mype,intpl_comm,time,gstep,nphi_total)

  integer :: mype, intpl_mype, intpl_comm, time, gstep, nphi_total
                



  end subroutine update_sim_param
end module sim_param
