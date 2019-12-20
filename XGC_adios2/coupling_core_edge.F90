#ifndef GENE_SIDE
#include "adios_macro.h"
#endif
module coupling_core_edge
#ifdef ADIOS2
  use adios2
  use adios2_comm_module
#endif
  use adios_comm_module
  use adios_read_mod
  use perf_monitor

  IMPLICIT NONE

  integer :: cce_side ! 0:core, 1: edge
  integer :: cce_density_model
  
  integer, dimension(:),allocatable:: cce_surface_first_node,cce_surface_last_node
  integer :: cce_first_surface,cce_last_surface,cce_surface_number,cce_all_surface_number
  integer :: cce_first_node,cce_last_node,cce_node_number
  integer :: cce_step,cce_field_step,cce_field_model
  integer :: cce_comm_density_mode,cce_comm_field_mode
  integer :: cce_field_first_node,cce_field_last_node,cce_field_node_number
  integer :: cce_first_surface_field,cce_last_surface_field
  
  integer :: cce_first_surface_coupling,cce_last_surface_coupling
  integer :: cce_first_surface_coupling_axis,cce_last_surface_coupling_axis
  
  real, dimension(:), allocatable :: cce_density,cce_pot0,cce_dpot0,cce_dpot1
  real*8 :: cce_alpha
  
  character(256) :: cce_folder  
  character(5) :: cce_my_side,cce_other_side  
  logical ::   cce_bcast_dpot,cce_bcast_pot0,cce_nzeromode,cce_dpot_index0
  
  integer :: cce_npsi
  real (8), dimension(:),allocatable :: cce_varpi,cce_psi
#ifdef ADIOS2
  type(adios2_io),save :: send_io, read_io
  type(adios2_engine),save :: send_engine, read_engine
  type(adios2_variable),save:: field_id
#endif
  real*8,dimension(:),allocatable,save::field_out

contains

  subroutine cce_initialize()

    use sml_module

    !TODO clean this and make the minimum common 
    integer :: ipsi

    namelist /coupling/ cce_side,cce_density_model,cce_folder,         &
       cce_first_surface, cce_last_surface,           &
       cce_first_node,cce_last_node,                  &
       cce_comm_density_mode,cce_all_surface_number,  &
       cce_alpha,cce_field_model,cce_comm_field_mode, &
       cce_first_surface_coupling,cce_last_surface_coupling, &
       cce_first_surface_field,cce_last_surface_field, &
       cce_npsi,cce_first_surface_coupling_axis,      &
       cce_last_surface_coupling_axis,cce_nzeromode,  &
       cce_dpot_index0

    namelist /surfaces/ cce_surface_first_node,cce_surface_last_node
    
    namelist /varpi/ cce_varpi,cce_psi
    
    if(.not.allocated(cce_density))then
     
       cce_alpha=0.5D0
       cce_density_model=0
       cce_step=1 !In case of restart, one might want to change this
       cce_field_step=1
       cce_first_node=10
       cce_last_node=0
       cce_comm_density_mode=2
       cce_all_surface_number=1
       cce_first_surface_coupling=-1
       cce_last_surface_coupling=-1
       cce_field_model=0
       cce_comm_field_mode=0
       cce_npsi=-1
       cce_nzeromode=.false.
       cce_dpot_index0=.false.
       
       open(unit=20,file='coupling.in', status='old',action='read')
       READ(20, NML=coupling)
       close(unit=20)

       if(cce_side==0)then
          cce_my_side='core'
          cce_other_side='edge'
       else
          cce_my_side='edge'
          cce_other_side='core'
       endif
       
       cce_surface_number=cce_last_surface-cce_first_surface+1
       
       cce_node_number=cce_last_node-cce_first_node+1
  
       cce_field_node_number=cce_node_number
       cce_field_first_node=cce_first_node
       cce_field_last_node =cce_last_node
       
       if(cce_node_number.LE.0)then
          allocate(cce_surface_first_node(cce_all_surface_number))
          allocate(cce_surface_last_node(cce_all_surface_number))
          
          open(unit=20,file='coupling.in', status='old',action='read')
          READ(20, NML=surfaces)
          close(unit=20)
          
          cce_first_node=cce_surface_first_node(cce_first_surface)
          cce_last_node=cce_surface_last_node(cce_last_surface)
          cce_node_number=cce_last_node-cce_first_node+1
          
          cce_field_first_node=cce_surface_first_node(1)
          cce_field_last_node =cce_surface_last_node(cce_all_surface_number)
          cce_field_node_number=cce_field_last_node-cce_field_first_node+1
       endif
       
       allocate(cce_density(cce_first_node:cce_last_node))
       
       if(cce_comm_field_mode.NE.0)then
#ifdef XGC_COUPLING_CORE_EDGE_FIELD
          allocate(cce_dpot0(cce_field_first_node:cce_field_last_node))
          allocate(cce_dpot1(cce_field_first_node:cce_field_last_node))
          allocate(cce_pot0(cce_field_first_node:cce_field_last_node))
#else
          print *,'You need to compile with -DXGC_COUPLING_CORE_EDGE_FIELD'
#endif
       endif
!   print *,sml_intpl_mype,'first-last',cce_first_node,cce_last_node
       
       if(cce_first_surface_coupling.lt.1)cce_first_surface_coupling=cce_first_surface
       if(cce_last_surface_coupling.lt.1)cce_last_surface_coupling=cce_last_surface
   
 !  if(sml_restart)then
 !   cce_step=2*sml_istep
 !   cce_field_step=cce_step
 !  endif

       if(cce_npsi>0)then
          allocate(cce_varpi(cce_npsi))
          allocate(cce_psi(cce_npsi))
          open(unit=20,file='coupling.in', status='old',action='read')
          READ(20, NML=varpi)
          close(unit=20)
          cce_first_surface_coupling=-1
          cce_last_surface_coupling=-1
          do ipsi=1,cce_npsi
             if(cce_first_surface_coupling.eq.-1.and.cce_varpi(ipsi)>0.0001)then
                cce_first_surface_coupling=ipsi
             endif
             if(cce_last_surface_coupling.eq.-1.and.cce_varpi(ipsi)>0.9999)then
                cce_last_surface_coupling=ipsi
             endif
          enddo
          print *,"psi_min=",cce_first_surface_coupling,"psi_max=",cce_last_surface_coupling
          if(cce_side.EQ.0)cce_varpi(:)=1D0-cce_varpi(:)
       endif
       
       
    endif
 
!    call adios2_declare_io(read_io, adios2obj, "density_coupling", err)
!    call adios2_open(read_engine, read_io, trim(cce_folder) // &
!      &'/' // 'density.bp', adios2_mode_read, sml_intpl_comm, err)


  end subroutine cce_initialize


  subroutine cce_destroy()
    
    if(allocated(cce_surface_first_node)) deallocate(cce_surface_first_node)
    if(allocated(cce_surface_last_node)) deallocate(cce_surface_last_node)
    if(allocated(cce_density)) deallocate(cce_density)
    if(allocated(cce_dpot0)) deallocate(cce_dpot0)
    if(allocated(cce_dpot1)) deallocate(cce_dpot1)
    if(allocated(cce_pot0)) deallocate(cce_pot0)
    
  end subroutine cce_destroy
  

#ifndef GENE_SIDE
  subroutine cce_receive_density()
    use adios_read_mod
    use sml_module
    implicit none
    include 'mpif.h'
    
    logical :: ex
    character(5)::cce_stepstr,planestr
    character(512) :: cce_filename
    integer*8 :: buf_id, buf_size, total_size
    integer :: err
    integer*8 :: sel2, start2(2), count2(2)
    real*8, dimension(:),allocatable :: arrtmp
#ifdef ADIOS2
    type(adios2_variable) :: varid
    character(4):: fld_name='data'
#endif
    call t_startf("CCE_RECEIVE_DENSITY")
    
    cce_density=0D0

    !if(cce_side.eq.1.and.cce_comm_density_mode.GT.1)then
    if(cce_comm_density_mode.eq.2.or.cce_comm_density_mode.eq.3)then      
       allocate(arrtmp(cce_node_number))
       arrtmp=0D0
       start2(1) = 0
       count2(1) = cce_node_number
       start2(2) = sml_intpl_mype
       count2(2) = 1

#ifdef ADIOS2
       if (cce_step .eq. 1) then
          call adios2_declare_io(read_io, adios2obj, 'density_coupling',&
               & err)
          call adios2_open(read_engine,read_io,trim(cce_folder)//'/density.bp',&
               &adios2_mode_read, sml_intpl_comm, err)
       endif
       call adios2_begin_step(read_engine, adios2_step_mode_read, err)
       call adios2_inquire_variable(varid,read_io,fld_name, err)
       if (.not.varid%valid)then
          if (sml_mype.eq.0) print *, fld_name, ' variable not found'
       else
          if (sml_mype.eq.0) print *, fld_name, ' found'
       endif
       call adios2_set_selection(varid, 2, start2, count2, err)
       call adios2_get(read_engine,varid,arrtmp,adios2_mode_deferred,err)
       call adios2_end_step(read_engine, err)

#else
       ex=.false.
       
       write(cce_stepstr,'(I0.5)') cce_step
       write(planestr,'(I0.5)') sml_intpl_mype
       
#ifndef CCE_TEST_1

       cce_filename=trim(cce_folder)//'/density_'//trim(cce_other_side)//&
            &'_'//trim(cce_stepstr)//'.unlock'
       !    print *,'Wait unlock filename',cce_filename
       do while((.not.ex).and.(staging_read_method.eq.0))
          inquire(file=cce_filename,EXIST=ex)
       end do
       
       cce_filename=trim(cce_folder)//'/density_'//trim(cce_other_side)//&
            &'_'//trim(cce_stepstr)//'.bp'
       if (staging_read_method.gt.0) then
          cce_filename='density_'//trim(cce_other_side)//'.bp'
#ifdef CCE_DEBUG
          if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
          call adios_read_open (buf_id, cce_filename, staging_read_method,&
               &sml_intpl_comm, ADIOS_LOCKMODE_ALL, -1.0, err)
       else
#ifdef CCE_DEBUG
          if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
          call adios_read_open_file (buf_id, cce_filename, staging_read_method, sml_intpl_comm, err)
       end if
       if(err/=0) then
          print *, 'coupling receive error: could not open file', cce_filename
          !stop
       endif

       call adios_selection_boundingbox(sel2, 2, start2, count2)
       call adios_schedule_read (buf_id, sel2, 'data', 0, 1, arrtmp, err)
       call adios_perform_reads (buf_id, err)
       call adios_read_close (buf_id, err)
       
#else
       
       allocate(arrtmp(cce_node_number))
       arrtmp=0D0
       if(sml_intpl_mype.eq.0) arrtmp(90000)=1D0

       cce_filename=trim(cce_folder)//'/density_test_'&
            &//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.bp'
       ADIOS_OPEN(buf_id,'coupling',cce_filename,'w',sml_intpl_comm,err)  !sml_comm_null,err)
       buf_size= 4*8 + 8 + 8*cce_node_number  + 100 !last 100 is buffer
       
       ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
       ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
       ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
       ADIOS_WRITE_LBL(buf_id,'first_node',cce_first_node,err)
       ADIOS_WRITE_LBL(buf_id,'last_node',cce_last_node,err)
       ADIOS_WRITE_LBL(buf_id,'node_number',cce_node_number,err)
       ADIOS_WRITE_LBL(buf_id,'cce_side',cce_side,err)
       ADIOS_WRITE_LBL(buf_id,'cce_model',cce_density_model,err)
       ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
       ADIOS_WRITE_LBL(buf_id,'step',sml_gstep,err)
       !actual data
       ADIOS_WRITE_LBL(buf_id,'data', arrtmp,err)
       
       ADIOS_CLOSE(buf_id,err)
       
#endif

#endif
       
       cce_density(cce_first_node:cce_last_node)=arrtmp(1:cce_node_number)
       
       !print *,trim(cce_my_side)//'density_'//trim(planestr)//'_'//trim(cce_stepstr)//'_R',arrtmp(1),arrtmp(cce_node_number)
       deallocate(arrtmp)
       
    endif
    call t_stopf("CCE_RECEIVE_DENSITY")
  end subroutine cce_receive_density
  
  subroutine cce_process_density(density)
    use sml_module
    implicit none
    include 'mpif.h'

    real*8, dimension(:), intent(inout) :: density
    real*8 :: alpha

    integer :: ipsi,ipsi0,ipsi1

    ! ADIOS
    integer*8 :: buf_id, buf_size, total_size
    integer :: err

    character(5)::cce_stepstr,planestr
    character(512) :: cce_filename
    call t_startf("CCE_PROCESS_DENSITY")

    write(cce_stepstr,'(I0.5)') cce_step
    write(planestr,'(I0.5)') sml_intpl_mype

    !cce_step=cce_step+1

    select case (cce_density_model)
    case (-1)
       density(cce_first_node:cce_last_node)=0D0
    case (0)
       !Do nothing
    case(1)
       !Linear coupling
       if((cce_side.EQ.1).AND.(cce_first_surface.LT.cce_first_surface_coupling))then
          ipsi0=cce_surface_first_node(cce_first_surface)
          ipsi1=cce_surface_last_node(cce_first_surface_coupling)
          density(ipsi0:ipsi1)=cce_density(ipsi0:ipsi1)
       endif
       do ipsi=cce_first_surface_coupling,cce_last_surface_coupling
          alpha=dble(ipsi-cce_first_surface_coupling)/dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
          ipsi0=cce_surface_first_node(ipsi)
          ipsi1=cce_surface_last_node(ipsi)
          if(cce_side.EQ.1)alpha=1D0-alpha
          density(ipsi0:ipsi1)=(1D0-alpha)*density(ipsi0:ipsi1)+alpha*cce_density(ipsi0:ipsi1)
       enddo
       if((cce_side.EQ.0).AND.(cce_last_surface.GT.cce_last_surface_coupling))then
          ipsi0=cce_surface_first_node(cce_last_surface_coupling)
          ipsi1=cce_surface_last_node(cce_last_surface)
          density(ipsi0:ipsi1)=cce_density(ipsi0:ipsi1)
       endif
    case(2)
       !Average
       alpha=cce_alpha
       ipsi0=cce_first_node
       ipsi1=cce_last_node
       density(ipsi0:ipsi1)=(1D0-alpha)*density(ipsi0:ipsi1)+alpha*cce_density(ipsi0:ipsi1)
    case(3)
       ipsi0=cce_first_node
       ipsi1=cce_last_node
       density(ipsi0:ipsi1)=1D0*density(ipsi0:ipsi1)
    case(4)
       alpha=cce_alpha
       ipsi0=cce_first_node
       ipsi1=cce_last_node
       density(ipsi0:ipsi1)=(1D0-alpha)*density(ipsi0:ipsi1)
    case(5)
       if((cce_side.EQ.1))then
          ipsi0=cce_surface_first_node(cce_first_surface)
          ipsi1=cce_surface_last_node(cce_last_surface_coupling)
          density(ipsi0:ipsi1)=cce_density(ipsi0:ipsi1)
       endif
       if((cce_side.EQ.0))then
          ipsi0=cce_surface_first_node(cce_first_surface_coupling)
          ipsi1=cce_surface_last_node(cce_last_surface)
          density(ipsi0:ipsi1)=cce_density(ipsi0:ipsi1)
       endif
#ifdef XGC_COUPLING_CORE_EDGE_VARPI2
    case(6)
       ipsi0=cce_surface_first_node(cce_first_surface) !cce_first_node
       ipsi1=cce_surface_last_node(cce_last_surface)   !cce_last_node
       !print *,'case(6) density(ipsi0:ipsi1)',ipsi0,ipsi1
       density(ipsi0:ipsi1)=density(ipsi0:ipsi1)+cce_density(ipsi0:ipsi1)
#endif
    case default
       print *,'Unknown coupling model'
       stop
    end select
    if(.false.)then
       cce_filename=trim(cce_folder)//'/'//trim(cce_my_side)//'after_'//trim(planestr)//'_'//trim(cce_stepstr)//'.bp'
#ifdef CCE_DEBUG
       if (sml_mype.eq.0) print *, 'writing ', trim(cce_filename)
#endif
       ADIOS_OPEN(buf_id,'coupling',cce_filename,'w',MPI_COMM_SELF,err)  !sml_comm_null,err)
       buf_size= 4*8 + 8 + 8*cce_node_number  + 100 !last 100 is buffer

       ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
       ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
       ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
       ADIOS_WRITE_LBL(buf_id,'first_node',cce_first_node,err)
       ADIOS_WRITE_LBL(buf_id,'last_node',cce_last_node,err)
       ADIOS_WRITE_LBL(buf_id,'node_number',cce_node_number,err)
       ADIOS_WRITE_LBL(buf_id,'cce_side',cce_side,err)
       ADIOS_WRITE_LBL(buf_id,'cce_model',cce_density_model,err)
       ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
       ADIOS_WRITE_LBL(buf_id,'step',sml_gstep,err)
       !actual data
       ADIOS_WRITE_LBL(buf_id,'data', density(cce_first_node:cce_last_node),err)

       ADIOS_CLOSE(buf_id,err)
    endif

    cce_step=cce_step+1

    call t_stopf("CCE_PROCESS_DENSITY")
  end subroutine
#endif

#ifndef GENE_SIDE
#ifdef XGC_COUPLING_CORE_EDGE_FIELD
  !Send the density to the other side
  subroutine cce_send_field(dpot0,dpot1,pot0,flag_pot0)
    use sml_module
    implicit none
    include 'mpif.h'

    real*8, dimension(:), intent(in) :: dpot0,dpot1,pot0
    !real*8, dimension(:,:), intent(in) :: dpot
    integer, intent(in) :: flag_pot0

    character(5)::fld_name
    integer :: myrank, mysize, ierr
    character(5)::cce_stepstr,planestr
    character(512) :: cce_filename
    
    ! ADIOS
    integer*8 :: buf_id, buf_size, total_size
    integer :: err,mode
#ifdef ADIOS2
!    type(adios2_variable) :: varid
    integer(8), dimension(2) :: gdims, goffset, ldims
#endif
    real*8, dimension(:),allocatable :: arrtmp
    call t_startf("CCE_SEND_FIELD")

    if(cce_side.eq.1.and.cce_comm_field_mode.GT.0)then

       write(cce_stepstr,'(I0.5)') cce_field_step
       write(planestr,'(I0.5)') sml_intpl_mype


#ifdef ADIOS2
       if (cce_field_step == 1) then
          allocate(field_out(cce_field_node_number))
          gdims(1) = cce_field_node_number
          gdims(2) = sml_nphi_total
          goffset(1) = 0
          goffset(2) = sml_intpl_mype
          ldims(1) = cce_field_node_number
          ldims(2) = 1
          call adios2_declare_io(send_io, adios2obj, 'field_coupling', err)
          call adios2_define_variable(field_id, send_io, 'dadat',  &
               & adios2_type_dp, 2, gdims, goffset,&
               & ldims, adios2_constant_dims, err)
          call adios2_open(send_engine, send_io, trim(cce_folder) // '/'&
               & //'field.bp', adios2_mode_write,&
               & sml_intpl_comm, err)
          call MPI_Comm_rank(sml_intpl_comm, myrank, ierr);
          call MPI_Comm_size(sml_intpl_comm, mysize, ierr);
          print *, sml_mype, 'opens as rank ', myrank,' of size', mysize

       endif
#endif
       if (.not.allocated(field_out)) allocate(field_out(cce_field_node_number))

       if (cce_dpot_index0) then
          fld_name='dpot0'
          field_out(1:cce_field_node_number)=dpot0(cce_field_first_node:cce_field_last_node)
          if(cce_nzeromode) field_out(1:cce_field_node_number)=field_out(1:cce_field_node_number)+&
               &pot0(cce_field_first_node:cce_field_last_node)
       else  !else of "if(cce_dpot_index0)then"
          fld_name='dpot1'
          field_out(1:cce_field_node_number)=dpot1(cce_field_first_node:cce_field_last_node)
          if(cce_nzeromode) field_out(1:cce_field_node_number)=field_out(1:cce_field_node_number)+&
               pot0(cce_field_first_node:cce_field_last_node)
       endif


#ifdef ADIOS2
       if (sml_mype.eq.0) print *, 'sending field as ', fld_name, kind(field_out)
       call adios2_begin_step(send_engine, adios2_step_mode_append, err) 
       call adios2_put(send_engine, field_id, field_out, err)
       call adios2_end_step(send_engine, err)
#else
       if (staging_read_method.gt.0) then
          cce_filename= trim(fld_name)//'_'//trim(cce_my_side)//'.bp'
       else
          cce_filename=trim(cce_folder)//'/'//trim(fld_name)//&
               &'_'//trim(cce_my_side)//'_'//trim(cce_stepstr)//'.bp'
       endif
#ifdef CCE_DEBUG
       if (sml_mype.eq.0) print *, 'writing ', trim(cce_filename)
#endif
       ADIOS_OPEN(buf_id,'coupling',cce_filename,'w',sml_intpl_comm,err)  !sml_comm_null,err)
       buf_size= 4*8 + 8 + 8*cce_field_node_number  + 100 !last 100 is buffer
       ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
       ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
       ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
       ADIOS_WRITE_LBL(buf_id,'first_node',cce_field_first_node,err)
       ADIOS_WRITE_LBL(buf_id,'last_node',cce_field_last_node,err)
       ADIOS_WRITE_LBL(buf_id,'node_number',cce_field_node_number,err)
       ADIOS_WRITE_LBL(buf_id,'cce_side',cce_side,err)
       ADIOS_WRITE_LBL(buf_id,'cce_model',cce_field_model,err)
       ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
       ADIOS_WRITE_LBL(buf_id,'step',sml_gstep,err)
       ADIOS_WRITE_LBL(buf_id,'data', field_out,err)
       ADIOS_CLOSE(buf_id,err)

       !Create an unlock file: one file per step
       call mpi_barrier(sml_intpl_comm, err)
       if((sml_intpl_mype.eq.0).and.(staging_read_method.eq.0))then
          cce_filename=trim(cce_folder)//'/field_'//trim(cce_my_side)//&
               &'_'//trim(cce_stepstr)//'.unlock'
          open(20, file=cce_filename, status="new", action="write")
          close(20)
       endif
#endif
    endif


    call t_stopf("CCE_SEND_FIELD")
  end subroutine cce_send_field



subroutine cce_process_field(dpot0,dpot1,pot0, flag_pot0)
  use sml_module
  implicit none
  include 'mpif.h'

  real*8, dimension(:)  , intent(inout) :: pot0,dpot0,dpot1
  !real*8, dimension(:,:), intent(inout) :: dpot
  integer, intent(in) :: flag_pot0
  !integer :: ipsi,ipsi0,ipsi1

  !real*8 :: alpha
  call t_startf("CCE_PROCESS_FIELD")

  !cce_bcast_dpot=.false.
  !cce_bcast_pot0=.false.

  select case (cce_field_model)
    case (0)
      !Do nothing
    case(1)
      cce_bcast_dpot=.true.
      dpot0(:)=cce_dpot0(:)
      dpot1(:)=cce_dpot1(:)
      if(flag_pot0.eq.0)then
        cce_bcast_pot0=.true.
        pot0(:)=cce_pot0(:)
      endif
    case default
      print *,'Unknown coupling model'
      stop
  end select

  cce_field_step=cce_field_step+1
  call t_stopf("CCE_PROCESS_FIELD")
end subroutine
#endif
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! useless 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cce_receive_field(flag_pot0)
  use sml_module
  implicit none
  include 'mpif.h'

  integer, intent(in) :: flag_pot0

  logical :: ex
  character(5)::cce_stepstr,planestr
  character(512) :: cce_filename
  integer*8 :: buf_id
  integer :: err
  integer*8 :: sel2, start2(2), count2(2)
  real*8, dimension(:),allocatable :: arrtmp
  call t_startf("CCE_RECEIVE_FIELD")
  
  cce_dpot0=0D0
  cce_dpot1=0D0
  cce_pot0=0D0

  if(cce_side.eq.0.and.cce_comm_field_mode.GT.1)then

    ex=.false.

    write(cce_stepstr,'(I0.5)') cce_field_step
    write(planestr,'(I0.5)') sml_intpl_mype

    cce_filename=trim(cce_folder)//'/field_'//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.unlock'
    do while((.not.ex).and.(staging_read_method.eq.0))
      inquire(file=cce_filename,EXIST=ex)
    end do

    cce_filename=trim(cce_folder)//'/dpot0_'//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.bp'
    if (staging_read_method.gt.0) then
      cce_filename='dpot0_'//trim(cce_other_side)//'.bp'
#ifdef CCE_DEBUG
      if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
      call adios_read_open (buf_id, cce_filename, staging_read_method, sml_intpl_comm, ADIOS_LOCKMODE_ALL, -1.0, err)
    else
#ifdef CCE_DEBUG
      if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
      call adios_read_open_file (buf_id, cce_filename, staging_read_method, sml_intpl_comm, err)
    end if
    if(err/=0) then
      print *, 'coupling receive error: could not open file', cce_filename
      stop
    endif
    allocate(arrtmp(cce_field_node_number))
    arrtmp=0D0

    start2(1) = 0
    count2(1) = cce_field_node_number
    start2(2) = sml_intpl_mype
    count2(2) = 1
    call adios_selection_boundingbox(sel2, 2, start2, count2)
    call adios_schedule_read (buf_id, sel2, 'data', 0, 1, arrtmp, err)
    call adios_perform_reads (buf_id, err)
    call adios_read_close (buf_id, err)
    cce_dpot0(cce_field_first_node:cce_field_last_node)=arrtmp(1:cce_field_node_number)

    cce_filename=trim(cce_folder)//'/dpot1_'//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.bp'
    if (staging_read_method.gt.0) then
      cce_filename='dpot1_'//trim(cce_other_side)//'.bp'
#ifdef CCE_DEBUG
      if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
      call adios_read_open (buf_id, cce_filename, staging_read_method, sml_intpl_comm, ADIOS_LOCKMODE_ALL, -1.0, err)
    else
#ifdef CCE_DEBUG
      if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
      call adios_read_open_file (buf_id, cce_filename, staging_read_method, sml_intpl_comm, err)
    end if
    if(err/=0) then
      print *, 'coupling receive error: could not open file', cce_filename
      stop
    endif
    arrtmp=0D0
    call adios_schedule_read (buf_id, sel2, 'data', 0, 1, arrtmp, err)
    call adios_perform_reads (buf_id, err)
    call adios_read_close (buf_id, err)
    cce_dpot1(cce_field_first_node:cce_field_last_node)=arrtmp(1:cce_field_node_number)

    if(flag_pot0.eq.0)then
      arrtmp=0D0
      cce_filename=trim(cce_folder)//'/pot0_'//trim(cce_other_side)//'_'//trim(cce_stepstr)//'.bp'
      if (staging_read_method.gt.0) then
        cce_filename='pot0_'//trim(cce_other_side)//'.bp'
#ifdef CCE_DEBUG
        if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
        call adios_read_open (buf_id, cce_filename, staging_read_method, sml_intpl_comm, ADIOS_LOCKMODE_ALL, -1.0, err)
      else 
#ifdef CCE_DEBUG
        if (sml_mype.eq.0) print *, 'reading ', trim(cce_filename)
#endif
        call adios_read_open_file (buf_id, cce_filename, staging_read_method, sml_intpl_comm, err)
      end if
      if(err/=0) then
        print *, 'coupling receive error: could not open file', cce_filename
        stop
      endif
      call adios_schedule_read (buf_id, sel2, 'data', 0, 1, arrtmp, err)
      call adios_perform_reads (buf_id, err)
      call adios_read_close (buf_id, err)
      cce_pot0(cce_field_first_node:cce_field_last_node)=arrtmp(1:cce_field_node_number)
    endif

    deallocate(arrtmp)
  endif

  call t_stopf("CCE_RECEIVE_FIELD")
end subroutine

subroutine cce_send_density(density)
  use sml_module
  implicit none
  include 'mpif.h'

  real*8, dimension(:), intent(in) :: density
  character(5)::cce_stepstr,planestr
  character(512) :: cce_filename

  ! ADIOS
  integer*8 :: buf_id, buf_size, total_size
  integer :: err

  real*8, dimension(:),allocatable :: arrtmp
  call t_startf("CCE_SEND_DENSITY")

  cce_density=0D0

  if(cce_comm_density_mode.eq.1.or.cce_comm_density_mode.eq.2)then

    write(cce_stepstr,'(I0.5)') cce_step
    write(planestr,'(I0.5)') sml_intpl_mype

    allocate(arrtmp(cce_node_number))
    arrtmp(1:cce_node_number)=density(cce_first_node:cce_last_node)
    !cce_density(cce_first_node:cce_last_node)=density(cce_first_node:cce_last_node)
    if (staging_read_method.gt.0) then
      cce_filename='density_'//trim(cce_my_side)//'.bp'
    else
      cce_filename=trim(cce_folder)//'/density_'//trim(cce_my_side)//'_'//trim(cce_stepstr)//'.bp'
    endif

  !print *,'Send filename',cce_filename

    !TODO: Open the file and write the density for each node point
#ifdef CCE_DEBUG
    if (sml_mype.eq.0) print *, 'writing ', trim(cce_filename)
#endif
    ADIOS_OPEN(buf_id,'coupling',cce_filename,'w',sml_intpl_comm,err)  !sml_comm_null,err)
    buf_size= 4*8 + 8 + 8*cce_node_number  + 100 !last 100 is buffer

    ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
    ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
    ADIOS_WRITE_LBL(buf_id,'first_node',cce_first_node,err)
    ADIOS_WRITE_LBL(buf_id,'last_node',cce_last_node,err)
    ADIOS_WRITE_LBL(buf_id,'node_number',cce_node_number,err)
    ADIOS_WRITE_LBL(buf_id,'cce_side',cce_side,err)
    ADIOS_WRITE_LBL(buf_id,'cce_model',cce_density_model,err)
    ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
    ADIOS_WRITE_LBL(buf_id,'step',sml_gstep,err)
    !actual data
    ADIOS_WRITE_LBL(buf_id,'data', arrtmp,err)

    ADIOS_CLOSE(buf_id,err)

!    print *,trim(cce_my_side)//'_'//trim(planestr)//'_'//trim(cce_stepstr)//'_W',arrtmp(2),arrtmp(cce_node_number)

    deallocate(arrtmp)

    !Create an unlock file: one file per step
    if((sml_intpl_mype.eq.0).and.(staging_read_method.eq.0))then
      cce_filename=trim(cce_folder)//'/density_'//trim(cce_my_side)//'_'//trim(cce_stepstr)//'.unlock'
      open(20, file=cce_filename, status="new", action="write")
      close(20)
    endif

  endif
  call t_stopf("CCE_SEND_DENSITY")
end subroutine


!#ifdef XGC_COUPLING_CORE_EDGE_VARPI2
#ifndef GENE_SIDE
    subroutine cce_varpi_grid(rho_ff)
      real (8), dimension(:), intent(inout) :: rho_ff
      integer :: ipsi,ipsi0,ipsi1
      real (8) :: varpi

      call cce_initialize()
      if(cce_density_model.eq.6)then

        if(cce_npsi>0)then
          if((cce_side.EQ.1).AND.(cce_first_surface.LT.cce_first_surface_coupling))then
            ipsi0=cce_surface_first_node(cce_first_surface)
            ipsi1=cce_surface_last_node(cce_first_surface_coupling-1)
            rho_ff(ipsi0:ipsi1)=0D0!rho_ff(ipsi0:ipsi1)
          endif
          do ipsi=cce_first_surface_coupling,cce_last_surface_coupling
            ipsi0=cce_surface_first_node(ipsi)
            ipsi1=cce_surface_last_node(ipsi)
            rho_ff(ipsi0:ipsi1)=cce_varpi(ipsi)*rho_ff(ipsi0:ipsi1)
          enddo
          if((cce_side.EQ.0).AND.(cce_last_surface.GT.cce_last_surface_coupling))then
            ipsi0=cce_surface_first_node(cce_last_surface_coupling+1)
            ipsi1=cce_surface_last_node(cce_last_surface)
            rho_ff(ipsi0:ipsi1)=0D0 !rho_ff(ipsi0:ipsi1)
          endif

        else
           !if(ipsi.le.cce_first_surface_coupling_axis)then
           !   varpi=0D0
           !elseif(ipsi.le.cce_last_surface_coupling_axis)then
           !   varpi=1D0-dble(ipsi-cce_first_surface_coupling_axis)/dble(cce_last_surface_coupling_axis-cce_first_surface_coupling_axis+1)
           !elseif (ipsi.le.cce_first_surface_coupling) then
           !   varpi=1D0
           !elseif (ipsi.gt.cce_last_surface_coupling)then
           !   varpi=0D0
           !else
           !   varpi=1D0-dble(ipsi-cce_first_surface_coupling)/dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
           !endif

          !Linear weight
          if(cce_side.EQ.0)then !core
            !1. Zero near axis
            ipsi0=cce_surface_first_node(cce_first_surface)
            ipsi1=cce_surface_last_node(cce_first_surface_coupling_axis)
            rho_ff(ipsi0:ipsi1)=0D0 !rho_ff(ipsi0:ipsi1)
            !2. Linear increasing weight 
            do ipsi=cce_first_surface_coupling_axis,cce_last_surface_coupling_axis
              varpi=dble(ipsi-cce_first_surface_coupling_axis)/&
                   &dble(cce_last_surface_coupling_axis-cce_first_surface_coupling_axis+1)
              ipsi0=cce_surface_first_node(ipsi)
              ipsi1=cce_surface_last_node(ipsi)
              rho_ff(ipsi0:ipsi1)=varpi*rho_ff(ipsi0:ipsi1)
            enddo
            !3. Unity in the middle
            
            !4. Linear decreasinging weight
            do ipsi=cce_first_surface_coupling,cce_last_surface_coupling
              varpi=1D0-dble(ipsi-cce_first_surface_coupling)/dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
              ipsi0=cce_surface_first_node(ipsi)
              ipsi1=cce_surface_last_node(ipsi)
              !varpi=1D0-varpi
              rho_ff(ipsi0:ipsi1)=varpi*rho_ff(ipsi0:ipsi1)
            enddo
            !5. Zero near edge
            ipsi0=cce_surface_first_node(cce_last_surface_coupling)
            ipsi1=cce_surface_last_node(cce_last_surface)
            rho_ff(ipsi0:ipsi1)=0D0 !rho_ff(ipsi0:ipsi1)
          elseif(cce_side.EQ.1)then !edge
            ipsi0=cce_surface_first_node(cce_last_surface_coupling_axis)
            ipsi1=cce_surface_last_node(cce_first_surface_coupling)
            rho_ff(ipsi0:ipsi1)=0D0!rho_ff(ipsi0:ipsi1)

            do ipsi=cce_first_surface_coupling_axis,cce_last_surface_coupling_axis
              varpi=dble(ipsi-cce_first_surface_coupling_axis)/&
                   &dble(cce_last_surface_coupling_axis-cce_first_surface_coupling_axis+1)
              ipsi0=cce_surface_first_node(ipsi)
              ipsi1=cce_surface_last_node(ipsi)
              varpi=1D0-varpi
              rho_ff(ipsi0:ipsi1)=varpi*rho_ff(ipsi0:ipsi1)
            enddo
            do ipsi=cce_first_surface_coupling,cce_last_surface_coupling
              varpi=dble(ipsi-cce_first_surface_coupling)/&
                   &dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
              ipsi0=cce_surface_first_node(ipsi)
              ipsi1=cce_surface_last_node(ipsi)
              rho_ff(ipsi0:ipsi1)=varpi*rho_ff(ipsi0:ipsi1)
            enddo
          endif
        endif
      endif
    end subroutine
#else
  !GENE SIDE
  function cce_varpi_grid(ipsi) result(varpi)
    real :: varpi
    integer :: ipsi

!      call cce_initialize()
    if(cce_density_model.eq.6)then
      if(cce_npsi>0)then
          stop('Not implemented');
      else
        !Linear weight
          if (cce_side.EQ.0)then
             if(ipsi.le.cce_first_surface_coupling_axis)then
                varpi=0D0
             elseif(ipsi.le.cce_last_surface_coupling_axis)then
                varpi=dble(ipsi-cce_first_surface_coupling_axis)/dble(cce_last_surface_coupling_axis-cce_first_surface_coupling_axis+1)
             elseif (ipsi.le.cce_first_surface_coupling) then
                varpi=1D0
             elseif (ipsi.gt.cce_last_surface_coupling)then
                varpi=0D0
             else
                varpi=1D0-dble(ipsi-cce_first_surface_coupling)/dble(cce_last_surface_coupling-cce_first_surface_coupling+1)
             endif
           else
             stop('GENE is for the core, put cce_side=0');
          endif
        endif
      endif
    end function
#endif

end module



!<adios-group name="coupling">
!  <var name="nphi" type="integer"/>
!  <var name="iphi" type="integer"/>
!  <var name="first_node" type="integer"/>
!  <var name="last_node" type="integer"/>
!  <var name="nodes_number" type="integer"/>
!  <var name="cce_side" type="integer"/>
!  <var name="cce_density_model" type="integer"/>
!  <var name="time" type="real*8"/>
!  <var name="step" type="integer"/>
!
!!  <global-bounds dimensions="nphi,nodes_number" offsets="iphi,0">
!!     <var name="density" type="real*8" dimensions="1,nodes_number"/>
!!  </global-bounds>
!  <global-bounds dimensions="node_number" offsets="0">
!     <var name="density" type="real*8" dimensions="node_number"/>
!  </global-bounds>
!</adios-group>
!
!<method priority="3" method="MPI" iterations="100" group="coupling"/>
