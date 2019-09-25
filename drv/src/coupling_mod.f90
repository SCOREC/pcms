! from xgc
module coupling_module
  use mpi
  use adios_comm_module
  implicit none
  logical :: sml_coupling_on
  integer :: sml_coupling_mstep
  integer :: sml_gstep0
  integer :: sml_nstep_xgc1
  integer :: sml_coupling_firstfile, sml_coupling_lastfile
  integer :: sml_coupling_nskip
  character (len=65) :: sml_coupling_xgc1_dir
  character (len=65) :: sml_coupling_xgca_dir
  real (8), allocatable :: turb_pot(:,:,:), turb_pot0(:,:)
contains
  !------------------------------------------------------------------
  subroutine coupling_init()
    implicit none
    call check_point('coupling_init')
    sml_coupling_on = .false.
    sml_coupling_mstep = 0
    sml_gstep0 = 0
    sml_nstep_xgc1 = 0
    sml_coupling_firstfile = 0
    sml_coupling_lastfile = 0
    sml_coupling_nskip = 0
    write (sml_coupling_xgc1_dir, "(A)") "."
    write (sml_coupling_xgca_dir, "(A)") "."
  end subroutine coupling_init
  !------------------------------------------------------------------
  subroutine coupling_setup(n)
    implicit none
    integer :: n ! number of nodes
    call check_point('coupling_setup')
    allocate(turb_pot0(n,2))
    allocate(turb_pot(n,-1:2,2))
  end subroutine coupling_setup
  !------------------------------------------------------------------
  subroutine coupling_particle_write(spall)
    use f0_module
    use sml_module
    use ptl_module
    use perf_monitor
    implicit none
    type(species_type):: spall(0:ptl_nsp_max)
    integer :: i,j
    character (len=50) :: filename, dirname
    integer*8 :: buf_id, buf_size, total_size
    integer :: err
    real*8 start_time, end_time !SAK added this for timers
    integer :: np
    real (8), allocatable :: phase(:,:), ephase(:,:)
    integer (8), allocatable :: gid(:), egid(:)
    integer, parameter :: ict1=ptl_nphase+1
    integer, parameter :: ict2=ptl_nphase+ptl_nconst
    ! f0 grid
    integer :: ndata, mudata, vpdata
    integer :: inum ! ion particle number to be saved per proc
    integer*8 :: ioffset, itotal
    integer :: enum ! electron particle number to be saved per proc
    integer*8 :: eoffset, etotal
    integer   :: istep0, ierr
    real (8) :: phasesum(ict2)
    call check_point('coupling_init')
    phasesum=0.D0
    np=spall(1)%num
    allocate(phase(ict2,np),gid(np))
    do i=1, np
       phase(1:ptl_nphase,i)=spall(1)%ptl(i)%ph
       phase(ict1:ict2,i)=spall(1)%ptl(i)%ct
       gid(i)=spall(1)%ptl(i)%gid
    enddo
    if(sml_electron_on)then
       allocate(ephase(ict2,np),egid(np))
       ! copy to temp. memory
       do i=1, np
          ephase(1:ptl_nphase,i)=spall(0)%ptl(i)%ph
          ephase(ict1:ict2   ,i)=spall(0)%ptl(i)%ct
          egid(i)=spall(0)%ptl(i)%gid
       enddo
    endif
    !! Sketch (Jong)
    !! Save in interplane major order. I.e.,
    !! | (interplane #1) | ... | (interplance #N) |
    !! where N = PEs per plane (= sml_plane_totalpe)
    inum = spall(1)%num
    call get_offset(inum,ioffset,itotal)
    buf_size = 4*8 + 8*4 + 8*spall(1)%num + 8*spall(1)%num*ict2
    if (sml_electron_on) then
       buf_size = buf_size + 4*3 + 8*3 + 8*spall(0)%num + 8*spall(0)%num*ict2
    endif
    !! call t_startf("COUPLINGP_WRITE")
    write(filename,'("xgc.couplingp.xgc1.",i5.5,".bp")') sml_gstep
    if (sml_mype==0) print *, 'Writing', filename
    ADIOS_OPEN(buf_id,'couplingp',filename,'w',sml_comm,err)
    ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
    istep0=max(sml_gstep-sml_istep+1,1)
    ADIOS_WRITE_LBL(buf_id,'timestep0',istep0,err)
    ADIOS_WRITE_LBL(buf_id,'timestep',sml_gstep,err)
    ADIOS_WRITE_LBL(buf_id,'run_count',sml_run_count,err)
    ADIOS_WRITE_LBL(buf_id,'time',sml_time,err)
    ADIOS_WRITE_LBL(buf_id,'totalpe',sml_totalpe,err)
    ADIOS_WRITE_LBL(buf_id,'mype',sml_mype,err)
    ADIOS_WRITE_LBL(buf_id,'inum_arr',spall(1)%num,err)
    ADIOS_WRITE_LBL(buf_id,'ioff_arr',ioffset,err)
    ADIOS_WRITE_LBL(buf_id,'inum_total',itotal,err)
    ADIOS_WRITE_LBL(buf_id,'inum',spall(1)%num,err)
    ADIOS_WRITE_LBL(buf_id,'ioff',ioffset,err)
    ADIOS_WRITE_LBL(buf_id,'inphase',ict2,err)
    ADIOS_WRITE_LBL(buf_id,'igid',gid,err)
    ADIOS_WRITE_LBL(buf_id,'iphase',phase,err)
    !! Check values. Copy and paste to Matlab. In Matlab, use adiosread function
    call mpi_allreduce(sum(phase, 2),phasesum,ict2,MPI_DOUBLE,MPI_SUM,sml_comm,ierr)
    if (sml_mype == 0) then
       do i=1, ict2
          print *, 'phasesum(', i, ') =', phasesum(i), ';'
       end do
    end if
    if (sml_electron_on) then
       enum = spall(0)%num
          call get_offset(enum, eoffset,etotal)
          ADIOS_WRITE_LBL(buf_id,'enum_arr',spall(0)%num,err)
          ADIOS_WRITE_LBL(buf_id,'eoff_arr',eoffset,err)
          ADIOS_WRITE_LBL(buf_id,'enum_total',etotal,err)
          ADIOS_WRITE_LBL(buf_id,'enum',spall(0)%num,err)
          ADIOS_WRITE_LBL(buf_id,'eoff',eoffset,err)
          ADIOS_WRITE_LBL(buf_id,'enphase',ict2,err)
          ADIOS_WRITE_LBL(buf_id,'egid',egid,err)
          ADIOS_WRITE_LBL(buf_id,'ephase',ephase,err)
       endif
       ADIOS_CLOSE(buf_id,err)
       call write_unlock(filename)
       call t_stopf("COUPLINGP_WRITE")
       deallocate(phase,gid)
       if(sml_electron_on)   deallocate(ephase,egid)
       ADIOS_OPEN(buf_id,'coupling.info','xgc.coupling-info.bp','w',sml_comm,err)
       ADIOS_WRITE_LBL(buf_id,'sml_gstep',sml_gstep,err)
       ADIOS_WRITE_LBL(buf_id,'sml_run_count',sml_run_count,err)
       ADIOS_CLOSE(buf_id,err)
       call write_unlock('xgc.coupling-info.bp')
       return
     contains
       subroutine get_offset(num,offset,total)
         integer, intent(in) :: num
         integer (8), intent(out) :: offset, total
         integer*8 :: num8, inum_all(sml_intpl_totalpe) ! inum array per plane
         integer*8 :: inumsum, inumsum_all(sml_plane_totalpe) ! sum of inums in each plane
         ! two step mpi_allgather to avoid mpi_allgather in com_world
         num8 = num
         call mpi_allgather(num8,1,MPI_INTEGER8,inum_all,1,MPI_INTEGER8,sml_intpl_comm,err)
         inumsum = sum(inum_all)
         call mpi_allgather(inumsum,1,MPI_INTEGER8,inumsum_all,1,MPI_INTEGER8,sml_plane_comm,err)
         offset = sum(inumsum_all(1:sml_plane_mype)) + sum(inum_all(1:sml_intpl_mype))  !mype has zero base
         total = sum(inumsum_all)
       end subroutine get_offset
     end subroutine coupling_particle_write
     !------------------------------------------------------------------
     subroutine coupling_particle_read(spall)
       use sml_module
       use ptl_module
       use adios_read_mod
       use grid_class
       implicit none
       type(species_type) :: spall(0:ptl_nsp_max)
       integer, parameter :: ict1=ptl_nphase+1
       integer, parameter :: ict2=ptl_nphase+ptl_nconst
       integer*8 :: fh, sel=0
       integer*8 :: sel1, start1(1), count1(1)
       integer*8 :: sel2, start2(2), count2(2)
       integer :: method = 0
       integer :: i, j, i_beg, i_end, ierr
       integer :: istep_restart, istep0_restart
       integer :: inum, enum
       integer*8 :: ioff, eoff
       integer, parameter :: proc_coll=1
       integer :: inum_tmp(proc_coll)
       integer*8 :: ioff_tmp(proc_coll)
       real (8), allocatable :: phase(:,:)
       integer (8), allocatable :: gid(:)
       real (8), allocatable :: ephase(:,:)
       integer (8), allocatable :: egid(:)
       real (8) :: phasesum(ict2)
       integer (8) :: maxgid=0
       character (len=50) :: filename
       call check_point('coupling_particle_read')
       if (staging_read_method==0) then
          write(filename,'(A,"/","xgc.couplingp.xgc1.",i5.5,".bp")') trim(sml_coupling_xgc1_dir), sml_gstep+sml_gstep0
       else
          write(filename,'("xgc.couplingp.xgc1.",i5.5,".bp")') sml_gstep+sml_gstep0
       endif
       if (sml_mype == 0) print *, 'Loading particle information ... ', trim(filename)
       method = staging_read_method
       if (method==0) then
          call check_unlock(filename)
          call adios_read_open_file (fh, trim(filename), method, sml_comm, ierr)
       else
          call adios_read_open (fh, trim(filename), method, sml_comm, ADIOS_LOCKMODE_ALL, -1.0, ierr)
       endif
       call assert(ierr.eq.0, 'file open error', ierr)
       call adios_schedule_read (fh, sel, 'timestep0', 0, 1, sml_gstep0, ierr)
       call adios_schedule_read (fh, sel, 'timestep', 0, 1, istep_restart, ierr)
       call adios_schedule_read (fh, sel, 'run_count', 0, 1, sml_run_count, ierr)
       call adios_schedule_read (fh, sel, 'time', 0, 1, sml_time, ierr)
       call adios_perform_reads (fh, ierr)
       ! sml_gstep0 is the global time step index of the previous XGC1 run (read from restart file)
       !rh sml_gstep0 = istep_restart - istep0_restart !! time index at start of simulation
       call assert(istep_restart.eq.sml_gstep, 'sml_gstep is not consistent. Check timestep.dat', istep_restart, sml_gstep)
       sml_gstep = istep_restart
       sml_run_count=sml_run_count+1
       sml_nstep_xgc1 = sml_gstep-sml_gstep0+1
       ! Try to merge particles on a few cores to one.
       if(proc_coll > 1 .and. sml_mype==0) print *,'Note: merging particles on ',proc_coll,' cores. Hoping for the best.'
       start1(1) = proc_coll*sml_mype
       count1(1) = proc_coll
       call adios_selection_boundingbox(sel1, 1, start1, count1)
       call adios_schedule_read (fh, sel1, 'inum_arr', 0, 1, inum_tmp, ierr)
       call adios_schedule_read (fh, sel1, 'ioff_arr', 0, 1, ioff_tmp, ierr)
       call adios_perform_reads (fh, ierr)
       call adios_selection_delete (sel1)
       !   temporary!!!
       allocate(phase(ict2,maxval(inum_tmp)), gid(maxval(inum_tmp)))
       i_end=0
       do j=1, proc_coll
          i_beg=i_end+1
          inum=inum_tmp(j)
          i_end=i_end+inum
          ioff=ioff_tmp(j)
          start1(1) = ioff
          count1(1) = inum
          call adios_selection_boundingbox(sel1, 1, start1, count1)
          start2(1) = 0
          count2(1) = ict2
          start2(2) = ioff
          count2(2) = inum
          call adios_selection_boundingbox(sel2, 2, start2, count2)
          call adios_schedule_read (fh, sel1, 'igid', 0, 1, gid, ierr)
          call adios_schedule_read (fh, sel2, 'iphase', 0, 1, phase, ierr)
          call adios_perform_reads (fh, ierr)
          call adios_selection_delete (sel1)
          call adios_selection_delete (sel2)
          !! Check values. Copy and paste to Matlab. In Matlab, use adiosread function
          call mpi_allreduce(sum(phase, 2),phasesum,ict2,MPI_DOUBLE,MPI_SUM,sml_comm,ierr)
          if (sml_mype == 0) then
             do i=1, ict2
                print *, 'phasesum(', i, ') =', phasesum(i), ';'
             end do
          end if
          if(spall(1)%maxnum < i_end) print *,'Too many, ',i_end
          do i=1, inum
             spall(1)%ptl(i_beg-1+i)%ph=phase(1:ptl_nphase,i)
             spall(1)%ptl(i_beg-1+i)%ct=phase(ict1:ict2   ,i)
             spall(1)%ptl(i_beg-1+i)%gid=gid(i)
          enddo
          maxgid=maxval(gid)
       enddo
       spall(1)%num = i_end
       spall(1)%maxgid = maxgid
       call mpi_allreduce(MPI_IN_PLACE,spall(1)%maxgid,1,MPI_INTEGER8,MPI_MAX,sml_comm,ierr)
       deallocate(phase,gid)
       if(sml_electron_on) then
          start1(1) = sml_mype
          count1(1) = 1
          call adios_selection_boundingbox(sel1, 1, start1, count1)
          call adios_schedule_read (fh, sel1, 'enum_arr', 0, 1, enum, ierr)
          call adios_schedule_read (fh, sel1, 'eoff_arr', 0, 1, eoff, ierr)
          call adios_perform_reads (fh, ierr)
          call adios_selection_delete (sel1)
          allocate(ephase(ict2,enum), egid(enum))
          start1(1) = eoff
          count1(1) = enum
          call adios_selection_boundingbox(sel1, 1, start1, count1)
          start2(1) = 0
          count2(1) = ict2
          start2(2) = eoff
          count2(2) = enum
          call adios_selection_boundingbox(sel2, 2, start2, count2)
          call adios_schedule_read (fh, sel1, 'egid', 0, 1, egid, ierr)
          call adios_schedule_read (fh, sel2, 'ephase', 0, 1, ephase, ierr)
          call adios_perform_reads (fh, ierr)
          call adios_selection_delete (sel1)
          call adios_selection_delete (sel2)
          do i=1, enum
             spall(0)%ptl(i)%ph=ephase(1:ptl_nphase,i)
             spall(0)%ptl(i)%ct=ephase(ict1:ict2   ,i)
             spall(0)%ptl(i)%gid=egid(i)
          enddo
          spall(0)%num = enum
          spall(0)%maxgid = maxval(egid)
          call mpi_allreduce(MPI_IN_PLACE,spall(0)%maxgid,1,MPI_INTEGER8,MPI_MAX,sml_comm,ierr)
          deallocate(ephase,egid)
       endif
       call adios_read_close (fh, ierr)
     end subroutine coupling_particle_read
     !------------------------------------------------------------------
#ifndef XGCA
     subroutine coupling_turb_write(istep,irk,grid,psn)
       use sml_module
       use diag_module
       use grid_class
       use psn_class
       implicit none
       type(grid_type) :: grid
       type(psn_type) :: psn
       integer   :: istep,irk
       integer   :: nnode,nphi
       integer :: i,iphi
       integer*8, save :: buf_id
       integer*8 :: buf_size, total_size
       integer :: err
       real (kind=8), save, allocatable :: dpot(:,:), pot0(:,:)
       character (len=30) :: filename
       if(.not.allocated(dpot)) allocate(dpot(grid%nnode,1:2),pot0(grid%nnode,1:2))
       nnode=grid%nnode
       nphi=grid%nphi
       if (sml_plane_mype==0) then
          if (sml_iter_solver) then
             dpot(:,irk) = psn%dpot(:,0) -(psn%pot0m-psn%pot0)
             pot0(:,irk) = psn%pot0m
             !  Testing only: save the whole shebang
             !        dpot = psn%dpot(:,1)+psn%pot0
             buf_size = 4*4 + nnode*sml_nrk*8 + 1000
             if(sml_mype==0) buf_size = buf_size + nnode*sml_nrk*8
             !jc build dpot step by step and write in the last step
             if (irk == sml_nrk) then
                write(filename,'("xgc.couplingt",".",i5.5,".bp")') sml_gstep
                ADIOS_OPEN(buf_id,'couplingt',filename,'w',sml_intpl_comm,err)
                ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
                ADIOS_WRITE_LBL(buf_id,'nnode',nnode,err)
                ADIOS_WRITE_LBL(buf_id,'nphi',sml_nphi_total,err)
                ADIOS_WRITE_LBL(buf_id,'iphi',sml_intpl_mype,err)
                ADIOS_WRITE_LBL(buf_id,'nrk',sml_nrk,err)
                ADIOS_WRITE_LBL(buf_id,'dpot',dpot,err)
                if(sml_mype==0) ADIOS_WRITE_LBL(buf_id,'pot0',pot0,err)
                ADIOS_CLOSE(buf_id,err)
                call write_unlock(filename)
             endif
          endif
       end if
       ! deallocate on the last step for good practice
       if(sml_mstep==istep.and.irk==2) deallocate(dpot,pot0)
     end subroutine coupling_turb_write
#endif
     !------------------------------------------------------------------
     subroutine coupling_turb_read(grid,psn,istep)
       use grid_class
       use sml_module
       use psn_class
       use adios_read_mod
       implicit none
       type(grid_type),intent(in) :: grid
       type(psn_type), intent(inout) :: psn
       integer, intent(in) :: istep
       integer :: nnode=0, nphi=0, nrk
       integer, parameter :: numdim=3
       integer*8 :: fh
       integer*8 :: sel=0, start(numdim), count(numdim)
       integer :: method=0
       integer :: ierr
       real (8),allocatable :: tmp(:,:,:), tmp2(:,:,:)
       integer :: vartype, nsteps, ndim
       integer*8 :: dims(numdim)
       character (len=256) :: filename
       integer :: fnum, iphi, ipc
       integer :: iphi_first, iphi_last, i, n_pass, is(2), ie(2), is2, ie2
       logical, save :: first = .true.
       integer, save :: l_fnum, l_nnode, l_nphi, l_nrk
       integer, save :: l_vartype, l_nsteps, l_ndim
       integer*8, save :: l_dims(numdim)
       real (8), allocatable, save :: l_turb_pot0(:,:)
       real (8), allocatable, save :: l_tmp(:,:,:)
       call check_point('coupling_turb_read')
       if(sml_plane_mype == 0) then
          !allocate(tmp(sml_nphi_total,grid%nnode))
          allocate(tmp(-1:2,grid%nnode,sml_nrk))
          tmp = 0D0
          if (sml_coupling_firstfile .gt. 0 .and. sml_coupling_lastfile .gt. 0  &
               .and. sml_coupling_firstfile .le. sml_coupling_lastfile) then
             fnum = sml_coupling_firstfile + mod(istep,sml_coupling_lastfile-sml_coupling_firstfile+1)
          else
             fnum = sml_gstep0 + mod(istep,sml_nstep_xgc1) + sml_coupling_nskip  !! sml_gstep0...sml_gstep0 + (sml_nstep_xgc1-1), periodic
          endif
          if (.not. first .and. l_fnum==fnum) then
             print *, 'Use saved values'
             fnum = l_fnum
             nnode = l_nnode
             nphi = l_nphi
             nrk = l_nrk
             vartype = l_vartype
             nsteps = l_nsteps
             ndim = l_ndim
             dims = l_dims
             turb_pot0 = l_turb_pot0
             tmp = l_tmp
          else
             if (staging_read_method==0) then
                write(filename,'(A,"/","xgc.couplingt",".",i5.5,".bp")') trim(sml_coupling_xgc1_dir), fnum
             else
                write(filename,'("xgc.couplingt",".",i5.5,".bp")') fnum
             endif
             if (sml_mype == 0) print *, 'Loading turbulence information ... ', trim(filename)
             method = staging_read_method
             if (method==0) then
                call check_unlock(filename)
                call adios_read_open_file (fh, trim(filename), method, sml_intpl_comm, ierr)
             else
                call adios_read_open (fh, trim(filename), method, sml_intpl_comm, ADIOS_LOCKMODE_ALL, -1.0, ierr)
             endif
             !call assert(ierr.eq.0, 'file open error', ierr, trim(filename))
             call adios_schedule_read (fh, sel, 'nnode', 0, 1, nnode, ierr)
             call adios_schedule_read (fh, sel, 'nphi', 0, 1, nphi, ierr)
             call adios_schedule_read (fh, sel, 'nrk', 0, 1, nrk, ierr)
             call adios_schedule_read (fh, sel, 'pot0', 0, 1, turb_pot0, ierr)
             call adios_perform_reads (fh, ierr)
             call adios_inq_var (fh, 'dpot', vartype, nsteps, ndim, dims, ierr)
             call assert(nnode == grid%nnode, 'the number of node mismatch', nnode, grid%nnode)
             call assert(nphi == sml_nphi_total, 'total nphi mismatch', nphi, sml_nphi_total)
             call assert(nrk == sml_nrk, 'total nrk mismatch', nrk, sml_nrk)
             ! I assume that plane indices in the ADIOS file run from 0 to sml_nphi_total-1
             ! 1) Determine indices of first and last plane to be read
             if (sml_plane_index==0) then
                iphi_first=sml_nphi_total-1
             else
                iphi_first=sml_plane_index-1
             endif
             iphi_last=mod(sml_plane_index+2,sml_nphi_total)
             ! 2) Find out whether we have to read in one or two steps
             !    and determine first and last plane to be read in each step
             if (iphi_first .lt. iphi_last) then
                n_pass=1
                is(1)=iphi_first
                ie(1)=iphi_last
             else
                n_pass=2
                is(1)=iphi_first
                ie(1)=sml_nphi_total-1
                is(2)=0
                ie(2)=iphi_last
             endif
             start(2) = 0
             count(2) = nnode
             start(3) = 0
             count(3) = sml_nrk
             ! 3) Read the data in one or two steps (if ADIOS knows what it is supposed
             !    to do if plane indices are <0 or >nphi-1, we do not need the loop)
             do i=1,n_pass
                ! 3a) Determine start(1), count(1) and the indices of tmp to which
                !     the data will be written
                start(1) = is(i) !0
                count(1) = ie(i)-is(i)+1 !nphi
                is2=-1 + (i-1) * (ie(1)-is(1)+1)
                ie2=-1 + (ie(1)-is(1)) + (i-1) * (ie(2)-is(2)+1)
                ! 3b) Read the data and store into tmp
                !allocate(tmp2(sml_nrk, count(1), count(2), sml_nrk))
                allocate(tmp2(count(1), count(2), count(3)))
                call adios_selection_boundingbox (sel, numdim, start, count)
                ! This is a delay read. Adios don't support auto created pointer.
                !call adios_schedule_read (fh, sel, 'dpot', 0, 1, tmp(is2:ie2,:), ierr)
                call adios_schedule_read (fh, sel, 'dpot', 0, 1, tmp2, ierr)
                call adios_perform_reads (fh, ierr)
                call adios_selection_delete (sel)
                call assert(ie2-is2+1.eq.count(1), 'index error', ierr)
                tmp(is2:ie2,:,:) = tmp2
                deallocate(tmp2)
             enddo
             call adios_read_close (fh, ierr)
          endif
          if (first) then
             l_fnum = fnum
             l_nnode = nnode
             l_nphi = nphi
             l_nrk = nrk
             l_vartype = vartype
             l_nsteps = nsteps
             l_ndim = ndim
             l_dims = dims
             allocate(l_turb_pot0(nnode,2))
             l_turb_pot0 = turb_pot0
             allocate(l_tmp(-1:2,grid%nnode,sml_nrk))
             l_tmp = tmp
             first = .false.
          endif
          ! 4) Copy turbulence data from tmp to psn%turb_pot
          do iphi=-1,2 !1,sml_nphi_total
             !size(turb_pot)=(nnode,4,2)
             do nrk=1,sml_nrk
                turb_pot(:,iphi,nrk)=tmp(iphi,:,nrk)
             enddo
          enddo
       endif ! sml_plane_mype==0
       if(sml_mype==0) print *,'Broadcasting turbulence data.'
       call MPI_bcast(turb_pot,4*grid%nnode*sml_nrk,MPI_REAL8,0,sml_plane_comm,ierr)
       call MPI_bcast(turb_pot0,grid%nnode*sml_nrk,MPI_REAL8,0,sml_plane_comm,ierr)
       call check_point('end of diag_turb_read')
       if(sml_plane_mype==0) deallocate(tmp)
     end subroutine coupling_turb_read
     !------------------------------------------------------------------
     subroutine write_unlock(fname)
       use sml_module
       implicit none
       character(len=*), intent(in) :: fname
       character(512) :: lname
       character(5) :: mypestr
       write(mypestr,'(I0.5)') sml_mype
       lname=trim(fname)//'_'//trim(mypestr)//'.unlock'
       print *,'Write unlock filename ',trim(lname)
       open(20, file=lname, status="new", action="write")
       close(20)
     end subroutine write_unlock
     subroutine check_unlock(fname)
       use sml_module
       implicit none
       character(len=*), intent(in) :: fname
       character(512) :: lname
       character(5) :: mypestr
       logical :: ex
       write(mypestr,'(I0.5)') sml_mype
       lname=trim(fname)//'_'//trim(mypestr)//'.unlock'
       print *,'Wait unlock filename ',trim(lname)
       ex=.false.
       do while(.NOT.ex)
          call sleep(1)
          inquire(file=lname,EXIST=ex)
       end do
     end subroutine check_unlock
   end module
