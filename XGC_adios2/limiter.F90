!***********************************************
! Limiter check routine 
! first created 2003/11/03
! last modifed 2003/03/07
! bug fix - 2003/03/27
! bug fix - 2003/11/03
!************************************************

subroutine limiter_setup
  use lim_module
  use eq_module, only : is_rgn12,eq_x_psi,eq_axis_r,eq_axis_z, eq_x_z
  use sml_module, only : sml_2pi,&
       sml_bd_min_r, sml_bd_max_r, sml_bd_min_z, sml_bd_max_z, sml_mype
  implicit none
  real (kind=8) :: zmin, zmax,r,z,psi,psi_interpol,dr,old_r
  real (kind=8), allocatable :: tmp_r(:,:), tmp_z(:,:)
  integer :: i,j,izmin,izmax,offset,zindex(2),rl,size,index
  integer, parameter :: unit=92

integer, parameter :: outpe=0
  if(sml_mype==outpe) print *, 'enter limiter_setup'
!  if(lim_filename/='') then
     if(sml_mype==outpe) print *, 'limiter read call'
     call limiter_read
     call check_point('limiter read finished:')
!  endif
  allocate(tmp_r(lim_mdata,2), tmp_z(lim_mdata,2))
  tmp_r=eq_axis_r
  tmp_z=eq_axis_z
  zmin=0D0
  zmax=0D0

  do  i=1, lim_mdata
     if(zmin > lim_org_z(i)) then
        zmin=lim_org_z(i)
        izmin=i
     elseif(zmax < lim_org_z(i)) then
        zmax=lim_org_z(i)
        izmax=i
     endif
  enddo

  lim_zmin=zmin
  lim_zmax=zmax
  lim_dz=(lim_zmax-lim_zmin)/real(lim_store_mz)
  lim_r0_up=lim_org_r(izmax)
  lim_r0_down=lim_org_r(izmin)

  if( izmax<izmin ) then
     offset=lim_mdata
  else
     offset=0
  endif
 if(sml_mype==outpe) print *,'izmax==', izmax,izmin,lim_zmin,lim_zmax
! stop
!------- One direction ---------------------------
  rl=1
  
  tmp_z(1,rl)=lim_org_z(izmin)
  tmp_r(1,rl)=lim_org_r(izmin)

  j=2
  do i= izmin+1, izmax+offset
     index=mod(i-1,lim_mdata)+1
     if(lim_org_z(index) > tmp_z(j-1,rl) ) then
        tmp_z(j,rl)=lim_org_z(index)
        tmp_r(j,rl)=lim_org_r(index)
        j=j+1
     else if ( lim_org_z(index) == tmp_z(j-1,rl) ) then
        tmp_z(j,rl)=lim_org_z(index) + 1E-5
        tmp_r(j,rl)=lim_org_r(index)
        j=j+1
     endif
  enddo
  zindex(rl)=j-1

!------- the other direction ---------------------------
  rl=2
  
  tmp_z(1,rl)=lim_org_z(izmin)
  tmp_r(1,rl)=lim_org_r(izmin)

  j=2
  do i= izmin-1, izmax+offset-lim_mdata , -1
     index=mod(i+lim_mdata-1,lim_mdata)+1
     if(lim_org_z(index) > tmp_z(j-1,rl) ) then
        tmp_z(j,rl)=lim_org_z(index)
        tmp_r(j,rl)=lim_org_r(index)
        j=j+1
     else if ( lim_org_z(index) == tmp_z(j-1,rl) ) then
        tmp_z(j,rl)=lim_org_z(index) + 1E-5
        tmp_r(j,rl)=lim_org_r(index)
        j=j+1
     endif
  enddo
  zindex(rl)=j-1

!----------- setting up ---------------------------------
  size=maxval(zindex)
  allocate( lim_r(size,2), lim_z(size,2) )
  allocate( lim_weight(lim_store_mz,2), lim_en(lim_store_mz,2),lim_ds(lim_store_mz,2))
  lim_weight=0D0
  lim_en=0D0

  if( maxval(tmp_r(:,1)) > maxval(tmp_r(:,2)) ) then
     lim_r=tmp_r(1:size,:)
     lim_z=tmp_z(1:size,:)
     lim_zindex=zindex
  else
     lim_r(:,1)=tmp_r(1:size,2)
     lim_r(:,2)=tmp_r(1:size,1)
     lim_z(:,1)=tmp_z(1:size,2)
     lim_z(:,2)=tmp_z(1:size,1)
     lim_zindex(1)=zindex(2)
     lim_zindex(2)=zindex(1)
  endif
  deallocate(tmp_r,tmp_z)
  ! finding smallest psi value -----------------------------
  lim_psi_min=1D30
  do rl=1, 2
     do i=1,lim_zindex(rl)
        r=lim_r(i,rl)
        z=lim_z(i,rl)
        !if inside of interpolation range
        if(r > sml_bd_min_r .and. r < sml_bd_max_r &
             .and. z > sml_bd_min_z .and. z < sml_bd_max_z ) then
           psi=psi_interpol(r,z,0,0)
!           if(z>eq_x_z .AND. psi < lim_psi_min) then
           if(is_rgn12(r,z,psi) .and. z>eq_x_z .and. psi < lim_psi_min) then
              lim_psi_min=psi
           endif
        endif
     enddo
  enddo

  !writing out the limiter data for plot and debug------------
  if(sml_mype==0) then
  write(unit, *)'# limiter setup output'
  write(unit, *)'# num =', lim_zindex
  write(unit, *)'# psi_min = ', lim_psi_min/eq_x_psi
  write(unit, *)'# (r,z) values of limiter position'
  do rl=1,2
     do i=1, lim_zindex(rl)
        write(unit, *) lim_r(i,rl),lim_z(i,rl)
     enddo
     write(unit,*) ' '
  enddo
  endif

  !store the ds values

  do rl=1, 2
     j=1 
     old_r=lim_r0_down
     do i=1, lim_store_mz
        z=lim_zmin + real(i)*lim_dz
        z=min(lim_zmax,z)
        ! find proper j value
        do while( z<lim_z(j,rl) .OR. lim_z(j+1,rl) <z ) 
           j=j+1
        enddo
        r=  (lim_r(j+1,rl)-lim_r(j,rl)) * (z-lim_z(j,rl)) / (lim_z(j+1,rl)-lim_z(j,rl))&
             + lim_r(j,rl)
        dr=r-old_r
        lim_ds(i,rl)=sqrt(lim_dz**2+dr**2)*sml_2pi*( 0.5*(r+old_r))
        old_r=r

     enddo
  enddo

end subroutine limiter_setup

subroutine limiter_read
  use lim_module
  use neu_module, only : neu_sep_mtheta_file,neu_sep_r_file, neu_sep_z_file, neu_limfile, neu_sepfile
  use sml_module, only : sml_mype
  implicit none
  integer, parameter :: fileunit=114
  integer, parameter :: fileunit_sep=115
  integer :: i
  character :: char(10)

  open(unit=fileunit,file=neu_limfile,status='old',action='read')
!  read(fileunit,100) char
  read(fileunit,*) lim_mdata
  allocate(lim_org_r(lim_mdata),lim_org_z(lim_mdata))
  call check_point('start reading limiter file')
  if(sml_mype==0) print *, 'lim_mdata=', lim_mdata
  do i=1, lim_mdata
     !read(fileunit,102) lim_org_r(i),lim_org_z(i)
     read(fileunit,*) lim_org_r(i),lim_org_z(i)
     !if(sml_mype==0)print *, 'lim_org=', lim_org_r(i), lim_org_z(i)
  enddo
  call check_point('after reading lim.dat')
  close(fileunit)

  open(unit=fileunit_sep,file=neu_sepfile,status='old',action='read')
  read(fileunit_sep,*) neu_sep_mtheta_file
  allocate(neu_sep_r_file(neu_sep_mtheta_file),neu_sep_z_file(neu_sep_mtheta_file))
  if(sml_mype==0) print *, 'neu_sep_mtheta_file=', neu_sep_mtheta_file
  do i=1, neu_sep_mtheta_file
     !read(fileunit_sep,102) neu_sep_r_file(i),neu_sep_z_file(i)
     read(fileunit_sep,*) neu_sep_r_file(i),neu_sep_z_file(i)
!     print *, 'neu_sep=', neu_sep_r_file(i), neu_sep_z_file(i)
  enddo
  close(fileunit_sep)
!  call check_point('after reading sep.dat')
100 format (a2)
101 format (i8)
102 format (e19.13,1x,e19.13)
end subroutine

subroutine calbdpoints(itheta,xlim,ylim,limitr,xb,yb)
  use sml_module
  use neu_module
  use eq_module
  implicit none
  integer, intent(in) :: itheta, limitr
  real (kind=8), intent(in)  :: xlim(limitr+1),ylim(limitr+1)
  real (kind=8), intent(out) :: xb(itheta),yb(itheta)
  real (kind=8) :: wb(itheta)
  real (kind=8) :: wa(limitr),wc(limitr),wd(limitr),we(limitr),wf(limitr),wg(limitr)
  real (kind=8) :: dp1,dp2,dtheta
  integer, external :: modn1
  integer :: i,j,i0,jj

 if(sml_mype==0) print *, 'ylim=', ylim(1),xlim(1),eq_axis_r,eq_axis_z
  dp1=atan2((ylim(1)-eq_axis_z),(xlim(1)-eq_axis_r)); if(dp1<0D0) dp1=dp1+sml_2pi
  dp2=atan2((ylim(2)-eq_axis_z),(xlim(2)-eq_axis_r)); if(dp2<0D0) dp2=dp2+sml_2pi

! reverse order
  if(dp2<dp1) then
     do i = 1, limitr
        j = limitr - i + 1
        wc(j) = xlim(i)
        wd(j) = ylim(i)
     enddo
  else
     do i = 1, limitr
        wc(i) = xlim(i)
        wd(i) = ylim(i)
     enddo
  endif

  do j = 1, limitr
     wa(j)=atan2((wd(j)-eq_axis_z),(wc(j)-eq_axis_r)); if(wa(j)<0D0) wa(j)=wa(j)+sml_2pi
  enddo

  i0 = 0
  do j = 1, limitr-1
     if((wa(j+1)-wa(j)) < -sml_pi) i0 = j
  enddo

  do j = 1, limitr
     jj = modn1(j + i0,limitr)
     we(j) = wc(jj)
     wf(j) = wd(jj)
     wg(j) = wa(jj)
  enddo
  do j = 1, limitr
     wc(j) = we(j)
     wd(j) = wf(j)
     wa(j) = wg(j)
  enddo

! interpolate to new boundary points with uniform distribution in poloidal angle
  dtheta = sml_2pi/real(itheta)
  do j = 1, itheta
     wb(j) = (j - 1D0)*dtheta
  enddo

  call interpsl(wb,xb,yb,itheta,wa,wc,wd,limitr)

end subroutine calbdpoints

integer function modn1(n,l)
  implicit none
  integer, intent(in) :: n, l
  integer m

  m=n
  do while ((m>l).or.(m<1))
     if(m>l) then
        m=m-l
     elseif(m<1) then
        m=m+l
     endif
  enddo
  modn1=m
end function modn1

subroutine interpsl(rr,pr,pz,num0,r1,p1,pz1,num)
  use sml_module, only: sml_2pi,sml_pi
  use eq_module,  only:eq_axis_r,eq_axis_z
  implicit none
  real (kind=8), intent(in)  :: rr(num0),r1(num),p1(num),pz1(num)
  real (kind=8), intent(out) :: pr(num0),pz(num0)
  real (kind=8) :: rri,delr,denom,r1j,r1jp
  integer, external :: modn1
  integer, intent(in) :: num0,num
  integer :: i,j,jp

  do i=1,num0
     rri=rr(i)
     do j=1,num
        jp=modn1(j+1,num)
        r1j=r1(j)
        r1jp=r1(jp)
        if(r1jp<(r1j-sml_pi)) then
           r1j=r1j-sml_2pi
           if(rri>sml_pi) rri=rri-sml_2pi
        endif
        denom=r1jp-r1j
        delr=2D0
        if(denom/=0D0) delr=(rri-r1j)/denom
        if((delr>=0D0).and.(delr<1D0)) then
           pr(i)=((pz1(j)-eq_axis_z)*(p1(jp)-eq_axis_r)-(p1(j)-eq_axis_r)*(pz1(jp)-eq_axis_z))/ &
                 ((p1(jp)-p1(j))*dtan(rri)-(pz1(jp)-pz1(j)))
           pz(i)=pr(i)*dtan(rri)+eq_axis_z
           pr(i)=pr(i)+eq_axis_r
           goto 10
        endif
     enddo
     10 continue
  enddo
end subroutine
