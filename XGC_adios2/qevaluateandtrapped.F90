subroutine q_evaluation(grid,R,Z,nw,nh,npsi,psi_in)

   use EQ_module
   use grid_class
   use sml_module, only : sml_mype, sml_pi, sml_2pi
   implicit none

   ! We want to have the trapped particle fraction, qsafety, minor/major radius and aspect ratio
   ! in the grid file xgc.mesh.bp -->
   type(grid_type), intent(inout) :: grid

   ! nw and nh are the dimensions of the equilibrium data from the eqd file
   integer, intent(in) :: nw,nh, npsi
   real(kind=8) ,intent(in) :: R(nw),Z(nh), psi_in(npsi)
   real(kind=8) :: Bpol(nw,nh),BR(nw,nh),BZ(nw,nh)
   real(kind=8) :: Btot(nw,nh)
   real(kind=8) :: Bphi(nw,nh)
   real(kind=8) , external :: I_interpol,psi_interpol

   ! These parameters define the size of the working grid in (Psi,theta) -->
   integer, parameter :: mtheta=512,gridm=200
   real(kind=8) :: r_con(mtheta,gridm),z_con(mtheta,gridm),psi_con(mtheta,gridm)
   real(kind=8) :: r_tmp(gridm),z_tmp(gridm),psi_tmp(gridm)
   real(kind=8) :: distUP,distDOWN,distRIGHT,distLEFT,length1,length2,&
                   length3,length4
   real(kind=8) :: lcheck,rcheck,zcheck
   real(kind=8) :: delr,delz,distanceR,distanceZ
   real(kind=8) :: r_psitmp(mtheta),z_psitmp(mtheta),psipsi
   real(kind=8) :: dtheta,theta,dpsi
   real(kind=8) ::  mu0
   real(kind=8) :: dx,dy,R_R,dl

   ! In original version mpsi=eq_mpsi
   integer, parameter :: mpsi=129
   real(kind=8) :: r_psi(mpsi,mtheta),z_psi(mpsi,mtheta)
   real(kind=8) :: psigrid(mpsi)
   real(kind=8) :: Bp_psi(mpsi,mtheta),Bp0
   real(kind=8) :: Btot_psi(mpsi,mtheta),Btot2_psi(mpsi,mtheta)
   real(kind=8) :: Bphi_psi(mpsi,mtheta)
   real(kind=8) :: Bmax(mpsi),psi,Bmax_tmp
   real(kind=8) :: Btot2_fluxavg(mpsi), dbdlpar2_fluxavg(mpsi), norm(mpsi)
   real(kind=8) :: lamda,del_lamda
   real(kind=8), dimension(mpsi,mtheta) :: dbdlpar, dbdlpar2
   integer , parameter :: n_lamda=50
   real(kind=8) :: func_psi(mpsi,mtheta)
   real(kind=8) :: Func_fluxavgtmp(mpsi)
   real(kind=8) :: Func_fluxavg(n_lamda,mpsi)
   real(kind=8) :: trapped(mpsi), qpsi(mpsi)
   real(kind=8) :: trapped_FM(mpsi),sqrtoneb(mpsi,mtheta),fluxavg(mpsi)
   real(kind=8) :: trapped2(mpsi),Rmax,Rmin,epspar(mpsi,3)
   integer :: i, j, idx, ip, jp, jm
   real (kind=8) :: wp(2), dl1, dl2

   mu0=4D0*sml_pi*1D-7

   allocate(grid%epspar(grid%npsi00,3),grid%qsafety(grid%npsi00),grid%trapped(grid%npsi00))
   !rh For mode filtering in em code
   allocate(grid%qsafety2(npsi))


!######################################################################
!##                      evaluate Bpol field                         ##
!######################################################################
!   allocate(eq_Bmax(eq_mpsi))

!   print *,'nw,nh,eq_mr,eq_mz:',nw,nh,eq_mr,eq_mz

   do j=1,eq_mr !eq_mz
      do i=1,eq_mz !eq_mr
         call bvec_interpol(r(j),z(i),0.D0,BR(j,i),BZ(j,i),Bphi(j,i))

!         if(i .eq. 1) then
!            BR(j,i)=(EQ_psi_rz(j,i+1)-EQ_psi_rz(j,i))/(Z(i+1)-Z(i))
!            BR(j,i)=-1D0*BR(j,i)/R(j)
!         elseif (i .eq. nh) then
!            BR(j,i)=(EQ_psi_rz(j,i)-EQ_psi_rz(j,i-1))/(Z(i)-Z(i-1))
!            BR(j,i)=-1D0*BR(j,i)/R(j)
!         else
!            BR(j,i)=(EQ_psi_rz(j,i+1)-EQ_psi_rz(j,i-1))/(Z(i+1)-Z(i-1))
!            BR(j,i)=-1D0*BR(j,i)/R(j)
!         endif
!       enddo
!   enddo
!
!   do j=1, eq_mz
!      do i=1, eq_mr
!         if ( i .eq. 1) then
!            BZ(i,j)=(EQ_psi_rz(i+1,j)-EQ_psi_rz(i,j))/(R(i+1)-R(i))
!            BZ(i,j)=BZ(i,j)/R(i)
!        elseif( i .eq. nw) then
!            BZ(i,j)=(EQ_psi_rz(i,j)-EQ_psi_rz(i-1,j))/(R(i)-R(i-1))
!            BZ(i,j)=BZ(i,j)/R(i)
!         else
!            BZ(i,j)=(EQ_psi_rz(i+1,j)-EQ_psi_rz(i-1,j))/(R(i+1)-R(i-1))
!            BZ(i,j)=BZ(i,j)/R(i)
!         endif
!      enddo
!   enddo
!   do j=1, eq_mz
!      do i=1, eq_mr
!         psi=psi_interpol(R(i),Z(j),0,0)
!         if(psi<eq_x_psi .and. z(j)<eq_x_z) then
!           Bphi(i,j)=I_interpol(psi,0,3)/R(i)
!         else
!           Bphi(i,j)=I_interpol(psi,0,1)/R(i)
!         endif
      enddo
   enddo
   Bpol(:,:)=sqrt(BR(:,:)**2 + BZ(:,:)**2)
!   if (sml_s_alpha .and. sml_s_alpha_bapprox) then
     ! NEO assumes B~Bphi in s-alpha geometry
!     Btot(:,:)=dabs(Bphi(:,:))
!   else
     Btot(:,:)=sqrt(BR(:,:)**2 + BZ(:,:)**2 + Bphi(:,:)**2)
!   endif

!###########################################################################
!##                 calculate R and Z on constant Psi                     ##
!###########################################################################
     distUP=abs(z(eq_mz)-EQ_axis_z)
     distDOWN=abs(z(1)-EQ_axis_z)
     distRIGHT=abs(r(eq_mr)-EQ_axis_r)
     distLEFT=abs(r(1)-EQ_axis_r)
     dtheta=sml_2pi/real(mtheta)
     length1=sqrt((r(eq_mr)-EQ_axis_r)**2+(z(eq_mz)-EQ_axis_z)**2)
     length2=sqrt((r(1)-EQ_axis_r)**2+(z(eq_mz)-EQ_axis_z)**2)
     length3=sqrt((r(1)-EQ_axis_r)**2+(z(1)-EQ_axis_z)**2)
     length4=sqrt((r(eq_mr)-EQ_axis_r)**2+(z(1)-EQ_axis_z)**2)

     do i=1,mtheta
        theta=dtheta*real(i-1)
        if ((theta .ge. 0d0) .and.(theta .lt. sml_pi/2d0)) then
           lcheck=abs(distRIGHT/cos(theta))
           if (lcheck .le. length1) then
               rcheck=r(nw)
               zcheck=EQ_axis_z+distRIGHT*tan(theta)
           elseif(lcheck .gt. length1) then
               rcheck=EQ_axis_r+distUP*tan(sml_pi/2d0-theta)
               zcheck=z(eq_mz)
           endif
        elseif((theta .ge. sml_pi/2D0) .and. (theta .lt. sml_pi)) then
           lcheck=abs(distLEFT/cos(sml_pi-theta))
           if (lcheck .le. length2) then
              rcheck=r(1)
              zcheck=EQ_axis_z+abs(distLEFT*tan(sml_pi-theta))
           elseif (lcheck .gt. length2) then
              rcheck=EQ_axis_r-abs(distUP*tan(theta-sml_pi/2d0))
              zcheck=z(eq_mz)
           endif
        elseif((theta .ge. sml_pi) .and. (theta .lt.3D0*sml_pi/2D0)) then
           lcheck=abs(distLEFT/cos(theta-sml_pi))
           if (lcheck .le. length3) then
              rcheck=r(1)
              zcheck=EQ_axis_z-abs(distLEFT*tan(theta-sml_pi))
           elseif(lcheck .gt. length3) then
              rcheck=EQ_axis_r-abs(distDOWN*tan(3D0*sml_pi/2D0-theta))
              zcheck=z(1)
           endif
        elseif((theta .ge. 3D0*sml_pi/2D0) .and. (theta .lt. sml_2pi)) then
           lcheck=abs(distRIGHT/cos(sml_2pi-theta))
           if (lcheck .le. length4) then
              rcheck=r(eq_mr)
              zcheck=EQ_axis_z-abs(distRIGHT*tan(sml_2pi-theta))
           elseif(lcheck .gt. length4) then
              rcheck=EQ_axis_r+abs(distDOWN*tan(theta-3D0*sml_pi/2D0))
              zcheck=z(1)
           endif
        endif
!        write(888,100) rcheck,zcheck

        distanceR=rcheck-EQ_axis_r
        distanceZ=zcheck-EQ_axis_z
        delr=distanceR/real(gridm)
        delz=distanceZ/real(gridm)
!        print *,rcheck,zcheck

        do j=1,gridm
           r_con(i,j)=EQ_axis_r+delr*real(j-1)
           z_con(i,j)=EQ_axis_z+delz*real(j-1)
!           print *,r_con(i,j),z_con(i,j)
           r_tmp(j)=r_con(i,j)
           z_tmp(j)=z_con(i,j)
        enddo
        call findpsi(r_tmp,z_tmp,gridm,r,z,eq_mr,eq_mz,EQ_psi_rz,psi_tmp)
        psi_con(i,:)=psi_tmp(:)
    enddo

    dpsi=(EQ_psi_grid(eq_mpsi)-EQ_psi_grid(1))/real(mpsi-1)
    psigrid(1)=EQ_psi_grid(1)
    do i=2, mpsi
       psigrid(i)=psigrid(1)+dpsi*real(i-1,8)
       psipsi=EQ_psi_grid(1)+dpsi*real(i-1,8)
       call constpsi_contour(psipsi,r_con,z_con,psi_con,mtheta,gridm,r_psitmp,z_psitmp)
       r_psi(i,:)=r_psitmp(:)
       z_psi(i,:)=z_psitmp(:)
    enddo
!rh ifort recommends 20.13 instead of 19.13
100 format(2(e20.13," "))
!############################################################################
!##      find maximun B field on flux surface                              ##
!############################################################################


!############################################################################

!###########################################################################
!##          calculate Function values on constant Psi                    ##
!###########################################################################
     call findFonconstPsi(R,Z,Bpol,eq_mr,eq_mz,r_psi,z_psi,Bp_psi,mpsi,mtheta)
     call findFonconstPsi(R,Z,Bphi,eq_mr,eq_mz,r_psi,z_psi,Bphi_psi,mpsi,mtheta)
     call findFonconstPsi(R,Z,Btot,eq_mr,eq_mz,r_psi,z_psi,Btot_psi,mpsi,mtheta)

     do i=2,mpsi
       do j=1,mtheta
         jp=mod(j,mtheta)+1
         jm=j-1
         if (jm .eq. 0) jm=mtheta
         dl1=sqrt((r_psi(i,jp)-r_psi(i,j))**2+(z_psi(i,jp)-z_psi(i,j))**2)
         dl2=sqrt((r_psi(i,j)-r_psi(i,jm))**2+(z_psi(i,j)-z_psi(i,jm))**2)
         dbdlpar(i,j)=Btot_psi(i,j)*(dl1/dl2-dl2/dl1)+Btot_psi(i,jp)*dl2/dl1-Btot_psi(i,jm)*dl1/dl2
         dbdlpar(i,j)=dbdlpar(i,j)/(dl1+dl2)*(Bp_psi(i,j)/Btot_psi(i,j))
         dbdlpar2(i,j)=dbdlpar(i,j)**2
       enddo
     enddo
     dbdlpar(1,:)=0D0
     dbdlpar2(1,:)=0D0

     dbdlpar2_fluxavg(:)=0D0
     norm(:)=0D0
     do i=2,mpsi
        qpsi(i)=0D0
        do j=2,mtheta-1
           bp0=(Bp_psi(i,j-1)+Bp_psi(i,j))/2D0
           dx=r_psi(i,j)-r_psi(i,j-1)
           dy=z_psi(i,j)-z_psi(i,j-1)
           R_R=(r_psi(i,j)+r_psi(i,j+1))/2D0
!           if (sml_s_alpha) then
             ! We need the cylindrical q here (neglect finite epsilon corrections except for B-dependence)
!             R_R=R_R*EQ_axis_r
!           else
             R_R=R_R**2
!           endif
           dl=sqrt(dx**2+dy**2)
           qpsi(i)=qpsi(i)+dl/bp0/R_R

           dbdlpar2_fluxavg(i)=dbdlpar2_fluxavg(i)+(dl/bp0)*(dbdlpar2(i,j-1)+dbdlpar2(i,j))/2D0
           norm(i)=norm(i)+(dl/bp0)
        enddo
        qpsi(i)=abs(i_interpol(psigrid(i),0,1))*qpsi(i)/sml_2pi
        !if (sml_mype==0) print *,i,dbdlpar2_fluxavg(i),norm(i)
        dbdlpar2_fluxavg(i)=dbdlpar2_fluxavg(i)/norm(i)
     enddo
     qpsi(1)=qpsi(2)+(qpsi(3)-qpsi(2))*dpsi*psigrid(1)
     norm(1)=1D0
     dbdlpar2_fluxavg(1)=0D0

!#############################################################################
!##     evaluate trapped particle fraction                                  ##
!#############################################################################
     call findFonconstPsi(R,Z,Btot,eq_mr,eq_mz,r_psi,z_psi,Btot_psi,mpsi,mtheta)
     Btot2_psi(:,:)=Btot_psi(:,:)**2

     do i=1,mpsi
        Bmax_tmp=-1D10
        do j=1,mtheta
!           if(Btot_psi(i,j) .ge. Bmax_tmp) then
             Bmax_tmp=max(Bmax_tmp,Btot_psi(i,j))
!           endif
        enddo
        Bmax(i)=Bmax_tmp
!        eq_bmax(i)=Bmax_tmp
     enddo
     do i=1,mpsi
        Btot2_psi(i,:)=Btot2_psi(i,:)/(Bmax(i)**2)
        Btot_psi(i,:)=Btot_psi(i,:)/Bmax(i)
        sqrtoneb(i,:)=sqrt(1D0-Btot_psi(i,:))
     enddo
     call integrate(r_psi,z_psi,Bp_psi,Btot2_psi,mpsi,mtheta,Btot2_fluxavg)
     call integrate(r_psi,z_psi,Bp_psi,sqrtoneb,mpsi,mtheta,fluxavg)
     trapped_FM=fluxavg

     del_lamda=1D0/real(n_lamda)

     do i=1, n_lamda
        lamda=del_lamda*real(i-0.5D0)
        func_psi(:,:)=sqrt(1D0-lamda*Btot_psi(:,:))
        call integrate(r_psi,z_psi,Bp_psi,func_psi,mpsi,mtheta,func_fluxavgtmp)
        func_fluxavg(i,:)=func_fluxavgtmp(:)
     enddo

     trapped=0D0
     do j=2, mpsi-1
        do i=1, n_lamda
           lamda=del_lamda*real(i-0.5D0)
           trapped(j)=trapped(j)+lamda*del_lamda/func_fluxavg(i,j)
        enddo
     enddo

     do i=2, mpsi-1
        trapped(i)=1D0-3D0/4D0*Btot2_fluxavg(i)*trapped(i)
!        eq_ft(i)=trapped(i)
     enddo

     trapped(1)=trapped(2)+(trapped(3)-trapped(2))/(psigrid(3)-psigrid(2))&
*(psigrid(1)-psigrid(2))
     trapped(mpsi)=trapped(mpsi-1)+(trapped(mpsi-2)-trapped(mpsi-3))&
    /(psigrid(mpsi-2)-psigrid(mpsi-3))*(psigrid(mpsi)-psigrid(mpsi-1))

!     eq_ft(1)=trapped(1)
!     eq_ft(eq_mpsi)=trapped(eq_mpsi)
!##########################################################################
     if (sml_mype==0) then
       open(unit=125,file='trapped_fraction.dat',status='replace',action='write')
       ! IBM compiler on BGQ seems to need double quotes because of the triple question mark
       write(125,*) "# psi_normal -- trapped fr. (est.) -- trapped fr. -- 1.46*sqrt(eps) -- ???"
       open(unit=126,file='aspect_ratio.dat',status='replace',action='write')
       write(126,*) '# psi_normal -- rmin=(Rmax-Rmin)/2 -- rmaj=(Rmax+Rmin)/2 -- rmin/rmaj'
       open(unit=127,file='magnetic_field_data.dat',status='replace',action='write')
       write(127,*) '# psi_normal -- q_safety -- <(Bdl_par)^2> -- <B^2>'
     endif
     do i=1,mpsi
        Rmax=0D0
        Rmin=10D3
        do j=1,mtheta
           Rmax=max(Rmax,r_psi(i,j))
           Rmin=min(Rmin,r_psi(i,j))
        enddo
        if (i .gt. 1) then
          epspar(i,1)=(Rmax-Rmin)/2D0  ! Neoclassical minor radius
          epspar(i,2)=(Rmax+Rmin)/2d0  ! Neoclassical major radius
          epspar(i,3)=(Rmax-Rmin)/(Rmax+Rmin)  ! neoclassical inverse aspect ratio
        else
          epspar(i,1)=0.D0
          epspar(i,2)=EQ_axis_r
          epspar(i,3)=0.D0
        endif
        trapped2(i)=(1.46D0*sqrt(epspar(i,1))+2.40D0*epspar(i,1))/(1.46D0*sqrt(epspar(i,1))+  &
                     2.40D0*epspar(i,1)+((sqrt(1-epspar(i,1)))**3))
!#ifdef DEBUG
        if(sml_mype==0) then
            write(126,112) psigrid(i)/eq_x_psi,epspar(i,1),epspar(i,2),epspar(i,3)
        endif
!#endif
     enddo
!#ifdef DEBUG
     if(sml_mype==0) then

        do i=1,mpsi
            write(125,111) psigrid(i)/eq_x_psi, trapped2(i),trapped(i),1.46*sqrt(epspar(i,1)),trapped_FM(i)
            write(127,112) psigrid(i)/eq_x_psi, qpsi(i), dbdlpar2_fluxavg(i), btot2_fluxavg(i)
        enddo
        close(125)
        close(126)
        close(127)
      endif
!      eq_ft2=trapped_FM
!rh ifort recommends 20.13 instead of 19.13
111 format(5(e20.13,' '))
112 format(4(e20.13,' '))
113 format(3(e20.13,' '))
!#endif

  ! Interpolate epspar and trapped to final psi-grid and save in grid
  do i=1,grid%npsi00
    psipsi=grid%psi00min+grid%dpsi00*real(i-1,8)
    ip=min(floor(psipsi/dpsi)+1,mpsi-1)
    wp(2)=psipsi/dpsi-real(ip-1,8)
    wp(1)=1.D0-wp(2)
    grid%qsafety(i)=wp(1)*qpsi(ip)+wp(2)*qpsi(ip+1)
    grid%trapped(i)=wp(1)*trapped(ip)+wp(2)*trapped(ip+1)
    do j=1,3
      grid%epspar(i,j)=wp(1)*epspar(ip,j)+wp(2)*epspar(ip+1,j)
    enddo
  enddo

  ! Interpolate qsafety to psi_in
  do i=1,npsi
    psipsi=psi_in(i)
    ip=min(floor(psipsi/dpsi)+1,mpsi-1)
    wp(2)=psipsi/dpsi-real(ip-1,8)
    wp(1)=1.D0-wp(2)
    grid%qsafety2(i)=wp(1)*qpsi(ip)+wp(2)*qpsi(ip+1)
  enddo

end subroutine q_evaluation

subroutine findpsi(r_r,z_z,m,r,z,nw,nh,psi_rz,psi_temp)

    implicit none
    integer,intent(in) :: m,nw,nh
    real(kind=8),intent(in) :: r_r(m),z_z(m),r(nw),z(nh),psi_rz(nw,nh)
    real(kind=8),intent(out) :: psi_temp(m)
    integer :: i,j,k,p,q
    real(kind=8) :: delr,delz,psi1,psi2
    real(kind=8), external :: psi_interpol

    do k=1,m
      psi_temp(k)=psi_interpol(r_r(k),z_z(k),0,0)
!rh!       print *,r_r(m),z_z(m)
!rh       do i=1,nw-1
!rh          delr=(r_r(k)-r(i))/(r(i+1)-r(i))
!rh          if((delr .ge. 0D0) .and.(delr .lt. 1D0)) then
!rh              p=i
!rh              do j=1,nh-1
!rh                 delz=(z_z(k)-z(j))/(z(j+1)-z(j))
!rh                 if((delz .ge. 0D0) .and. (delz .lt.1D0)) then
!rh                     q=j
!rh                     !print *,psi_rz(i,j)
!rh                     psi1=psi_rz(p,q)+(psi_rz(p,q+1)-psi_rz(p,q))*delz
!rh                     psi2=psi_rz(p+1,q)+(psi_rz(p+1,q+1)-psi_rz(p+1,q))*delz
!rh
!rh                     psi_temp(k)=psi1+(psi2-psi1)*delr
!rh!                     print *,psi_temp(k)
!rh                     goto 10
!rh                 endif
!rh              enddo
!rh          endif
!rh       enddo
!rh       10 continue
     enddo

end subroutine findpsi

subroutine constpsi_contour(psipsi,r_con,z_con,psi_con,mtheta,gridm,r_psitmp,z_psitmp)

    implicit none
    integer ,intent(in) :: mtheta,gridm
    real(kind=8),intent(in) :: r_con(mtheta,gridm),z_con(mtheta,gridm),psi_con(mtheta,gridm)
    real(kind=8),intent(in) :: psipsi
    real(kind=8),intent(out) :: r_psitmp(mtheta),z_psitmp(mtheta)

    integer :: i,j,k
    real(kind=8) :: dpsi

    do i=1,mtheta
       do j=1,gridm-1
          dpsi=(psipsi-psi_con(i,j))/(psi_con(i,j+1)-psi_con(i,j))
          if((dpsi .ge. 0D0).and.(dpsi .lt. 1D0)) then
             r_psitmp(i)=r_con(i,j)+(r_con(i,j+1)-r_con(i,j))*dpsi
             z_psitmp(i)=z_con(i,j)+(z_con(i,j+1)-z_con(i,j))*dpsi
             goto 100
          endif
       enddo
       100 continue
    enddo

end subroutine  constpsi_contour

subroutine findFonconstPsi(r,z,f,nw,nh,r_psi,z_psi,fun_psi,np,mtheta)

    implicit none
    integer,intent(in) :: nw,nh,np,mtheta
    real(kind=8),intent(in) :: r(nw),z(nh),f(nw,nh)
    real(kind=8),intent(in) :: r_psi(np,mtheta),z_psi(np,mtheta)
    real(kind=8),intent(out) :: fun_psi(np,mtheta)

    integer :: i,j,k,l
    real(kind=8) :: delr,delz
    real(kind=8) :: fun1,fun2

    do  i=1,mtheta
        do j=2,np
           do k=1,nw-1
              delr=(r_psi(j,i)-r(k))/(r(k+1)-r(k))
              if ((delr .ge. 0D0) .and. (delr .lt. 1D0)) then
                 do l=1,nh-1
                    delz=(z_psi(j,i)-z(l))/(z(l+1)-z(l))
                    if ((delz .ge.0D0) .and. (delz .lt. 1D0)) then
                        fun1=f(k,l)+(f(k,l+1)-f(k,l))*delz
                        fun2=f(k+1,l)+(f(k+1,l+1)-f(k+1,l))*delz

                        fun_psi(j,i)=fun1+(fun2-fun1)*delr
                        goto 100
                     endif
                 enddo
              endif
           enddo
           100 continue
        enddo
    enddo
end subroutine findFonconstPsi

subroutine integrate(r_psi,z_psi,Bp_psi,fun_psi,np,mtheta,F_psi)

    implicit none
    integer,intent(in) :: np,mtheta
    real(kind=8),intent(in) :: r_psi(np,mtheta),z_psi(np,mtheta),&
                          Bp_psi(np,mtheta),fun_psi(np,mtheta)
    real(kind=8),intent(out) :: F_psi(np)

    integer :: i,j,k
    real(kind=8) :: fnorm,fave,f,dl,bpol,dlbpol,slobp
    real(kind=8) :: dx,dy

    do i=2,np-1
       fave=0D0
       fnorm=0D0
       slobp=0d0
       do j=2,mtheta
          f=(fun_psi(i,j-1)+fun_psi(i,j))/2D0
          bpol=(Bp_psi(i,j-1)+Bp_psi(i,j))/2D0
          dx=r_psi(i,j)-r_psi(i,j-1)
          dy=z_psi(i,j)-z_psi(i,j-1)
          dl=sqrt(dx**2+dy**2)
          dlbpol=dl/bpol
          fnorm=fnorm+dlbpol
          fave=fave+dlbpol*f
       enddo
       fave=fave/fnorm
       slobp=fnorm
       F_psi(i)=fave
    enddo
end subroutine integrate

