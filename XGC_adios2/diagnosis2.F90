#include "adios_macro.h"

subroutine output_bfield
  use eq_module
  use sml_module
  implicit none
  real (kind=8) :: psi,br,bz,bphi
  real (kind=8) :: r,z,rstart,rend,dr
  real (kind=8), external:: psi_interpol
  integer :: i,N=1000
#ifdef ADIOS
  integer*8 :: buf_id, buf_size, total_size
  integer :: err
  real (kind=8), allocatable :: rvec(:),psi_eq_x_psi(:),dpdr(:),abs_bphi(:)
  allocate(rvec(N),psi_eq_x_psi(N),dpdr(N),abs_bphi(N))
#endif
  z=0
  rstart=sml_bd_min_r ! general r start --2002/02/01
  rend=sml_bd_max_r ! general r end -- 2002/02/01
  dr=(rend-rstart)/real(N)
#ifndef ADIOS
  write(60,*) '#R   Psi/eq_x_psi  dpdr  dpdz   d2pdr2 d2pdrdz d2pdz2'
#endif
  do i=1, N
     r=rstart+dr*real(i)
     psi=psi_interpol(r,z,0,0)
     call bvec_interpol(r,z,3.141592/2.,br,bz,bphi)
#ifdef ADIOS
     rvec(i)=r
     psi_eq_x_psi(i)=psi/eq_x_psi;
     dpdr(i)=dsqrt(br**2+bz**2)
     abs_bphi(i)=abs(bphi)
#else
     write(60,1000) (r),psi/eq_x_psi, dsqrt(br**2+bz**2), abs(bphi), r
!          if(psi > NC_x_psi*1.1) then
!             return
!          endif
#endif 
  enddo
#ifdef ADIOS
  buf_size = 4 + 4 * N * 8
  ADIOS_OPEN(buf_id,'output.bfield','xgc.bfieldm.bp','a',sml_comm_null,err)
  ADIOS_SET_PATH(buf_id,'/bfield',err)
  ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
  ADIOS_WRITE(buf_id,N,err)
  ADIOS_WRITE(buf_id,rvec,err)
  ADIOS_WRITE(buf_id,psi_eq_x_psi,err)
  ADIOS_WRITE(buf_id,dpdr,err)
  ADIOS_WRITE(buf_id,abs_bphi,err)
  ADIOS_CLOSE(buf_id,err)
  deallocate(rvec,psi_eq_x_psi,dpdr,abs_bphi)
#endif
1000 format (5(e19.13,' '))
end subroutine output_bfield

subroutine output_volumes(grid)
  use grid_class
  use diag_module
  use sml_module
  implicit none
  type(grid_type)::grid
#ifdef ADIOS
  integer*8 :: buf_id, buf_size, total_size
  integer :: err

  buf_size = (4 + 8 * diag_1d_npsi) + ( 0 )
  ADIOS_OPEN(buf_id,'output.volumes','xgc.volumes.bp','a',sml_comm_null,err)
  ADIOS_GROUP_SIZE(buf_id,buf_size,total_size,err)
  ADIOS_WRITE_LBL(buf_id,'diag_1d_npsi',diag_1d_npsi,err) 
  ADIOS_WRITE_LBL(buf_id,'diag_1d_vol',diag_1d_vol,err) 

  ADIOS_CLOSE(buf_id,err)
#endif
end subroutine

! routines to check interpolation validity
subroutine output_psi_der
  use eq_module
  use sml_module
  implicit none
  real (kind=8) :: psi,br,bz,bphi
  real (kind=8) :: r,z,rstart,rend,dr,zstart,zend,dz
  real (kind=8), external:: psi_interpol
  integer :: i,Nr,Nz,j
  Nz=10
  Nr=1000
  zstart=sml_bd_min_z !
  zend=sml_bd_max_z 
  dz=(zend-zstart)/real(Nz)

  rstart=sml_bd_min_r! general r start --2002/02/01
  rend=sml_bd_max_r ! general r end -- 2002/02/01
  dr=(rend-rstart)/real(Nr)
  write(60,*) '#minor_r   Psi/eq_x_psi  bp  bt   R'

  do j=1, Nz

  do i=1, Nr
     r=rstart+dr*real(i)
     z=zstart+dz*real(j)
     psi=psi_interpol(r,z,0,0)
     call bvec_interpol(r,z,3.141592/2.,br,bz,bphi)
     write(1060,1000) (r),psi/eq_x_psi, psi_interpol(r,z,1,0),psi_interpol(r,z,0,1),&
          psi_interpol(r,z,2,0),psi_interpol(r,z,1,1),psi_interpol(r,z,0,2)
     !     if(psi > NC_x_psi*1.1) then
     !        return
     !     endif
  enddo
     write(1060,*) ' '
  enddo



1000 format (7(e19.13,' '))
end subroutine output_psi_der


