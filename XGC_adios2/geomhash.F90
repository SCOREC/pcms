        subroutine geomhash(n,rz,nr,nz,                                    &
     &        rmin,dr_inv,zmin,dz_inv,xstart,list)
        implicit none
        integer, intent(in) :: n, nr,nz
        real(kind=8), dimension(1:2,n), intent(in) :: rz
        real(kind=8), intent(in) :: rmin,zmin,dr_inv,dz_inv
        integer, dimension(nr*nz+1) :: xstart
        integer, dimension(n) :: list

        integer :: ij, i,j,k, ifree, ipos, ierr
        integer :: istart,iend,isize
        integer, parameter :: nb = 1024

        real(kind=8), dimension(nb) :: veci, vecj
        integer, dimension(nr,nz) :: icount
        real(kind=8), parameter :: one = 1.0d0
        real(kind=8) :: roff,zoff
        real(kind=8) :: veci_min, veci_max, vecj_min, vecj_max
        real(kind=8) :: veci_k, vecj_k, rk, zk


        do j=1,nz
        do i=1,nr
          icount(i,j) = 0
        enddo
        enddo

!       ------------------------------------
!       1st pass to count how many particles
!       land in the box
!       ------------------------------------
        roff = (one - rmin*dr_inv)
        zoff = (one - zmin*dz_inv)
        veci_min = dble(1)
        veci_max = dble(nr)
        vecj_min = dble(1)
        vecj_max = dble(nz)

        do istart=1,n,nb
          iend = min(n,istart+nb-1)
          isize = iend - istart + 1

          do k=1,isize
            rk = rz(1, (istart-1) + k)
            zk = rz(2, (istart-1) + k)

!           veci(k) = one + (rk - rmin) * dr_inv
!           vecj(k) = one + (zk - zmin) * dz_inv

            veci_k  = rk*dr_inv + roff
            vecj_k  = zk*dz_inv + zoff

            veci(k) = max(veci_min,min(veci_max, veci_k  ) )
            vecj(k) = max(vecj_min,min(vecj_max, vecj_k  ) )
    
          enddo

          do k=1,isize
!            i = max(1,min(nr, int(veci(k)) ) )
!            j = max(1,min(nz, int(vecj(k)) ) )

            i = int(veci(k))
            j = int(vecj(k))
            icount(i,j) = icount(i,j) + 1
          enddo
        enddo

!       --------------------
!       setup data structure
!       --------------------
        ifree = 1
        do j=1,nz
        do i=1,nr
          ij = i + (j-1)*nr
          xstart(ij) = ifree
          ifree = ifree + icount(i,j)
        enddo
        enddo
        xstart(nz*nr+1) = ifree


!       -------------------------------
!       second pass to deposit particle
!       -------------------------------

        do j=1,nz
        do i=1,nr
          ij = i + (j-1)*nr
          icount(i,j) = xstart(ij)
        enddo
        enddo

        do istart=1,n,nb
          iend = min(n,istart+nb-1)
          isize = iend - istart + 1

          do k=1,isize

            rk = rz(1, (istart-1) + k)
            zk = rz(2, (istart-1) + k)

!           veci(k) = one + (rk - rmin) * dr_inv
!           vecj(k) = one + (zk - zmin) * dz_inv

            veci_k  = rk*dr_inv + roff
            vecj_k  = zk*dz_inv + zoff

            veci(k) = max(veci_min,min(veci_max, veci_k  ) )
            vecj(k) = max(vecj_min,min(vecj_max, vecj_k  ) )
          enddo

          do k=1,isize
!            i = max(1,min(nr, int(veci(k)) ) )
!            j = max(1,min(nz, int(vecj(k)) ) )

             i = int( veci(k) )
             j = int( vecj(k) )
            ipos = icount(i,j)
            list(ipos) = (istart-1) + k
            icount(i,j) = icount(i,j) + 1
          enddo
        enddo



        return
        end subroutine 
