#ifdef USE_GPU
       attributes(global)                                                &
#endif
       subroutine gen_perm_gpu_pass2(nx,ny,xmin,ymin,                    &
     &         inv_dx,inv_dy, n,gid,xy,xydim,ncopy,lcount,iperm,xstart)
       use cudafor
       use precision_mod_gpu
       implicit none

       integer, value :: nx,ny,n,ncopy
       real(kind=work_p),value :: xmin, ymin, inv_dx, inv_dy
       integer*8 :: gid(n)
       integer, value :: xydim
       real(kind=work_p) :: xy(xydim,2)
       integer(kind=4) :: lcount(nx*ny, ncopy)
       integer :: iperm(n)
       integer :: xstart(nx*ny+1)

#ifdef USE_GPU
       attributes(device) :: gid, xy, lcount, xstart, iperm
#endif

       integer :: ith,nthread,iblock, icopy
       integer :: nblock,gnthread,gith
       integer :: i,ix,iy,k,ip

       ith = threadIdx%x + (threadIdx%y - 1)*blockDim%x
       nthread = blockDim%x * blockDim%y
       iblock = blockIdx%x + (blockIdx%y-1)*gridDim%x

       nblock = gridDim%x * gridDim%y
       gnthread = nblock * nthread
       gith = ith + (iblock-1)*nthread

       icopy = 1 + mod(iblock-1,ncopy)

!      -------------------------
!      perform geometric hashing
!      -------------------------
       do i=gith,n,gnthread
        if (gid(i) <= 0) then
         ix = nx
         iy = ny
        else
         ix = 1+int( (xy(i,1) - xmin) * inv_dx )
         iy = 1+int( (xy(i,2) - ymin) * inv_dy )


         if ((iy < 0).or.(iy > ny)) then
!           print*,'i,iy,xy(i,2) ',i,iy,xy(i,2)
           iy = max(1,min(ny,iy))
         endif

         if ((ix < 0).or.(ix > ny)) then
!           print*,'i,ix,xy(i,1) ',i,ix,xy(i,1)
           ix = max(1,min(ny,ix))
         endif

         ix = max(1, min(nx, ix) )
         iy = max(1, min(ny, iy) )
        endif

         k = ix + (iy-1)*nx

         ip = atomicadd( lcount(k,icopy), 1 )
         iperm(ip) =  i
       enddo

       return
       end subroutine  gen_perm_gpu_pass2
