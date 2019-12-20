#ifdef USE_GPU
       attributes(global)                                                  &
#endif
       subroutine gen_perm_gpu_pass1( nx,ny,xmin,ymin,               &
     &               inv_dx,inv_dy, n,gid,xy,xydim, ncopy,lcount)
       use cudafor
       use precision_mod_gpu
       integer, value :: nx,ny,ncopy
       real(kind=work_p), value :: xmin, ymin, inv_dx, inv_dy

       integer, value :: n, xydim
       integer*8 :: gid(n)
       real(kind=work_p) :: xy(xydim,2)
       integer(kind=4) :: lcount( nx*ny, ncopy)

#ifdef USE_GPU
       attributes(device) :: gid,xy,lcount
#endif

       integer :: ith, nthread, iblock, nblock
       integer :: gith, gnthread
       integer :: i, k, ix,iy , idummy,icopy


!      ------------------------------------
!      assume lcount(nx*ny,:) to be zero
!      ------------------------------------

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

         ix = max(1,min(nx,ix))
         iy = max(1,min(ny,iy))

        endif


         k = ix + (iy-1)*nx
         idummy = atomicadd( lcount(k,icopy), 1 )
       enddo

       return
       end subroutine gen_perm_gpu_pass1
