!        -------------------------------
!        solver using lapack band solver
!        -------------------------------
         subroutine bsolver(n,nnz,nrhs,values,rowind,colptr,b,ldb,info)
           use perf_monitor
           implicit none
           integer :: n,nnz,nrhs,rowind(nnz),colptr(n+1),ldb,info
           real*8 :: values(nnz), b(ldb,nrhs)

           integer :: i,j,k,kstart,kend
           real*8 :: aij

           integer :: kl,ku, ldAB,istat
           integer :: ipiv(n)
           real*8, allocatable, dimension(:,:) :: AB

           integer, parameter :: idebug = 0
           integer, parameter :: outdev = 17

!          ------------------
!          lapack band solver
!          ------------------
           interface
             subroutine dgbsv(n,kl,ku,nrhs,AB,ldAB,ipiv,b,ldb,info)
             implicit none
             integer :: n,kl,ku,nrhs,ldAB,ldb,info
             integer :: ipiv(n)
             real*8 AB(ldAB,n)
             real*8 b(ldb,nrhs)
             end subroutine dgbsv
           end interface

           call t_startf("COL_F_SOLVER_PREP")

!          -------------------------------------
!          estimate  lower/upper bandwidth kl/ku
!          -------------------------------------
           kl = 0
           ku = 0
!$omp      parallel do default(none)                                       &
!$omp&     private(j,k,kstart,kend,i)                                      &
!$omp&     reduction(max:kl,ku)                                            &
!$omp&     shared(n,colptr,rowind)                                         &
!$omp&     num_threads(col_f_nthreads)
           do j=1,n
             kstart = colptr(j)
             kend = colptr(j+1)-1
             do k=kstart,kend
               i = rowind(k)
               kl = max(kl,i-j)
               ku = max(ku,j-i)
             enddo
           enddo

           if (idebug .ge.1) then
             write(*,*) 'bsolver: kl,ku ', kl,ku
           endif

           if (idebug .ge. 2) then
             open(outdev,file='A.dat',form='formatted',                    &
                  access='sequential')
             rewind(outdev)
             write(outdev,*) '% i j aij '
             do j=1,n
               kstart = colptr(j)
               kend = colptr(j+1)-1
               do k=kstart,kend
                 i = rowind(k)
                 aij = values(k)
                 write(outdev,*) i, j, aij
               enddo
             enddo
             close(outdev)
           endif

!          -------------------------------------
!          need extra kl entries due to pivoting
!          -------------------------------------
           ldAB = kl + (kl+ku) + 1
! !$omp      critical (alloc1)
           allocate( AB(ldAB,n), stat=istat)
! !$omp      end critical (alloc1)
           if (istat.ne.0) then
             write(*,*) 'bsolver: allocate(AB(',                           &
     &                  ldb,',',n,') return ',istat
             info = -1
             return
           endif
!$omp      parallel do default(none)                                       &
!$omp&     private(j)                                                      &
!$omp&     shared(n,AB)                                                    &
!$omp&     num_threads(col_f_nthreads)
           do j=1,n
             AB(:,j) = 0
           enddo

!$omp      parallel do default(none)                                       &
!$omp&     private(j,k,kstart,kend,i,aij)                                  &
!$omp&     shared(n,kl,ku,colptr,rowind,values,AB)                         &
!$omp&     num_threads(col_f_nthreads)
           do j=1,n
             kstart = colptr(j)
             kend = colptr(j+1)-1
             do k=kstart,kend
               i = rowind(k)
               aij = values(k)
               AB(kl+ku+1+i-j,j) = aij
             enddo
           enddo
           call t_stopf("COL_F_SOLVER_PREP")

           call t_startf("COL_F_SOLVER_DGBSV")
           call dgbsv(n,kl,ku,nrhs,AB,ldAB, ipiv, b,ldb,info)
           call t_stopf("COL_F_SOLVER_DGBSV")
           if (info.ne.0) then
               write(*,*) 'bsolver: dgbsv return info=',info
           endif
! !$omp      critical (alloc1)
           deallocate( AB, stat=istat )
! !$omp      end critical (alloc1)
           if (istat.ne.0) then
             write(*,*) 'bsolver: deallocate(AB) return ',istat
           endif
           return
         end subroutine
