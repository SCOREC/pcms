      subroutine gen_v_node(nnodes, ntrig, max_t, max_v,                    &
     &                      nd, num_t_node, tr_node,                        &
     &                      num_v_node,v_node)
      implicit none
!   ------------------------------
!   generate the node to node list
!   num_v_node(1:nnodes)
!   v_node(1:max_v,1:nnodes)
!   ------------------------------

      integer :: nnodes, ntrig, max_t, max_v
      integer, dimension(1:3,1:ntrig) :: nd
      integer, dimension(1:nnodes) :: num_t_node
      integer, dimension(1:max_t,1:nnodes) :: tr_node

      integer, dimension(1:nnodes) :: num_v_node
      integer, dimension(1:max_v,1:nnodes) :: v_node

      intent(in) :: nd, num_t_node, tr_node
      intent(out) :: num_v_node, v_node
!     ---------------
!     local variables
!     ---------------
      integer, dimension(nnodes) :: isseen
      integer :: i,j,k, vj, vk, itr, ndeg_i, max_deg
      logical :: is_new, isvalid

      isseen(1:nnodes) = 0

      max_deg = 0
      do i=1,nnodes
        isseen(i) = i
        ndeg_i = 1
        v_node(1,i) = i

        do j=1,num_t_node(i)
           itr  = tr_node(j,i)
           do k=1,3
             vk = nd(k,itr)
             is_new = (isseen(vk) .ne. i)
             if (is_new) then
                ndeg_i = ndeg_i + 1
                isseen(vk) = i
                v_node(ndeg_i,i) = vk
             endif
           enddo
        enddo
        num_v_node(i) = ndeg_i

        isvalid = all( (1 <= v_node( 1:num_v_node(i),i) ) .and.          &
     &                 (v_node(1:num_v_node(i),i) <= nnodes ) )
        if (.not.isvalid) then
              write(*,*) 'i, num_v_node(i) ',i,num_v_node(i)
              write(*,*) 'v_node(j,i) ', v_node(:,i)
              stop '** error **'
        endif

       enddo

       max_deg = maxval( num_v_node(1:nnodes) )
       isvalid = (max_deg <= max_v)
       if (.not.isvalid) then
          write(*,*) 'increase max_v to ',max_deg
          stop '** error ** '
       endif
!      -------------
!      double check
!      -------------
       do i=1,nnodes
         do j=1,num_v_node(i)
            vj = v_node(j,i)
            isvalid = (1 <= vj).and.(vj <= nnodes)
            if (.not.isvalid) then
              write(*,*) 'gen_v_node:'
              write(*,*) 'i, num_v_node(i) ',i,num_v_node(i)
              write(*,*) 'v_node(j,i) ', v_node(:,i)
              stop '** error **'
            endif
          enddo
        enddo
       return
       end subroutine gen_v_node
