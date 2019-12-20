module boundary_class_gpu
use dimensions_mod_gpu
implicit none


  type range_type
     integer :: start,end     
  end type range_type

  type boundary2_type
     integer :: nseg
     type(range_type) :: iseg(iseg_dim)
  end type boundary2_type


end module boundary_class_gpu
