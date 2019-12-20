      attributes(host) &
      subroutine init_push_mod_gpu( grid )

      use grid_class, only :  grid_type


      use grid_class_gpu, only : &
          update_device_grid_type

      implicit none
      type(grid_type) :: grid

      ! --------------------------------------------------------
      ! perform one-time initialization related to push() on GPU
      ! --------------------------------------------------------
      call push_update_device_gpu()

      call update_device_grid_type(grid )



      return
      end subroutine init_push_mod_gpu
