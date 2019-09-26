  subroutine SUBNAME(label,arr)
    character(len=*) :: label
    REAL_COMPLEX,dimension(DB_OUT_DIM) :: arr
    real :: local_sum, global_sum

    call calculate_test_sum(arr,local_sum, global_sum)
    IF (mype.EQ.0) WRITE(*,"(A,ES20.12)") trim(label),global_sum
  end subroutine SUBNAME
#undef SUBNAME
#undef DB_OUT_DIM
