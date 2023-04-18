
  subroutine ci_acc_evaluate_c()

  use dimensions
  use runtime_choices
  use physics
  use mpi_tasks
  use mpi_utilities
  use fciqmc_module

    implicit none
    ! temp variables
    integer, dimension(:,:), allocatable :: occ_num_1
    ! temp indices
    integer :: i_state, j_state, a_state, b_state, i_spin
    integer :: i_ci
    integer :: errnum
    !character*128 :: info_str
    logical :: exit_flag
    real*8  :: configuration_eigenvalue
    Double precision :: ci_acc_rtmp
    integer :: index_start(4)
    integer :: wnum,iii

    
    if (.not.allocated(occ_num_1)) then
        allocate(occ_num_1(n_states,2), stat=errnum)
        call check_allocation(errnum, 'occ_num_1 (1D) in CI')
    endif

    c_0    = c_0 + (w_0 + (E_HF - E_ci) * c_0) / (E_ci-E_HF)

    do i_ci = 1, n_configuration_table(myid+1)
      ci_acc_rtmp = E_ci - ci_acc_Hdiag(i_ci)
      c_vect(i_ci) = c_vect(i_ci) + (w_vect(i_ci) + (ci_acc_Es(i_ci) &
                     + V_00 - E_ci) * c_vect(i_ci)) / ci_acc_rtmp

    end do

    call normalization_vector(c_vect,n_configuration_table(myid+1),c_0,norm_A,1)

end subroutine ci_acc_evaluate_c

