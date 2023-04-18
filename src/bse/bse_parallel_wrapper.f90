subroutine bse_parallel_wrapper( )
  use dimensions
  use physics
  use prodbas
  use constants
  use runtime_choices
  use synchronize_mpi
  use mpi_tasks
  use timing
  use hartree_fock
  implicit none
! variables for BSE
  real*8, dimension(:), allocatable :: qpe
  real*8, dimension(:, :, :, :), allocatable :: coulomb_4ks_vcvc, coulomb_4ks_vccv
  real*8, dimension(:, :, :, :), allocatable :: screened_coulomb_4ks_vvcc, screened_coulomb_4ks_vccv
  integer :: val, cond, v1, v2, c1, c2, counter, col, row
! varaibles for ovlp_3ks_bse
   real*8, dimension(:, :), allocatable :: ovlp_3fn_bse
   real*8, dimension(:, :, :, :), allocatable :: ovlp_3ks_bse
!  real*8, dimension(n_basis_pairs,n_basbas) :: ovlp_3fn_bse
!  real*8, dimension(n_basbas,n_states,n_states,n_spin) :: ovlp_3ks_bse
  integer :: i, j, k, l
  integer :: info
! variables for screened coulomb
  real*8, dimension(:,:), allocatable :: screened_coulomb_bse
  integer :: n_blk
!  integer, allocatable :: n_homo(:)
  integer :: n_first(n_spin)
  integer :: i_spin, i_state
  real*8 :: omega_full_i_freq
  real*8  :: occ_numbers_bse(n_states,n_spin)
  integer :: n_high_state
  real*8, dimension(:,:), allocatable :: polar_freq
! varaibles for calculating 4ks integrals
  integer :: errnum
  character(*), parameter :: deffmt = '2X'
!  integer :: n_pair_states, occ_states, unocc_states, num_prow, num_pcol
!  integer :: my_bse_ctxt
!  real*8 :: total_mem, proc_mem
!  integer :: a_state, k_cnt
!  integer :: na_rows, na_cols, my_prow, my_pcol, num_blocks
!  integer :: nprow, npcol, myrow, mycol, errnum
!  integer, external :: numroc
!  integer :: bse_desc(9)
!  character*128 :: info_str
  character(*), parameter :: func = 'bse_parallel_wrapper'

  call get_timestamps(time_total_bse, clock_time_total_bse)
  call get_timestamps(time_bse_ovlp_3ks, clock_time_bse_ovlp_3ks)
  call get_timestamps(time_bse_construct_mat, clock_time_bse_construct_mat)
  if(myid == 0) then
    print*, 'n_states', n_states
    print*, 'n_states_k(1)', n_states_k(1)
    print*, 'n_electrons', n_electrons
    print*, 'n_tasks', n_tasks
    print*, "n_basis_pairs", n_basis_pairs
    print*, "n_loc_prodbas", n_loc_prodbas
    print*, "n_basbas", n_basbas
  end if
  allocate(qpe(n_states_k(1)))
  qpe = 0
!  call qpe_calculation()
!  call qpe_calculation_bse(qpe)
  if(myid == 0) print*, 'read qpe from energy_qp file'
  if(myid == 0) then
    open(unit = 99, file = "energy_qp")
    read(99, *) qpe
    close(99)
  end if

  call mpi_bcast(qpe, n_states_k(1), mpi_double_precision, 0, mpi_comm_world, errnum)
!  call aims_stop()
  allocate(ovlp_3fn_bse(n_basis_pairs,n_loc_prodbas))
  call get_times(time_bse_ovlp_3fn, clock_time_bse_ovlp_3fn)
  ovlp_3fn_bse = ovlp_3fn
  call get_times(time_bse_ovlp_3fn, clock_time_bse_ovlp_3fn)
!
!  ndim1_o3ks = (n_states_k(1)-1)/np1_o3ks + 1
!  ndim2_o3ks = (n_states_k(1)-1)/np2_o3ks + 1
  ndim1_o3ks = (n_states-1)/np1_o3ks + 1
  ndim2_o3ks = (n_states-1)/np2_o3ks + 1

!  if(myid == 0) then
!      print*, 'ndim1_o3ks', ndim1_o3ks
!      print*, 'ndim2_o3ks', ndim2_o3ks
!  end if
  allocate(ovlp_3ks_bse(n_basbas,ndim1_o3ks,ndim2_o3ks,n_spin))

!  if (use_2d_corr) then
    call transform_ovlp3fn_2(n_states,KS_eigenvector, ovlp_3fn_bse, ovlp_3ks_bse)
!  else
!    call transform_ovlp3fn(n_states,KS_eigenvector, ovlp_3fn_bse, ovlp_3ks_bse)
!  end if

  call get_times(time_bse_ovlp_3ks, clock_time_bse_ovlp_3ks)
!  call aims_stop()
  call bse_2d(ndim1_o3ks,ndim2_o3ks, ovlp_3ks_bse, qpe)
  call get_times(time_total_bse, clock_time_total_bse)
!  call output_timeheader(deffmt, 'Detailed time accounting')
!            call output_times(deffmt, 'Time for get_coeff_3fn_v_2d', &
!                  time_bse_ovlp_3fn, clock_time_bse_ovlp_3fn)
!            call output_times(deffmt, 'Time for get ovlp_3ks', &
!                  time_bse_ovlp_3ks, clock_time_bse_ovlp_3ks)
!            call output_times(deffmt, 'Total time for BSE', &
!                  time_total_bse, clock_time_total_bse)
!            call output_times(deffmt, 'Time for BSE intiliazaiton', &
!                  time_bse_init, clock_time_bse_init)
!            call output_times(deffmt, 'Time for building <ia|V|jb>', &
!                  time_bse_build_V, clock_time_bse_build_V)
!            call output_times(deffmt, 'Time for calculating W0 abf', &
!                  time_bse_w0, clock_time_bse_w0)
!            call output_times(deffmt, 'Time for calculating W * ovlp_3ks', &
!                  time_bse_w_ovlp_3ks, clock_time_bse_w_ovlp_3ks)
!            call output_times(deffmt, 'Time for building <ij|W|ab>', &
!                  time_bse_build_W, clock_time_bse_build_W)
!            call output_times(deffmt, 'Time for solving BSE matrix', &
!                  time_bse_solve_mat, clock_time_bse_solve_mat)
end subroutine bse_parallel_wrapper
