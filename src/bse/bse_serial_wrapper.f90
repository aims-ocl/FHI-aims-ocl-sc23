subroutine bse_serial_wrapper( )
  use dimensions
  use physics
  use prodbas
  use constants
  use runtime_choices
  use synchronize_mpi
  implicit none
! variables for BSE
  real*8, dimension(:), allocatable :: qpe
  real*8, dimension(:, :, :, :), allocatable :: coulomb_4ks_vcvc, coulomb_4ks_vccv
  real*8, dimension(:, :, :, :), allocatable :: screened_coulomb_4ks_vvcc, screened_coulomb_4ks_vccv
  integer :: val, cond, v1, v2, c1, c2, counter, col, row
! varaibles for ovlp_3ks_bse
  real*8, dimension(n_basis_pairs,n_basbas) :: ovlp_3fn_bse
  real*8, dimension(n_basbas,n_states,n_states,n_spin) :: ovlp_3ks_bse
  integer :: i, j, k, l
  real*8, dimension(:,:), allocatable :: screened_coulomb_bse
  integer :: info
  character(*), parameter :: func = 'bse_serial_wrapper'

  print*, '----------------------------------------------------------------------'
  print*, 'BSE calculation starts...'

  print*, 'calculating quasi-particle energy'
  allocate(qpe(n_states))
  if(read_qpe) then
    print*, 'read qpe from energy_qp file'
    open(unit = 99, file = "energy_qp")
    read(99, *) qpe
    close(99)
  else
    call qpe_calculation_bse(qpe)
  end if
! bse_reduce_matrix
  if (bse_reduce_matrix) then
    do i = 1, n_states
      if (qpe(i) > bse_reduce_occ) then
        exit
      end if
    end do
    bse_lower_limit = i
    do i = n_states, 1, -1
      if (qpe(i) < bse_reduce_unocc) then
        exit
      end if
    end do
    bse_upper_limit = i
    if(myid == 0) then
      print*, 'bse_reduce_occ', bse_reduce_occ
      print*, 'bse_lower_limit', bse_lower_limit
      print*, 'bse_reduce_unocc', bse_reduce_unocc
      print*, 'bse_upper_limit', bse_upper_limit
    end if
  end if
 
  !print*, 'qpe', qpe

  print*, 'calculating ovlp_3ks_bse'
  call get_coeff_3fn_v_1d(ovlp_3fn_bse)
!  open(unit = 99, file = 'ovlp_3fn_bse_serial')
!  do i = 1, n_basis_pairs
!    do j = 1, n_loc_prodbas
!      write(99, *) i, j, ovlp_3fn_bse(i, j)
!    enddo
!  enddo
!  close(99)
  call transform_ovlp3fn(n_states, KS_eigenvector, ovlp_3fn_bse, ovlp_3ks_bse)
!<<<<<<< HEAD:bse/bse_serial_wrapper.f90
!  open(unit = 99, file = 'ovlp_3ks_bse')
!  do i = 1, n_loc_prodbas, 1
!    do j = 1, n_states, 1
!      do k = 1, n_states, 1
!=======
!!  open(unit = 99, file = 'ovlp_3ks_bse_serial')
!!  do i = 1, n_loc_prodbas, 1
!>>>>>>> 73e3766ec8026f2b4a56c471a608bbe6abf8cab6:bse/bse_wrapper.f90
!!    do j = 1, n_states, 1
!        do l = 1, n_spin, 1
!          write(99, *) ovlp_3ks_bse(i, j, k, l)
!!          write(99, *) i, j, k, l, ovlp_3ks_bse(i, j, k, l)
!        enddo
!      enddo
!    enddo
!  enddo
!  close(99)
! calculate coulomb_4ks
  print*, 'calculating coulomb_4ks'
  val = int(n_electrons) / 2
  cond = n_states - val
  print*, 'n_states, n_electrons', n_states, n_electrons
  allocate(coulomb_4ks_vcvc(val, cond, val, cond))
  allocate(coulomb_4ks_vccv(val, cond, cond, val))

  call coulomb_4ks(val, cond, ovlp_3ks_bse, coulomb_4ks_vcvc, coulomb_4ks_vccv)

  print*, 'calculating screened_coulomb_abf'
  allocate(screened_coulomb_bse( n_basbas, n_loc_prodbas ), stat=info)
  call check_allocation(info, 'screened_coulomb', func)

  call screened_coulomb_abf(ovlp_3ks_bse, screened_coulomb_bse)
!  open(unit = 99, file = 'screened_coulomb_bse')
!  do i = 1, n_loc_prodbas, 1
!    do j = 1, n_loc_prodbas, 1
!      if(abs(screened_coulomb_bse(i, j)) > 1e-9) write(99, *) i, j, screened_coulomb_bse(i, j)
!    enddo
!  enddo
!  close(99)

  print*, 'calculating screened_coulomb_4ks'
  allocate(screened_coulomb_4ks_vvcc(val, val, cond, cond))
  allocate(screened_coulomb_4ks_vccv(val, cond, cond, val))
  call screened_coulomb_4ks(val, cond, screened_coulomb_bse, ovlp_3ks_bse, screened_coulomb_4ks_vvcc, screened_coulomb_4ks_vccv)
  print*, 'construct and solve BSE matrix'
  call construct_solve_bse_mat(val, cond, qpe, coulomb_4ks_vcvc, coulomb_4ks_vccv, screened_coulomb_4ks_vvcc, screened_coulomb_4ks_vccv)

  deallocate(qpe)
  deallocate(screened_coulomb_4ks_vvcc)
  deallocate(screened_coulomb_4ks_vccv)
  deallocate(coulomb_4ks_vcvc)
  deallocate(coulomb_4ks_vccv)
end subroutine

