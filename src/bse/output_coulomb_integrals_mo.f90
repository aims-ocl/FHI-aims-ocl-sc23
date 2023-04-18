subroutine output_coulomb_integrals_mo ()
  use dimensions
  use physics
  use prodbas
  use constants
  use runtime_choices
  use synchronize_mpi
  implicit none

  real*8, dimension(:, :), allocatable :: ovlp_3fn_bse
  real*8, dimension(:, :, :, :), allocatable :: ovlp_3ks_bse
  integer :: i, j, k, l, counter
  integer :: info
  real*8 :: temp_4ks
  real*8, dimension(n_states, n_states, n_states, n_states) :: coulomb_integrals_mo

  character(*), parameter :: func = 'output_coulomb_integrals_mo'
  
!  print*, 'output_coulomb_integrals_mo'
!  call initialize_prodbas()
  if(n_tasks > 1) then
    call output_coulomb_integrals_mo_2d()
  else
!    print*, 'calculating ovlp_3ks_bse'
    allocate(ovlp_3fn_bse(n_basis_pairs,n_loc_prodbas), stat=info)
    call check_allocation(info, 'ovlp_3fn_bse', func)
    call get_coeff_3fn_v_1d(ovlp_3fn_bse)

    allocate(ovlp_3ks_bse(n_loc_prodbas, n_states, n_states, n_spin),stat=info)
    call check_allocation(info, 'ovlp_3KS (1D)', func)
!!    print*, 'after get_coeff_3fn_v_1d'
    call transform_ovlp3fn(n_states, KS_eigenvector, ovlp_3fn_bse, ovlp_3ks_bse)
!    print*, 'after transform_ovlp3fn'
    if (coulomb_integral_format == 'qtmp') then
    open(unit = 99, file = "bielec_mo")
    do l = 1, n_states
      do k = 1, n_states
        do j = l, n_states
          do i = max(j, k), n_states
            if (i == j .and. k > l) continue
            temp_4ks = 0
            do counter = 1, n_basbas
              temp_4ks = temp_4ks + ovlp_3ks_bse(counter, i, k, n_spin) *& 
              ovlp_3ks_bse(counter, j, l, n_spin)
            end do
            if (abs(temp_4ks) > 1E-16) write(99, *) i, j, k, l, temp_4ks
          end do
        end do
      end do
    end do
    else if (coulomb_integral_format == 'full') then
      open(unit = 99, file = "coulomb_integrals_mo.out")
      do i = 1, n_states
        do j = 1, n_states
          do k = 1, n_states
            do l = 1, n_states
              temp_4ks = 0
              do counter = 1, n_basbas
                temp_4ks = temp_4ks + ovlp_3ks_bse(counter, i, k, n_spin) *& 
                ovlp_3ks_bse(counter, j, l, n_spin)
              end do
  
              write(99, *) i, j, k, l, temp_4ks
            end do
          end do
        end do
      end do
    end if
    close(99)

!    do i = 1, n_states
!      do j = 1, n_states
!        do k = 1, n_states
!          do l = 1, n_states
!            temp_4ks = 0
!            do counter = 1, n_basbas
!              temp_4ks = temp_4ks + ovlp_3ks_bse(counter, i, j, n_spin) *& 
!              ovlp_3ks_bse(counter, k, l, n_spin)
!            end do
!            write(99, *) i, j, k, l, temp_4ks
!  !          write(99, *) temp_4ks
!            coulomb_integrals_mo(i, j, k, l) = temp_4ks
!          end do
!        end do
!      end do
!    end do
    close(99)
   deallocate(ovlp_3ks_bse)
   deallocate(ovlp_3fn_bse)
  endif
end subroutine output_coulomb_integrals_mo
