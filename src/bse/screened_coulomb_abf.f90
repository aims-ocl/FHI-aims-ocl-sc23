subroutine screened_coulomb_abf(ovlp_3ks_bse, screened_coulomb_bse)
  use dimensions
  use physics
  use prodbas
  use constants
  use runtime_choices
  use synchronize_mpi
  use evaluate_polarisability_freq, only: evaluate_polarisability_freq_1
  implicit none
! varaibles for screened_coulomb_bse
  integer :: n_low_state
  integer :: n_high_state
  integer, allocatable :: n_homo(:)
  real*8 :: omega_full_i_freq
  integer :: n_first(n_spin)
  real*8, dimension(n_basbas,n_states,n_states,n_spin) :: ovlp_3ks_bse
  real*8, dimension(n_basbas, n_loc_prodbas) :: screened_coulomb_bse
  real*8, dimension(:,:), allocatable :: polar_freq
  real*8, dimension(:,:), allocatable :: aux_ovlp3KS
  real*8, dimension(:,:), allocatable :: screened_coulomb_f
  real*8, dimension(:,:,:,:), allocatable :: ovlp_3KS_redist(:,:,:,:)
  integer :: i, j
  integer :: i_state, j_state, k_state, l_state, i_basbas, j_basbas, k_basbas, i_index, i_freq, i_freq_1, i_spin, i_task, n_states_loc, i_state_loc
  integer :: info
  character(*), parameter :: func = 'screened_coulomb_abf'
  print*, 'calculate screened_coulomb'
         n_low_state=1
         n_high_state=n_states
      n_first(:) = 1
      do i_spin = 1, n_spin
       do i_state = 1, n_states
        if (abs(occ_numbers(i_state, i_spin, 1)-dble(2/n_spin)) &
                        .lt.1.d-8) then
         n_first(i_spin)= i_state + 1
        endif
       enddo
       if(n_first(i_spin) .gt. n_states) then
         n_first(i_spin) = n_states
       endif
      enddo
!
!      temp_matr(:) = 0.d0
!
!      if(myid.eq.0) then
!        write(use_unit,'(2X, A,A,4I5)') &
!                   "HOMO and first non-fully-occupied", &
!                   " orbitals:", n_homo(:), n_first(:)
!        write(use_unit,*)
!      endif
! allocation
      n_states_loc = (n_high_state-1)/n_tasks + 1

      allocate(polar_freq(n_basbas, n_loc_prodbas), stat=info)
      call check_allocation(info, 'polar_freq', func)
!      allocate(screened_coulomb_bse( n_basbas, n_loc_prodbas ), stat=info)
!      call check_allocation(info, 'screened_coulomb', func)

      allocate(aux_ovlp3KS(n_basbas, n_states ), stat=info)
      call check_allocation(info, 'aux_ovlp3KS', func)
!      allocate(aux_coulomb_matr &
!             (n_states_loc,n_states,n_full_freq,n_spin), stat=info)
!      call check_allocation(info, 'aux_coulomb_matr', func)

!      allocate(screened_exchange(n_high_state), stat=info)
!      call check_allocation(info, 'screened_exchange', func)
!      allocate(coulomb_hole(n_high_state,n_spin), stat=info)
!      call check_allocation(info, 'coulomb_hole', func)

      allocate(ovlp_3KS_redist(n_basbas, n_states, n_states_loc, n_spin), stat=info)
      call check_allocation(info, 'ovlp_3KS_redist', func)
      allocate(screened_coulomb_f( n_basbas, n_basbas ), stat=info)
      call check_allocation(info, 'screened_coulomb_f', func)
! Redistribute ovlp_3KS, the parallelization is over the 3rd parameter in
! ovlp_3KS_redist
      ovlp_3KS_redist = 0
      do i_spin = 1, n_spin
        do i_state = 1,n_high_state
          aux_ovlp3KS = 0.
          do j_state = 1, n_states
            do j_basbas = 1, n_loc_prodbas
              i_index = map_prodbas(j_basbas,myid+1)
              aux_ovlp3KS(i_index,j_state) = ovlp_3ks_bse(j_basbas,j_state,i_state,i_spin)
            enddo
          enddo
          call sync_matrix(aux_ovlp3KS, n_basbas, n_states)
          if(MOD(i_state-1,n_tasks) == myid) then
            i_state_loc = (i_state-1)/n_tasks + 1
            ovlp_3KS_redist(:,:,i_state_loc,i_spin) = aux_ovlp3KS(:,:)
          endif
        enddo
      enddo
! i_freq
          i_freq = 1
!    evaluate the polarisability at frequency point i_freq
  print*, 'evaluate the polarisability at frequency point i_freq'
  if(.not.allocated(n_homo)) then
     allocate(n_homo(n_spin))
     n_homo(:) = 0
     do i_spin = 1, n_spin, 1
       do i_state = 1, n_states, 1
         if(occ_numbers(i_state,i_spin,1).gt.1.e-6) then
            n_homo(i_spin) = i_state
         endif
       enddo
      enddo
  endif
  omega_full_i_freq = 0
!  omega_full_i_freq = 0.000441
          call  evaluate_polarisability_freq_1 &
             ( 1, n_homo, n_first, n_high_state, occ_numbers, &
               omega_full_i_freq, &
               KS_eigenvalue, ovlp_3KS_redist, screened_coulomb_f &
             )
          do j_basbas = 1, n_loc_prodbas, 1
            i_index = map_prodbas(j_basbas,myid+1)
            do i_basbas = 1, n_basbas, 1
               polar_freq(i_basbas, j_basbas) = screened_coulomb_f(i_basbas, i_index)
             enddo
          enddo
!    now calculate the screened Coulomb interaction W
          call screened_coulomb_interaction &
             ( polar_freq, screened_coulomb_bse &
             )
!          if (write_screened_coulomb) then
!            open(unit = 99, file = 'screened_coulomb_abf')
!            do i = 1, n_basbas, 1
!              do j = 1, n_loc_prodbas, 1
!                write(99, *) screened_coulomb_bse(i, j)
!!                write(99, *) i, j, screened_coulomb_bse(i, j)
!!                if(abs(screened_coulomb_bse(i, j)) > 1e-9) write(99, *) i, j, screened_coulomb_bse(i, j)
!              enddo
!            enddo
!            close(99)
!          end if

!            open(unit = 99, file = 'screened_coulomb_abf')
!            do i = 1, n_loc_prodbas, 1
!              do j = 1, n_loc_prodbas, 1
!                if(abs(screened_coulomb_bse(i, j)) > 1e-9) write(99, *) i, j, screened_coulomb_bse(i, j)
!              enddo
!            enddo
!            close(99)
end subroutine
