!****s* FHI-aims/evaluate_Born_effective_charges_phonon_reduce_memory
!  NAME
!    evaluate_Born_effective_charges_phonon_reduce_memory
!  SYNOPSIS

subroutine evaluate_Born_effective_charges_phonon_reduce_memory( & 
           occ_numbers, Omega_MO, first_order_U_complex, & 
           Born_effective_charges)


!  PURPOSE
!    calculate Born_effective_charges.

!  Born_effective_charges= 
!
!         4                       
!   =   ----- * 4*(sum_(i,j,k) Omega_ij U_ji  
!         Nk                                 
!  USES

  use dimensions
  use pbc_lists , only : k_weights
  use mpi_tasks
  use synchronize_mpi
  use localorb_io


!  ARGUMENTS

  implicit none
 
  real*8, dimension(n_states, n_spin,n_k_points), intent(IN) :: occ_numbers 
  complex*16, dimension(n_basis, n_basis, n_k_points_task,3), intent(IN) :: Omega_MO
  complex*16, dimension(n_basis, n_basis, n_k_points_task), intent(IN) :: first_order_U_complex
 
  real*8, dimension(3), intent(out) :: Born_effective_charges

!  INPUTS
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  dielectric_constant

  integer :: i_state, j_state,i_coord, i_spin
  integer :: n_occ_states(n_spin)
  integer :: i_k_task, i_k_point
  integer :: i_basis
  character*1000 :: info_str
 
   Born_effective_charges(1:3) = 0.0d0

   i_k_task = 0
   do i_k_point = 1,n_k_points, 1
   if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
        i_k_task = i_k_task + 1

       do i_spin = 1, n_spin, 1
            n_occ_states(i_spin) = 0
            do i_state = n_states, 1, -1
            !if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.0.d0) then
            ! The following line implies the system does not have any fractional occupation number
            if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.1.e-6) then
               n_occ_states(i_spin) = i_state
            exit
            endif
            enddo
       enddo

      do i_coord = 1, 3
         do i_state = 1, n_occ_states(1)
         do j_state = n_occ_states(1)+1, n_states
 
          ! 2 for occ_number, 2 for exchange.
        Born_effective_charges(i_coord) = Born_effective_charges(i_coord) & 
          + dble(k_weights(i_k_point) * 4.0d0 * Omega_MO(i_state,j_state,  i_k_task, i_coord) &  
          * first_order_U_complex(j_state,i_state, i_k_task))

         enddo  ! i_state
         enddo  ! j_state
      enddo     ! i_coord 
         
   endif ! i_k_task   
   enddo ! i_k_point


  !-------shanghui begin parallel------
   call  sync_vector(Born_effective_charges, 3)
  !-------shanghui end parallel------

   
end subroutine evaluate_Born_effective_charges_phonon_reduce_memory
!******
