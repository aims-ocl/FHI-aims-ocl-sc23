!****s* FHI-aims/evaluate_Omega_MO_diag
!  NAME
!    evaluate_Omega_MO_v1
!  SYNOPSIS

subroutine evaluate_Omega_MO_diag( & 
           r_complex, &
           Omega_MO_diag)
               

!  PURPOSE
!    Omega_MO_diag(i)  = < i(k)|-r|i(k) > 
!


!  USES

  use dimensions
  use runtime_choices, only : real_eigenvectors
  use mpi_tasks
  use physics, only : occ_numbers, KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue  

!  ARGUMENTS

  implicit none
  

  complex*16, dimension( n_basis, n_basis, n_k_points_task), intent(IN) :: r_complex
  complex*16, dimension( n_basis , n_k_points_task), intent(OUT) :: Omega_MO_diag

!  INPUTS
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  Omega_MO


  integer :: i_basis, j_basis, i_spin
  integer :: i_k_task, i_k_point
  integer :: max_occ_number(n_spin)
  integer :: n_occ_states,i_state

  Omega_MO_diag = (0.0d0,0.0d0) 

  i_k_task = 0
  do i_k_point = 1,n_k_points, 1
  if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
        i_k_task = i_k_task + 1

     do i_spin = 1, n_spin, 1
         max_occ_number(i_spin) = 0
         do i_state = n_states, 1, -1
          !if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.0.d0) then
          ! The following line implies the system does not have any fractional occupation number
          if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.1.e-6) then
           max_occ_number(i_spin) = i_state
           exit
          endif
         enddo
     enddo
     n_occ_states=max_occ_number(1) 
  ! shanghui need to change this in the future for n_spin


     do i_state=1,n_occ_states
        do i_basis=1,n_basis
        do j_basis=1,n_basis
   
         Omega_MO_diag( i_state, i_k_task)=         &
         Omega_MO_diag( i_state, i_k_task) +        &
         dconjg(KS_eigenvector_complex(i_basis,i_state,1,i_k_task )) * & 
         r_complex(i_basis,j_basis, i_k_task) * & 
         KS_eigenvector_complex(j_basis,i_state,1,i_k_task )
        enddo
        enddo 

     enddo

  endif ! i_k_task   
  enddo ! i_k_point



end subroutine evaluate_Omega_MO_diag
!*****
