!****s* FHI-aims/evaluate_Omega_MO_v1
!  NAME
!    evaluate_Omega_MO_v1
!  SYNOPSIS

subroutine evaluate_Omega_MO_v1( & 
           momentum_matrix, &
           Omega_MO)
               

!  PURPOSE
!    Omega_MO  = < i(k)|-r|j(k) > 
!              = (C^+ momentum_matrix_complex C)_ij/( Ei(k) - Ej(k) ) 
!
!    This is the first version for Omega_MO by using momentum_matrix_complex. 
! 
!                                                        shanghui,2017.01.09


!  USES

  use dimensions
  use runtime_choices, only : real_eigenvectors
  use mpi_tasks
  use physics, only : occ_numbers, KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue  

!  ARGUMENTS

  implicit none
  

  complex*16, dimension( n_basis, n_basis, n_k_points_task), intent(IN) :: momentum_matrix
  complex*16, dimension( n_states, n_states, n_k_points_task), intent(OUT) :: Omega_MO

!  INPUTS
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  Omega_MO


  integer :: i_state, j_state, i_spin
  integer :: i_k_task, i_k_point
  integer :: max_occ_number(n_spin)
  integer :: n_occ_states

  complex*16, allocatable     ::  temp_first_order(:,:)
  complex*16, allocatable     ::  temp_1(:,:)
  complex*16, allocatable     ::  temp_H(:,:)
  complex*16, allocatable     ::  temp_M(:,:)
  complex*16, allocatable     ::  temp_eigenvector(:,:)
  complex*16     ::  zero, one

  zero = (0.d0,0.d0)
  one = (1.d0,0.d0)


  allocate(temp_first_order(n_basis,n_basis))
  allocate(temp_1(n_basis,n_basis))
  allocate(temp_M(n_basis,n_basis))
  allocate(temp_eigenvector(n_basis,n_basis))

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


   if(real_eigenvectors) then
          temp_eigenvector(:,:)=          &
          cmplx(KS_eigenvector(:,:,1,i_k_task))
   else
          temp_eigenvector(:,:)=          &
          !dconjg(KS_eigenvector_complex(:,:,1,i_k_task))
          KS_eigenvector_complex(:,:,1,i_k_task)
   endif

 !--------- Ct*momentum*C ----------------------------------------
    temp_first_order(1:n_basis,1:n_basis)=          &
    momentum_matrix(1:n_basis,1:n_basis,i_k_task)

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
 
    CALL zgemm("C","N",n_states,n_basis,n_basis,&
               one, temp_eigenvector, n_basis,temp_first_order, n_basis,&
               zero,temp_1,n_states)

    temp_M(1:n_basis,1:n_basis)=(0.0d0,0.0d0)

    CALL zgemm("N","N",n_states,n_states,n_basis,&
                one,temp_1, n_states,temp_eigenvector, n_basis,&
                zero,temp_M, n_states)

    ! Nath to Hui: I don't think this comment is needed anymore, right ? I leave it here as a memory :D
    !if(i_k_point.eq.1) then
    !write(use_unit,*) ' ' ! only when I add this output the complie -O3 at thnec could get right. 
    !endif 

   do i_state=1,n_occ_states
      do j_state=n_occ_states+1, n_basis
!----------------------we do not need this term, we calculate here just for debug--------
     Omega_MO( i_state,j_state, i_k_task)=         &
      temp_M(i_state,j_state)/(KS_eigenvalue(i_state,1,i_k_point)-KS_eigenvalue(j_state,1,i_k_point)) 

     Omega_MO( j_state,i_state, i_k_task)= & 
      temp_M(j_state,i_state)/(KS_eigenvalue(j_state,1,i_k_point)-KS_eigenvalue(i_state,1,i_k_point))  

      enddo
   enddo

  endif ! i_k_task   
  enddo ! i_k_point


  deallocate(temp_first_order)
  deallocate(temp_1)
  deallocate(temp_M)
  deallocate(temp_eigenvector)

end subroutine evaluate_Omega_MO_v1
!*****
