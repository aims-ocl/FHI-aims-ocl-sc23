!****s* FHI-aims/evaluate_Omega_MO_v2
!  NAME
!    evaluate_Omega_MO_v2
!  SYNOPSIS

subroutine evaluate_Omega_MO_v2( & 
           Omega_part_1_complex, &
           first_order_Q_k_complex, &
           S_complex,  &          
           Omega_MO)
               

!  PURPOSE
!    Omega_MO  = < i(k)|-r|j(k) >
!              = Omega_MO_part_1 + Omega_MO_part_2
!              = C^+ Omega_part_1_complex C) + C^+ S_complex C^{1}
!              = C^+ Omega_part_1_complex C) + C^+ S_complex C Q
!
!    This is the second version for Omega_MO by using k space method. 
! 
!                                                        shanghui,2017.01.09


!  USES

  use dimensions
  use runtime_choices, only : real_eigenvectors
  use mpi_tasks
  use physics, only : occ_numbers, KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue  
  use localorb_io, only: use_unit

!  ARGUMENTS

  implicit none
  

  complex*16, dimension( n_basis, n_basis, n_k_points_task), intent(IN) ::  Omega_part_1_complex
  complex*16, dimension( n_basis, n_basis, n_k_points_task), intent(IN) ::  first_order_Q_k_complex
  complex*16, dimension( n_basis, n_basis, n_k_points_task), intent(IN) ::  S_complex
  complex*16, dimension( n_basis, n_basis, n_k_points_task), intent(OUT) :: Omega_MO

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

  complex*16, allocatable     ::  temp_first_order(:,:)
  complex*16, allocatable     ::  temp_1(:,:)
  complex*16, allocatable     ::  temp_part_1(:,:)
  complex*16, allocatable     ::  temp_part_2(:,:)
  complex*16, allocatable     ::  temp_eigenvector(:,:)
  complex*16     ::  zero, one

  zero = (0.d0,0.d0)
  one = (1.d0,0.d0)


  allocate(temp_first_order(n_basis,n_basis))
  allocate(temp_1(n_basis,n_basis))
  allocate(temp_part_1(n_basis,n_basis))
  allocate(temp_part_2(n_basis,n_basis))
  allocate(temp_eigenvector(n_basis,n_basis))

  i_k_task = 0
  do i_k_point = 1,n_k_points, 1
  if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
        i_k_task = i_k_task + 1

!     if(i_k_point.eq.1) then
!        write(use_unit,*) 'Omega_part_1_complex', Omega_part_1_complex(:,:,1)
!     endif

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

!--------- part_2 = Ct* Omega_part_1_complex *C ----------------------------------------
    temp_first_order(1:n_basis,1:n_basis)=          &
    Omega_part_1_complex(1:n_basis,1:n_basis,i_k_task)

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
 
    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               one, temp_eigenvector, n_basis,temp_first_order, n_basis,&
               zero,temp_1,n_basis)

    temp_part_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)

    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                one,temp_1, n_basis,temp_eigenvector, n_basis,&
                zero,temp_part_1, n_basis)

!--------- part_2 = Ct* S_complex *C Q = Q ----------------------------------------
!  temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
!
!  CALL zgemm("C","N",n_basis,n_basis,n_basis,&
!             1.0d0, temp_eigenvector, n_basis, S_complex(:,:,i_k_task), n_basis,&
!             0.0d0,temp_1,n_basis)
!
!  temp_part_2(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
!
!  CALL zgemm("N","N",n_basis,n_basis,n_basis,&
!              1.0d0,temp_1, n_basis,temp_eigenvector, n_basis,&
!              0.0d0,temp_part_2, n_basis)
!
!     !  write(use_unit,*) 'S ', i_k_point,  S_complex(:,:,i_k_task)
!     !  write(use_unit,*) 'C ', i_k_point,  KS_eigenvector_complex(:,:,1,i_k_task) 
!     !  write(use_unit,*) 'Ct * S * C', i_k_point, temp_part_2
!

   if(i_k_point.eq.1) then
    write(use_unit,*) ' ' ! only when I add this output the complie -O3 at thnec could get right. 
    endif 

   do i_basis=1,n_occ_states
      do j_basis=n_occ_states+1, n_basis

      Omega_MO( i_basis,j_basis, i_k_task)=         &
       -temp_part_1(i_basis,j_basis) + first_order_Q_k_complex(i_basis,j_basis,i_k_task) 

      Omega_MO( j_basis,i_basis, i_k_task)= & 
       -temp_part_1(j_basis,i_basis) + first_order_Q_k_complex(j_basis,i_basis,i_k_task) 

      enddo
   enddo


   !write(use_unit,*) 'par1_1:',-temp_part_1
   !write(use_unit,*) 'part_2:',-first_order_Q_k_complex(:,:,i_k_task)

  endif ! i_k_task   
  enddo ! i_k_point



  deallocate(temp_first_order)
  deallocate(temp_1)
  deallocate(temp_part_1)
  deallocate(temp_part_2)
  deallocate(temp_eigenvector)

end subroutine evaluate_Omega_MO_v2
!*****
