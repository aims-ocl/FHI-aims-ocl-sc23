!****s* FHI-aims/evaluate_first_order_DM_polarizability
!  NAME
!    evaluate_first_order_DM_polarizability
!  SYNOPSIS

subroutine evaluate_first_order_DM_polarizability(   &
           KS_eigenvector, KS_eigenvalue, occ_numbers, max_occ_number,  &
           first_order_U,first_order_density_matrix)
               

!  PURPOSE
!    calculate the first-order DM for polarizability 

!  shanghui,2013.12.12
!  USES

  use dimensions
  use runtime_choices
  use synchronize_mpi
  use scalapack_wrapper

!  ARGUMENTS

  implicit none
  

  real*8, dimension(n_basis, n_states, n_spin,n_k_points), intent(IN) :: KS_eigenvector
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers
  integer, dimension(n_spin), intent(IN) :: max_occ_number
  real*8, dimension(3, n_states, n_states, n_spin), intent(IN) :: first_order_U
  real*8, dimension(3, n_basis, n_basis, n_spin),intent(INOUT) :: first_order_density_matrix

!  INPUTS
!   o KS_eigenvector -- Kohn-Sham eigenvectors (real format)
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_density_matrix


  integer :: i_basis,j_basis ,i_coord, i_spin
  integer :: i_state,j_state,n_occ_states

  real*8, allocatable :: temp_first_order(:,:)
  real*8, allocatable :: temp_1(:,:)
  real*8, allocatable :: temp_2(:,:)
  real*8, allocatable :: KS_eigen_times_occnumber(:,:)

  allocate(temp_first_order(n_states,n_states))
  allocate(temp_1(n_basis,n_states))
  allocate(temp_2(n_basis,n_basis))
  allocate(KS_eigen_times_occnumber(n_basis,n_states))


if (use_scalapack) then

   do i_spin = 1, n_spin

      ! i)  Construct the first-order density matrix (dense, scalapack format)
      call construct_first_order_dm_polar_scalapack(occ_numbers, i_spin) 

     ! ii) Put the density matrix back in natural form (n_basis dim) so that it could be used in the following
      call get_first_order_dm_polar_scalapack( first_order_density_matrix, i_spin )

   end do ! i_spin

   ! iii) Synchronize
   call sync_vector(first_order_density_matrix,3*n_basis*n_basis*n_spin)

else ! lapack version

do i_spin = 1, n_spin

  n_occ_states=max_occ_number(i_spin) 

  do i_state=1, n_states
     KS_eigen_times_occnumber(1:n_basis,i_state)=   &
     KS_eigenvector(1:n_basis,i_state,i_spin,1)*occ_numbers(i_state,i_spin,1)
  enddo



  do i_coord=1,3 
    
     first_order_density_matrix(i_coord,:,:,i_spin)=0.0d0

!--------------------DM_ov+ DM_vo--------------------------
  temp_first_order = 0.0d0
  do i_state=1, n_occ_states
  do j_state=n_occ_states+1,n_states
     temp_first_order(j_state,i_state) = first_order_U(i_coord, j_state,i_state, i_spin) 
  enddo 
  enddo  

  temp_1=0.0d0
  temp_2=0.0d0 
  CALL dgemm("N","N",n_basis,n_states,n_states,&
             1.0d0, KS_eigenvector(:,:,i_spin,1), n_basis,temp_first_order, n_states,&
             0.0d0,temp_1,n_basis)
  CALL dgemm("N","T",n_basis,n_basis,n_states,&
              1.0d0,temp_1, n_basis,KS_eigen_times_occnumber, n_basis,&
              0.0d0,temp_2, n_basis)
  
  first_order_density_matrix(i_coord,:,:,i_spin) = & 
  first_order_density_matrix(i_coord,:,:,i_spin) + temp_2(:,:) 

  temp_1=0.0d0
  temp_2=0.0d0 
  CALL dgemm("N","T",n_basis,n_states,n_states,&
             1.0d0, KS_eigen_times_occnumber, n_basis,temp_first_order, n_states,&
             0.0d0,temp_1,n_basis)
  CALL dgemm("N","T",n_basis,n_basis,n_states,&
              1.0d0,temp_1, n_basis,KS_eigenvector(:,:,i_spin,1), n_basis,&
              0.0d0,temp_2, n_basis)

  first_order_density_matrix(i_coord,:,:,i_spin) = & 
  first_order_density_matrix(i_coord,:,:,i_spin) + temp_2(:,:) 

!-----------------------the loop version (slow)------------------------------
! do i_basis=1,n_basis
!    do j_basis=1,n_basis
!     first_order_density_matrix(i_coord,i_basis,j_basis)=0.0d0
!
!     do i_state=1, n_occ_states
!     do j_state=n_occ_states+1,n_states
!     first_order_density_matrix(i_coord,i_basis,j_basis) = &
!     first_order_density_matrix(i_coord,i_basis,j_basis) + &
!         KS_eigen_times_occnumber(i_basis,i_state) &
!       * first_order_U(i_coord, j_state,i_state )                  &
!       * KS_eigenvector(j_basis,j_state,1,1) + KS_eigenvector(i_basis,j_state,1,1) &
!       * first_order_U(i_coord, j_state,i_state )                  &
!       *  KS_eigen_times_occnumber(j_basis,i_state) 
!     enddo
!     enddo 
!
!     enddo
! enddo
!-----------------------end loop version (slow)----------------------------


   enddo ! i_coord
enddo ! i_spin

end if ! lapack/scalapack

   deallocate(temp_first_order)
   deallocate(temp_1)
   deallocate(temp_2)
   deallocate(KS_eigen_times_occnumber)

end subroutine evaluate_first_order_DM_polarizability
!******
