!****s* FHI-aims/evaluate_first_order_DM_supercell_p1
!  NAME
!    evaluate_first_order_DM_supercell_p1
!  SYNOPSIS

subroutine evaluate_first_order_DM_supercell_p1(  &
           first_order_S_sparse,first_order_H_sparse,  & 
           first_order_density_matrix_sparse)
               

!  PURPOSE
!    calculate the first-order DM 

!  DM(1)=
!  shanghui,2012.05.02
!  USES

  use dimensions
  use mpi_tasks
  use runtime_choices, only : use_local_index, use_scalapack_DFPT_phonon
  use synchronize_mpi, only : sync_vector
  use physics,         only : n_electrons
  use pbc_lists,       only : n_cells_in_sc_DFPT
  use scalapack_wrapper, only : construct_first_order_overlap_supercell_scalapack, & 
                                construct_first_order_hamiltonian_supercell_scalapack, & 
                                construct_first_order_dm_supercell_scalapack, & 
                                evaluate_first_order_U_supercell_scalapack, & 
                                get_first_order_dm_sparse_matrix_from_supercell_scalapack 
  use DFPT_phonon_supercell, only : KS_eigenvalue_supercell, KS_eigenvector_supercell 
  use localorb_io, only: use_unit

!  ARGUMENTS

  implicit none
  

  real*8, dimension(3, n_centers_in_sc_DFPT, n_hamiltonian_matrix_size), intent(IN) :: first_order_S_sparse
  real*8, dimension(3, n_centers_in_sc_DFPT, n_hamiltonian_matrix_size), intent(IN) :: first_order_H_sparse
  real*8, dimension(3, n_centers_in_sc_DFPT, n_hamiltonian_matrix_size),intent(INOUT) :: first_order_density_matrix_sparse

!  INPUTS

!
!  OUTPUT
!  first_order_density_matrix


  integer :: i_basis,j_basis ,k_coord, k_atom
  integer :: i_state,j_state,n_occ_states

  real*8, allocatable :: temp_first_order(:,:) 
  real*8, allocatable :: temp_1(:,:) 
  real*8, allocatable :: temp_S(:,:) 
  real*8, allocatable :: temp_H(:,:) 
  real*8, allocatable :: first_order_U(:,:) 
  real*8, allocatable :: KS_eigen_times_occnumber(:,:) 




 !--------------initialize----------------------
  first_order_density_matrix_sparse(1:3,1:n_centers_in_sc_DFPT,1:n_hamiltonian_matrix_size) = 0.0d0

  n_occ_states = NINT(n_electrons * n_cells_in_sc_DFPT / 2.0d0)  ! n_basis_sc_DFPT / 2 ! todo : extend
 
   if(mod(nint(n_electrons),2).ne.0) then
      call aims_stop('n_electrons is odd, please check', 'evaluate_first_order_DM_supercell_p1')
   endif

  if(myid.eq.0) then 
  write(use_unit,*) 'n_occ_states in DM1:', n_occ_states 
  endif

 if(use_scalapack_DFPT_phonon)then

     do k_coord=1, 3
     do k_atom = 1,n_atoms, 1

        !-----(1) get S1(k_coord,k_atom, n_basis_sc_DFPT,n_basis_sc_DFPT)------------
        call construct_first_order_overlap_supercell_scalapack( k_coord,k_atom,first_order_S_sparse )    
 
        !-----(2) get H1(k_coord,k_atom, n_basis_sc_DFPT,n_basis_sc_DFPT)------------
        call construct_first_order_hamiltonian_supercell_scalapack( k_coord,k_atom,first_order_H_sparse )
  
        !-----(3) get U1(k_coord,k_atom, n_basis_sc_DFPT,n_basis_sc_DFPT)------------
        call evaluate_first_order_U_supercell_scalapack() 

        !-----(4) calcualte first_order_dm (k_coord,k_atom,n_basis_sc_DFPT,n_basis_sc_DFPT)-----------
        ! DM1 =      C(C^+ S1 C)C^+ occ_number  
        !          + (CU1) C^+ occ_number 
        !          +  C   (CU1)^+ occ_number   
        call construct_first_order_dm_supercell_scalapack()

        !-----(6) calcualte first_order_density_matrix_sparse(k_coord,k_centers,n_hamiltonian_matrix_size)---
        call get_first_order_dm_sparse_matrix_from_supercell_scalapack & 
             (k_coord, k_atom, first_order_density_matrix_sparse)
            

     enddo !k_atom 

       !-----(7) sync ------------------------------   
        call sync_vector(first_order_density_matrix_sparse(k_coord,:,:),n_centers_in_sc_DFPT*n_hamiltonian_matrix_size)


     enddo !k_coord 


 else

  allocate(temp_first_order(n_basis_sc_DFPT,n_basis_sc_DFPT))
  allocate(temp_1(n_basis_sc_DFPT,n_basis_sc_DFPT))
  allocate(temp_S(n_basis_sc_DFPT,n_basis_sc_DFPT))
  allocate(temp_H(n_basis_sc_DFPT,n_basis_sc_DFPT))
  allocate(first_order_U(n_basis_sc_DFPT,n_basis_sc_DFPT))
  allocate(KS_eigen_times_occnumber(n_basis_sc_DFPT,n_basis_sc_DFPT))


!-----------shanghui test DM0 supercell-------------------------
! KS_eigen_times_occnumber = 0.0d0
! do i_state=1, n_occ_states
!    KS_eigen_times_occnumber(1:n_basis_sc_DFPT,i_state)=   &
!    KS_eigenvector(1:n_basis_sc_DFPT,i_state)*dsqrt(2.0d0)
! enddo
!
! temp_1 = 0.0d0
! call dsyrk('U', 'N', n_basis_sc_DFPT, n_occ_states, 1.d0, KS_eigen_times_occnumber, n_basis_sc_DFPT, &
!                    0.d0, temp_1, n_basis_sc_DFPT)
!
! temp_1 = temp_1 + transpose(temp_1)
! do i_basis = 1, n_basis_sc_DFPT
!    temp_1(i_basis, i_basis) = temp_1(i_basis, i_basis)*0.5d0
! end do
!
!
! call add_matrix_to_sparse(1, temp_1, first_order_density_matrix_sparse(1,:,:))
!
! if(myid.eq.0) then
!    write(use_unit,*) '===============(2) DM0 from KS_eigen =================' 
!   !-----matrix----------- 
!   do i_basis=1,n_basis_sc_DFPT
!    write(use_unit,'(90f9.5)') (temp_1(i_basis,j_basis),j_basis=1,n_basis_sc_DFPT)
!   enddo 
!
!   !----sparse-----------
!   !write(use_unit,'(418f20.15)') (first_order_density_matrix_sparse(1,1,i_basis), & 
!   !                        i_basis=1,n_hamiltonian_matrix_size-1)
! endif
!-----------shanghui end test DM0 supercell-------------------------

  do k_coord=1, 3 
  do k_atom = 1,n_atoms, 1

!----------------the first term-------------------------------
  !temp_S=Ct S(1) C----->

  call trans_first_order_sparse_to_matrix(k_coord,k_atom,first_order_S_sparse,temp_first_order) 

  temp_1 = 0.0d0
  CALL dgemm("T","N",n_basis_sc_DFPT,n_basis_sc_DFPT,n_basis_sc_DFPT,&
             1.0d0, KS_eigenvector_supercell, n_basis_sc_DFPT,temp_first_order, n_basis_sc_DFPT,&
             0.0d0,temp_1,n_basis_sc_DFPT)
  temp_S = 0.0d0
  CALL dgemm("N","N",n_basis_sc_DFPT,n_basis_sc_DFPT,n_basis_sc_DFPT,&
              1.0d0,temp_1, n_basis_sc_DFPT,KS_eigenvector_supercell, n_basis_sc_DFPT,&
              0.0d0,temp_S, n_basis_sc_DFPT)


  KS_eigen_times_occnumber = 0.0d0
  do i_state=1, n_occ_states
     KS_eigen_times_occnumber(1:n_basis_sc_DFPT,i_state)=   &
     KS_eigenvector_supercell(1:n_basis_sc_DFPT,i_state)*2.0d0  !*occ_numbers(i_state,1,1)
  enddo

  temp_first_order = 0.0d0
  do i_state=1, n_occ_states
  do j_state=1, n_occ_states
     temp_first_order(i_state,j_state) = - temp_S(i_state,j_state)
  enddo
  enddo


  temp_1 = 0.0d0
  CALL dgemm("N","N",n_basis_sc_DFPT,n_basis_sc_DFPT,n_basis_sc_DFPT,&
             1.0d0, KS_eigen_times_occnumber, n_basis_sc_DFPT,temp_first_order, n_basis_sc_DFPT,&
             0.0d0,temp_1,n_basis_sc_DFPT)
  temp_first_order = 0.0d0
  CALL dgemm("N","T",n_basis_sc_DFPT,n_basis_sc_DFPT,n_basis_sc_DFPT,&
              1.0d0,temp_1, n_basis_sc_DFPT,KS_eigenvector_supercell, n_basis_sc_DFPT,&
              0.0d0,temp_first_order, n_basis_sc_DFPT)



  call add_matrix_to_sparse(k_atom, temp_first_order, first_order_density_matrix_sparse(k_coord,:,:)) 


  !-----------------begin old loop code------------------------------------- 
  ! temp_1 = Cui * ( temp_S)ij * Cvj----->
  !do i_basis=1,n_basis_sc_DFPT
  !do j_basis=1,n_basis_sc_DFPT
  !   temp_1(i_basis,j_basis)=0.0d0
  !   do i_state=1, n_occ_states
  !   do j_state=1, n_occ_states
 
  !     temp_1(i_basis,j_basis)= temp_1(i_basis,j_basis) &
  !     - 2.0d0*KS_eigenvector(i_basis,i_state) &
  !     * temp_S(i_state,j_state)              &
  !     * KS_eigenvector(j_basis,j_state) 
 
  !   enddo
  !   enddo

  !enddo
  !enddo
  !call add_matrix_to_sparse(temp_1, first_order_density_matrix_sparse(i_coord,i_center,:)) 
  !-----------------end old loop code------------------------------------- 
!-----------------------end the first term----------------------------

  
!-----------------------the second term-------------------------------
  !(1) temp_S = temp_S *E= Ct *S1 *C*E=  Cij=Aij*Bjj
  do i_basis=1, n_basis_sc_DFPT
     temp_S(1:n_basis_sc_DFPT,i_basis)=    &
     temp_S(1:n_basis_sc_DFPT,i_basis)*KS_eigenvalue_supercell(i_basis)
  enddo


  !(2) temp_H = Ct*H1*C 
  call trans_first_order_sparse_to_matrix(k_coord,k_atom,first_order_H_sparse,temp_first_order) 

  temp_1 = 0.0d0
  CALL dgemm("T","N",n_basis_sc_DFPT,n_basis_sc_DFPT,n_basis_sc_DFPT,&
             1.0d0, KS_eigenvector_supercell, n_basis_sc_DFPT,temp_first_order, n_basis_sc_DFPT,&
             0.0d0,temp_1,n_basis_sc_DFPT)
  temp_H = 0.0d0
  CALL dgemm("N","N",n_basis_sc_DFPT,n_basis_sc_DFPT,n_basis_sc_DFPT,&
              1.0d0,temp_1, n_basis_sc_DFPT,KS_eigenvector_supercell, n_basis_sc_DFPT,&
              0.0d0,temp_H, n_basis_sc_DFPT)


  !(3)  first_order_U(a,i) = (Ct *S1 *C*E - Ct*H1*C)ai/(Ea-Ei) 
   first_order_U = 0.0d0 
   do i_basis = 1,n_occ_states
      do j_basis = n_occ_states+1, n_basis_sc_DFPT
!----------------in fact, we only neet U(unocc,occ)=U(a,i) like this:----------
     first_order_U(j_basis,i_basis)= &
    (temp_S(j_basis,i_basis)-temp_H(j_basis,i_basis))/ &
    (KS_eigenvalue_supercell(j_basis)-KS_eigenvalue_supercell(i_basis))
      enddo
   enddo

  temp_1=0.0d0
  CALL dgemm("N","N",n_basis_sc_DFPT,n_basis_sc_DFPT,n_basis_sc_DFPT,&
             1.0d0, KS_eigenvector_supercell, n_basis_sc_DFPT,first_order_U, n_basis_sc_DFPT,&
             0.0d0,temp_1,n_basis_sc_DFPT)
  temp_first_order=0.0d0
  CALL dgemm("N","T",n_basis_sc_DFPT,n_basis_sc_DFPT,n_basis_sc_DFPT,&
              1.0d0,temp_1, n_basis_sc_DFPT,KS_eigen_times_occnumber, n_basis_sc_DFPT,&
              0.0d0,temp_first_order, n_basis_sc_DFPT)
  call add_matrix_to_sparse(k_atom,temp_first_order, first_order_density_matrix_sparse(k_coord,:,:)) 

  temp_1=0.0d0
  CALL dgemm("N","T",n_basis_sc_DFPT,n_basis_sc_DFPT,n_basis_sc_DFPT,&
             1.0d0, KS_eigen_times_occnumber, n_basis_sc_DFPT,first_order_U, n_basis_sc_DFPT,&
             0.0d0,temp_1,n_basis_sc_DFPT)
  temp_first_order=0.0d0
  CALL dgemm("N","T",n_basis_sc_DFPT,n_basis_sc_DFPT,n_basis_sc_DFPT,&
              1.0d0,temp_1, n_basis_sc_DFPT,KS_eigenvector_supercell, n_basis_sc_DFPT,&
              0.0d0,temp_first_order, n_basis_sc_DFPT)
  call add_matrix_to_sparse(k_atom,temp_first_order, first_order_density_matrix_sparse(k_coord,:,:)) 



!-----------------begin old loop code------------------------------------- 
! do i_basis=1,n_basis_sc_DFPT
! do j_basis=1,n_basis_sc_DFPT
!
!     temp_first_order(i_basis,j_basis)=0.0d0
!
!     do i_state = 1, n_occ_states
!     do j_state = n_occ_states+1,n_basis_sc_DFPT
!       temp_first_order(i_basis,j_basis)= temp_first_order(i_basis,j_basis) &
!       + 2.0d0*KS_eigenvector(i_basis,i_state) &
!       * first_order_U( j_state,i_state )                  &
!       * KS_eigenvector(j_basis,j_state)   &  
!       + 2.0d0*KS_eigenvector(i_basis,j_state) &
!       * first_order_U(j_state,i_state )                  &
!       *  KS_eigenvector(j_basis,i_state) 
!     enddo
!     enddo 
!
! enddo
! enddo
!
! call add_matrix_to_sparse(temp_first_order, first_order_density_matrix_sparse(i_coord,i_center,:)) 
!-----------------end old loop code------------------------------------- 
!-----------------------end the second term----------------------------


   enddo ! k_atom

   enddo ! k_coord


  deallocate(temp_first_order)
  deallocate(temp_1)
  deallocate(temp_S)
  deallocate(temp_H)
  deallocate(first_order_U)
  deallocate(KS_eigen_times_occnumber)

 endif ! scalapack 

end subroutine evaluate_first_order_DM_supercell_p1
!******
