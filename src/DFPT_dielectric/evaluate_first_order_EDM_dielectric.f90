!****s* FHI-aims/evaluate_first_order_EDM_dielectric
!  NAME
!    evaluate_first_order_EDM_dielectric
!  SYNOPSIS

subroutine evaluate_first_order_EDM_dielectric( &
           KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers,  &
           Omega_MO_diag, &
           first_order_H_complex,   &
           first_order_U_complex,first_order_energy_density_matrix_sparse)
               

!  PURPOSE
!    calculate the first-order DM sparse for PBC case 
!    shanghui, 2015.07.30

!  USES

  use dimensions
  use localorb_io, only: use_unit
  use pbc_lists
  use synchronize_mpi

!  ARGUMENTS

  implicit none
  
  real*8, dimension(n_basis, n_states, n_spin,n_k_points_task),  intent(IN) :: KS_eigenvector
  complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task),  intent(IN) :: KS_eigenvector_complex
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers

  complex*16, dimension( n_basis, n_k_points_task), intent(IN) :: Omega_MO_diag
  complex*16, dimension( n_basis, n_basis, n_k_points_task), intent(IN) :: first_order_H_complex
  complex*16, dimension( n_basis, n_basis, n_k_points_task), intent(IN) :: first_order_U_complex
  real*8, dimension(n_hamiltonian_matrix_size),intent(INOUT) :: first_order_energy_density_matrix_sparse

!  INPUTS
!   o KS_eigenvector -- Kohn-Sham eigenvectors (real format)
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates
!   o KS_eigenvalue -- Kohn-Sham eigenvalues

!
!  OUTPUT
!  first_order_density_matrix_sparse


  integer :: i_basis,j_basis ,i_spin , i_cell, i_index
  integer :: i_state,j_state,i_k_point,  i_k_task
  integer, dimension(n_spin, n_k_points_task)   :: max_occ_number

  complex*16, allocatable     ::  temp_first_order(:,:)
  complex*16,  allocatable    ::  temp_1(:,:)
  complex*16, allocatable     ::  temp_complex(:,:) 
  complex*16, allocatable     ::  temp_H(:,:)
  complex*16, allocatable     ::  temp_eigenvector(:,:)
  complex*16     ::  zero, one

  zero = (0.d0,0.d0)
  one = (1.d0,0.d0)

  allocate(temp_first_order(n_basis,n_basis)) 
  allocate(temp_1(n_basis,n_basis)) 
  allocate(temp_complex(n_basis,n_basis))
  allocate(temp_H(n_basis,n_basis))
  allocate(temp_eigenvector(n_basis,n_basis))

!-----------------initialize----------------------
  first_order_energy_density_matrix_sparse(1:n_hamiltonian_matrix_size)= 0.0d0

! find the max_occ_number
  i_k_task = 0
  do i_k_point = 1, n_k_points,1

     if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
        i_k_task = i_k_task + 1
        do i_spin = 1, n_spin, 1
           max_occ_number(i_spin,i_k_task) = 0

          do i_state = n_states, 1, -1
             !if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.0.d0) then
             ! The following line implies the system does not have any fractional occupation number
             if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.1.e-6) then
              max_occ_number(i_spin,i_k_task) = i_state
             exit
             endif
          enddo

        enddo ! i_spin
     endif
  enddo ! i_k_point


!-----here EDM^{0}_uv = Ei Ciu_aims Cvi_aims^(*)-------------- 
!-----here EDM^{1}_uv = Ei^{1} Ciu_aims Cvi_aims^(*)-------------- 
!                     + Ei Ciu_aims^{1} Cvi_aims^(*)
!                     + Ei Ciu_aims Cvi_aims^{1*}
!==============================(2) complex_eigenvectors==============================
   i_k_task = 0
   do i_k_point = 1, n_k_points,1
   if (myid.eq.MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
       i_k_task = i_k_task + 1

 !---------begin  Ct*H(1)*C ----------------------------------------
    temp_eigenvector(1:n_basis,1:n_basis)= & 
    KS_eigenvector_complex(1:n_basis,1:n_basis,1,i_k_task)

    temp_first_order(1:n_basis,1:n_basis)=          &
    first_order_H_complex(1:n_basis,1:n_basis, i_k_task)

    temp_1(1:n_basis,1:n_basis)=(0.0d0,0.0d0)
 
    CALL zgemm("C","N",n_basis,n_basis,n_basis,&
               one, temp_eigenvector, n_basis,temp_first_order, n_basis,&
               zero,temp_1,n_basis)

    temp_H(1:n_basis,1:n_basis)=(0.0d0,0.0d0)

    CALL zgemm("N","N",n_basis,n_basis,n_basis,&
                one,temp_1, n_basis,temp_eigenvector, n_basis,&
                zero,temp_H, n_basis)


    !shanghui note:  only when I add this output the complie -O3 
    !                at thnec/draco could get rid of NA error. 
    !                The same as evaluate_first_order_U_dielectric.f90
    if(i_k_point.eq.1) then
    write(use_unit,*) ' '  
    endif
 !---------end  Ct*H(1)*C ----------------------------------------

     do i_basis=1,n_basis
     do j_basis=1,n_basis
         temp_complex(i_basis,j_basis)=(0.0d0,0.0d0)

         !Ei^{1} Ciu_aims Cvi_aims^(*)
         do i_state=1, max_occ_number(1,i_k_task)
          temp_complex(i_basis,j_basis)= temp_complex(i_basis,j_basis) &
           + occ_numbers(i_state,1,i_k_point) & 
         !shanghui note 1: currently, Omega_MO_diag is set to zero. If we want this,
         !               term, the Berry phase need to be calcualted.
         !shanghui note 2: I use r_sparse to get Omega_MO_diag, but the result
         !               is not as good as set to zero.   
         !*(Omega_MO_diag(i_state,i_k_task)+temp_H(i_state,i_state)) &
           *( temp_H(i_state,i_state) ) &
           * KS_eigenvector_complex(i_basis,i_state,1,i_k_task) & 
           * dconjg(KS_eigenvector_complex(j_basis,i_state,1,i_k_task )) 
         enddo 


         !Ei Ciu_aims^{1} Cvi_aims^(*) + Ei Ciu_aims Cvi_aims^{1*}
         do i_state=1, max_occ_number(1,i_k_task)
         do j_state=max_occ_number(1,i_k_task)+1,n_states

          temp_complex(i_basis,j_basis)= temp_complex(i_basis,j_basis) &
           + occ_numbers(i_state,1,i_k_point)*KS_eigenvalue(i_state,1,i_k_point) &
           * KS_eigenvector_complex(i_basis,j_state,1,i_k_task)    &
           * first_order_U_complex( j_state,i_state,i_k_task)       &
           * dconjg(KS_eigenvector_complex(j_basis,i_state,1,i_k_task ))           &
           + occ_numbers(i_state,1,i_k_point)*KS_eigenvalue(i_state,1,i_k_point)  &
           * KS_eigenvector_complex(i_basis,i_state,1,i_k_task)   &
           * dconjg(KS_eigenvector_complex(j_basis,j_state,1,i_k_task ))           &
           * dconjg(first_order_U_complex( j_state,i_state,i_k_task))

         enddo
         enddo

     enddo ! i_basis
     enddo ! j_basis


     !-----------------add k_phase ----------------
     do i_basis = 1, n_basis
     do i_cell = 1,n_cells_in_hamiltonian-1

             if (index_hamiltonian(1,i_cell, i_basis) > 0) then
             do i_index = index_hamiltonian(1, i_cell, i_basis), &
                &         index_hamiltonian(2, i_cell, i_basis)
                   j_basis =  column_index_hamiltonian(i_index)

                  first_order_energy_density_matrix_sparse(i_index)= &
                  first_order_energy_density_matrix_sparse(i_index)  &
                          + dble( (temp_complex(j_basis,i_basis) ) &
                           * dconjg(k_phase(i_cell,i_k_point) ) ) &  
                           * k_weights(i_k_point)

             end do
             end if

     end do ! i_cell 
     end do ! i_basis

    endif ! i_k_task
    enddo ! i_k_point




  !-------shanghui begin parallel------
   call sync_sparse_matrix(first_order_energy_density_matrix_sparse) 
  !-------shanghui end parallel------
 
  deallocate(temp_first_order) 
  deallocate(temp_1) 
  deallocate(temp_complex)
  deallocate(temp_H)
  deallocate(temp_eigenvector)

end subroutine evaluate_first_order_EDM_dielectric
!******
