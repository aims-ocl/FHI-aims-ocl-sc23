!****s* FHI-aims/evaluate_zero_order_DM_phonon_reduce_memory
!  NAME
!    evaluate_zero_order_DM_phonon_reduce_memory
!  SYNOPSIS

subroutine evaluate_zero_order_DM_phonon_reduce_memory &
           (KS_eigenvector, KS_eigenvector_complex, occ_numbers,density_matrix_sparse)
               

!  PURPOSE
!    calculate the zero-order DM for PBC case 

!  DM(uR1, vR2)(0)=sum{k} sum{i} Cui*(k) Cvi(k) exp(-ik(R2-R1))
!  shanghui,2013.12.30

!  move to phonon_reduce_memory
!  shanghui, 2015.07.30
!  USES

  use dimensions
  use pbc_lists
  use synchronize_mpi
  use scalapack_wrapper
  use physics, only : overlap_matrix, hamiltonian, n_electrons, &
                      KS_eigenvalue, &
                      chemical_potential, chemical_potential_spin 
  use runtime_choices

!  ARGUMENTS

  implicit none
  
  real*8,     dimension(n_basis, n_states, n_spin,n_k_points_task),  intent(IN) :: KS_eigenvector
  complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task),  intent(IN) :: KS_eigenvector_complex
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers

  real*8, dimension(n_hamiltonian_matrix_size),intent(INOUT) :: density_matrix_sparse

!  INPUTS
!   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
!   o occ_numbers -- occupations of eigenstates

!
!  OUTPUT
!  density_matrix


  integer :: i_basis,j_basis ,i_coord, i_atom,i_spin,i_index,i_cell
  integer :: i_state,j_state,i_k_point, i_k_task
  integer, dimension(n_spin, n_k_points_task)   :: max_occ_number

  real*8 :: occu
  integer :: n_max_occupied,info

  real*8, dimension(n_basis,n_basis) :: temp
  real*8,  dimension(n_basis, n_states) :: KS_scaled
  complex*16 , dimension(n_basis, n_basis)  :: temp_cmplx
  complex*16,  dimension(n_basis, n_states) :: KS_scaled_cmplx

 

  select case(packed_matrix_format)
  case(PM_none)   
    call aims_stop('DFPT_phonon_reduce_memory only use sparse matrix', 'evaluate_zero_order_DM_phonon_reduce_memory')
   
  
  case(PM_index)

    density_matrix_sparse(1:n_hamiltonian_matrix_size) = 0.0d0

    if(use_scalapack)then  
      !---------------------start scalapack version for DM0-----------------------
      do i_spin = 1, n_spin 
         !-------note, this occ_numbers have been with k_weights at cpscf_solver.f90
         !-----(1) calcualte first_order_dm at every k point------------
         ! DM0(my_k_point) = CC^{+} occ_number  
         call construct_dm_scalapack(occ_numbers, i_spin)
    
         !-----(2) add the k point phase ------------------------------
         ! DM0_sparse =  [ DM0(my_k_point) * exp(-ikR) ]
         call get_sparse_matrix_scalapack(density_matrix_sparse,i_spin)

         !-----(3) sum over {all k points and all matrix} at every core------
         call sync_density_matrix_sparse(density_matrix_sparse)
      enddo 

    else  
    !---------------------start lapack version for DM0-----------------------
    i_k_task = 0
    do i_k_point = 1, n_k_points,1

      if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
       i_k_task = i_k_task + 1

        do i_spin = 1, n_spin, 1
            max_occ_number(i_spin,i_k_task) = 0
            do i_state = n_states, 1, -1
             if (dabs(occ_numbers(i_state,i_spin,i_k_point)).gt.0.d0) then
              max_occ_number(i_spin,i_k_task) = i_state
              exit
             endif
            enddo
        enddo

      
       if (real_eigenvectors) then
          temp = 0.0d0
       else 
          temp_cmplx = (0.d0, 0.d0)
       endif


!---------------------------------(1) use lapacek------------------------------------
!      n_max_occupied = 1
!      do i_state = 1, n_states, 1
!         occu =  occ_numbers(i_state, 1, i_k_point)
!         if (occu .gt. 0.d0) then
!            n_max_occupied = i_state
!            if (real_eigenvectors) then
!               KS_scaled(:, i_state) = KS_eigenvector(:, i_state, 1, i_k_task) * sqrt(occu)
!            else 
!               KS_scaled_cmplx(:, i_state) = KS_eigenvector_complex(:, i_state, 1, i_k_task) &
!                                            * sqrt(occu)
!            endif
!         endif
!      end do

!      !-----here DM = Cui^(*) Cvi = Ci_aims * Ci_aims^(diagger)  the same as PM_none-----  
!       if (real_eigenvectors) then
!          call dsyrk('U', 'N', n_basis, n_max_occupied, 1.d0, KS_scaled, n_basis, &
!          &          0.d0, temp, n_basis)
!       else 
!          call zherk('U', 'N', n_basis, n_max_occupied, 1.d0, KS_scaled_cmplx, n_basis, &
!          &          0.d0, temp_cmplx, n_basis)
!       endif
! 
!      !----zherk only give upder part, so we need add the down park
!       if (real_eigenvectors) then
!          temp = temp + transpose(temp)
!          do i_basis = 1, n_basis
!             temp(i_basis, i_basis) = temp(i_basis, i_basis)*0.5d0
!          end do
!       else 
!          temp_cmplx = temp_cmplx + transpose(dconjg(temp_cmplx))
!          do i_basis = 1, n_basis
!             temp_cmplx(i_basis, i_basis) = temp_cmplx(i_basis, i_basis)*0.5d0
!          end do
!       endif
!---------------------------------(1) end use lapacek------------------------------------


!---------------------------------(2) use loop------------------------------------
!------------------  DM = Cui^(*) Cvi= Cui_aims * Cvi_aims^(*)---------------- 
     do i_basis=1,n_basis
       do j_basis=1,n_basis
        do i_state=1,  max_occ_number(1,i_k_task)

         if (real_eigenvectors) then
            temp(i_basis,j_basis)=temp(i_basis,j_basis)       &
          + KS_eigenvector(i_basis,i_state,1,i_k_task) & 
          * occ_numbers(i_state,1,i_k_point)                            &
          * KS_eigenvector(j_basis,i_state,1,i_k_task )            
         else    
            temp_cmplx(i_basis,j_basis)=temp_cmplx(i_basis,j_basis)       &
          + KS_eigenvector_complex(i_basis,i_state,1,i_k_task) & 
          * occ_numbers(i_state,1,i_k_point)                            &
          * dconjg(KS_eigenvector_complex(j_basis,i_state,1,i_k_task ))              
         endif

        enddo
       enddo
     enddo  
!----------------------------(2) end loop---------------------------------------
  

      do i_basis = 1, n_basis
          do i_cell = 1,n_cells_in_hamiltonian-1
             if (index_hamiltonian(1,i_cell, i_basis) > 0) then
                do i_index = index_hamiltonian(1, i_cell, i_basis), &
                &            index_hamiltonian(2, i_cell, i_basis)
                   j_basis =  column_index_hamiltonian(i_index)

                if (real_eigenvectors) then
                  density_matrix_sparse(i_index) &
                  & = density_matrix_sparse(i_index) &
                  & + temp(i_basis, j_basis) *dble(k_phase(i_cell,i_k_point)) &
                    * k_weights(i_k_point)
                else
                  density_matrix_sparse(i_index) &
                  & = density_matrix_sparse(i_index) &
                  & + dble(temp_cmplx(i_basis, j_basis)*(k_phase(i_cell,i_k_point))) &
                    * k_weights(i_k_point)
                  ! DM(uR,v) = temp(u,v)* exp(-ik(-R))
                endif

                

                end do
             end if
          end do
       end do


      endif 
    enddo ! i_k_point

    call sync_sparse_matrix(density_matrix_sparse) 

    endif !   use_scalapack

   
   ! if(myid.eq.0) then 
   ! write(use_unit,*) '================================(04) write DM in evaluate_zero_order_DM_phonon_reduce_memory '
   !  write(use_unit,*) density_matrix_sparse(:) 
   ! endif 
   !  do i_cell = 1,n_cells_in_hamiltonian-1
   !  do i_basis = 1, n_basis
   !  do j_basis = 1,n_basis

   !  !write(use_unit,*) i_cell,i_basis,j_basis, density_matrix_sparse_full(i_cell,i_basis,j_basis)

   !   end do
   !   end do
   !   end do


    case default

       call aims_stop('Invalid packing type', 'evaluate_zero_order_DM_phonon_reduce_memory')

    end select


end subroutine evaluate_zero_order_DM_phonon_reduce_memory
!******
