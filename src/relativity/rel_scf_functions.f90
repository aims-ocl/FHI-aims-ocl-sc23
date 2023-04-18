
 subroutine evaluate_free_atom_sums_rel( i_r, dir_tab, pce, pce_gradient, &
   n_atom_list, atom_list, n_compute_atoms, i_compute2i_atom, temp_free_rho_rel) 

!  PURPOSE
!  Tabulate the superposition of free-atom densities PCE correction on the entire integration grid.
  
   use pbc_lists
   use species_data
   use spline
   use grids, only: n_grid
   use dimensions, only : n_spin, use_initial_rho, use_density_gradient
   use constants, only : pi4_inv
   use localorb_io, only : use_unit
   use rel_x2c_mod, only : free_drho_dr_diff_spl
   implicit none

   integer :: n_atom_list
   integer :: atom_list(n_atom_list)
   integer:: n_compute_atoms
   integer:: i_compute2i_atom(n_compute_atoms)
   real*8 :: temp_free_rho_rel(n_compute_atoms)
   real*8, dimension(n_atom_list) :: i_r
   real*8, dimension(3,n_atom_list) :: dir_tab
   real*8, dimension(n_spin) :: pce
   real*8, dimension(3,n_spin) :: pce_gradient

!  INPUTS
!    o n_atom_list -- number of atoms
!    o atom_list -- list of atoms
!    o n_compute_atoms -- number of relevant atoms
!    o i_compute2i_atom -- this list contains the indices of relevant atom centers in atom_list
!    o temp_free_rho_rel -- this data is added to free_rho_superpos
!    o i_r  the distance from current integration point to all atoms in units of the logarithmic grid, i(r)
!    o dir_tab -- direction to atoms
!   
!  OUTPUT
!    o pce -- picture change error correction for X2C and Q4C; initialized to be free_rho_superpos/(4*pi)
!    o pce_gradient -- gradient of the pce term
!
   real*8 :: free_rho_superpos
   real*8 :: free_rho_gradient_superpos(3)
   real*8 :: atomic_gradient
   integer :: i_compute_atom, i_atom, i_coord, i_center

   pce = 0.d0; pce_gradient = 0.d0
   free_rho_superpos = 0.d0
   free_rho_gradient_superpos = 0.d0

   ! We only count occupied sites here, no ghost atoms
   do i_compute_atom = 1, n_compute_atoms, 1
     if (species_pseudoized(species_center(atom_list(i_compute_atom))))then
        write(use_unit,*)'WARNING: Pseudopotential was used for fully-relativistic calculations.'
        write(use_unit,*)'WARNING: Results may be incorrect.'
        write(use_unit,*)'WARNING: X2C does not support pseudopotential!'
        cycle
     end if

     free_rho_superpos = free_rho_superpos + temp_free_rho_rel(i_compute_atom)

     i_atom = i_compute2i_atom(i_compute_atom)
     i_center = atom_list(i_atom)

     if (use_density_gradient) then
        ! density gradient contribution
        ! factor pi4 already divided out!! different to free_rho_superpos!!!

        atomic_gradient = &
             val_spline( i_r(i_atom), &
             free_drho_dr_diff_spl(1,1,species_center(i_center)), &
             n_grid(species_center(i_center)) )

        atomic_gradient = atomic_gradient * pi4_inv

        do i_coord = 1, 3
           free_rho_gradient_superpos(i_coord) = &
                free_rho_gradient_superpos(i_coord) + &
                atomic_gradient * dir_tab(i_coord, i_atom)
        enddo
     end if

   enddo

   if (.not.use_initial_rho) then
      ! non-polarized case
      ! (Rundong) For fully-relativistic cases, we currently use this branch.
      pce(1)      = pi4_inv * free_rho_superpos
      pce(n_spin) = pi4_inv * free_rho_superpos
   else
      write(use_unit,*)'Error in subroutine evaluate_free_atom_sums_rel:'
      write(use_unit,*)'use_initial_rho = .true.'
      write(use_unit,*)'X2C does not support this setting, please change it.'
      stop
   end if

  ! if required, obtain density gradient also
  if (use_density_gradient) then
     if (use_initial_rho) then
        write(use_unit,*)'Error in subroutine evaluate_free_atom_sums_rel:'
        write(use_unit,*)'use_initial_rho = .true.'
        write(use_unit,*)'X2C/Q4C does not support this setting, please set it to .false.'
        stop
     end if
     pce_gradient(:,1)      = free_rho_gradient_superpos(:)
     pce_gradient(:,n_spin) = free_rho_gradient_superpos(:)
  end if

 end subroutine evaluate_free_atom_sums_rel


 subroutine evaluate_k_densmat_rel(kdm_complex, occ_numbers, KS_eigenvector_complex, i_spin, i_k_point, i_k)
   !  PURPOSE
   !    Calculate the k-dependend density matrix for one k and one spin only for real_eigenvectors.
   !  USES
   use aims_memory_tracking, only : aims_allocate, aims_deallocate
   use mpi_tasks, only: myid
   use dimensions, only : n_basis, n_states, n_spin, n_k_points, n_k_points_task
   use rel_x2c_mod, only : dim_matrix_rel, scf_iteration
   use localorb_io
   implicit none

   !  ARGUMENTS
   complex*16, intent(OUT) :: kdm_complex(2*dim_matrix_rel, 2*dim_matrix_rel)
   real*8, intent(IN) :: occ_numbers(n_states, n_spin, n_k_points)
   ! NOTE:  For this routine to give the correct results, the occ_numbers variable passed in should have already been
   !        properly k-weighted (and are thus not the "true" occupation numbers.)
   !        In scf_solver, this is done by the kweight_occs function (and undone at the end by de_kweight_occs)
   complex*16, intent(IN) :: KS_eigenvector_complex(2*dim_matrix_rel, n_states, n_spin, n_k_points_task)
   integer, intent(IN) :: i_spin, i_k_point, i_k

   !  INPUTS
   !    o occ_numbers -- occupation of states, with k-weighting already applied
   !    o KS_eigenvector_complex -- eigencoefficients
   !    o i_spin -- spin component
   !    o i_k_point -- k-point index
   !    o i_k -- node-local k-point index
   !  OUTPUTS
   !    o kdm_complex -- density matrix of this k-point

   integer :: i_state, n_max_occupied, info, i_bas, i,j,k
   real*8 :: occu
   complex*16, allocatable :: KS_scaled_cmplx(:,:), denmat_spinor(:,:)

   character*50 :: file_name
   character*500 :: info_str

   call aims_allocate( denmat_spinor, 2*dim_matrix_rel, 2*dim_matrix_rel, "denmat_spinor" )
   denmat_spinor = (0.d0, 0.d0)

   if(i_spin.eq.1)then
   ! write(info_str,'(2X,A,I2,A)') "i_spin=", i_spin, "   occ_numbers:"
   ! call localorb_info(info_str,use_unit,'(A)',OL_norm)
   !!write(use_unit,"(20f8.5)")occ_numbers(1:n_states, i_spin, i_k_point)
   ! write(info_str,'(2X,40(F8.5))') occ_numbers(1:max(n_states,40),i_spin,i_k_point)
   ! call localorb_info(info_str,use_unit,'(A)',OL_norm)
   endif

   call aims_allocate( KS_scaled_cmplx, 2*dim_matrix_rel, n_states, "KS_scaled_cmplx" )

   n_max_occupied = 1
   do i_state = 1, n_states, 1
      occu =  occ_numbers(i_state, i_spin, i_k_point)
      if (occu .gt. 0.d0) then
         n_max_occupied = i_state
      endif
      KS_scaled_cmplx(:, i_state) = KS_eigenvector_complex(:, i_state, i_spin, i_k) * dsqrt(occu)
   end do

   call zherk('U', 'N', 2*dim_matrix_rel, n_max_occupied, 1.d0, KS_scaled_cmplx, 2*dim_matrix_rel, 0.d0, denmat_spinor, 2*dim_matrix_rel)


   denmat_spinor = denmat_spinor + transpose(dconjg(denmat_spinor))
   do i_bas = 1, 2*dim_matrix_rel
      denmat_spinor(i_bas, i_bas) = denmat_spinor(i_bas, i_bas) * 0.5d0
   end do

   kdm_complex = denmat_spinor

  ! Allocatable arrays that are tracked
  if(allocated(KS_scaled_cmplx)) call aims_deallocate( KS_scaled_cmplx, "KS_scaled_cmplx" )
  if(allocated(denmat_spinor)) call aims_deallocate( denmat_spinor, "denmat_spinor" )

 end subroutine evaluate_k_densmat_rel


 subroutine accumulate_k_densmat_rel( n_k_points, kdm_complex, density_matrix )
  !  PURPOSE
  !    Accumulate the k-dependent density matrices to the real-space representation.
  !  USES
  use dimensions, only: n_basis, n_centers_basis_I
  use pbc_lists
  use rel_x2c_mod, only: dim_matrix_rel, scf_iteration
  use mpi_tasks,   only: myid,n_tasks
  use synchronize_mpi_basic, only: sync_vector
  implicit none
  integer,intent(in) :: n_k_points
  complex*16, intent(in) :: kdm_complex(2*dim_matrix_rel,2*dim_matrix_rel,n_k_points)
  real*8, intent(inout) :: density_matrix(n_centers_basis_I,n_centers_basis_I)
  !    o kdm{_cmplx} -- k-dependent density matrix for i_k_point
  !    o density_matrix{_sparse} -- real-space density matrix
  complex*16 :: add, temp_comp
  character*50 :: file_name
  integer :: i_k_point, i_bas1, i_bas2, i_basT1, i_basT2
  integer :: i,j, i_basis_1, i_basis_2, temp, i_cell_1, i_cell_2, i_cell, i_full_cell_1, i_full_cell_2
  integer :: n_useful_cells ! Not every cell in n_cells (see pbc_lists) is useful, we need to filter the useless ones out,
                            ! or the memory will explode.
  integer,allocatable :: full_cell_index(:), cell_index_tmp(:), cbasis_to_cell(:), full_bas_index(:,:)
  complex*16,allocatable :: mat_spinor(:,:) ! spinor density matrix: mat_spinor(2*dim_matrix_rel,2*dim_matrix_rel)
  complex*16,allocatable :: mat_scalar(:,:) ! temporary scalar density matrix: mat_scalar(n_basis,n_basis)

  allocate( cell_index_tmp(n_cells) )

  cell_index_tmp = 0
  n_useful_cells = 0
 ! first, I need to get the value of n_useful_cells
  do i_basis_1=1, n_centers_basis_I

     temp = center_to_cell( Cbasis_to_center(i_basis_1) )
     do i=1, n_useful_cells
        if( temp.eq.cell_index_tmp(i) ) goto 1001 ! this cell_index is already contained
     enddo

     n_useful_cells = n_useful_cells + 1
     cell_index_tmp(n_useful_cells) = center_to_cell( Cbasis_to_center(i_basis_1) )

 1001 continue
  enddo

  allocate( full_cell_index(n_useful_cells), cbasis_to_cell(n_centers_basis_I), &
            full_bas_index(n_basis,n_useful_cells) )
 ! then, I calculate full_cell_index and cbasis_to_cell
  full_cell_index = 0; cbasis_to_cell = 0; full_bas_index = 0
  i_cell_1 = 0
  do i_basis_1=1, n_centers_basis_I

     temp = center_to_cell( Cbasis_to_center(i_basis_1) )
     do i=1, i_cell_1
        if( temp.eq.full_cell_index(i) )then ! this cell_index is already contained
            cbasis_to_cell(i_basis_1) = i
            i_bas1 = Cbasis_to_basis(i_basis_1) ! i_bas1 go from 1 to n_basis
            full_bas_index(i_bas1,i) = i_basis_1
            goto 1002 
        endif
     enddo

     i_cell_1 = i_cell_1 + 1
     cbasis_to_cell(i_basis_1) = i_cell_1   ! i_basis_1 belongs to the i_cell_1 useful cell
     full_cell_index(i_cell_1) = center_to_cell( Cbasis_to_center(i_basis_1) ) ! i_cell_1 points to the full_cell_index(i_cell_1) cell in the space

     i_bas1 = Cbasis_to_basis(i_basis_1)  ! i_bas1 go from 1 to n_basis
     full_bas_index(i_bas1,i_cell_1) = i_basis_1

 1002 continue
  enddo

  deallocate(cell_index_tmp)

  allocate ( mat_scalar(n_basis,n_basis) )

 ! Sum over k:

 ! (Rundong) The serial code here used to be: 
 ! do i_cell_2 = 1, n_useful_cells
 ! do i_cell_1 = 1, n_useful_cells
 ! For MPI parallelization, we switch it to a single loop like now:
  do i_cell = 1, n_useful_cells*n_useful_cells
  if( (myid == mod(i_cell,n_tasks) .and. myid <= n_useful_cells*n_useful_cells) ) then
     i_cell_2 = int((i_cell-1.d-5)/n_useful_cells)+1
     i_cell_1 = mod(i_cell,n_useful_cells)
     if(i_cell_1.eq.0)i_cell_1=n_useful_cells

     allocate ( mat_spinor(2*dim_matrix_rel,2*dim_matrix_rel) )
     mat_spinor = (0.d0,0.d0)

     i_full_cell_1 = full_cell_index(i_cell_1);   i_full_cell_2 = full_cell_index(i_cell_2)

     do i_k_point = 1, n_k_points
       
       temp_comp = k_phase(i_full_cell_1,i_k_point) * dconjg(k_phase(i_full_cell_2,i_k_point))

       do i=1, 2*dim_matrix_rel
       do j=1, 2*dim_matrix_rel

          add = kdm_complex(j,i,i_k_point) * temp_comp
          mat_spinor(j,i) = mat_spinor(j,i) + add

       enddo
       enddo

     enddo ! end of i_k_point

     call transform_spinor2scalar_denmat(mat_spinor,mat_scalar)

     do i_bas2=1, n_basis
        i_basT2 = full_bas_index(i_bas2,i_cell_2)
        if(i_basT2.eq.0)cycle
        do i_bas1=1, n_basis
           i_basT1 = full_bas_index(i_bas1,i_cell_1)
           if(i_basT1.eq.0)cycle
           density_matrix(i_basT1,i_basT2) = dble(mat_scalar(i_bas1,i_bas2))
        enddo
     enddo

     deallocate ( mat_spinor )

  endif
  enddo ! end of i_cell

  call sync_vector(density_matrix,n_centers_basis_I*n_centers_basis_I)

  deallocate( full_cell_index, cbasis_to_cell, full_bas_index, mat_scalar )

 end subroutine accumulate_k_densmat_rel
 


