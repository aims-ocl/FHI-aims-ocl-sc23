!****s* FHI-aims/evaluate_rho_from_densmat
!  NAME
!    evaluate_rho_from_densmat
!  SYNOPSIS

subroutine evaluate_rho_from_densmat(n_points, points, is_sparse, densmat, rho)

  !  PURPOSE
  !
  !    Evaluate real function defined by densmat and basis.f90.
  !
  !  USES

  use dimensions
  use pbc_lists
  use basis
  use grids
  use mpi_tasks, only: check_allocation
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_points
  real*8, intent(IN) :: points(3, n_points)
  logical, intent(IN) :: is_sparse
  real*8, intent(IN) :: densmat(*)
  real*8, intent(OUT) :: rho(n_points)

  !  INPUTS
  !    o n_points -- Number of grid points
  !    o points -- Grid points
  !    o is_sparse -- Sparse storage of density matrix?
  !    o densmat -- Full density matrix
  !                 size: n_hamiltonian_matrix size for i_sparse,
  !                       (n_centers_basis_T, n_centers_basis_T) else.
  !  OUTPUTS
  !    o rho -- Density on grid points
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  integer :: n_max_comp, n_comp
  integer, allocatable :: center2basis_off(:)
  integer, allocatable :: comp2center(:), comp2basis_sp(:)
  integer, allocatable :: comp2fn(:), comp2l(:), comp2m(:)
  integer, allocatable :: comp2Cbasis(:)
  real*8, allocatable :: wave(:,:), local_dm(:,:), work(:,:)
  integer :: i_comp, j_comp, i_center, i_atom, i_bas, i_point, i_func
  integer :: i_Cbasis, i_basis, i_cell
  integer :: info
  real*8 :: dummy(1)
  character(*), parameter :: func = 'evaluate_rho_from_densmat'

  integer :: i_comp_2

  logical :: founderror = .false.
  real*8, allocatable :: dummy_radius(:)
  integer :: outer_radius_dimension

  !Create a dummy array for outer radius with large entries
  !This will prevent kinky cube files, without
  !explicitely changing outer_radius
  !
  !However, this will also change the list of basis functions that
  !will be included as "relevant" for the density at the present set
  !of grid points.
  !
  outer_radius_dimension=size(outer_radius)
  allocate(dummy_radius(outer_radius_dimension), stat=info)
  call check_allocation(info, 'dummy_radius', func)
  !dummy_radius=maxval(outer_radius)
  dummy_radius=outer_radius

  ! --- Prune (get comp2...)

  ! Worst case: All functions on all centers have maximum angular momentum.
  ! Improbable.  But we are talking only about some indexing array sizes.
  n_max_comp = n_centers * max_n_basis_fnLsp * (max_basis_L+1) * (2*max_basis_L+1)

  allocate(comp2center(n_max_comp), comp2basis_sp(n_max_comp), stat=info)
  call check_allocation(info, 'comp2center, comp2basis_sp', func)
  allocate(comp2fn(n_max_comp), stat=info)
  call check_allocation(info, 'comp2fn', func)
  allocate(comp2l(n_max_comp), comp2m(n_max_comp), stat=info)
  call check_allocation(info, 'comp2l, comp2m', func)

  
  call prune_general_basis(n_points, points, &
  &                    n_centers, coords_center, species_center, &
  &                    n_species, max_basis_L, max_n_basis_fnLsp, n_max_comp, &
  &                    Lsp2n_basis_fnLsp, Lsp2basis_fn, Lsp2basis_sp, &
  &                    n_basis_fns, dummy_radius, &
  &                    n_comp, comp2center, comp2basis_sp, &
  &                    comp2fn, comp2l, comp2m)

  ! --- comp2Cbasis

  allocate(center2basis_off(n_centers), stat=info)
  call check_allocation(info, 'center2basis_off', func)
  do i_center = 1, n_centers
     i_atom = center_to_atom(i_center)
     center2basis_off(i_center) = atom2basis_off(i_atom)
  end do

  allocate(comp2Cbasis(n_comp), stat=info)
  call check_allocation(info, 'comp2Cbasis', func)
  comp2Cbasis = 0   ! we rely on the assumption that unfound basis functions lead to a zero in comp2Cbasis below ...
  do i_comp = 1, n_comp
     i_center = comp2center(i_comp)
     i_cell = center_to_cell(i_center)
     i_basis = center2basis_off(i_center) + comp2basis_sp(i_comp)
     do i_Cbasis = 1, n_centers_basis_I
        if (i_basis == Cbasis_to_basis(i_Cbasis) .and. &
        &   i_cell == center_to_cell(Cbasis_to_center(i_Cbasis))) then
           comp2Cbasis(i_comp) = i_Cbasis
           exit
        end if
     end do
     if (i_Cbasis > n_centers_basis_I) then
!      Degbug testing
        comp2Cbasis(i_comp)=0 !Make sure this function is skipped
!        write(use_unit,*) '****Cbasis not found'
!        write(use_unit,*) 'i_center: ', i_center
!        write(use_unit,*) 'Koordinates of this center: '
!        write(use_unit,*) coords_center(1:3,i_center)
!        write(use_unit,*) 'i_cell: ', i_cell
!        write(use_unit,*) 'i_basis: ' , i_basis
!        founderror = .true.
     end if
  end do

! VB: In our normal density update, comp2Cbasis(i_compute) could never be zero - that
!     would mean that we counted a basis function as active that does not exist.
!
!     In the counting of active basis functions for the "cube" file output of the density
!     (threedimensional output of density on cartesian grid for plotting purposes),
!     however, comp2Cbasis(i_compute) = 0 can happen.
!
!     Well, it must not happen.
!
!     Since at this point, n_comp has been counted and only comp2Cbasis(i_comp) depends
!     on it, we re-compact the arrays associated with it and throw out the zeros.
!
!     This is, in principle, a very bad idea, as it will cost time. 
!
!     It also is a very bad idea as it does not address the problem. If a basis function
!     is zero, this fact should already be established by prune_basis_* .... there should
!     not be any basis functions that are _later_ not found.
!
!     However, it's the only safe fix
!     that I (VB) can think of right now, after the fact. At least it will tell us where the whole 
!     problem can strike in the first place.

  i_comp = 1
  do while (i_comp.le.n_comp)
      if (comp2Cbasis(i_comp).eq.0) then
          ! this basis function is a problem. We did not find it in
          ! the list of tabulated basis functions up to n_centers_basis_I .
          ! We throw it out.

          ! remove offending basis function from all arrays that index something to do with i_comp
          do i_comp_2 = i_comp+1,n_comp
             comp2Cbasis(i_comp_2-1)   = comp2Cbasis(i_comp_2)
             comp2center(i_comp_2-1)   = comp2center(i_comp_2)  
             comp2basis_sp(i_comp_2-1) = comp2basis_sp(i_comp_2)
             comp2fn(i_comp_2-1)       = comp2fn(i_comp_2)      
             comp2l(i_comp_2-1)        = comp2l(i_comp_2)       
             comp2m(i_comp_2-1)        = comp2m(i_comp_2)       
          enddo
          n_comp = n_comp-1
          founderror = .true.
      else 
          i_comp = i_comp+1
      end if
  enddo
  ! End recompacting. After this point, the density evaluation proceeds without any zeroes
  ! in the arrays that index the non-zero basis functions.

  !if (founderror) then
  !        write(use_unit,*) 
  !        write(use_unit,*) "** Warning: Number of active basis functions reduced in evaluate_rho_from_densmat.f90."
  !        ! ... but as detailed above, this should never have happened in the first place.
  !        write(use_unit,*) 
  !end if

  ! if (founderror)  call aims_stop('Cbasis not found (forgot map_to_center_cell?)', func)

  ! --- Evaluate atomic orbitals

  allocate(wave(n_points, n_comp), stat=info)
  call check_allocation(info, 'wave', func)

  call evaluate_waves_mult_point_center_fn(&
  &                  n_points, points, &
  &                  n_centers, coords_center, &
  &                  n_comp, comp2center, comp2fn, comp2l, comp2m, &
  &                  n_species, species_center, &
  &                  4*n_max_grid, n_grid, r_grid_min, r_grid_inc, &
  &                  n_basis_fns, basis_wave_spl, dummy_radius, &
  &                  .false., .false., wave, .false., dummy)

  ! --- Prune density matrix

  allocate(local_dm(n_comp, n_comp), work(n_comp, n_points), stat=info)
  call check_allocation(info, 'local_dm, work', func)

  if (is_sparse) then
     call prune_density_matrix_sparse(densmat, local_dm, n_comp, comp2Cbasis)
  else
     call prune_density_matrix(densmat, local_dm, n_comp, comp2Cbasis)
  end if

  ! --- Evaluate rho

  call evaluate_KS_density_densmat(n_points, transpose(wave), n_comp, rho,  &
  & n_comp, n_centers_basis_T, local_dm, work)

  ! --- Tidy up

  deallocate(work, local_dm, wave, center2basis_off)
  deallocate(comp2center, comp2basis_sp, comp2fn, comp2l, comp2m)
  deallocate(comp2Cbasis)
  deallocate(dummy_radius)

end subroutine evaluate_rho_from_densmat
!******
