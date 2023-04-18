!****s* FHI-aims/evaluate_rho_from_densmat
!  NAME
!    evaluate_rho_from_densmat
!  SYNOPSIS

subroutine evaluate_rho_gradient_from_densmat(n_points, points, is_sparse, densmat, rho_gradient)

  !  PURPOSE
  !
  !    Evaluate real function defined by densmat and basis.f90.
  !    WARNING: This function determines the density gradient numerically
  !             It is both slow and relatively inaccurate. It is intended for
  !             testing purposes onlym so do not use it for production runs. 
  !
  !  USES

  use dimensions
  use pbc_lists
  use basis
  use grids
  use localorb_io, only: use_unit
  use mpi_tasks, only: check_allocation, aims_stop
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_points
  real*8, intent(IN) :: points(3, n_points)
  logical, intent(IN) :: is_sparse
  real*8, intent(IN) :: densmat(*)
  real*8, intent(OUT) :: rho_gradient(n_points,3)

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
  character(*), parameter :: func = 'evaluate_rho_gradient_from_densmat'
  real*8, dimension(3,7) :: localpoints !Dimensions are 3 spatial dimensions on 7 local points.
  real*8 :: rho(7) !Density on local and neighbouring points
  integer :: ilocalpoint,i_coord

 logical :: founderror = .false.
 real*8 :: dummy(1)
  real, parameter :: stepsize = 1e-3 !Stepsize for finite differences

!Wrapper
  do i_point=1,n_points,1
     !Localpoints: current point and a set of close points
     do i_coord=1,3,1
      do ilocalpoint=1,7,1
       localpoints(i_coord,ilocalpoint)=points(i_coord,i_point)
      enddo
     enddo
     !offset in x: point 2 and 3
     localpoints(1,2)=points(1,i_point)-stepsize
     localpoints(1,3)=points(1,i_point)+stepsize
     !offset in y: point 4 and 5
     localpoints(2,4)=points(2,i_point)-stepsize
     localpoints(2,5)=points(2,i_point)+stepsize
     !offset in x: point 6 and 7
     localpoints(3,6)=points(3,i_point)-stepsize
     localpoints(3,7)=points(3,i_point)+stepsize
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
  &                    n_basis_fns, outer_radius, &
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
        write(use_unit,*) 'Cbasis not found'
        write(use_unit,*) 'i_center: ', i_center
        write(use_unit,*) 'Coordinates of this center: '
        write(use_unit,*) coords_center(1:3,i_center)
        write(use_unit,*) 'i_cell: ', i_cell
        write(use_unit,*) 'i_basis: ' , i_basis
        founderror = .true.
     end if
  end do

  if (founderror)  call aims_stop('Cbasis not found (forgot map_to_center_cell?)', func)

  ! --- Evaluate atomic orbitals

  allocate(wave(7, n_comp), stat=info)
  call check_allocation(info, 'wave', func)

   call evaluate_waves_mult_point_center_fn(&
   &                  7, localpoints, &
   &                  n_centers, coords_center, &
   &                  n_comp, comp2center, comp2fn, comp2l, comp2m, &
   &                  n_species, species_center, &
   &                  4*n_max_grid, n_grid, r_grid_min, r_grid_inc, &
   &                  n_basis_fns, basis_wave_spl, outer_radius, &
   &                  .false., .false., wave,.false., dummy)

  ! --- Prune density matrix

  allocate(local_dm(n_comp, n_comp), work(n_comp, 7), stat=info)
  call check_allocation(info, 'local_dm, work', func)

  if (is_sparse) then
     call prune_density_matrix_sparse(densmat, local_dm, n_comp, comp2Cbasis)
  else
     call prune_density_matrix(densmat, local_dm, n_comp, comp2Cbasis)
  end if

  ! --- Evaluate rho on the 7 points

  call evaluate_KS_density_densmat(7, transpose(wave), n_comp, rho,  &
  & n_comp, n_centers_basis_T, local_dm, work)

  !Evaluate the gradient
  rho_gradient(i_point,1) = 0.5*((rho(3)-rho(1))+(rho(1)-rho(2)))/stepsize
  rho_gradient(i_point,2) = 0.5*((rho(5)-rho(1))+(rho(1)-rho(4)))/stepsize
  rho_gradient(i_point,3) = 0.5*((rho(7)-rho(1))+(rho(1)-rho(6)))/stepsize

  ! --- Tidy up

  deallocate(work, local_dm, wave, center2basis_off)
  deallocate(comp2center, comp2basis_sp, comp2fn, comp2l, comp2m, comp2Cbasis)

  enddo
end subroutine evaluate_rho_gradient_from_densmat
!******
