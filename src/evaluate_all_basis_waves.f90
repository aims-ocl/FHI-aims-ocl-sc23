!****s* FHI-aims/evaluate_all_basis_waves
!  NAME
!    evaluate_all_basis_waves
!  SYNOPSIS

subroutine evaluate_all_basis_waves(n_points, points, &
&                                   n_cent, cent2coords, cent2species, &
&                                   n_comp, extrapolate_to_zero, &
&                                   wave, want_dwave, dwave)

  !  PURPOSE
  !
  !     Wrapper around evaluate_waves_mult_point_center_fn() for general
  !     geometry but basis functions from basis.f90.
  !
  !     Evaluate the basis functions (defined in basis.f90) at the specified
  !     centers.  The order of wave functions in wave(:,:) is the "standard
  !     order" for the given species within each center (that means, e.g. that
  !     there are n_basis_sp(cent2species(i_cent)) functions).  The blocks for
  !     the different centers are then concetenated similar to what is done
  !     during basis set construction for all atoms (see
  !     generate_full_bas.f90).
  !
  !     Please note that this routine is not optimized for performance and
  !     mainly meant for debugging purposes.
  !
  !  USES

  use dimensions
  use basis
  use prodbas
  use grids
  use mpi_tasks
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_points
  real*8, intent(IN) :: points(3, n_points)
  integer, intent(IN) :: n_cent
  real*8, intent(IN) :: cent2coords(3, n_cent)
  integer, intent(IN) :: cent2species(n_cent)
  integer, intent(IN) :: n_comp
  logical, intent(IN) :: extrapolate_to_zero
  real*8, intent(OUT) :: wave(n_points, n_comp)
  logical, intent(IN) :: want_dwave
  real*8, intent(OUT) :: dwave(3, n_points, n_comp)

  !  INPUTS
  !    o n_points, points -- grid points
  !    o n_cent, cent2... -- Atomic centers from which to take bas funcs.
  !    o n_comp -- number of all basis funcs (sum over n_basis_sp; is checked)
  !    o extrapolate_to_zero -- extrapolate loggrid spline to zero instead of
  !                             cutting; compatible code should use .false.
  !  OUTPUTS
  !    o wave -- Resulting function values
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

  integer, allocatable :: comp2center(:), comp2fn(:), comp2l(:), comp2m(:)
  integer, allocatable :: cent2bas_off(:)
  integer :: i_cent, i_species, i_fnLsp, i_bas_fn, i_bas_sp, L, M
  integer :: n_comp_uptonow, i_comp
  integer :: info
  real*8 :: dummy(1)
  character(*), parameter :: func = 'evaluate_all_basis_waves'

  ! --- Construct cent2bas_off

  allocate(cent2bas_off(n_cent), stat=info)
  call check_allocation(info, 'cent2bas_off', func)
  n_comp_uptonow = 0
  do i_cent = 1, n_cent
     cent2bas_off(i_cent) = n_comp_uptonow
     i_species = cent2species(i_cent)
     n_comp_uptonow = n_comp_uptonow + sp2n_basis_sp(i_species)
     if (n_comp_uptonow > n_comp) then
        call aims_stop('Consistenty error: Too many comp funcs.', func)
     end if
  end do
  if (n_comp_uptonow /= n_comp) then
     call aims_stop('Consistenty error: Too few comp funcs.', func)
  endif

  ! --- Construct comp2...

  allocate(comp2center(n_comp), comp2fn(n_comp), stat=info)
  call check_allocation(info, 'comp2center, comp2fn', func)
  allocate(comp2l(n_comp), comp2m(n_comp), stat=info)
  call check_allocation(info, 'comp2l, comp2m', func)

  do i_cent = 1, n_cent
     i_species = cent2species(i_cent)
     do L = 0, max_basis_L
        do i_fnLsp = 1, Lsp2n_basis_fnLsp(L, i_species)
           i_bas_fn = Lsp2basis_fn(i_fnLsp, L, i_species)
           i_bas_sp = Lsp2basis_sp(i_fnLsp, L, i_species)
           do M = -L, L
              i_comp = cent2bas_off(i_cent) + i_bas_sp + M
              comp2center(i_comp) = i_cent
              comp2fn(i_comp) = i_bas_fn
              comp2l(i_comp) = L
              comp2m(i_comp) = M
           end do
        end do
     end do
  end do
  
  call evaluate_waves_mult_point_center_fn(&
  &                  n_points, points, &
  &                  n_cent, cent2coords, &
  &                  n_comp, comp2center, comp2fn, comp2l, comp2m, &
  &                  n_species, cent2species, &
  &                  4*n_max_grid, n_grid, r_grid_min, r_grid_inc, &
  &                  n_basis_fns, basis_wave_spl, outer_radius, .false., &
  &                  extrapolate_to_zero, wave, want_dwave, dwave)

  deallocate(comp2center, comp2fn, comp2l, comp2m)
  deallocate(cent2bas_off)

end subroutine evaluate_all_basis_waves
!******
