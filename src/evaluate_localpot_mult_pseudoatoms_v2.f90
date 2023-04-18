!****s* FHI-aims/evaluate_localpot_mult_pseudoatoms_v2
!  NAME
!    evaluate_localpot_mult_pseudoatoms_v2
!  SYNOPSIS

subroutine evaluate_localpot_mult_pseudoatoms_v2(&
&                  n_points, points, &
&                  n_centers, center2coords, &
&                  n_species, center2species, &
&                  n_grid, r_grid_min, r_grid_inc, &
&                  fn2radius, &
&                  wave)

  !  PURPOSE
  !
  !    Evaluate the local pseudopot compontent from all the pseudoatoms
  !    and put it on the global integration grid.
  !
  !    Everything is done (copied from) like in evaluate_waves_mult_point_center_fn.f90
  !
  !    Only the multiplication with Y_lm is left out, since we only need
  !    spherical symmetric local potentials.
  !
  !  USES

  use spline
  use grids, only: invert_log_grid
  use pseudodata
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_points
  real*8, intent(IN) :: points(3, n_points)
  integer, intent(IN) :: n_centers
  real*8, intent(IN) :: center2coords(3, n_centers)
  integer, intent(IN) :: n_species
  integer, intent(IN) :: center2species(n_centers)
  integer, intent(IN) :: n_grid(n_species)
  real*8, intent(IN) :: r_grid_min(n_species), r_grid_inc(n_species)
  real*8, intent(IN) :: fn2radius(n_centers)
  real*8, intent(OUT) :: wave(n_points, n_centers)

  !  INPUTS
  !   Grid points
  !    o n_points, points -- Integration points at which to evaluate waves
  !   Basis specs
  !    o n_centers, center2coords -- Centers of centeric functions
  !   Grid related quantities
  !    o n_species, center2species -- Center to (grid-)species map
  !    o n_grid, r_grid_... -- Loggrids on which fns are defined
  !   Other radial part related quantities
  !    o fn2wave -- Radial splines (B-)spline
  !    o fn2radius -- Outer radius of this radial function.
  !              Evaluates to zero outside; set to large value to avoid this.
  !  OUTPUTS
  !    o wave -- Resulting function values
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  integer :: i_center, i_point, i_species
  real*8 :: relvec(3), dist, ilog
  integer :: i1
  real*8 :: rad_wave
  integer :: info
  character(*), parameter :: func = 'evaluate_local_pseudopot_mult_pseudoatom'


  ! --- center2Lmax, Lmax
  wave = 0.d0

  ! --- Main loop
  i1 =1

  do i_point = 1, n_points
     do i_center = 1, n_centers
           i_species = center2species(i_center)

           relvec = points(:, i_point) - center2coords(:, i_center)
           dist = sqrt(sum(relvec**2))
           ilog = invert_log_grid(dist, &
           &                    r_grid_min(i_species), r_grid_inc(i_species))


           if (dist < fn2radius(i_center)) then
              rad_wave = val_spline ( ilog, local_pseudopot_spl(:,:,i_center), n_grid(i_species) )
           else
              ! TODO: here we should add a interpolation for the farfield
              rad_wave = 0.d0
           end if
          wave(i_point, i_center) = rad_wave
     end do
  end do


end subroutine evaluate_localpot_mult_pseudoatoms_v2
!******
