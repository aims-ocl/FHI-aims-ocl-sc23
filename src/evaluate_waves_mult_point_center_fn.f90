!****s* FHI-aims/evaluate_waves_mult_point_center_fn
!  NAME
!    evaluate_waves_mult_point_center_fn
!  SYNOPSIS

subroutine evaluate_waves_mult_point_center_fn(&
&                  n_points, points, &
&                  n_centers, center2coords, &
&                  n_comp, comp2center, comp2fn, comp2l, comp2m, &
&                  n_species, center2species, &
&                  n_wave_size, n_grid, r_grid_min, r_grid_inc, &
&                  n_fn, fn2wave, fn2radius, use_bwave, &
&                  extrapolate_to_zero, &
&                  wave, want_dwave, dwave)

  !  PURPOSE
  !
  !    Evaluate multiple numeric atomic orbitals (NAOs) at multiple centers on
  !    a set of multiple points.
  !
  !    Optimized for subsequent entries in comp2center (and comp2fn).
  !
  !    This routine is mainly meant to contain the difficult parts of a whole
  !    bunch of wrappers.  Most probably, there is already some fitting
  !    wrapper and you do not need to bother with all the arguments.
  !
  !  USES

  use spline
  use bspline
  use grids, only: invert_log_grid
  use mpi_tasks, only: check_allocation
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_points
  real*8, intent(IN) :: points(3, n_points)
  integer, intent(IN) :: n_centers
  real*8, intent(IN) :: center2coords(3, n_centers)
  integer, intent(IN) :: n_comp
  integer, intent(IN) :: comp2center(n_comp), comp2fn(n_comp)
  integer, intent(IN) :: comp2l(n_comp), comp2m(n_comp)
  integer, intent(IN) :: n_species
  integer, intent(IN) :: center2species(n_centers)
  integer, intent(IN) :: n_wave_size
  integer, intent(IN) :: n_grid(n_species)
  real*8, intent(IN) :: r_grid_min(n_species), r_grid_inc(n_species)
  integer, intent(IN) :: n_fn
  real*8, intent(IN) :: fn2wave(n_wave_size, n_fn)
  real*8, intent(IN) :: fn2radius(n_fn)
  logical, intent(IN) :: use_bwave
  logical, intent(IN) :: extrapolate_to_zero
  real*8, intent(OUT) :: wave(n_points, n_comp)
  logical, intent(IN) :: want_dwave
  real*8, intent(OUT) :: dwave(3, n_points, n_comp)

  !  INPUTS
  !   Grid points
  !    o n_points, points -- Integration points at which to evaluate waves
  !   Basis specs
  !    o n_centers, center2coords -- Centers of centeric functions
  !    o n_comp, comp2... -- Functions to compute
  !   Grid related quantities
  !    o n_species, center2species -- Center to (grid-)species map
  !    o n_wave_size -- if (use_bwave): n_max_grid+2; else: 4*n_max_grid
  !    o n_grid, r_grid_... -- Loggrids on which fns are defined
  !   Other radial part related quantities
  !    o n_fn -- Number of radial functions
  !    o fn2wave -- Radial splines (B-)spline
  !    o fn2radius -- Outer radius of this radial function.
  !              Evaluates to zero outside; set to large value to avoid this.
  !    o use_bwave -- Is it B- or ordinary cubic spline.
  !    o want_dwave -- If .true., calculate dwave.
  !    o extrapolate_to_zero -- extrapolate loggrid spline to zero instead of
  !                             cutting; compatible code should use .false.
  !  OUTPUTS
  !    o wave -- Resulting function values
  !    o dwave -- if (want_dwave) gradients of function values,
  !               not referenced otherwise.
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

  integer :: i_comp, i_center, i_point, i_fn, i_l
  integer :: i_fn_last, i_l_last, i_center_last, i_species
  integer :: center2Lmax(n_centers), center_Lmax, Lmax
  integer :: n_lm, n_max_lm, i_lm
  real*8 :: relvec(3), dist, ilog
  integer :: i1
  real*8, allocatable :: blender(:,:), dblender(:,:)
  real*8 :: dr_drelvec(3)
  real*8, allocatable :: rlylm(:), drlylm_drelvec(:,:)
  real*8 :: f_r, df_r
  integer :: info
  character(*), parameter :: func = 'evaluate_waves_mult_point_center_fn'

  ! --- center2Lmax, Lmax

  wave = 0.d0
  center2Lmax = -1
  do i_comp = 1, n_comp
     i_center = comp2center(i_comp)
     center2Lmax(i_center) = max(center2Lmax(i_center), comp2l(i_comp))
  end do
  Lmax = maxval(center2Lmax)
  n_max_lm = (Lmax+1)**2
  allocate(rlylm(n_max_lm), blender(4, 0:Lmax), stat=info)
  call check_allocation(info, 'rlylm, blender', func)
  if (want_dwave) then
     allocate(drlylm_drelvec(3, n_max_lm), dblender(4, 0:Lmax), stat=info)
     call check_allocation(info, 'drlylm_drelvec, dblender', func)
  else 
     allocate(drlylm_drelvec(1, 1), dblender(1, 1), stat=info)
     call check_allocation(info, 'drlylm_drelvec, dblender', func)
  end if

  ! JW: In principle, we could preprocess comp2{center,fn,l,m} right now,
  ! hoping that most 'comp2m's are subsequent.  Stuff like i_lm could be
  ! precalculated as comp2lm.  This should make indexing less indirect.
  ! loops: (point,) i_center, i_fn, m


  ! --- Main loop

  do i_point = 1, n_points
     i_center_last = 0
     i_fn_last = 0
     i_l_last = -1
     do i_comp = 1, n_comp
        ! New center -> New spherical coords, new Ylms
        i_center = comp2center(i_comp)
        if (i_center /= i_center_last) then
           i_fn_last = 0
           i_l_last = -1
           i_center_last = i_center
           i_species = center2species(i_center)
           center_Lmax = center2Lmax(i_center)
           n_lm = (center_Lmax+1)**2
           relvec = points(:, i_point) - center2coords(:, i_center)
           dist = sqrt(sum(relvec**2))
           call log_spline_blender(dist, &
           & n_grid(i_species), r_grid_min(i_species), r_grid_inc(i_species), &
           & center_Lmax, use_bwave, want_dwave, extrapolate_to_zero, &
           & i1, blender, dblender)
           if (want_dwave) then
              call get_rlylm_and_derivs(relvec, center_Lmax, n_lm, &
              &                         rlylm, drlylm_drelvec)
              if (dist /= 0.d0) then
                 dr_drelvec = relvec / dist
              else
                 dr_drelvec = 0.d0
              end if
           else
              call get_rlylm(relvec, center_Lmax, n_lm, rlylm)
           end if
        end if
        ! New fn -> New f_r
        i_fn = comp2fn(i_comp)
        i_l = comp2l(i_comp)
        if (i_fn /= i_fn_last .or. i_l /= i_l_last) then
           i_fn_last = i_fn
           i_l_last = i_l
           if (dist < fn2radius(i_fn)) then
              f_r = dot_product(fn2wave(i1:i1+3, i_fn), blender(:, i_l))
              if (want_dwave) then
                 df_r = dot_product(fn2wave(i1:i1+3, i_fn), dblender(:, i_l))
              end if
           else
              ! JW: To avoid that choice of clustering has impact on final
              ! numbers, always set wave function to zero outside
              ! outer_radius.
              f_r = 0.d0
              df_r = 0.d0
           end if
        end if
        i_lm = i_l**2 + i_l + comp2m(i_comp) + 1
        wave(i_point, i_comp) = f_r * rlylm(i_lm)
        if (want_dwave) then
           dwave(:, i_point, i_comp) = f_r * drlylm_drelvec(:, i_lm) &
           &                         + df_r * dr_drelvec * rlylm(i_lm)
        end if
     end do
  end do

  deallocate(rlylm, blender)
  if (allocated(drlylm_drelvec)) deallocate(drlylm_drelvec)
  if (allocated(dblender)) deallocate(dblender)

end subroutine evaluate_waves_mult_point_center_fn
!******
