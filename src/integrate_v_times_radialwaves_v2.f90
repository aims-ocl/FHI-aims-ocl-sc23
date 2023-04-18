!****s* FHI-aims/integrate_v_times_radialwaves_v2
!  NAME
!   integrate_v_times_radialwaves_v2
!  SYNOPSIS

subroutine integrate_v_times_radialwaves_v2 &
( v_times_radialwaves_spl,fn_l_basbas,species_basbas_fns, multipole)

  !  PURPOSE
  !  
  !    This subroutine computes the integration over the (possibly screened)
  !    Coulomb potential and the radial wave function.  The difference to
  !    integrate_v_times_radialwaves() is that this is done for all radial
  !    parts (n_basbas_fn) and not for all local product basis functions
  !    (n_loc_prodbas).
  !
  !  USES

  use dimensions
  use basis
  use runtime_choices
  use grids
  use geometry
  use species_data
  use spline
  use prodbas
  use constants
  use mpi_tasks
  use basbas_fn_coulomb
  implicit none

  !  ARGUMENTS
  real*8, dimension(n_max_spline, n_max_grid, n_basbas_fns), intent(OUT) &
  ::  v_times_radialwaves_spl
  integer, dimension(n_basbas_fns), intent(IN) :: fn_l_basbas
  integer, dimension(n_basbas_fns), intent(IN) :: species_basbas_fns
  real*8, intent(OUT) :: multipole(n_basbas_fns)

  !  INPUTS 
  !   o fn_l_basbas -- i_basbas_fn -> i_l (angular momentum)
  !   o species_basbas_fn -- i_basbas_fn -> i_species
  !  OUTPUTS
  !   o v_times_radialwaves_spl -- the integration over Coulomb interaction v
  !         and the radial part of the auxiliary basis function
  !   o multipole -- multipole moments
  !
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE

  !  local variables
  !  variables for hse_correction_term & testvariables
  real*8, dimension(:,:,:,:), allocatable :: hse_matrix
  integer :: max_prodbas_L
  real*8, dimension(:), allocatable :: v_times_radial_waves

  !  counters
  integer i_basbas_fn
  integer i_species
  integer i_grid
  integer i_l

  integer :: info
  character(*), parameter :: func = 'integrate_v_times_radialwaves_v2'
  !  begin work

  allocate(v_times_radial_waves(n_max_grid))
  !
  ! start allocate variables for hse_correction_term
  !
  if ((use_hse.and.hse_omega_hf.ne.0.0d0.and..not. use_gw_and_hse)&
     .or. lrc_pt2_started) then
     if (.not. use_logsbt_for_radial_hse_integration) then
        max_prodbas_L = maxval(fn_l_basbas)
        allocate(hse_matrix(n_max_grid,n_max_grid, &
        &                   max_prodbas_L+1,n_species), stat=info)
        call check_allocation(info, 'hse_matrix', func)

        do i_species = 1, n_species
           call integrate_errorfunction(max_prodbas_l, &
           & n_grid(i_species), n_max_grid, r_grid(:, i_species), &
           & hse_omega_hf, hse_matrix(:,:,:, i_species))
        end do
     end if
  endif
  ! end allocating errorfunction expansion


  do i_basbas_fn = 1, n_basbas_fns

     i_species = species_basbas_fns(i_basbas_fn)
     i_l = fn_l_basbas(i_basbas_fn)
     v_times_radial_waves = 0.d0

     if ((.not.(use_hse.and.(hse_omega_hf.ne.0.0d0)) .or. &
         use_gw_and_hse) .and. .not. lrc_pt2_started) then
        if(Adams_Moulton_integrator)then

           call adams_moulton_wave_integrator( &
           & i_l, basbas_wave_spl(1, :, i_basbas_fn), &
           & n_grid(i_species), r_grid(:, i_species), &
           & v_times_radial_waves, multipole(i_basbas_fn))

        else

           call trapezoidal_wave_integrator( &
           & i_l, basbas_wave_spl(1, :, i_basbas_fn), &
           & n_grid(i_species), r_grid(:, i_species), &
           & v_times_radial_waves, multipole(i_basbas_fn))

        endif

     else

        if (use_logsbt_for_radial_hse_integration) then

           call hse_logsbt_integrator_grid(i_l, &
           & basbas_wave_spl(:,:,i_basbas_fn), n_grid(i_species), &
           & r_grid_min(i_species), r_grid_inc(i_species), &
           & hse_omega_hf, v_times_radial_waves)

        else

           call hse_wave_integrator(i_l, &
           & basbas_wave_spl(1, :, i_basbas_fn), n_grid(i_species), &
           & n_hartree_grid, r_grid(:, i_species), &
           & hse_matrix(:,:, i_l+1, i_species), &
           & v_times_radial_waves)

        end if

        ! No multipole-like far field.
        multipole(i_basbas_fn) = 0.d0

     endif

     !  first cubic spline the v_times_radialwave integration on the log grid,
     !  then convert it to the radial grid

     call cubic_spline &
     ( v_times_radial_waves, n_grid(i_species), &
     v_times_radialwaves_spl( :, :,  i_basbas_fn))


  enddo

  deallocate (v_times_radial_waves)
  if(allocated(hse_matrix))then
     deallocate(hse_matrix)
  endif
  return

end subroutine integrate_v_times_radialwaves_v2
!******
