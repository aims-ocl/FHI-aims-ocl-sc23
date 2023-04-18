!----------------------------------------------------------------------------
!****s* FHI-aims/generate_basbasfn
!  NAME
!    generate_basbasfn
!  SYNOPSIS

subroutine generate_basbasfn(basbas_wave, fn_to_l, fn_to_species)

  !  PURPOSE
  !
  !    Given the n_basbas_fns radial parts in basbas_wave, prepare other
  !    arrays (basbas_wave_spl, field_radius_basbas_fn,
  !    charge_radius_basbas_fn, multipole_basbas_fn, v_times_radialbasbas_spl).
  !
  !  USES

  use prodbas
  use basbas_fn_coulomb
  use localorb_io
  use grids
  use species_data
  use runtime_choices
  use spline
  use debug_output
  use dimensions, only: use_hse, n_max_grid, n_basbas_fns, n_species, &
      use_gw_and_hse
  use mpi_tasks, only: aims_stop, check_allocation, myid
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: basbas_wave(n_max_grid, n_basbas_fns)
  integer, intent(IN) :: fn_to_l(n_basbas_fns)
  integer, intent(IN) :: fn_to_species(n_basbas_fns)

  !  INPUTS
  !    o basbas_wave -- radial part of the product basis fns [u(r)=r*P(r)]
  !    o fn_to_l -- i_fn -> angular momentum quantum number
  !    o fn_to_species -- i_fn -> atom types
  !  OUTPUTS
  !    - sets product basis spline arrays in prodbas.f90.
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
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  integer :: i_basbas_fn, i_l, i_species, max_L, info
  integer, allocatable :: n_in_channel(:,:)
  character*150 :: info_str, filename
  character*50 :: charge_str, field_str, mp_str
  character(*), parameter :: func = 'generate_basbasfn'

  if (n_basbas_fns <= 0) then
     call aims_stop('No product basis function radial parts given.', func)
  end if
  if (.not. allocated(basbasfn_l) .or. &
  &   .not. allocated(field_radius_basbas_fn)) then
     call aims_stop('Missing allocations.', func)
  end if

  call localorb_info('Product basis:', use_unit, '(2X,A)', OL_norm)
  call localorb_info('charge radius: extent of product basis function', &
  &                  use_unit, "(2X,'| ',A)", OL_norm)
  call localorb_info('field radius: extent of its Coulomb potential', &
  &                  use_unit, "(2X,'| ',A)", OL_norm)
  if ((use_hse.and.(hse_omega_hf.ne.0.0d0).and. .not. use_gw_and_hse) &
      .or. lrc_pt2_started) then
     write(info_str, "(A7,' ',A3,' ',A14,2X,A14)") &
     & 'Species', '  l', 'charge radius', ' field radius'
  else
     write(info_str, "(A7,' ',A3,' ',A14,2X,A14,2X,A14)") &
     & 'Species', '  l', 'charge radius', ' field radius', 'multipol moment'
  end if
  call localorb_info(info_str, use_unit, "(2X,'| ',A)", OL_norm)

  ! Get basbas_wave_spl
  do i_basbas_fn = 1, n_basbas_fns, 1
     i_l = fn_to_l(i_basbas_fn)
     i_species = fn_to_species(i_basbas_fn)
     call cubic_spline &
     ( basbas_wave(1, i_basbas_fn), n_grid(i_species), &
     basbas_wave_spl(1, 1,i_basbas_fn) )
  enddo

  ! Get v_times_radialbasbas_spl
  call integrate_v_times_radialwaves_v2(v_times_radialbasbas_spl, &
  &                   fn_to_l, fn_to_species, multipole_basbas_fn)

  ! Get extents & stuff.
  do i_basbas_fn = 1, n_basbas_fns, 1
     i_l = fn_to_l(i_basbas_fn)
     i_species = fn_to_species(i_basbas_fn)
     call inspect_wave_field(i_l, basbas_wave(:, i_basbas_fn), &
     & v_times_radialbasbas_spl(1,:, i_basbas_fn), &
     & multipole_basbas_fn(i_basbas_fn), &
     & n_grid(i_species), r_grid(:, i_species), wave_threshold, &
     & charge_radius_basbas_fn(i_basbas_fn), &
     & field_radius_basbas_fn(i_basbas_fn))

     if (use_hse.and.(hse_omega_hf.ne.0.0d0).and. .not. use_gw_and_hse &
         .or. lrc_pt2_started) then
        write(info_str, &
        & "(3X,A4,' ',I3,' ',F12.6,' A',2X,F12.6,' A')") &
        & species_name(i_species), i_l, &
        & charge_radius_basbas_fn(i_basbas_fn)*bohr, &
        & field_radius_basbas_fn(i_basbas_fn)*bohr
     else if (abs(multipole_basbas_fn(i_basbas_fn)) > 1d-12) then
        write(info_str, &
        & "(3X,A4,' ',I3,' ',F12.6,' A',2X,ES12.6,' A',2X,F18.6,' a.u.')") &
        & species_name(i_species), i_l, &
        & charge_radius_basbas_fn(i_basbas_fn)*bohr, &
        & field_radius_basbas_fn(i_basbas_fn)*bohr, &
        & multipole_basbas_fn(i_basbas_fn)
     else
        write(info_str, &
        & "(3X,A4,' ',I3,' ',F12.6,' A',2X,F12.6,' A',2X,ES18.6,' a.u.')") &
        & species_name(i_species), i_l, &
        & charge_radius_basbas_fn(i_basbas_fn)*bohr, &
        & field_radius_basbas_fn(i_basbas_fn)*bohr, &
        & multipole_basbas_fn(i_basbas_fn)
     end if
     call localorb_info(info_str, use_unit, "(2X,'| ',A)", OL_norm)
  enddo
  call localorb_info('', use_unit, "(A)", OL_norm)

  ! --- output

  if (out_basis) then
     max_L = maxval(fn_to_l(1:n_basbas_fns))
     allocate(n_in_channel(0:max_L, n_species), stat=info)
     call check_allocation(info, 'n_in_channel', func)
     n_in_channel = 0
     do i_basbas_fn = 1, n_basbas_fns
        i_l = fn_to_l(i_basbas_fn)
        i_species = fn_to_species(i_basbas_fn)
        n_in_channel(i_l, i_species) = n_in_channel(i_l, i_species) + 1
        write(filename, "('prodbas-',A,'-',I1,'-',I2.2,'.dat')") &
        & trim(species_name(i_species)), i_l, n_in_channel(i_l, i_species)
        call debug_plot_data(filename, &
        &      r_grid(1:n_grid(i_species), i_species), &
        &      basbas_wave(1:n_grid(i_species), i_basbas_fn), &
        &      v_times_radialbasbas_spl(1, 1:n_grid(i_species), i_basbas_fn))
     end do
     deallocate(n_in_channel)
  end if

  ! --- Check

  do i_basbas_fn = 1, n_basbas_fns, 1
     i_species = fn_to_species(i_basbas_fn)
     if (all(abs(basbas_wave(1:n_grid(i_species), i_basbas_fn)) < wave_threshold)) then
        write(use_unit,'(1X,A,I7,A,A)') &
        & "* Task", myid, &
        & ": Warning - an auxiliary basis function is lower than the ", &
        & "wave_threshold value for integrations everywhere." 
     end if
  end do

end subroutine generate_basbasfn
!******
