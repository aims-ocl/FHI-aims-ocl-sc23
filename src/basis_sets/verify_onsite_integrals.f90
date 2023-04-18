!****s* FHI-aims/verify_onsite_integrals
!  NAME
!   verify_onsite_integrals
!  SYNOPSIS
subroutine verify_onsite_integrals ( )
!  PURPOSE
!   Subroutine verify_onsite_integrals:
      ! Test onsite integration of all basis functions now, after the generation of
      ! all radial functions is complete.
!
!     This routine tests only the radial integration grids (not the angular part, which is
!     essentially un-testable by onsite integrals, by construction of the Lebedev grids).
!
      ! It must be ensured that the (sparse) radial grid in the 3D overlapping
      ! atom-centered integrations later gives the same result as the much
      ! denser, one-dimensional logarithmic grid used to generate the basis
      ! functions - or, indeed, that an analytic integration would give.
      ! This applies especially to Gaussian-type basis functions with very high 
      ! exponents.

      ! We only test the diagonal terms. Importantly, this includes the <1s| H |1s>
      ! term.
      ! Other (off-diagonal) onsite terms would be trivial to add, if it ever becomes
      ! necessary.
      ! 
!  USES
  use constants,       only : bohr, hartree, light_speed_sq
  use mpi_tasks,       only : aims_stop
  use dimensions,      only : n_max_grid, n_max_radial, n_basis_fns, &
                              n_max_spline, n_species
  use runtime_choices, only : flag_rel, REL_zora, onsite_accuracy_threshold, &
                              out_onsite, override_integration_accuracy
  use localorb_io,     only : use_unit, localorb_info
  use grids,           only : r_grid, n_grid, r_radial, w_radial, &
                              n_radial, r_grid_inc, invert_log_grid_p2
  use basis,           only : basis_wave_spl, basis_kinetic_spl, &
                              basis_deriv_spl, basisfn_species, &
                              basis_wave_spl
  use free_atoms,      only : free_potential, free_potential_spl
  use spline,          only : val_spline, val_spline_deriv
  use species_data,    only : species_name
  use timing,          only : warn_integrals
  implicit none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180, 2175-2196 (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Created by VB. Release version, FHI-aims (2013).
!  SOURCE

  integer :: i_species, i_basis_fn

  ! parts of the on-site integrals on the logarithmic grid

  integer :: i_grid
  real*8 :: alpha
  real*8 :: wave, kinetic, potential
  real*8 :: potential_deriv, wave_deriv
  real*8, dimension(n_max_grid) :: h_phi_log
  real*8, dimension(n_max_grid) :: phi_h_phi_log
  real*8, dimension(n_max_grid) :: r2_phi_h_phi_log
  real*8 :: integral_log

  ! parts of the on-site integrals on the sparser 'radial' grid

  integer :: i_radial
  real*8 :: i_r
  real*8, dimension(n_max_radial) :: h_phi_rad
  real*8, dimension(n_max_radial) :: phi_h_phi_rad
  real*8, dimension(n_max_radial) :: r2_phi_h_phi_rad
  real*8 :: integral_rad

  logical :: tolerance_exceeded = .false.

  character*150 :: info_str

  character*50 :: file_base
  character*50 :: file_name

  character(*), parameter :: func = 'verify_onsite_integrals'

  call localorb_info("Testing on-site integration grid accuracy.",use_unit,'(2X,A)')
  call localorb_info("|  Species  Function  <phi|h_atom|phi> (log., in eV)  <phi|h_atom|phi> (rad., in eV)",& 
                     use_unit,'(2X,A)')

  do i_species = 1, n_species
     ! Run through splined radial function array and pull out the
     ! radial functions from this species

     ! logarithmic mesh increment for this species
     alpha = log(r_grid(2,i_species)/r_grid(1,i_species))

     do i_basis_fn = 1, n_basis_fns, 1
        if ( i_species .eq. basisfn_species(i_basis_fn) ) then

           ! Evaluate ( T - v_free(r) ) * phi_i(r) on logarithmic grid first.
           ! The potential here is the DFT potential
           ! (electrons plus nuclear potential) 
           ! of the spherically symmetric free atom

           write(info_str,'(2X,I10,1X,I8)') i_species, i_basis_fn
           call localorb_info(info_str,use_unit,'(A,$)')

           integral_log = 0.d0
           do i_grid = 1,n_grid(i_species),1

              wave = basis_wave_spl(1,i_grid,i_basis_fn)
              kinetic = basis_kinetic_spl(1,i_grid,i_basis_fn)
              potential = free_potential(i_grid,i_species)

              ! For scaled ZORA, the "kinetic" part above only contains
              ! the non-relativistic parts, which are not yet consistent 
              ! with the potential. This leads to seemingly large integration 
              ! errors, where the real errors are much smaller.
              ! 
              ! For scaled ZORA, we therefore need to add the relativistic
              ! kinetic energy correction by hand.
              !
              ! For atomic ZORA, the exact same terms were already included
              ! in kinetic, subroutine shrink_fixed_basis_phi_thres() .
              if (flag_rel .eq. REL_zora) then
                 ! Add the missing terms to kinetic right here.

                 potential_deriv = & 
                   val_spline_deriv(dble(i_grid),free_potential_spl(1,1,i_species),&
                                    n_grid(i_species) ) &
                   / ( log(r_grid_inc(i_species)) * r_grid(i_grid,i_species) )

                 wave_deriv = basis_deriv_spl(1,i_grid,i_basis_fn)

                 kinetic = &
                   kinetic &
                   * 2 * light_speed_sq / (  2 * light_speed_sq - potential ) &
                   !
                   - light_speed_sq / (  2*light_speed_sq - potential)**2 &
                   !
                   * potential_deriv * ( wave_deriv - wave / r_grid(i_grid,i_species) )

              end if

              ! h phi_i(r) 
              ! Our definition is phi_i(r) = [ u(r)/r ] * Y_lm(Omega)

              h_phi_log(i_grid) = &
                ( potential & 
                    * wave & 
                  + kinetic ) & 
                * ( 1.d0/r_grid(i_grid,i_species) )

              ! phi_i(r) * h phi_i(r) 
              phi_h_phi_log(i_grid) = &
                wave &
                * h_phi_log(i_grid) & 
                * ( 1.d0/r_grid(i_grid,i_species) )

              ! r^2 * phi_i(r) * h phi_i(r) 
              ! and its integral
              r2_phi_h_phi_log(i_grid) = &
                r_grid(i_grid,i_species)**2.d0 * phi_h_phi_log(i_grid)

              integral_log = integral_log + alpha * r_grid(i_grid,i_species) * r2_phi_h_phi_log(i_grid)

           enddo 

           write(info_str,'(10X,F20.10)') integral_log*hartree
           call localorb_info(info_str,use_unit,'(A,$)')

           ! write logarithmic functions here
           if (out_onsite) then
              if (i_basis_fn.lt.10) then
                 write (file_base,'(A,A,I1)') trim(species_name(i_species)),"_",i_basis_fn
              else if (i_basis_fn.lt.100) then
                 write (file_base,'(A,A,I2)') trim(species_name(i_species)),"_",i_basis_fn
              else if (i_basis_fn.lt.1000) then
                 write (file_base,'(A,A,I3)') trim(species_name(i_species)),"_",i_basis_fn
              else if (i_basis_fn.lt.10000) then
                 write (file_base,'(A,A,I4)') trim(species_name(i_species)),"_",i_basis_fn
              else if (i_basis_fn.lt.100000) then
                 write (file_base,'(A,A,I5)') trim(species_name(i_species)),"_",i_basis_fn
              else
                 write (file_base,'(A,A,A)') trim(species_name(i_species)),"_",'nonumber'
              end if
              file_base=trim(file_base)

              write(file_name,'(A,A,A)') "Onsite_h_phi_log.",trim(file_base),".dat"
              file_name=trim(file_name)
              open(50,file=file_name)
                do i_grid = 1,n_grid(i_species),1
                   write(50,'(F20.10,1X,F20.10)') r_grid(i_grid,i_species)*bohr, h_phi_log(i_grid)
                enddo
              close(50)

              write(file_name,'(A,A,A)') "Onsite_phi_h_phi_log.",trim(file_base),".dat"
              file_name=trim(file_name)
              open(50,file=file_name)
                do i_grid = 1,n_grid(i_species),1
                   write(50,'(F20.10,1X,F20.10)') r_grid(i_grid,i_species)*bohr, phi_h_phi_log(i_grid)
                enddo
              close(50)

              write(file_name,'(A,A,A)') "Onsite_r2_phi_h_phi_log.",trim(file_base),".dat"
              file_name=trim(file_name)
              open(50,file=file_name)
                do i_grid = 1,n_grid(i_species),1
                   write(50,'(F20.10,1X,F20.10)') r_grid(i_grid,i_species)*bohr, r2_phi_h_phi_log(i_grid)
                enddo
              close(50)

           end if

           ! Next, evaluate ( T - v_free(r) ) * phi_i(r) on the sparser 'radial' integration grid
           ! which is used to set up the three-dimensional, overlapping atom-centered grids.

           integral_rad = 0.d0
           do i_radial = 1, n_radial(i_species), 1

              ! The radial functions are tabulated as splines on the logarithmic grid.
              ! To get their value on the 'radial' grid points, we need to evaluate
              ! these splines.

              ! get the coordinate of the current radius in units of the logarithmic grid
              i_r = invert_log_grid_p2 ( r_radial(i_radial,i_species), i_species )

              ! Spline-evaluate radial function and kinetic energy part. 
              ! The basis_wave array is shaped for use with spline_vector though.
              ! To use it with val_spline, we need to reshape the array implicitly as done below.
              ! Only done here because it is not time critical. (Else, use spline_vector.)
              wave = val_spline ( i_r, basis_wave_spl(1:n_max_spline,1:n_max_grid,i_basis_fn), &
                                  n_max_grid )
              kinetic = val_spline ( i_r, basis_kinetic_spl(1:n_max_spline,1:n_max_grid,i_basis_fn), &
                                  n_max_grid )
              potential = val_spline ( i_r, free_potential_spl(1,1,i_species), n_grid(i_species) )

              ! Require correction to kinetic part for scaled ZORA only - as above!
              if (flag_rel .eq. REL_zora) then
                 ! Add the missing terms to kinetic right here.

                 potential_deriv = & 
                   val_spline_deriv(i_r,free_potential_spl(1,1,i_species),&
                                    n_grid(i_species) ) &
                   / ( log(r_grid_inc(i_species)) * r_radial(i_radial,i_species) )

                 wave_deriv = val_spline ( i_r, & 
                                  basis_deriv_spl(1:n_max_spline,1:n_max_grid,i_basis_fn), &
                                  n_max_grid )

                 kinetic = &
                   kinetic &
                   * 2 * light_speed_sq / (  2 * light_speed_sq - potential ) &
                   !
                   - light_speed_sq / (  2*light_speed_sq - potential)**2 &
                   !
                   * potential_deriv * ( wave_deriv - wave / r_radial(i_radial,i_species) )

              end if

              ! h phi_i(r) [without zora, however]
              h_phi_rad(i_radial) = ( potential * wave + kinetic ) * 1.d0/r_radial(i_radial,i_species)

              ! phi_i(r) * h phi_i(r) 
              phi_h_phi_rad(i_radial) = wave * h_phi_rad (i_radial) * 1.d0/r_radial(i_radial,i_species)

              ! r^2 * phi_i(r) * h phi_i(r) 
              r2_phi_h_phi_rad(i_radial) = r_radial(i_radial,i_species)**2.d0 * phi_h_phi_rad(i_radial)
              integral_rad = integral_rad + w_radial(i_radial,i_species)*r2_phi_h_phi_rad(i_radial)

           enddo

           write(info_str,'(10X,F20.10)') integral_rad*hartree
           call localorb_info(info_str,use_unit,'(A,$)')

           ! write functions tabulated on 'radial; grid here
           if (out_onsite) then
              ! file_base was already set above

              write(file_name,'(A,A,A)') "Onsite_h_phi_rad.",trim(file_base),".dat"
              file_name=trim(file_name)
              open(50,file=file_name)
                do i_radial = 1,n_radial(i_species),1
                   write(50,'(F20.10,1X,F20.10)') r_radial(i_radial,i_species)*bohr, h_phi_rad(i_radial)
                enddo
              close(50)

              write(file_name,'(A,A,A)') "Onsite_phi_h_phi_rad.",trim(file_base),".dat"
              file_name=trim(file_name)
              open(50,file=file_name)
                do i_radial = 1,n_radial(i_species),1
                   write(50,'(F20.10,1X,F20.10)') r_radial(i_radial,i_species)*bohr, phi_h_phi_rad(i_radial)
                enddo
              close(50)

              write(file_name,'(A,A,A)') "Onsite_r2_phi_h_phi_rad.",trim(file_base),".dat"
              file_name=trim(file_name)
              open(50,file=file_name)
                do i_radial = 1,n_radial(i_species),1
                   write(50,'(F20.10,1X,F20.10)') r_radial(i_radial,i_species)*bohr, r2_phi_h_phi_rad(i_radial)
                enddo
              close(50)

           end if

           if ( abs(integral_log-integral_rad).gt.onsite_accuracy_threshold ) then

              tolerance_exceeded = .true.

              write(info_str,'(1X,A)') &
                '***'
              call localorb_info (info_str,use_unit,'(A)')

              if (integral_log.lt.0.d0) then
                 ! This radial function describes a bound electron.
                 ! Here, we should give a warning because this could impact accuracy.
                 warn_integrals = .true.
              end if

              if (.not.override_integration_accuracy) then
                 write(info_str,'(1X,A)') &
                   '* Stopping the calculation. Use the "override_integration_accuracy" flag to suppress this stop.'
                 call localorb_info (info_str,use_unit,'(A)')
                 call aims_stop('Radial integration grid not accurate enough.',func)
              end if

           else

              write(info_str,'(1X,A)') &
                ' '
              call localorb_info (info_str,use_unit,'(A)')

           end if

        end if
     end do
  end do

           if (tolerance_exceeded) then
              write(info_str,'(1X,A)') &
                '* Note: Onsite integrals marked "***" above are less accurate than'
              call localorb_info (info_str,use_unit,'(A)')
              write(info_str,'(1X,A,F10.5,A,A)') &
                '* onsite_accuracy_threshold = ', onsite_accuracy_threshold*hartree, ' eV.', &
                ' Usually, this is harmless.'
              call localorb_info (info_str,use_unit,'(A)')
              if (warn_integrals) then
                 write(info_str,'(1X,A)') &
                   '* However, one integral did have a negative value. This integral could be important.'
                 call localorb_info (info_str,use_unit,'(A)')
              end if
              write(info_str,'(1X,A)') &
                '* When in doubt, tighten the "radial" and/or "radial_multiplier" flags to check.'
              call localorb_info (info_str,use_unit,'(A)')
           end if

              write(info_str,'(2X,A)') &
                ' '
              call localorb_info (info_str,use_unit,'(A)')

  return
end subroutine verify_onsite_integrals
!******
