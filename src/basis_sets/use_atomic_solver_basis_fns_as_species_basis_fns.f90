!****s* FHI-aims/use_atomic_solver_basis_fns_as_species_basis_fns
!  NAME
!   use_atomic_solver_basis_fns_as_species_basis_fns
!  SYNOPSIS
subroutine use_atomic_solver_basis_fns_as_species_basis_fns ( ) 
!  PURPOSE
!    Subroutine use_atomic_solver_basis_fns_as_species_basis_fns reuses the wavefunctions
!    and derivative quantities that came from the chosen atomic solver as the minimal 
!    basis.  Essentially a stripped-down version of get_species_basis_fns().
!  USES
  use constants,       only : bohr, hartree
  use dimensions,      only : n_species, use_basis_gradients
  use runtime_choices, only : wave_threshold
  use grids,           only : n_grid, r_grid
  use species_data,    only : atomic_wave, species_name, atomic_n, atomic_l, atomic_eigenval, &
                              atomic_outer_radius, n_atomic, r_cutoff, atomic_cutoff, &
                              atomic_eigenval, atomic_kinetic, atomic_wave_deriv
  use free_atoms,      only : free_wave_eigenval, free_wave, free_kinetic, free_wave_deriv
  use localorb_io,     only : use_unit, localorb_info
  implicit none
!  INPUT
!   none
!  OUTPUT
!   none
!  AUTHOR
!    William Huhn
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
!    Part of atom_sphere integration, October 2016
!  SOURCE
  character*100 :: info_str
  integer i_species, i_fn, i_grid

!  begin work

  write (info_str,'(2X,A,A)') &
       "Reusing the basis functions from the chosen atomic solver as the minimal/atomic basis."
  call localorb_info(info_str,use_unit,'(A)')

  do i_species = 1, n_species, 1
    do i_fn = 1, n_atomic(i_species), 1
      atomic_cutoff(i_species,i_fn) = r_cutoff(i_species)
      atomic_eigenval(i_species, i_fn) = free_wave_eigenval(i_species, i_fn)
      do i_grid = 1, n_grid(i_species),1
        atomic_wave(i_grid,i_species,i_fn) = free_wave(i_grid,i_species,i_fn)
        atomic_kinetic(i_grid, i_species, i_fn) = free_kinetic(i_grid,i_species,i_fn)
        if (use_basis_gradients) then
          atomic_wave_deriv(i_grid, i_species, i_fn) = free_wave_deriv(i_grid, i_species, i_fn)
        end if
      enddo
      ! atomic_n and atomic_l (and likely a number of other variables I'm
      ! forgetting) were set back in read_species_data

      i_grid = n_grid(i_species)
      do while (abs(atomic_wave(i_grid,i_species,i_fn)).lt.wave_threshold)
        i_grid = i_grid - 1 
      end do
      atomic_outer_radius(i_species,i_fn) = r_grid(i_grid,i_species)
    end do

!   output atomic basis data

    call localorb_info('',use_unit)

    write(info_str,'(2X,A,A,A)') &
         "Species ", species_name(i_species),":"
    call localorb_info(info_str,use_unit,'(A)')

    call localorb_info('',use_unit)

    write(info_str,'(2X,A)') &
         "List of atomic basis orbitals and eigenvalues: "
    call localorb_info(info_str,use_unit,'(A)')

    write(info_str,'(4X,A,4X,A,6X,A,4X,A,4X,A)') &
         "n", "l", "energy [Ha]", &
         "energy [eV]", "outer radius [A]"
    call localorb_info(info_str,use_unit,'(A)')

    do i_fn = 1, n_atomic(i_species), 1
      write(info_str,'(2X,I3,2X,I3,2X,F15.6,F15.4,F15.6)') &
           atomic_n(i_species, i_fn), atomic_l(i_species, i_fn), &
           atomic_eigenval(i_species, i_fn), &
           atomic_eigenval(i_species, i_fn)*hartree, &
           atomic_outer_radius(i_species,i_fn)*bohr
      call localorb_info(info_str,use_unit,'(A)')
    enddo
    call localorb_info('',use_unit)
  enddo

end subroutine use_atomic_solver_basis_fns_as_species_basis_fns
!******
