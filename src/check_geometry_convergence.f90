!****s* FHI-aims/check_geometry_convergence
!  NAME
!    check_geometry_convergence
!  SYNOPSIS

subroutine check_geometry_convergence( total_forces, forces_lv,  converged_geo )

!  PURPOSE
!  Check geometry convergence criterion using the present maximum 
!  force component
! USES

  use localorb_io
  use dimensions
  use constants
  use runtime_choices
  use geometry
  use species_data
  use relaxation
  implicit none

!  ARGUMENTS

  real*8, dimension(3, n_atoms)    :: total_forces
  real*8, dimension(3, n_periodic) :: forces_lv
  logical :: converged_geo

!  INPUTS
!    none ????UPDATE HERE?????????
!  OUTPUT
!    none ????UPDATE HERE?????????
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!     Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




  ! local variables

  real*8 :: max_force
  real*8 :: max_force_atom
  real*8 :: max_force_lv=0d0

  character*80 :: info_str

  real*8, external :: get_max_force

  ! begin work

  max_force_atom = get_max_force(total_forces)

  if (relax_unit_cell .ne.0) then
     max_force_lv= maxval(abs(forces_lv)) 
     ! FlK: summarize together with max force below
     ! write(info_str,'(2X,A,E14.6,A)') '|| Forces on lattice vectors || = '&
     !      , max_force_lv* hartree/bohr, ' eV/A.'
     ! call localorb_info(info_str)
  end if

  max_force  =  max (max_force_atom, max_force_lv)

  max_force = max_force * hartree/bohr

  converged_geo = (max_force.lt.relax_accuracy_forces)

  write(info_str,'(2X,A)') 'Net remaining forces (excluding translations, rotations) in present geometry:'
  call localorb_info(info_str)

  ! FlK:
  write(info_str,'(2X,A,E14.6,A)') '|| Forces on atoms   || = '&
        , max_force_atom* hartree/bohr, ' eV/A.'
    call localorb_info(info_str)
  if (relax_unit_cell .ne.0) then
    write(info_str,'(2X,A,E14.6,A)') '|| Forces on lattice || = '&
         , max_force_lv* hartree/bohr, ' eV/A^3.'
    call localorb_info(info_str)
  end if

  write(info_str,'(2X,A,E14.6,A)') 'Maximum force component is', max_force, ' eV/A.'
  call localorb_info(info_str)

  if (.not.converged_geo) then

    write(info_str,'(2X,A)') 'Present geometry is not yet converged.'
    call localorb_info(info_str)

  else

    write(info_str,'(2X,A)') 'Present geometry is converged.'
    call localorb_info(info_str)

  endif

  write(info_str,'(2X,A)') ''
  call localorb_info(info_str)

end subroutine check_geometry_convergence
!******
