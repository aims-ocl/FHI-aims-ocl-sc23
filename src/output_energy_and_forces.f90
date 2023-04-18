!****s* FHI-aims/output_energy_and_forces
!  NAME
!   output_energy_and_forces
!  SYNOPSIS

subroutine output_energy_and_forces()

!  PURPOSE
!  The subroutine prints out the energy and the forces with different components
!  The forces are printed only, if they are calculated.
!
!  USES

  use physics
  use dimensions
  use mpi_tasks
  use runtime_choices
  use constants
  use relaxation
  use geometry
  use localorb_io

  implicit none

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
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



  ! local variables
  real*8 :: conversion
  real*8 :: total_energy_corrected
  real*8 :: free_energy
  character(len=1024) :: info_str

  ! counter
  integer :: i_atom

  conversion = hartree / bohr

  total_energy_corrected = total_energy + entropy_correction
  free_energy = total_energy + 2 * entropy_correction
  if (occupation_type == 2) then
    ! methfessel-paxton
    total_energy_corrected = total_energy &
      + 2 * (dble(n_methfessel_paxton + 1) / dble(n_methfessel_paxton + 2)) &
      * entropy_correction
  end if

  ! While calculating numerical stress -> no output of forces
  if (.not. use_numerical_stress .or. num_stress_finished) then

    if (use_forces .and. final_forces_cleaned) then
      ! initialize the force cleaning
      if (remove_unitary_force_components == 2) then
        call initialize_cleaning_forces()
      end if
      !VA: IF one wants to have the cleaned forces on the lattice vectors
      !    here one should
      !
      !(1) clean the atomic forces - this gives us cleaned atomic forces
      !    and cleaned forces on some forces_lv (from previous relaxation step, or zeros)
      !    call clean_force_components(total_forces, forces_lv)
      !
      !(2) extract the forces on lattice vectors - based on the cleaned atomic forces
      !    and the current stress tensor
      !    call extract_lv_forces (forces_lv, stress_tensor ,atomic_forces)

      !(3) clean the lattice vector.
      !    here one also makes a full call clean_force_components(total_forces, forces_lv)
      !    because it should be a projector, and hence does not change the atomic subspace anymore
      call clean_force_components(total_forces, forces_lv)
    end if

    if (myid == 0) then
      call localorb_info("", use_unit)
      write (info_str,"(2X,A)") "Energy and forces in a compact form:"
      call localorb_info(info_str, use_unit)
      write (info_str,"(2X,A,1X,E30.15,A)") "| Total energy uncorrected      :", total_energy * hartree, " eV"
      call localorb_info(info_str, use_unit)
      write (info_str,"(2X,A,1X,E30.15,A)") "| Total energy corrected        :", total_energy_corrected * hartree, &
        " eV  <-- do not rely on this value for anything but (periodic) metals"
      call localorb_info(info_str, use_unit)
      write (info_str,"(2X,A,1X,E30.15,A)") "| Electronic free energy        :", free_energy * hartree, " eV"
      call localorb_info(info_str, use_unit)
      if(use_meta_gga_post) then
        write (info_str,"(2X,A,1X,E30.15,A)") "| Total energy Meta-GGA         :",meta_gga_total_energy * hartree, " eV"
        call localorb_info(info_str, use_unit)
      end if
      if (use_forces) then
        if (final_forces_cleaned) then
          if (use_relaxation_constraints) then
             write (info_str,"(2X,A)") &
               "Total atomic forces (unitary forces were cleaned, then relaxation constraints were applied) [eV/Ang]:"
             call localorb_info(info_str, use_unit)
          else
             write (info_str,"(2X,A)") "Total atomic forces (unitary forces cleaned) [eV/Ang]:"
             call localorb_info(info_str, use_unit)
          end if
        else
          write (info_str,"(2X,A)") "Total atomic forces (derivative of free energy) [eV/Ang]:"
          call localorb_info(info_str, use_unit)
        end if
        do i_atom = 1, n_atoms
          write (info_str,"(2X,A,1X,I4,1X,3(E30.15,1X))") "|",i_atom, total_forces(:,i_atom) * conversion
          call localorb_info(info_str, use_unit)
        end do

        if (use_qmmm) then
          write (info_str,"(2X,A)") "Total external charge forces (derivative of free energy) [eV/Ang]:"
          call localorb_info(info_str, use_unit)
          do i_atom = 1, n_multipoles
            if (multipole_order(i_atom) == 0) then
              write (info_str,"(2X,A,1X,I4,1X,3(E30.15,1X))") "|",i_atom, ext_charge_forces(:,i_atom) * conversion
              call localorb_info(info_str, use_unit)
            end if
          end do
        end if
      end if
      call localorb_info("", use_unit)
    end if
  end if ! use_numerical_stress

end subroutine output_energy_and_forces
!******
