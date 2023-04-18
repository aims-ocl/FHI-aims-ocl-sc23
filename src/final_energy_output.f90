!****s* FHI-aims/final_energy_output
!  NAME
!   final_energy_output
!  SYNOPSIS

  subroutine final_energy_output

!  PURPOSE
!
!  At the very end of am FHI-aims run, this routine provides some
!  final energy information for easier verification by the user.
!  No new energies should be computed here, and it must be absolutely
!  certain that the total energy values are not misleading. There may be
!  more than one possible "total energy" that a user may be interested in.
!
!  USES

  use dimensions
  use runtime_choices
  use constants
  use localorb_io
  use mpi_tasks
  use physics
  use xml_write
  use free_atoms, only: average_free_es_pot

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
!    Release version, FHI-aims (2012).
!  SOURCE

  implicit none

  character*150 info_str

  ! begin work

  call localorb_info("")

  if (calculate_atom_bsse) then
    ! corner case - no output

    write(info_str,"(2X,A)") &
      "Final total energy output will not be written in atomization BSSE run."
    call localorb_info(info_str)
    write(info_str,"(2X,A)") &
      "A listing of the relevant total energies for all subsystems can be found above."
    call localorb_info(info_str)

  else

    write(info_str,"(A)") "------------------------------------------------------------------------------"
    call localorb_info(info_str)
    write(info_str,"(2X,A)") "Final output of selected total energy values:"
    call localorb_info(info_str)

    call localorb_info("")

    write(info_str,"(2X,A)") "The following output summarizes some interesting total energy values"
    call localorb_info(info_str)
    write(info_str,"(2X,A)") "at the end of a run (AFTER all relaxation, molecular dynamics, etc.)."
    call localorb_info(info_str)

    write(info_str,"(2X,A)") ""
    call localorb_info(info_str)

    ! Actual total energy values:

    write(info_str,"(2X,A,F25.9,A)") &
      "| Total energy of the DFT / Hartree-Fock s.c.f. calculation      :", &
      (total_energy * hartree), " eV"
    call localorb_info(info_str)

    call xml_elem("energy", total_energy, name="Total energy", unit="Ha")
    call xml_elem("energy", en_vdw, name="van der Waals energy", unit="Ha")

    if (occupation_type == 2) then
      write(info_str,"(2X,A,F25.9,A)") &
        "| Final zero-broadening corrected energy (caution - metals only) :", &
        (total_energy + 2 * &
        (dble(n_methfessel_paxton + 1) / &
        dble(n_methfessel_paxton + 2)) * entropy_correction) &
        * hartree, " eV"
    else
      write(info_str,"(2X,A,F25.9,A)") &
        "| Final zero-broadening corrected energy (caution - metals only) :", &
        (total_energy + entropy_correction) * hartree, " eV"
    end if
    call localorb_info(info_str)

    ! VB: Please add further total energy values here with explanation above.
    !     Any of those further values should be placed in an "if" statement,
    !     i.e., only write it if it is actually computed (example MP2).

    if (post_scf_total_energy /= 0.d0) then
      write(info_str,"(2X,A,F25.9,A)") &
        "| Total energy after the post-s.c.f. correlation calculation     :", &
        post_scf_total_energy * hartree, " eV"
     call localorb_info(info_str)
    end if
    if (SOC_non_sc_total_energy /= 0.d0) then
      write(info_str,"(2X,A,F25.9,A)") &
        "| Total energy after non-SC SOC correction (DO NOT USE)          :", &
        SOC_non_sc_total_energy * hartree, " eV"
     call localorb_info(info_str)
    end if

    write(info_str,"(2X,A,F25.9,A)") &
      "| For reference only, the value of 1 Hartree used in FHI-aims is :", &
      hartree, " eV"
    call localorb_info(info_str)

    if (n_periodic /= 0) then
      write(info_str,"(2X,A)") &
        "| For reference only, the overall average (free atom contribution "
      call localorb_info(info_str)
      write(info_str,"(2X,A)") &
        "| + realspace contribution) of the electrostatic potential after  "
      call localorb_info(info_str)
      write(info_str,"(2X,A,F25.9,A)") &
        "| s.c.f. convergence is                                          :", &
        (average_free_es_pot+average_potential)*hartree, " eV"
      call localorb_info(info_str)
    end if

    call localorb_info("")

    ! Please add definitions of any given values below:

    write(info_str,"(2X,A)") "Before relying on these values, please be sure to understand exactly which"
    call localorb_info(info_str)
    write(info_str,"(2X,A)") "total energy value is referred to by a given number. Different objects may"
    call localorb_info(info_str)
    write(info_str,"(2X,A)") "all carry the same name 'total energy'. Definitions:"
    call localorb_info(info_str)
    write(info_str,"(2X,A)") ""
    call localorb_info(info_str)

    write(info_str,"(2X,A)") "Total energy of the DFT / Hartree-Fock s.c.f. calculation:"
    call localorb_info(info_str)
    write(info_str,"(2X,A)") "| Note that this energy does not include ANY quantities calculated after the"
    call localorb_info(info_str)
    write(info_str,"(2X,A)") "| s.c.f. cycle, in particular not ANY RPA, MP2, etc. many-body perturbation terms."
    call localorb_info(info_str)
    if (flag_rel == REL_zora) then
      write(info_str,"(2X,A)") "| Scalar relativistic effects have been included in the 'scaled ZORA' approximation."
      call localorb_info(info_str)
    end if

    call localorb_info("")

    write(info_str,"(2X,A)") "Final zero-broadening corrected energy: "
    call localorb_info(info_str)

    write(info_str,"(2X,A)") "| For metallic systems only, a broadening of the occupation numbers at the Fermi"
    call localorb_info(info_str)

    write(info_str,"(2X,A)") "| level can be extrapolated back to zero broadening by an electron-gas inspired"
    call localorb_info(info_str)

    write(info_str,"(2X,A)") "| formula. For all systems that are not real metals, this value can be"
    call localorb_info(info_str)

    write(info_str,"(2X,A)") "| meaningless and should be avoided."
    call localorb_info(info_str)

    write(info_str,"(2X,A)") ""
    call localorb_info(info_str)

    if (post_scf_total_energy /= 0.d0) then
      write(info_str,"(2X,A)") &
        "Total energy of the post-s.c.f. correlation calculation        :"
      call localorb_info(info_str)
      if(RI_type == RI_LVL) then
        write(info_str,"(2X,A)") &
          "| Resolution of identity used for the 2-electron Coulomb operator:          RI_LVL_fast"
      else if(RI_type == RI_V) then
        write(info_str,"(2X,A)") &
          "| Resolution of identity used for the 2-electron Coulomb operator:          RI_V"
      end if
      call localorb_info(info_str)
      write(info_str,"(2X,A)") "| The actual meaning of this energy depends on the setting for the 'total_energy_method' "
      call localorb_info(info_str)

      write(info_str,"(2X,A)") "| keyword (MP2, RPA, rPT2, etc). If not sure, compare to the results listed above."
      call localorb_info(info_str)

      write(info_str,"(2X,A)") ""
      call localorb_info(info_str)
    end if

    write(info_str,"(A)") "------------------------------------------------------------------------------"
    call localorb_info(info_str)

  end if

end subroutine final_energy_output
!******
