!****s* FHI-aims/mulliken_analysis
!  NAME
!   mulliken_analysis
!  SYNOPSIS
  subroutine dielectric_analysis ( converged_scf )
!  PURPOSE
!  Wrapper around actual Mulliken analysis that could also be
!  called from within scf_solver, e.g. per iteration for debugging purposes.
!  This also calls different versions of the Mulliken analysis lapack or scalapack,
!  depending do we have normal form of KS eigenvectors or only scalapack form.
!
!  USES
  use physics,               only : chemical_potential, chemical_potential_soc, &
                                    occ_numbers, occ_numbers_soc, &
                                    KS_eigenvalue, KS_eigenvalue_soc_perturbed, &
                                    KS_eigenvector, KS_eigenvector_complex, &
                                    KS_eigenvector_soc_perturbed, &
                                    partition_tab
  use dimensions,            only : n_basis, n_states, n_spin, spin_degeneracy
  use dimensions_soc,        only : n_basis_soc, n_saved_states_soc
  use localorb_io,           only : localorb_info
  use dimensions,            only : calculate_perturbative_soc
  use runtime_choices,       only : use_scalapack
  use load_balancing,        only : use_batch_permutation
  use species_data,          only : l_shell_max
  use KS_optical_properties, only : KS_dielectric_calculation, &
                                    tot_time_real_space_mommat, tot_time_bloch_mommat, &
                                    tot_time_plasma_interband, tot_time_sync_postprocess, &
                                    tot_clock_time_real_space_mommat, tot_clock_time_bloch_mommat, &
                                    tot_clock_time_plasma_interband, tot_clock_time_sync_postprocess
  use KS_optical_properties_tetrahedron, only : KS_dielectric_calculation_tetrahedron
  use timing,                only : tot_time_dielectric, tot_clock_time_dielectric
  use aims_memory_tracking,  only : aims_mem_current_output
  use scalapack_soc,         only : mxld_soc_vec, mxcol_soc_vec
  use scalapack_wrapper,     only : mxld, mxcol, eigenvec, eigenvec_complex
  implicit none
!  ARGUMENTS
  logical :: converged_scf
!  INPUTS
!   o converged_scf -- has the self consistant loop converged or not.
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
! Local variables
  character*100                                 :: info_str

  character(*), parameter :: func = 'dielectric_analysis'
  ! begin work

  write (info_str,'(2X, A)') "------------------------------------------------------------"
  call localorb_info(info_str)
  write (info_str,'(2X, A)') "             Starting Dielectric Calculations               "
  call localorb_info(info_str)
  write (info_str,'(2X, A)') "------------------------------------------------------------"
  call localorb_info(info_str)
  write(info_str,'(2X,A)')   ""
  call localorb_info(info_str)


  if (.not.converged_scf) then
    write(info_str,'(1X,A)') ' ***************************************************************************'
    call localorb_info(info_str)
    write(info_str,'(1X,A)') '* Warning: scf cycle was not converged!'
    call localorb_info(info_str)
    write(info_str,'(1X,A)') '* Conducting dielectric analysis from non-converged wave functions - '
    call localorb_info(info_str)
    write(info_str,'(1X,A)') '* the result could be complete rubbish.'
    call localorb_info(info_str)
    write(info_str,'(1X,A)') '* Do not trust it without further testing.'
    call localorb_info(info_str)
    write(info_str,'(1X,A)') ' ***************************************************************************'
    call localorb_info(info_str)
    write(info_str,'(2X,A)') ''
    call localorb_info(info_str)
  end if

  if (calculate_perturbative_soc) then
    if (use_scalapack) then
      ! KS_eigenvector is a dummy variable here
      call KS_dielectric_calculation( mxld_soc_vec, mxcol_soc_vec, 1, KS_eigenvector, KS_eigenvector_soc_perturbed, &
               n_saved_states_soc, KS_eigenvalue_soc_perturbed, occ_numbers_soc, chemical_potential_soc, 1.0d0, &
               partition_tab, l_shell_max )
    else
      ! KS_eigenvector is a dummy variable here
      call KS_dielectric_calculation( n_basis_soc, n_saved_states_soc, 1, KS_eigenvector, KS_eigenvector_soc_perturbed, &
               n_saved_states_soc, KS_eigenvalue_soc_perturbed, occ_numbers_soc, chemical_potential_soc, 1.0d0, &
               partition_tab, l_shell_max )
    end if
  else
    if (use_scalapack) then
      call KS_dielectric_calculation( mxld, mxcol, n_spin, eigenvec, eigenvec_complex, &
               n_states, KS_eigenvalue, occ_numbers, chemical_potential, spin_degeneracy, &
               partition_tab, l_shell_max )
    else
      call KS_dielectric_calculation( n_basis, n_states, n_spin, KS_eigenvector, KS_eigenvector_complex, &
               n_states, KS_eigenvalue, occ_numbers, chemical_potential, spin_degeneracy, &
               partition_tab, l_shell_max )
      ! YY: I don't want to change the behavior for now
      !call KS_dielectric_calculation_tetrahedron( n_basis, n_states, n_spin, KS_eigenvector, KS_eigenvector_complex, &
      !         n_states, KS_eigenvalue, occ_numbers, chemical_potential, spin_degeneracy, &
      !         partition_tab, l_shell_max )
    end if
  end if

  write(info_str, *)  " Dielectric calculations                                 :  max(cpu_time)    wall_clock(cpu1)"
  call localorb_info( info_str )
  write(info_str, "(2X, A,F15.3,F20.3)")  "| Time for real-space momentum matrix integration       :", &
      tot_time_real_space_mommat, tot_clock_time_real_space_mommat
  call localorb_info( info_str )
  write(info_str, "(2X, A,F15.3,F20.3)")  "| Time for calculating momentum matrix elements         :", &
      tot_time_bloch_mommat, tot_clock_time_bloch_mommat
  call localorb_info( info_str )
  write(info_str, "(2X, A,F15.3,F20.3)")  "| Time for plasma frequency and interband contributions :", &
      tot_time_plasma_interband, tot_clock_time_plasma_interband
  call localorb_info( info_str )
  write(info_str, "(2X, A,F15.3,F20.3)")  "| Time for syncing and postprocessing                   :", &
      tot_time_sync_postprocess, tot_clock_time_sync_postprocess
  call localorb_info( info_str )
  write(info_str, "(2X, A,F15.3,F20.3)")  "| Total time for dielectric calculations                :", &
      tot_time_dielectric, tot_clock_time_dielectric
  call localorb_info( info_str )
  write(info_str, *)
  call localorb_info( info_str )
  call aims_mem_current_output( )
  write(info_str, '(A)') "------------------------------------------------------------"
  call localorb_info( info_str )

end subroutine dielectric_analysis
!******
