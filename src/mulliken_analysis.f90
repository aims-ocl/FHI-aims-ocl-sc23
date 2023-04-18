!****s* FHI-aims/mulliken_analysis
!  NAME
!   mulliken_analysis
!  SYNOPSIS

      subroutine mulliken_analysis ( converged_scf )

!  PURPOSE
!  Wrapper around actual Mulliken analysis that could also be
!  called from within scf_solver, e.g. per iteration for debugging purposes.
!  This also calls different versions of the Mulliken analysis lapack or scalapack,
!  depending do we have normal form of KS eigenvectors or only scalapack form.
!
!  USES

      use physics,         only : chemical_potential, chemical_potential_soc, KS_eigenvalue, &
                                  KS_eigenvalue_soc_perturbed, KS_eigenvector, &
                                  KS_eigenvector_complex, KS_eigenvector_soc_perturbed, &
                                  n_electrons, occ_numbers, occ_numbers_soc, overlap_matrix
      use localorb_io,     only : localorb_info
      use mpi_tasks,       only : aims_stop_coll
      use runtime_choices, only : collect_eigenvectors, use_scalapack
      use dimensions,      only : calculate_perturbative_soc, n_basis, n_states, n_spin
      use dimensions_soc,  only : n_basis_soc, n_saved_states_soc
      use scalapack_soc,   only : mxld_soc_vec, mxcol_soc_vec
      use timing,          only : get_times, get_timestamps, tot_time_mulliken, &
                                  tot_clock_time_mulliken
      use mulliken,        only : output_mulliken, MULLIKEN_SR, MULLIKEN_SOC
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

      character*40  :: filename
      character*100 :: info_str
      real*8        :: time_mulliken_sr,  clock_time_mulliken_sr, &
                       time_mulliken_soc, clock_time_mulliken_soc
      real*8, dimension(1,1,1,1) :: dummy

      character(*), parameter :: func = 'mulliken_analysis'

      write (info_str,'(2X, A)') ""
      call localorb_info(info_str)
      write (info_str,'(2X, A)') "------------------------------------------------------------"
      call localorb_info(info_str)
      write (info_str,'(2X, A)') "                 Starting Mulliken Analysis                 "
      call localorb_info(info_str)
      write (info_str,'(2X, A)') "------------------------------------------------------------"
      call localorb_info(info_str)

      ! begin work

      if (.not.converged_scf) then

        write(info_str,'(2X,A)') &
        ''
        call localorb_info(info_str)

        write(info_str,'(1X,A)') &
        ' ***************************************************************************'
        call localorb_info(info_str)

        write(info_str,'(1X,A)') &
        '* Warning: scf cycle was not converged!'
        call localorb_info(info_str)

        write(info_str,'(1X,A)') &
        '* Conducting final Mulliken analysis from non-converged wave functions - '
        call localorb_info(info_str)

        write(info_str,'(1X,A)') &
        '* the result could be complete rubbish.'
        call localorb_info(info_str)

        write(info_str,'(1X,A)') &
        '* Do not trust it without further testing.'
        call localorb_info(info_str)

        write(info_str,'(1X,A)') &
        ' ***************************************************************************'
        call localorb_info(info_str)

        write(info_str,'(2X,A)') &
        ''
        call localorb_info(info_str)

      end if


      ! Safety only, should be enforced in read_control:
      if(.not. collect_eigenvectors .and. .not. use_scalapack) &
        call aims_stop_coll('Mulliken analysis needs collect_eigenvectors (unless Scalapack is used)',func)

      if (calculate_perturbative_soc) then 
        call get_timestamps( time_mulliken_sr, clock_time_mulliken_sr )
        write(filename,'(A19)') 'Mulliken.out.no_soc'
        call output_mulliken  & 
             ( n_basis, n_states, n_spin, KS_eigenvector, KS_eigenvector_complex, &
               overlap_matrix, n_states, KS_eigenvalue, occ_numbers, chemical_potential, n_electrons, &
               filename, MULLIKEN_SR )
        call get_times( time_mulliken_sr, clock_time_mulliken_sr, &
                             tot_time_mulliken, tot_clock_time_mulliken, .true. )
       
        call get_timestamps( time_mulliken_soc, clock_time_mulliken_soc )

        write(filename,'(A12)') 'Mulliken.out'
        if (use_scalapack) then
          call output_mulliken  & 
               ( mxld_soc_vec, mxcol_soc_vec, 1, dummy, KS_eigenvector_soc_perturbed, &
                 overlap_matrix, n_saved_states_soc, KS_eigenvalue_soc_perturbed, occ_numbers_soc, &
                 chemical_potential_soc, n_electrons, &
                 filename, MULLIKEN_SOC )
        else
          call output_mulliken  & 
               ( n_basis_soc, n_saved_states_soc, 1, dummy, KS_eigenvector_soc_perturbed, &
                 overlap_matrix, n_saved_states_soc, KS_eigenvalue_soc_perturbed, occ_numbers_soc, &
                 chemical_potential_soc, n_electrons, &
                 filename, MULLIKEN_SOC )
        end if
        call get_times( time_mulliken_soc, clock_time_mulliken_soc, &
                             tot_time_mulliken, tot_clock_time_mulliken, .true. )

      else
        call get_timestamps( time_mulliken_sr, clock_time_mulliken_sr )
        write(filename,'(A12)') 'Mulliken.out'
        call output_mulliken  & 
             ( n_basis, n_states, n_spin, KS_eigenvector, KS_eigenvector_complex, &
               overlap_matrix, n_states, KS_eigenvalue, occ_numbers, chemical_potential, n_electrons, &
               filename, MULLIKEN_SR )
        call get_times( time_mulliken_sr, clock_time_mulliken_sr, &
                             tot_time_mulliken, tot_clock_time_mulliken, .true. )
      end if

      write(info_str, *)  " Mulliken analysis                                       :  max(cpu_time)    wall_clock(cpu1)"
      call localorb_info( info_str )
      if (calculate_perturbative_soc) then
        write(info_str, "(2X, A,F15.3,F20.3)")  "| Time for Mulliken analysis (scalar relativistic)      :", &
            time_mulliken_sr, clock_time_mulliken_sr
        call localorb_info( info_str )
        write(info_str, "(2X, A,F15.3,F20.3)")  "| Time for Mulliken analysis (spin-orbit coupled)       :", &
            time_mulliken_soc, clock_time_mulliken_soc
        call localorb_info( info_str )
      else
        write(info_str, "(2X, A,F15.3,F20.3)")  "| Time for Mulliken analysis                            :", &
            time_mulliken_sr, clock_time_mulliken_sr
        call localorb_info( info_str )
      end if
      write(info_str, "(2X, A,F15.3,F20.3)")  "| Total time for Mulliken analysis                      :", &
          tot_time_mulliken, tot_clock_time_mulliken
      call localorb_info( info_str )
      write(info_str, '(A)') "------------------------------------------------------------"
      call localorb_info( info_str )

    end subroutine mulliken_analysis
!******
