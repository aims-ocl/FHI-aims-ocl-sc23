!****s* FHI-aims/lowdin_analysis
!  NAME
!   lowdin_analysis
!  SYNOPSIS

      subroutine lowdin_analysis ( converged_scf )

!  PURPOSE
!  Wrapper around actual Loewdin analysis that could also be
!  called from within scf_solver, e.g. per iteration for debugging purposes.
!  This also calls different versions of the Loewdin analysis lapack or scalapack,
!  depending do we have normal form of KS eigenvectors or only scalapack form.
!
!  USES

      use physics
      use localorb_io
      use runtime_choices
      use scalapack_wrapper
      use mpi_tasks, only: aims_stop_coll
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

      character*40 :: filename
      character*100 :: info_str

      character(*), parameter :: func = 'lowdin_analysis'


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
        '* Conducting final Loewdin analysis from non-converged wave functions - '
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
        call aims_stop_coll('Loewdin analysis needs collect_eigenvectors (unless Scalapack is used)',func)

      write(filename,'(A12)') 'Loewdin.out'

      call output_lowdin  & 
           ( KS_eigenvalue, KS_eigenvector, KS_eigenvector_complex, occ_numbers, overlap_matrix, filename, &
           chemical_potential, n_electrons)
         
    end subroutine lowdin_analysis
!******
