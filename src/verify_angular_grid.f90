!****s* FHI-aims/verify_angular_grid
!  NAME
!   verify_angular_grid
!  SYNOPSIS

      subroutine verify_angular_grid &
        ( n_angular, l_max, angular_new, flag_verify )

!  PURPOSE
!  The subroutine ensures the consistency of the
!  angular integration grid for the Hartree potential computation.
!
!  Notice that this is only a first line of defense - we still underestimate
!  the needed angular grids for accurate integrations of the Hamiltonian matrix
!  elements
!
!  USES
      use mpi_tasks, only : STDERR
      implicit none
!  ARGUMENTS

        integer, intent(in) :: n_angular
        integer, intent(in) :: l_max
        integer, intent(out) :: angular_new
        logical :: flag_verify

!  INPUTS
!  o n_angular -- Number of angular points in a grid
!  o l_max -- Desired integration accuracy for specified angular grids
!
!  OUTPUT
!  o n_angular_new --  Modified to higher value if required by l_max
!  o flag_verify --  set to true if a higher angular grid density was requested
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








!  begin work

        ! The actual grid is one size larger than needed
        ! to integrate l_max correctly. This is a heuristic fix to
        ! not only guarantee converged energies but forces as well ...

        flag_verify = .false.

        if (l_max.le.3) then
          if (n_angular.lt.14) then
            angular_new = 14
            flag_verify = .true.
          end if
        else if (l_max.le.5) then
          if (n_angular.lt.26) then
            angular_new = 26
            flag_verify = .true.
          end if
        else if (l_max.le.7) then
          if (n_angular.lt.50) then
            angular_new = 50
            flag_verify = .true.
          end if
        else if (l_max.le.11) then
          if (n_angular.lt.110) then
            angular_new = 110
            flag_verify = .true.
          end if
        else if (l_max.le.17) then
          if (n_angular.lt.194) then
            angular_new = 194
            flag_verify = .true.
          end if
        else if (l_max.le.23) then
          if (n_angular.lt.302) then
            angular_new = 302
            flag_verify = .true.
          end if
        else if (l_max.le.29) then
          if (n_angular.lt.434) then
            angular_new = 434
            flag_verify = .true.
          end if
        else if (l_max.le.35) then
          if (n_angular.lt.590) then
            angular_new = 590
            flag_verify = .true.
          end if
        else if (l_max.le.41) then
          if (n_angular.lt.770) then
            angular_new = 770
            flag_verify = .true.
          end if
        else if (l_max.le.47) then
          if (n_angular.lt.974) then
            angular_new = 974
            flag_verify = .true.
          end if
        else if (l_max.le.53) then
          if (n_angular.lt.1202) then
            angular_new = 1202
            flag_verify = .true.
          end if
        else if (l_max.le.59) then
          ! Highest Delley-supplied grid
          ! we should go one higher here, ideally ... but then,
          ! if this is not sufficient, what is?
          if (n_angular.lt.1202) then
            angular_new = 1202
            flag_verify = .true.
          end if
        else
          write(STDERR,'(1X,A)') &
          "* Verification of required angular grid density: "
          write(STDERR,'(1X,A)') &
          "* Your angular grid requirements are excessive. "
          write(STDERR,'(1X,A,A)') &
          "* This may be due to an unreasonably high expansion ", &
          "in the Hartree potential or in the RI_method 'V'."
          write(STDERR,'(1X,A,A)') &
          "* Please modify your input file control.in before ", &
          "proceeding - or disable this warning in the source code."
          stop
        end if

      end subroutine verify_angular_grid
!******		
