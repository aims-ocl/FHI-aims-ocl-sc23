!****s* FHI-aims/check_shell_convergence_p0
!  NAME
!    check_shell_convergence_p0
!  SYNOPSIS

      subroutine check_shell_convergence_p0 ( prev_shell, new_shell, &
             shell_accuracy, shell_converged  )

!  PURPOSE
!  Subroutine check_shell_convergence checks whether the difference 
!  between matrix integral contributions from different shells is
!  below a certain threshold
!
!  Note that we currently check the convergence for the entire basis set,
!  even though only a fraction of basis fns, n_compute, is relevant for each individual shell
!  in large systems.
!
!  The reason is that each full shell can be further subdivided into subdivisions,
!  each of which can account for different numbers of basis functions. On the whole,
!  however, we end up doing too much work ...
!
!  USES

      use dimensions
      implicit none

!  ARGUMENTS

      real*8, dimension( n_hamiltonian_matrix_size ) :: prev_shell
      real*8, dimension( n_hamiltonian_matrix_size ) :: new_shell
      real*8 :: shell_accuracy
      logical :: shell_converged

!  INPUTS
!   o prev_shell -- matrix from previous shell structure
!   o new_shell -- matrix from next shell structure
!   o shell_accuracy -- wanted accuracy 
!   o shell_converged -- ????????????
!
!  OUTPUT
!   o shell_converged -- did the shell converge of not.
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




!  local variables

      ! counters

      integer :: i_index

!  begin work

      i_index = 0

      do while ( (shell_converged) .and. &
                 (i_index < ( n_hamiltonian_matrix_size) ) )

        i_index = i_index+1

        shell_converged = &
        ( abs(new_shell(i_index)-prev_shell(i_index)) .lt. &
          shell_accuracy )

      enddo

!      if (i_index.gt.0) then
!        write(use_unit,*) shell_converged, i_index, 
!     +  abs(new_shell(i_index)-prev_shell(i_index))
!      else
!        write(use_unit,*) shell_converged, i_index
!      end if

      end subroutine check_shell_convergence_p0
!---------------------------------------------------------------------
!******
