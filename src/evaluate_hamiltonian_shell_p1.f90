!****s* FHI-aims/evaluate_hamiltonian_shell_p1
!  NAME
!    evaluate_hamiltonian_shell_p1
!  SYNOPSIS

subroutine evaluate_hamiltonian_shell_p1 &
     ( n_points, partition,  n_compute,  H_times_psi, &
     n_basis_list, wave, hamiltonian_shell &
     )

!  PURPOSE
!  Subroutine evaluate_hamiltonian_shell evaluates the Hamiltonian integral
!  contribution of several integration points.
!
!  USES

  use dimensions
  use runtime_choices
  implicit none

!  ARGUMENTS

!  imported variables

      integer :: n_points
      integer :: n_compute
      integer :: n_basis_list
      real*8  :: partition(n_points)
      real*8  :: wave(n_basis_list,n_points)
      real*8  :: H_times_psi(n_basis_list,n_points)
      real*8  :: hamiltonian_shell( n_compute,n_compute )

!  INPUTS
!   o  n_points     -- number of grid points in this grid batch
!   o  n_compute    -- number of non-zero basis functions in this grid batch
!   o  n_basis_list -- the total number  of basis functions
!   o  partition    -- values of partition function in this grid batch
!   o  wave         -- values of basis functions in this grid batch
!   o  H_times_psi  -- hamiltonian times basis functions in this grid batch
!
!  OUTPUT
!   o  hamiltonian_shell -- ( wave * hamiltonian * wave)
!                           the results of the integration.
!                            
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

      real*8 wave_compute(n_compute,n_points)
      
!     counters

      integer :: i_point, i_compute

!     begin work

!     Condense basis functions to only those that are used at present 
!     integration points
      do i_point = 1, n_points, 1
         wave_compute(1:n_compute, i_point) =  &
              partition(i_point)*wave(1:n_compute, i_point)
      enddo

      !do i_point = 1, n_points, 1
      !   do i_compute = 1, n_compute, 1
      !      print *, "Wave elem ", &
      !             i_compute + (i_point-1) * n_compute, &
      !             " : ", wave_compute(i_compute,i_point)
      !   end do
      !end do


!     Instead of calculating first 
!     hamiltonian_shell = H_times_psi * wave_compute**T
!     and then 
!     hamiltonian_shell = 0.5 * hamiltonian_shell * hamiltonian_shell**T
!     we can also compute:
!     hamiltonian_shell = 
!        0.5 * (H_times_psi * wave_compute**T + wave_compute*H_times_psi**T)
!     Using DSYRK for this operation results normally in a better performance
!     due to cache reuse etc.

      call dsyr2k('U', 'N', n_compute, n_points, &
            0.5d0, H_times_psi, n_basis_list, &
            wave_compute, n_compute, &
            0.0d0, hamiltonian_shell, n_compute)

    end subroutine evaluate_hamiltonian_shell_p1
!---------------------------------------------------------------------
!******
