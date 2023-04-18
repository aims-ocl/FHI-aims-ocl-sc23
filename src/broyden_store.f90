!****s* FHI-aims/broyden_store
!  NAME
!   broyden_store
!  SYNOPSIS

      subroutine broyden_store &
      ( rho, pulay_saved_iter, previous_rho &
      )

!  PURPOSE
!  Subroutine broyden_store stores the present density rho
!  in a packed array of previous densities for later use by the
!  Pulay mixer
!
!  USES

      use dimensions
      implicit none

!  ARGUMENTS

      real*8, dimension(n_full_points) :: rho
      integer :: pulay_saved_iter

      real*8, dimension(n_full_points, n_max_broyden) :: previous_rho

!  INPUTS
!  o rho -- electron density
!  o pulay_saved_iter -- number of iterations saved in previous_rho
!
!  OUTPUT
!  o previous_rho -- electron density from previous iterations
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





!  local variables

!     counters

      integer :: i_offset
      integer :: i_pulay_store

!  begin work

!     first, shift stored densities to make room for next density
      do i_pulay_store = pulay_saved_iter, 2, -1

         do i_offset = 1, n_full_points, 1
            previous_rho(i_offset, i_pulay_store) = &
                 previous_rho(i_offset, i_pulay_store-1)
         end do

      end do
!     next, store present density in first column of previous_rho
      do i_offset = 1, n_full_points, 1
         previous_rho(i_offset, 1) = rho(i_offset)
      end do
!  end work

      end subroutine broyden_store

!---------------------------------------------------------------------
!******	
