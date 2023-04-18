!****s* FHI-aims/pulay_grad_store
!  NAME
!   pulay_grad_store
!  SYNOPSIS

      subroutine pulay_grad_store &
      ( rho_gradient, pulay_saved_iter, &
        previous_rho_gradient &
      )

!  PURPOSE
!  Subroutine pulay_grad_store stores the present density gradient
!  rho_gradient in a packed array of previous density gradients for
!  later use by the Pulay mixer
!
!  Nore : This is a trivially modified variant of pulay_store; with a
!         slightly modified interface, can actually replace pulay_store, too ...
!
!  USES

      use dimensions
      implicit none

!  ARGUMENTS

      real*8, dimension(3, n_full_points) :: rho_gradient
      integer :: pulay_saved_iter

      real*8, dimension(3, n_full_points, n_max_pulay) :: previous_rho_gradient


!  INPUT
!  o rho_gradient -- gradient of electron density
!  o pulay_saved_iter -- number of previous gradient saved in previous_rho_gradient
!                
!  OUTPU
!  o previous_rho_gradient -- gradient of electron density from previous iterations
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
      integer :: i_coord
      integer :: i_pulay_store

!  begin work

!     first, shift stored density gradients to make room for next gradient
      do i_pulay_store = pulay_saved_iter, 2, -1

         do i_offset = 1, n_full_points, 1
            do i_coord = 1,3,1
               previous_rho_gradient(i_coord, i_offset, i_pulay_store) = &
                    previous_rho_gradient &
                    (i_coord, i_offset, i_pulay_store-1)
            enddo
         end do

      end do

!     next, store present density in first column of previous_rho
      do i_offset = 1, n_full_points, 1
         do i_coord = 1,3,1
            previous_rho_gradient(i_coord, i_offset, 1) = &
                 rho_gradient(i_coord, i_offset)
         enddo
      end do

!     end work

      end subroutine pulay_grad_store

!---------------------------------------------------------------------
!******	
