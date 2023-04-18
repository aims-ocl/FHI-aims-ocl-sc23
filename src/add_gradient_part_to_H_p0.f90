
!****s* FHI-aims/add_gradient_part_to_H_p0
!  NAME
!    add_gradient_part_to_H_p0
!  SYNOPSIS

      subroutine add_gradient_part_to_H_p0 &
           ( n_compute, gradient_basis_wave,  &
           gradient_deriv, H_times_psi )
 
!  PURPOSE
!   Calculates (gradient_deriv) dot_product ( gradient_basis_wave) for several
!   grid points at once.
!  
!  USES

      use dimensions
      implicit none

!  ARGUMENTS

      integer :: n_compute
      real*8, dimension(n_compute,3) :: gradient_basis_wave
      real*8, dimension(3) :: gradient_deriv
      real*8 H_times_psi(n_compute)

!  INPUTS
!   o n_compute -- number of nonzero basis functions in this grid batch.
!   o gradient_basis_wave -- gradient of the basis functions 
!   o gradient_deriv -- gradient which is added to H_psi.
!  
!  OUTPUT
!   o H_times_psi -- results of the products are added here.
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
!



!     local variables

!     counters

      integer :: i_compute
      integer :: i_coord

!     begin work

      do i_coord = 1, 3, 1

            H_times_psi(1:n_compute) = H_times_psi(1:n_compute) + &
                 gradient_deriv(i_coord) * &
                 gradient_basis_wave(1:n_compute,i_coord)

      enddo

    end subroutine add_gradient_part_to_H_p0
!---------------------------------------------------------------------
!******
