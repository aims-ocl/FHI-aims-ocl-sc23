!****s* FHI-aims/prune_basis_v1
!  NAME
!   prune_basis_v1
!  SYNOPSIS

subroutine prune_basis_v1 &
     ( dist_tab, n_compute, i_basis &
     )

!  PURPOSE
!
!     Reduces full basis to relevant functions only, for a given set of
!     integration points
!
!     The list of non-zero integration points is built incrementally by
!     calling this procedure subsequently for all integration points.
!
!  USES
!
  use dimensions
  use basis
  use grids
  use geometry
  implicit none

!  ARGUMENTS

  real*8 dist_tab(n_atoms)
  integer :: n_compute
  integer :: i_basis(n_basis)
  

!  INPUTS
!    o dist_tab -- distance to atoms
!
!  OUTPUT
!    o n_compute -- number of relevant atoms
!    o i_basis -- list of relevant atoms
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

      integer :: i_basis_1
      integer :: i_compute
      integer :: i_compute_1

      logical :: flag = .false.

!     begin work

!     tabulate total wave function value for each basis function

      if (n_compute.eq.0) then
!       this is the first point of the batch - simply list all non-zero wave functions

        i_compute = 0

        do i_basis_1 = 1, n_basis, 1

          if (dist_tab(basis_atom(i_basis_1)) .le. &
              outer_radius(basis_fn(i_basis_1)) ) then

!           nonzero basis function - use it ...
            i_compute = i_compute+1
            i_basis(i_compute) = i_basis_1

          end if

        enddo

        n_compute = i_compute

      else
!       this is not the first integration point of the batch - check whether non-zero
!       wave functions are already there, else add them to the list

         i_compute = 0
         do i_basis_1 = 1, n_basis, 1

            if (i_basis(i_compute+1).eq.i_basis_1) then
!     basis function already in list

               i_compute = i_compute+1

            else if (dist_tab(basis_atom(i_basis_1)) .le. &
                    outer_radius(basis_fn(i_basis_1)) ) then
!     basis function not in list - add it to the list of nonzero functions

               i_compute = i_compute+1

               do i_compute_1 = n_compute, i_compute, -1
                  i_basis(i_compute_1+1) = i_basis(i_compute_1)
               enddo

               i_basis(i_compute) = i_basis_1
               n_compute = n_compute+1

            end if

         enddo

      end if

      end subroutine prune_basis_v1
!---------------------------------------------------------------------
!******
