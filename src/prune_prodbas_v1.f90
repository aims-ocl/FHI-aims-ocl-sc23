!****s* FHI-aims/prune_prodbas_v1
!  NAME
!   prune_prodbas_v1
!  SYNOPSIS

      subroutine prune_prodbas_v1 &
      ( i_task, dist_tab, &
        n_compute, i_prodbas &
      )

!  PURPOSE
!     Subroutine prune_basis_v1
!     reduces full basis to relevant functions only, for a given set of integration 
!     points
!
!  USES

      use dimensions
      use basis
      use grids
      use geometry
      use prodbas

      implicit none

!  ARGUMENTS
      integer i_task
      real*8 dist_tab(n_atoms)
      integer :: n_compute
      integer :: i_prodbas(n_loc_prodbas)

!  INPUTS
!  o  i_task :: the current processor in cases of parallel calculations
!  o  dist_tab :: the distance of the current points with respect to all the atoms
!  o  n_compute :: the number of relevant basis functions for the part of the 
!          integraton batch which has already been covered up to now.
!  OUTPUT
!  o  n_compute :: the number of relevant basis functions now
!  o  i_prodbas :: integer array, maps the relevant basis functions to their actual
!  o        order in the original basis set.
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

      integer :: i_basis_1
      integer :: i_basbas
      integer :: i_compute
      integer :: i_compute_1

      logical :: flag = .false.

!     begin work

!     tabulate total wave function value for each basis function

      if (n_compute.eq.0) then
!       this is the first point of the batch - simply list all non-zero wave functions

        i_compute = 0

        do i_basis_1 = 1, n_loc_prodbas, 1

          i_basbas = map_prodbas(i_basis_1,i_task)

          if(i_basbas .gt.0 ) then
            if (dist_tab(basbas_atom(i_basbas)) .le. &
                outer_radius_prodbas(i_basis_1) ) then

!            nonzero basis function - use it ...
             i_compute = i_compute+1
             i_prodbas(i_compute) = i_basis_1

            end if
          end if

        enddo

        n_compute = i_compute

      else
!       this is not the first integration point of the batch - check whether non-zero
!       wave functions are already there, else add them to the list

         i_compute = 0
         do i_basis_1 = 1, n_loc_prodbas, 1

          i_basbas = map_prodbas(i_basis_1,i_task)
          if(i_basbas.gt.0) then
            if (i_prodbas(i_compute+1).eq.i_basis_1) then
!     basis function already in list

               i_compute = i_compute+1

            else if (dist_tab(basbas_atom(i_basbas)) .le. &
                    outer_radius_prodbas(i_basis_1) ) then
!     basis function not in list - add it to the list of nonzero functions

               i_compute = i_compute+1

               do i_compute_1 = n_compute, i_compute, -1
                  i_prodbas(i_compute_1+1) = i_prodbas(i_compute_1)
               enddo

               i_prodbas(i_compute) = i_basis_1
               n_compute = n_compute+1

            end if
           endif

         enddo

      end if

      end subroutine prune_prodbas_v1
!---------------------------------------------------------------------
!******
