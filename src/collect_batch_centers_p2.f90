!****s* FHI-aims/collect_batch_centers_p2
!  NAME
!    collect_batch_centers_p2
!  SYNOPSIS

      subroutine collect_batch_centers_p2 ( n_compute, i_basis, &
        n_basis_list, n_atom_list, inv_list, n_batch_centers, batch_center )


!  PURPOSE
!     For a given batch of integration points, we tabulate all integration centers
!     that are ever relevant in current batch, to remove any O(N^2) steps later
!     in prune_radial_basis ...
!
!  USES

      use dimensions
      use basis
      use grids
      use geometry
      use pbc_lists
      implicit none

!  ARGUMENTS

      integer :: n_basis_list, n_atom_list
      integer :: inv_list(n_centers)
      integer :: n_compute
      integer :: i_basis(n_compute)
      integer :: n_batch_centers
      integer :: batch_center(n_atom_list)

!  INPUTS
!   o n_basis_list -- total number of basis functions
!   o n_atom_list -- number of relevant atoms
!   o inv_list -- inverse citing list of atoms
!   o n_compute -- number of relevant basis functions
!   o i_basis -- list of relevant basis functions
! 
!  OUTPUT
!   o n_batch_centers -- number of batch centers 
!   o batch_center -- list of batch centers
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





!  local variables

!     counters

      integer :: i_compute

      integer :: i_basis_1
      integer :: i_center

      integer :: i_batch_center
      integer :: i_batch_center_2

      logical :: found

!     begin work

      n_batch_centers = 0
!      batch_center = 0

      do i_compute = 1, n_compute, 1
        i_basis_1 = i_basis(i_compute)

        ! number of center associated with current basis function
        i_center = inv_list(Cbasis_to_center(i_basis_1))

!test
!       write(use_unit,*) i_compute, i_center
!test end

        ! check if we already know of this center, else add it to the list

        found = .false.
        i_batch_center = 0
        do while ( (.not.found) .and. (i_batch_center.lt.n_batch_centers) )
          i_batch_center = i_batch_center+1

          if (batch_center(i_batch_center).eq.i_center) then
            ! atom already in list - no need to search further
            found = .true.

          else if (batch_center(i_batch_center).gt.i_center) then
            ! atom was apparently not in the list, since we create a list in
            ! order of increasing center index here ...

            n_batch_centers = n_batch_centers+1
            do i_batch_center_2 = n_batch_centers, i_batch_center+1, -1
              batch_center(i_batch_center_2) = batch_center(i_batch_center_2-1)
            enddo
            batch_center(i_batch_center) = i_center

            found = .true.

          end if

        enddo

        if (.not.found) then
          ! Atom was not yet in list and must have a higher index than any atom already in list

            n_batch_centers = n_batch_centers+1
            batch_center(n_batch_centers) = i_center
         
        end if
!test
!        write(use_unit,*) n_batch_centers
!        write(use_unit,*) batch_center(1:n_batch_centers)
!        write(use_unit,*)
!test end

      enddo

    end subroutine collect_batch_centers_p2
!---------------------------------------------------------------------
!******
