!****s* FHI-aims/reset_grid_partitions
!  NAME
!   reset_grid_partitions
!  SYNOPSIS

subroutine reset_grid_partitions ()

!  PURPOSE
! Subroutine reset_partitions resets everything to do with the batch partitioning of the grid.
! During relaxation, the geometry may change significantly, and mess up the batch partitioning
! from an earlier geometry.
!
!  USES

  use runtime_choices, only: prune_basis_once
  use dimensions
  use mpi_utilities
  use grids

!  ARGUMENTS
!  INPUTS
!    none
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




  implicit none

  ! local variables

  integer :: i_my_batch

  ! begin work

  ! VB: This is a good example why pointers are a realy bad idea for most of Fortran.
  ! We gain nothing from using the pointer, but deallocation it becomes a complicated operation.
  ! I would have liked to just deallocate the pointer, but that generates an error message for me.

  ! Presumably this pointer consists BOTH of allocatable memory and of non-allocatable memory (an integer).
  ! Thus, it lives both on the stack and in the allocated memory. Ugh.

  ! If I just nullified it, then the allocated memory inside would stay allocated. We'd have acompletely unobvious memory leak.
  ! So I need to deallocate the allocatable part of the pointer separately (MANY allocation statements) and then
  ! nullify the pointer iself. And even then I am not sure I did it all correctly.

  ! WPH: Volker's comment above was committed in June 2007 (and is likely older.)  Ville moved to
  !      using deallocate() in commit 887dc80f later in June 2007 to fix a bug with Compaq Fortran,
  !      which was then undone by Daniel's commit 7b04791e in March 2013 with no comment.
  !
  !      As Volker surmised, this does show up in valgrind as a memory leak and it WILL hurt us when
  !      we run aims as a library.  Deallocating batches does not generate a error message for me in
  !      January 2019, so I am assuming that Volker hit a compiler bug and reinstating Ville's change
  !      until someone provides evidence that a memory leak is acceptable relative to the alternatives.

  do i_my_batch = 1, n_my_batches, 1
!      if ( allocated(batches(i_my_batch) % points) ) then
!DB
!         nullify(batches(i_my_batch) % points)
        deallocate( batches(i_my_batch) % points )
!      end if
      if (prune_basis_once) then
!        if ( allocated(batches(i_my_batch) %batch_i_basis ) ) then
!DB
!           nullify(batches(i_my_batch)%batch_i_basis )
           deallocate( batches(i_my_batch)%batch_i_basis )
!        endif
      end if
  end do
!  nullify(batches)
  deallocate(batches)

  if (allocated(batch_task_list)) then
    deallocate(batch_task_list)
  end if

  ! in module grids: flag grids as not yet partitioned
  grid_partitioned = .false.

end subroutine reset_grid_partitions
!******	
