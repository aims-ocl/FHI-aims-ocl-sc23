!****h* FHI-aims/scalapack_utils
!  NAME
!    scalapack_utils
!  SYNOPSIS

module scalapack_utils

  !  PURPOSE
  !
  !    This module contains utilities for scalapack like block cyclic
  !    matrix distributions without actually resorting to the scalapack
  !    library.
  !
  !  USES

  implicit none

  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  integer, parameter :: LOCIND_BEFORE = -1
  integer, parameter :: LOCIND_ZERO   = 0
  integer, parameter :: LOCIND_AFTER  = 1

contains

  !----------------------------------------------------------------------------
  !****s* scalapack_utils/sclpck_loc_ind
  !  NAME
  !    sclpck_loc_ind
  !  SYNOPSIS

  integer function sclpck_loc_ind(idx, my_proc, num_procs, nblk, iflag)

    !  PURPOSE
    !
    !    Returns the local index for a given global index.
    !    If the global index has no local index on the
    !    processor my_proc behaviour is defined by iflag
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: idx, my_proc, num_procs, nblk, iflag

    !  INPUTS
    !
    !    idx     -- Global index
    !    my_proc -- Processor row/column for which to calculate local index
    !    num_procs -- Total number of processors along row/column
    !    nblk    -- Blocksize
    !    iflag -- Controls the behaviour if idx is not on local processor
    !                iflag< 0 : Return last local index before that row/col
    !                iflag==0 : Return 0
    !                iflag> 0 : Return next local index after that row/col
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !    Code is taken from elpa1.f90:local_index() from RJ.
    !  SOURCE

    character(*), parameter :: func = 'sclpck_loc_ind'

    integer iblk, tproc

    iblk = (idx-1)/nblk          ! global block number, 0 based
    tproc = mod(iblk, num_procs) ! responsible proc

    if(tproc == my_proc) then
       ! block is local, always return local row/col number
       sclpck_loc_ind = (iblk/num_procs)*nblk + mod(idx-1,nblk) + 1
    else
       if (iflag == LOCIND_ZERO) then
          sclpck_loc_ind = 0
       else
          sclpck_loc_ind = (iblk/num_procs)*nblk
          if(tproc > my_proc) sclpck_loc_ind = sclpck_loc_ind + nblk
          if (iflag == LOCIND_AFTER) sclpck_loc_ind = sclpck_loc_ind + 1
       endif
    endif

  end function sclpck_loc_ind
  !******
  !----------------------------------------------------------------------------
  !****s* scalapack_utils/sclpck_num_loc_block
  !  NAME
  !    sclpck_num_loc_block
  !  SYNOPSIS

  integer function sclpck_num_loc_block(N, NB, np, myp)

    !  PURPOSE
    !
    !     Figure out how many stripes of the global matrix are relevant for
    !     this given processor given a block-cyclic distribution.  A partial
    !     block is counted whole.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N, NB, np, myp

    !  INPUTS
    !    o N -- Global dimension along given direction
    !    o NB -- Blocking factor along given direction
    !    o np -- Number of processors along given direction
    !    o myp -- Own processor number along given direction
    !  OUTPUTS
    !    o sclpck_num_loc_block -- Number of distinct local blocks
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: num_global_blocks, min_loc_blocks
    character(*), parameter :: func = 'sclpck_num_loc_block'

    num_global_blocks = (N+NB-1) / NB  ! Include last (pot. partial)
    min_loc_blocks = num_global_blocks / np
    if (num_global_blocks - np * min_loc_blocks <= myp) then
       sclpck_num_loc_block = min_loc_blocks
    else
       sclpck_num_loc_block = min_loc_blocks + 1
    end if
  end function sclpck_num_loc_block
  !******
  !----------------------------------------------------------------------------
  !****s* scalapack_utils/sclpck_top_block
  !  NAME
  !    sclpck_top_block
  !  SYNOPSIS

  integer function sclpck_top_block(i_block, N, NB, np, myp)

    !  PURPOSE
    !     For a given local block number along some direction, figure out
    !     where we are in the global matrix.  "top" is the first index of
    !     the block.
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_block, N, NB, np, myp

    !  INPUTS
    !    o i_block -- Local block number
    !    o N -- Global dimension along given direction
    !    o NB -- Blocking factor along given direction
    !    o np -- Number of processors along given direction
    !    o myp -- Own processor number along given direction
    !  OUTPUTS
    !    o sclpck_top_block -- Global index of first block entry 
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_glob_block   ! 0 <= i <= sclpck_num_loc_block()-1
    character(*), parameter :: func = 'sclpck_top_block'

    i_glob_block = (i_block-1)*np + myp
    sclpck_top_block = i_glob_block * NB + 1
  end function sclpck_top_block
  !******
  !----------------------------------------------------------------------------
  !****s* scalapack_utils/sclpck_shp_block
  !  NAME
  !    sclpck_shp_block
  !  SYNOPSIS

  integer function sclpck_shp_block(i_block, N, NB, np, myp)


    !  PURPOSE
    !     For a given local block number along some direction, figure out
    !     the number of elements to be stored locally.  For every block but
    !     the last, this will be the blocking factor NB.
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_block, N, NB, np, myp

    !  SOURCE

    integer :: top, bot
    character(*), parameter :: func = 'sclpck_shp_block'
    top = sclpck_top_block(i_block, N, NB, np, myp)
    bot = min(top + NB, N+1)
    sclpck_shp_block = bot - top

  end function sclpck_shp_block
  !******

end module scalapack_utils
!******
