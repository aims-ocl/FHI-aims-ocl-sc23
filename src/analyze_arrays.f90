!****h* FHI-aims/analyze_arrays
!  NAME
!    analyze_arrays
!  SYNOPSIS

module analyze_arrays

!  PURPOSE
!
! This module analyzes sparsity of 1d-, 2d-, and 3d-arrays
! The user can request printout of the results or return the
! number of nonzeros in an array
! Relevant variables for the subroutines:
! array: the array itself
! dim1, dim2, dim3: dimensions of the array
! threshold: entries with absolute value less than threshold are considered zeros
! writeout: if true print out number of nonzeors and sparsity of the array
! nnz: if present contains the number of nonzeros in the array on return
!
!  USES

  implicit none


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
!******




contains

!****s* analyze_arrays/analyze_1d_array
!  NAME
!    analyze_1d_array
!  SYNOPSIS

  subroutine analyze_1d_array( array, dim1, threshold, writeout, nnz )

!  PURPOSE
!  Analyses and print outs how many number of the table is smaller than  threshold
!
!  USES
    use localorb_io, only: localorb_info, use_unit, OL_norm
    use mpi_tasks, only: myid
    use runtime_choices, only: use_local_index
    implicit none
!  ARGUMENTS


    integer :: dim1
    real*8 :: array(dim1)
    real*8 :: threshold
    logical :: writeout
    integer, optional :: nnz

!  INPUTS
!   o dim1 - dimension of the table
!   o array -- table
!   o threshold -- threshold value
!   o writeout -- is the date printed out?
!   
!  OUTPUT
!   o nnz -- number of numbers in the table smaller than threshold
! 
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




    ! locals
    integer :: i_index
    integer :: nonzeros
    character*100 :: info_str
    real*8 :: sparsity

    nonzeros = 0

    do i_index = 1, dim1, 1
       if ( abs(array(i_index)).gt.threshold ) then
          nonzeros = nonzeros + 1
       end if
    end do

    if (writeout) then
       write(info_str,'(2X,A,I8,A,I8,A)') "| Array has ", nonzeros, &
            " nonzero elements out of ", dim1, " elements"
       call localorb_info(info_str,use_unit,'(A)',OL_norm)
       if(use_local_index) print *,'Task: ',myid,info_str
       sparsity = 1.0 - dble(nonzeros)/dble(dim1)
       write(info_str,'(2X,A,F5.3)') "| Sparsity factor is ", sparsity
       call localorb_info(info_str,use_unit,'(A)',OL_norm)
       if(use_local_index) print *,'Task: ',myid,info_str
    end if

    if (present(nnz)) then
       nnz = nonzeros
    end if

  end subroutine analyze_1d_array

!******
!------------------------------------------------------------------------
!****s* analyze_arrays/analyze_2d_array
!  NAME
!    analyze_2d_array
!  SYNOPSIS

  subroutine analyze_2d_array( array, dim1, dim2, threshold, writeout, nnz )

!  PURPOSE
!    Analyses and print outs how many number of the table is smaller than  threshold
!  USES
     use localorb_io, only: localorb_info, use_unit
     implicit none
!  ARGUMENTS

    integer :: dim1, dim2
    real*8 :: array(dim1,dim2)
    real*8 :: threshold
    logical :: writeout
    integer, optional :: nnz

!  INPUTS
!   o dim1, dim2  - dimension of the table
!   o array -- table
!   o threshold -- threshold value
!   o writeout -- is the date printed out?
!   
!  OUTPUT
!   o nnz -- number of numbers in the table smaller than threshold
! 
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




    ! locals
    integer :: i_index
    integer :: nonzeros
    integer :: nnz_temp
    character*100 :: info_str
    real*8 :: sparsity
    integer :: dim_total

    dim_total = dim1*dim2

    nonzeros = 0

    ! analyze each subarray as a 1d-array
    do i_index = 1, dim2, 1

       call analyze_1d_array(array(:,i_index), dim1, threshold, .false., nnz_temp)
       nonzeros = nonzeros + nnz_temp

    end do

    if (writeout) then
       write(info_str,'(2X,A,I8,A,I8,A)') "| Array has ", nonzeros, &
            " nonzero elements out of ", dim_total, " elements"
       call localorb_info(info_str,use_unit,'(A)')
       sparsity = 1.0 - dble(nonzeros)/dble(dim_total)
       write(info_str,'(2X,A,F5.3)') "| Sparsity factor is ", sparsity
       call localorb_info(info_str,use_unit,'(A)')
    end if

    if (present(nnz)) then
       nnz = nonzeros
    end if

  end subroutine analyze_2d_array
!******
  !------------------------------------------------------------------------
!****s* analyze_arrays/analyze_3d_array
!  NAME
!    analyze_3d_array
!  SYNOPSIS

  subroutine analyze_3d_array( array, dim1, dim2, dim3, threshold, writeout, nnz )


!  PURPOSE
!    Analyses and print outs how many number of the table is smaller than threshold
!  USES
     use localorb_io, only: localorb_info, use_unit
     implicit none
!  ARGUMENTS

    integer :: dim1, dim2, dim3
    real*8 :: array(dim1,dim2,dim3)
    real*8 :: threshold
    logical :: writeout
    integer, optional :: nnz

!  INPUTS
!   o dim1, dim2, dim3  - dimensions of the table
!   o array -- table
!   o threshold -- threshold value
!   o writeout -- is the date printed out?
!  OUTPUT
!   o nnz -- number of numbers in the table smaller than threshold
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



    ! locals
    integer :: i_index
    integer :: nonzeros
    integer :: nnz_temp
    character*100 :: info_str
    real*8 :: sparsity
    integer :: dim_total

    dim_total = dim1*dim2*dim3

    nonzeros = 0

    ! analyze each subarray as a 2d-array
    do i_index = 1, dim3, 1

       call analyze_2d_array(array(:,:,i_index), dim1, dim2, threshold, .false., nnz_temp)
       nonzeros = nonzeros + nnz_temp

    end do

    if (writeout) then
       write(info_str,'(2X,A,I8,A,I8,A)') "| Array has ", nonzeros, &
            " nonzero elements out of ", dim_total, " elements"
       call localorb_info(info_str,use_unit,'(A)')
       sparsity = 1.0 - dble(nonzeros)/dble(dim_total)
       write(info_str,'(2X,A,F5.3)') "| Sparsity factor is ", sparsity
       call localorb_info(info_str,use_unit,'(A)')
    end if

    if (present(nnz)) then
       nnz = nonzeros
    end if

  end subroutine analyze_3d_array

end module analyze_arrays
!******
