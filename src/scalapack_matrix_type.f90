!****h* FHI-aims/scalapack_matrix_type
!  NAME
!    Routines for creating/setting scalapack matrices
!  SYNOPSIS
module scalapack_matrix_type

  use mpi_tasks

  implicit none
  
  private

  public :: matrix_structure_type
  public :: create_matrix_structure

! **************************************************************************************************
!> brief keeps the structure of a full matrix
!  o nrow_global/local number of rows of the matrix 
!  o ncol_global/local number of column of the matrix 
!  o nrow_block number of rows of a scalapack block
!  o ncol_block number of columns of a scalapack block
!  o icontxt blacs context
!  o first_p_pos position of the first processor (for scalapack)
!  o local_leading_dimension leading dimension of the data that is stored on this processor
!  o descriptor the scalapack descriptor of the matrices
! **************************************************************************************************
  type matrix_structure_type
    integer                 :: nrow_global, ncol_global
    integer                 :: nrow_local, ncol_local
    integer                 :: nrow_block, ncol_block
    integer                 :: icontxt      !blacs context
    integer                 :: local_leading_dimension
    integer, dimension(2)   :: first_p_pos
    integer, dimension(9)   :: descriptor
  end type matrix_structure_type

! **************************************************************************************************

contains

! **************************************************************************************************
!> brief creates the structure of a full (parallel) matrix
! **************************************************************************************************
   subroutine create_matrix_structure(fmstruct, nrow_global, ncol_global, nrow_local, ncol_local,&
                                     nrow_block, ncol_block, icontxt, first_p_pos)

      type(matrix_structure_type), intent(inout)       :: fmstruct
      integer, intent(in)                              :: nrow_global, ncol_global
      integer, intent(in)                              :: nrow_local, ncol_local
      integer, intent(in)                              :: nrow_block, ncol_block
      integer, intent(in)                              :: icontxt 
      integer, dimension(2), &
        intent(in), optional                           :: first_p_pos

      integer                                          :: stat     
 
      fmstruct%nrow_global = nrow_global
      fmstruct%ncol_global = ncol_global
      fmstruct%nrow_local =  nrow_local
      fmstruct%ncol_local =  ncol_local
      fmstruct%nrow_block  = nrow_block
      fmstruct%ncol_block  = ncol_block
      fmstruct%icontxt     = icontxt

      if(present(first_p_pos)) then
        fmstruct%first_p_pos(:) = first_p_pos 
      else
        fmstruct%first_p_pos(:) = (/0, 0/)     
      endif
  
      fmstruct%local_leading_dimension = MAX(1,fmstruct%nrow_local)

      call descinit(fmstruct%descriptor, fmstruct%nrow_global, &
                    fmstruct%ncol_global, fmstruct%nrow_block, &
                    fmstruct%ncol_block, fmstruct%first_p_pos(1), &
                    fmstruct%first_p_pos(2), fmstruct%icontxt, &
                    fmstruct%local_leading_dimension, stat)
      
      
   end subroutine create_matrix_structure

end module scalapack_matrix_type
