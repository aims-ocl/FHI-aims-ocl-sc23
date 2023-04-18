!****s* FHI-aims/RRS-PBC/get_rrs_pbc_cell_index()
!  NAME
!    get_rrs_pbc_cell_index() 
!  SYNOPSIS

    subroutine get_rrs_pbc_cell_index(i_center,i_equal,n_vec)

!  PURPOSE
!  This routine get the cell index of the equal atoms
!
!  USES

      use localorb_io
      use dimensions
      use geometry
      use runtime_choices
      use mpi_tasks
      implicit none

!  ARGUMENTS

!  INPUTS
!    none
!  OUTPUT
!    none
!  COPYRIGHT
!   
!   
!   
!  HISTORY
!   
!  SOURCE
!

      ! imported variables
      integer                :: i_center, i_equal
      real*8, dimension(3)   :: n_vec


      ! local variables
      character*132 :: rrs_pbc_info
      integer :: i, j
      real*8, dimension(3)   :: tmp_vec, tmp_vec_2

      tmp_vec(:) = &
          coords(:,rrs_pbc_center_atom(1,i_center)) - &
          coords(:,rrs_pbc_equal_atom(1,i_equal,i_center))

      !write(rrs_pbc_info,'(4X,A,3(F10.4))') &
      !    '| Cell vector                          :', (tmp_vec(i),i=1,3,1)
      !call localorb_info(rrs_pbc_info)

      do i=1,3,1
          n_vec(i)=0
          do j=1,3,1
              n_vec(i)=n_vec(i) + &
                  rrs_pbc_inv_lattice_vector(i,j) * tmp_vec(j)
          enddo
      enddo

      !write(rrs_pbc_info,'(4X,A,3(F8.1))') &
      !    '| Associated index                     :', (tmp_vec_2(i),i=1,3,1)
      !call localorb_info(rrs_pbc_info)

      !do i=1,3,1
      !    n_vec(i)=0
      !    do j=1,3,1
      !       n_vec(i) = n_vec(i) + &
      !           rrs_pbc_lattice_vector(i,j) * tmp_vec_2(j) 
      !    enddo
      !enddo
              
      end subroutine

