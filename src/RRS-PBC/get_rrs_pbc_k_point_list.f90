!****s* FHI-aims/RRS-PBC/get_rrs_pbc_k_point_list()
!  NAME
!    get_rrs_pbc_k_point_list() 
!  SYNOPSIS

    subroutine get_rrs_pbc_k_point_list()

!  PURPOSE
!  High-level wrapper around the post RRS-PBC projection depending on the
!  cluster calculation.
!  This routine can only be called after a succeful scf procedure. Necessary
!  conditions for running are:
!  * an converged overlap matrix
!  * an converged hamiltonian matrix,
!  and all of that with the correct array dimensions
!
!  USES

      use localorb_io
      use dimensions
      use physics
      use geometry
      use numerical_utilities
      use basis
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

      ! local variables
      character*132 :: rrs_pbc_info
      character*132 :: func = 'get_rrs_pbc_k_point_list()'
      integer :: i, j, o, p
      integer :: info
      integer :: i_k_point
      real*8  :: rrs_pbc_delta_k(3)
      real*8  :: one = 1.0d0

      ! allocate rrs_pbc_n_k_point_list
      rrs_pbc_n_k_points = 1
      do i = 1, 3, 1
          rrs_pbc_n_k_points = &
              rrs_pbc_n_k_points * rrs_pbc_n_k_points_xyz(i)
      enddo
      allocate(rrs_pbc_k_point_list(rrs_pbc_n_k_points,3),stat=info)
      call check_allocation(info, 'rrs_pbc_k_point_list                  ')
      allocate(rrs_pbc_KS_eigenvalue(rrs_pbc_n_center_basis,n_spin,&
          rrs_pbc_n_k_points),stat=info)
      call check_allocation(info, 'rrs_pbc_KS_eigenvalue                 ')
      allocate(rrs_pbc_band_info(rrs_pbc_n_center_basis,n_spin),stat=info)
      call check_allocation(info, 'rrs_pbc_band_info                     ')

      ! calculate the value of rrs_pbc_delta_k
      do i = 1, 3, 1
          rrs_pbc_delta_k(i) = one/real(rrs_pbc_n_k_points_xyz(i))
      enddo

      write(rrs_pbc_info,'(4X,A,I12)') &
          '| Now build the k_point_list           :', &
          rrs_pbc_n_k_points
      call localorb_info( rrs_pbc_info )
      i_k_point = 1
      do i = 1, rrs_pbc_n_k_points_xyz(1), 1
          do j = 1, rrs_pbc_n_k_points_xyz(2), 1
              do o = 1, rrs_pbc_n_k_points_xyz(3), 1
                  rrs_pbc_k_point_list(i_k_point,1) = &
                      rrs_pbc_delta_k(1)*(i-1)
                  rrs_pbc_k_point_list(i_k_point,2) = &
                      rrs_pbc_delta_k(2)*(j-1)
                  rrs_pbc_k_point_list(i_k_point,3) = &
                      rrs_pbc_delta_k(3)*(o-1)
                  write(rrs_pbc_info,'(4X,A,3(2X,F16.8))') &
                      '|',(rrs_pbc_k_point_list(i_k_point,p),p=1,3,1)
                  call localorb_info( rrs_pbc_info )
                  i_k_point = i_k_point + 1
              enddo
          enddo
      enddo

      end subroutine

