!!! YY interface to rlsy
!!! to make sure the datastructure of the batch the same
module batch_rlsy_interface
  use grids, only: batch_of_points, grid_point

  public batch_rlsy
  type batch_rlsy

      integer :: n_my_batches
      integer :: n_full_points

      type (batch_of_points), pointer :: batches(:)

      real*8, pointer :: partition_tab(:)
      real*8, pointer :: hartree_potential(:)
      real*8, pointer :: rho(:,:)
      real*8, pointer :: rho_gradient(:,:,:)
      real*8, pointer :: kinetic_density(:,:)
      integer, allocatable, dimension(:)              :: map_full_to_irr_icpu
      integer, allocatable, dimension(:)              :: map_full_to_irr_ipoint
      integer, allocatable, dimension(:)              :: map_full_to_irr_op
      integer, allocatable, dimension(:,:,:)          :: dict_at_rad_ang_to_icpu
      integer, allocatable, dimension(:,:,:)          :: dict_at_rad_ang_to_irr
      integer, allocatable, dimension(:,:,:)          :: dict_at_rad_ang_to_op

  end type

  type(batch_rlsy), public, target :: &
      batch_irr

  ! public routines
  public :: compute_batch_rlsy_interface,compute_batch_rlsy_full_to_irr,&
            compute_batch_rlsy_full_to_irr_points, &
            compute_batch_rlsy_full_to_irr_points_gradient , &
            compute_batch_rlsy_irr_to_full_points, &
            compute_batch_rlsy_irr_to_full_points_gradient

contains

  subroutine compute_batch_rlsy_interface
      use grids, only: batch_of_points, grid_point, batches
      use rlsy_interface, only: rlsy_h
      use dimensions, only:n_full_points, n_my_batches, n_max_batch_size
      use physics, only:partition_tab
      use mpi_tasks
      use synchronize_mpi_basic, only:sync_find_max, sync_vector_integer

      integer :: mpierr
      Character(len = 20) :: myid_str
      !get_local_points: block
        integer :: n_point_irr_local
      !mapping: block
      integer :: ipt,ip,ib,i,iop,icpu
      integer :: i_my_batch
      integer :: max_index_atom, max_index_radial, max_index_angular
      integer :: max_index_atom_global, max_index_radial_global, max_index_angular_global



  end subroutine compute_batch_rlsy_interface

  subroutine compute_batch_rlsy_full_to_irr(hartree_potential,rho,rho_gradient,kinetic_density)
      use grids, only: batch_of_points, grid_point, batches
      use rlsy_interface, only: rlsy_h
      use dimensions, only:n_full_points, n_my_batches
      real*8, intent(IN) :: hartree_potential(:)
      real*8, intent(IN) :: rho(:,:)
      real*8, intent(IN) :: rho_gradient(:,:,:)
      real*8, intent(IN) :: kinetic_density(:,:)
      !mapping: block
      integer :: ipt,ip,ib,i
      integer :: i_my_batch
      integer :: max_index_atom, max_index_radial, max_index_angular

      !mapping: block
      !end block mapping

  end subroutine compute_batch_rlsy_full_to_irr

  subroutine compute_batch_rlsy_full_to_irr_points(full_points, irr_points)
      use grids, only: batch_of_points, grid_point, batches
      use rlsy_interface, only: rlsy_h
      use dimensions, only:n_full_points, n_my_batches
      use mpi_tasks
      real*8, intent(IN) :: full_points(:)
      real*8, intent(OUT) :: irr_points(:)
      !mapping_mpi: block
      integer :: send_count(n_tasks), send_displ(n_tasks), i_count(n_tasks)
      integer :: recv_count
      integer :: i_my_batch, icpu_irr_points, i_cpu
      integer, allocatable :: send_irr_index(:), recv_irr_index(:)
      real*8, allocatable :: send_irr_values(:), recv_irr_values(:)


  end subroutine compute_batch_rlsy_full_to_irr_points

  subroutine compute_batch_rlsy_full_to_irr_points_gradient(full_points, irr_points)
      use grids, only: batch_of_points, grid_point, batches
      use rlsy_interface, only: rlsy_h
      use dimensions, only:n_full_points, n_my_batches
      real*8, intent(IN) :: full_points(:,:)
      real*8, intent(OUT) :: irr_points(:,:)
      !mapping: block


  end subroutine compute_batch_rlsy_full_to_irr_points_gradient

  subroutine compute_batch_rlsy_irr_to_full_points(irr_points, full_points)
      use grids, only: batch_of_points, grid_point, batches
      use rlsy_interface, only: rlsy_h
      use dimensions, only:n_full_points, n_my_batches
      use mpi_tasks
      real*8, intent(IN) :: irr_points(:)
      real*8, intent(OUT) :: full_points(:)
      !mapping_mpi: block
  end subroutine compute_batch_rlsy_irr_to_full_points

  subroutine compute_batch_rlsy_irr_to_full_points_gradient(irr_points, full_points)
      use grids, only: batch_of_points, grid_point, batches
      use rlsy_interface, only: rlsy_h
      use dimensions, only:n_full_points, n_my_batches
      real*8, intent(IN) :: irr_points(:,:)
      real*8, intent(OUT) :: full_points(:,:)
      !mapping: block
      !end block mapping

  end subroutine compute_batch_rlsy_irr_to_full_points_gradient

end module batch_rlsy_interface
