!!! YY interface to rlsy
!!! to make sure the datastructure of the batch the same
module batch_rlsy_interface
   use grids, only: batch_of_points, grid_point

   public batch_rlsy
   type batch_rlsy

      integer :: n_my_batches
      integer :: n_full_points

      type(batch_of_points), pointer :: batches(:)

      real*8, pointer :: partition_tab(:)
      real*8, pointer :: hartree_potential(:)
      real*8, pointer :: rho(:, :)
      real*8, pointer :: rho_gradient(:, :, :)
      real*8, pointer :: kinetic_density(:, :)
      integer, allocatable, dimension(:)              :: map_full_to_irr_icpu
      integer, allocatable, dimension(:)              :: map_full_to_irr_ipoint
      integer, allocatable, dimension(:)              :: map_full_to_irr_op
      integer, allocatable, dimension(:, :, :)          :: dict_at_rad_ang_to_icpu
      integer, allocatable, dimension(:, :, :)          :: dict_at_rad_ang_to_irr
      integer, allocatable, dimension(:, :, :)          :: dict_at_rad_ang_to_op

   end type

   type(batch_rlsy), public, target :: &
      batch_irr

   ! public routines
   public :: compute_batch_rlsy_interface, compute_batch_rlsy_full_to_irr, &
             compute_batch_rlsy_full_to_irr_points, &
             compute_batch_rlsy_full_to_irr_points_gradient, &
             compute_batch_rlsy_irr_to_full_points, &
             compute_batch_rlsy_irr_to_full_points_gradient

contains

   subroutine compute_batch_rlsy_interface
      use grids, only: batch_of_points, grid_point, batches
      use rlsy_interface, only: rlsy_h
      use dimensions, only: n_full_points, n_my_batches, n_max_batch_size
      use physics, only: partition_tab
      use mpi_tasks
      use synchronize_mpi_basic, only: sync_find_max, sync_vector_integer

      integer :: mpierr
      Character(len=20) :: myid_str
      !get_local_points: block
      integer :: n_point_irr_local
      !mapping: block
      integer :: ipt, ip, ib, i, iop, icpu
      integer :: i_my_batch
      integer :: max_index_atom, max_index_radial, max_index_angular
      integer :: max_index_atom_global, max_index_radial_global, max_index_angular_global

      batch_irr%n_my_batches = rlsy_h%grid%n_irr_batch

      !get_local_points: block
      n_point_irr_local = 0
      do ib = 1, rlsy_h%grid%n_irr_batch
         n_point_irr_local = n_point_irr_local + rlsy_h%grid%irr_batch(ib)%n_point
      end do
      batch_irr%n_full_points = n_point_irr_local
      !end block get_local_points

      allocate (batch_irr%batches(batch_irr%n_my_batches))
      allocate (batch_irr%partition_tab(batch_irr%n_full_points))
      allocate (batch_irr%hartree_potential(batch_irr%n_full_points))
      allocate (batch_irr%rho(1, batch_irr%n_full_points))
      allocate (batch_irr%rho_gradient(3, 1, batch_irr%n_full_points))
      allocate (batch_irr%kinetic_density(1, batch_irr%n_full_points))

      !mapping: block
      ipt = 0
      max_index_atom = 0
      max_index_radial = 0
      max_index_angular = 0
      bl2: do ib = 1, rlsy_h%grid%n_irr_batch
      do ip = 1, rlsy_h%grid%irr_batch(ib)%n_point
         ipt = ipt + 1
         do i = 1, rlsy_h%grid%irr_batch(ib)%unfold_ctr(ip)
            max_index_atom = MAX(rlsy_h%grid%irr_batch(ib)%unfold_atom(i, ip), max_index_atom)
            max_index_radial = MAX(rlsy_h%grid%irr_batch(ib)%unfold_index_radial(i, ip), max_index_radial)
            max_index_angular = MAX(rlsy_h%grid%irr_batch(ib)%unfold_index_angular(i, ip), max_index_angular)
         enddo
      enddo
      enddo bl2

      call sync_find_max(max_index_atom, max_index_atom_global)
      max_index_atom = max_index_atom_global
      call sync_find_max(max_index_radial, max_index_radial_global)
      max_index_radial = max_index_radial_global
      call sync_find_max(max_index_angular, max_index_angular_global)
      max_index_angular = max_index_angular_global

      allocate (batch_irr%dict_at_rad_ang_to_icpu(max_index_atom, max_index_radial, max_index_angular), stat=info)
      allocate (batch_irr%dict_at_rad_ang_to_irr(max_index_atom, max_index_radial, max_index_angular), stat=info)
      allocate (batch_irr%dict_at_rad_ang_to_op(max_index_atom, max_index_radial, max_index_angular), stat=info)
      batch_irr%dict_at_rad_ang_to_icpu = 0
      batch_irr%dict_at_rad_ang_to_irr = 0
      batch_irr%dict_at_rad_ang_to_op = 0

      ipt = 0

      bl3: do ib = 1, rlsy_h%grid%n_irr_batch
      do ip = 1, rlsy_h%grid%irr_batch(ib)%n_point
         ipt = ipt + 1
         do i = 1, rlsy_h%grid%irr_batch(ib)%unfold_ctr(ip)
            batch_irr%dict_at_rad_ang_to_icpu(rlsy_h%grid%irr_batch(ib)%unfold_atom(i, ip), &
                                              rlsy_h%grid%irr_batch(ib)%unfold_index_radial(i, ip), &
                                              rlsy_h%grid%irr_batch(ib)%unfold_index_angular(i, ip)) = myid
            batch_irr%dict_at_rad_ang_to_irr(rlsy_h%grid%irr_batch(ib)%unfold_atom(i, ip), &
                                             rlsy_h%grid%irr_batch(ib)%unfold_index_radial(i, ip), &
                                             rlsy_h%grid%irr_batch(ib)%unfold_index_angular(i, ip)) = ipt
            batch_irr%dict_at_rad_ang_to_op(rlsy_h%grid%irr_batch(ib)%unfold_atom(i, ip), &
                                            rlsy_h%grid%irr_batch(ib)%unfold_index_radial(i, ip), &
                                            rlsy_h%grid%irr_batch(ib)%unfold_index_angular(i, ip)) = &
               rlsy_h%grid%irr_batch(ib)%unfold_operation(i, ip)
         enddo
      enddo
      enddo bl3

      call sync_vector_integer(batch_irr%dict_at_rad_ang_to_icpu, &
                               max_index_atom*max_index_radial*max_index_angular, &
                               mpi_comm_global)
      call sync_vector_integer(batch_irr%dict_at_rad_ang_to_irr, &
                               max_index_atom*max_index_radial*max_index_angular, &
                               mpi_comm_global)
      call sync_vector_integer(batch_irr%dict_at_rad_ang_to_op, &
                               max_index_atom*max_index_radial*max_index_angular, &
                               mpi_comm_global)

      allocate (batch_irr%map_full_to_irr_icpu(n_full_points), stat=info)
      allocate (batch_irr%map_full_to_irr_ipoint(n_full_points), stat=info)
      allocate (batch_irr%map_full_to_irr_op(n_full_points), stat=info)
      batch_irr%map_full_to_irr_icpu = 0
      batch_irr%map_full_to_irr_ipoint = 0
      batch_irr%map_full_to_irr_op = 0
      i_full_points = 0
      bl4: do i_my_batch = 1, n_my_batches, 1
         ! loop over one batch
         do i_index = 1, batches(i_my_batch)%size, 1

            i_full_points = i_full_points + 1
            batch_irr%map_full_to_irr_icpu(i_full_points) = &
               batch_irr%dict_at_rad_ang_to_icpu(batches(i_my_batch)%points(i_index)%index_atom, &
                                                 batches(i_my_batch)%points(i_index)%index_radial, &
                                                 batches(i_my_batch)%points(i_index)%index_angular)
            batch_irr%map_full_to_irr_ipoint(i_full_points) = &
               batch_irr%dict_at_rad_ang_to_irr(batches(i_my_batch)%points(i_index)%index_atom, &
                                                batches(i_my_batch)%points(i_index)%index_radial, &
                                                batches(i_my_batch)%points(i_index)%index_angular)
            batch_irr%map_full_to_irr_op(i_full_points) = &
               batch_irr%dict_at_rad_ang_to_op(batches(i_my_batch)%points(i_index)%index_atom, &
                                               batches(i_my_batch)%points(i_index)%index_radial, &
                                               batches(i_my_batch)%points(i_index)%index_angular)
         end do
      end do bl4

      bl5: do ib = 1, rlsy_h%grid%n_irr_batch
         batch_irr%batches(ib)%size = rlsy_h%grid%irr_batch(ib)%n_point
         allocate (batch_irr%batches(ib)%points(batch_irr%batches(ib)%size), stat=info)
         do ip = 1, rlsy_h%grid%irr_batch(ib)%n_point
            batch_irr%batches(ib)%points(ip)%coords = rlsy_h%grid%irr_batch(ib)%folded_coordinate(:, ip)
            batch_irr%batches(ib)%points(ip)%index_atom = rlsy_h%grid%irr_batch(ib)%unfold_atom(1, ip)
            batch_irr%batches(ib)%points(ip)%index_radial = rlsy_h%grid%irr_batch(ib)%unfold_index_radial(1, ip)
            batch_irr%batches(ib)%points(ip)%index_angular = rlsy_h%grid%irr_batch(ib)%unfold_index_angular(1, ip)
         enddo
      enddo bl5

      i_full_points = 0
      batch_irr%partition_tab = 0.0d0
      call compute_batch_rlsy_full_to_irr_points(partition_tab, batch_irr%partition_tab)
      ipt = 0
      bl7: do ib = 1, rlsy_h%grid%n_irr_batch
      do ip = 1, rlsy_h%grid%irr_batch(ib)%n_point
         ipt = ipt + 1
         batch_irr%partition_tab(ipt) = batch_irr%partition_tab(ipt) &
                                        *rlsy_h%grid%irr_batch(ib)%unfold_ctr(ip) &
                                        /rlsy_h%spacegroup%n_operation

      enddo
      enddo bl7
      bl8: do ib = 1, rlsy_h%grid%n_irr_batch
         n_max_batch_size = max(n_max_batch_size, rlsy_h%grid%irr_batch(ib)%n_point)
      enddo bl8
      !end block mapping

   end subroutine compute_batch_rlsy_interface

   subroutine compute_batch_rlsy_full_to_irr(hartree_potential, rho, rho_gradient, kinetic_density)
      use grids, only: batch_of_points, grid_point, batches
      use rlsy_interface, only: rlsy_h
      use dimensions, only: n_full_points, n_my_batches
      real*8, intent(IN) :: hartree_potential(:)
      real*8, intent(IN) :: rho(:, :)
      real*8, intent(IN) :: rho_gradient(:, :, :)
      real*8, intent(IN) :: kinetic_density(:, :)
      !mapping: block
      integer :: ipt, ip, ib, i
      integer :: i_my_batch
      integer :: max_index_atom, max_index_radial, max_index_angular

      !mapping: block
      ipt = 0
      i_full_points = 0
      bl6: do i_my_batch = 1, n_my_batches, 1
         ! loop over one batch
         do i_index = 1, batches(i_my_batch)%size, 1
            i_full_points = i_full_points + 1
            i_irr_points = batch_irr%map_full_to_irr_ipoint(i_full_points)
            if (i_irr_points > 0) then
               batch_irr%hartree_potential(i_irr_points) = hartree_potential(i_full_points)
               batch_irr%rho(1, i_irr_points) = rho(1, i_full_points)
            end if

         end do
      end do bl6
      !end block mapping

   end subroutine compute_batch_rlsy_full_to_irr

   subroutine compute_batch_rlsy_full_to_irr_points(full_points, irr_points)
      use grids, only: batch_of_points, grid_point, batches
      use rlsy_interface, only: rlsy_h
      use dimensions, only: n_full_points, n_my_batches
      use mpi_tasks
      real*8, intent(IN) :: full_points(:)
      real*8, intent(OUT) :: irr_points(:)
      !mapping_mpi: block
      integer :: send_count(n_tasks), send_displ(n_tasks), i_count(n_tasks)
      integer :: recv_count
      integer :: i_my_batch, icpu_irr_points, i_cpu
      integer, allocatable :: send_irr_index(:), recv_irr_index(:)
      real*8, allocatable :: send_irr_values(:), recv_irr_values(:)

      !mapping: block
      !integer :: ipt,ip,ib,i,iop
      !integer :: i_my_batch
      !integer :: max_index_atom, max_index_radial, max_index_angular
      !ipt = 0
      !i_full_points = 0
      !bl6:do i_my_batch = 1, n_my_batches, 1
      !     ! loop over one batch
      !  do i_index = 1, batches(i_my_batch)%size, 1
      !     i_full_points = i_full_points + 1
      !     i_irr_points = batch_irr%map_full_to_irr_ipoint(i_full_points)
      !     iop = batch_irr%map_full_to_irr_op(i_full_points)
      !     !write(*,*) "YY i_full_points, i_irr_points", i_full_points, i_irr_points
      !     if (i_irr_points > 0 .and. iop .eq. 1) then
      !       irr_points(i_irr_points) = full_points(i_full_points)
      !     end if
      !
      !  end do
      !end do bl6
      !end block mapping
      !mapping_mpi: block

      i_full_points = 0
      send_count = 0
      send_displ = 0
      i_count = 0

      do i_my_batch = 1, n_my_batches, 1
         do i_index = 1, batches(i_my_batch)%size, 1
            i_full_points = i_full_points + 1
            icpu_irr_points = batch_irr%map_full_to_irr_icpu(i_full_points)
            i_irr_points = batch_irr%map_full_to_irr_ipoint(i_full_points)
            iop = batch_irr%map_full_to_irr_op(i_full_points)
            if ((i_irr_points > 0) .and. (iop .eq. 1)) then
               send_count(icpu_irr_points + 1) = send_count(icpu_irr_points + 1) + 1
            end if
         end do
      end do

      allocate (send_irr_index(sum(send_count)))
      allocate (send_irr_values(sum(send_count)))
      send_irr_index = 0
      send_irr_values = 0
      send_displ(1) = 0
      do i_cpu = 2, n_tasks
         send_displ(i_cpu) = send_displ(i_cpu - 1) + send_count(i_cpu - 1)
      end do
      !write(*,*) "YY send_displ myid", send_count, send_displ, myid
      i_full_points = 0
      do i_my_batch = 1, n_my_batches, 1
         do i_index = 1, batches(i_my_batch)%size, 1
            i_full_points = i_full_points + 1
            icpu_irr_points = batch_irr%map_full_to_irr_icpu(i_full_points)
            i_irr_points = batch_irr%map_full_to_irr_ipoint(i_full_points)
            iop = batch_irr%map_full_to_irr_op(i_full_points)
            if ((i_irr_points > 0) .and. (iop .eq. 1)) then
               i_count(icpu_irr_points + 1) = i_count(icpu_irr_points + 1) + 1
               send_irr_index(send_displ(icpu_irr_points + 1) + i_count(icpu_irr_points + 1)) = &
                  i_irr_points
               send_irr_values(send_displ(icpu_irr_points + 1) + i_count(icpu_irr_points + 1)) = &
                  full_points(i_full_points)
            end if
         end do
      end do
      do icpu = 1, n_tasks
         !write(*,*) "YY myid in ", myid
         !call MPI_Barrier(mpi_comm_global, mpierr)
         call MPI_scatter(send_count, 1, MPI_INTEGER, &
                          recv_count, 1, MPI_INTEGER, icpu - 1, mpi_comm_global, mpierr)
         !write(*,*) "YY myid out ", myid
         if (mpierr /= MPI_SUCCESS) call aims_stop('MPI_scatter error', '1')
         !write(*,*) "YY recv_count", recv_count
         allocate (recv_irr_index(recv_count))
         allocate (recv_irr_values(recv_count))
         recv_irr_index = 0
         recv_irr_values = 0
         !call MPI_Barrier(mpi_comm_global, mpierr)
         call MPI_scatterv(send_irr_index, send_count, send_displ, MPI_INTEGER, &
                           recv_irr_index, recv_count, MPI_INTEGER, &
                           icpu - 1, mpi_comm_global, mpierr)
         if (mpierr /= MPI_SUCCESS) call aims_stop('MPI_scatterv error', '2')
         !call MPI_Barrier(mpi_comm_global, mpierr)
         call MPI_scatterv(send_irr_values, send_count, send_displ, MPI_DOUBLE_PRECISION, &
                           recv_irr_values, recv_count, MPI_DOUBLE_PRECISION, &
                           icpu - 1, mpi_comm_global, mpierr)
         if (mpierr /= MPI_SUCCESS) call aims_stop('MPI_scatterv error', '3')
         do i_recv_point = 1, recv_count
            irr_points(recv_irr_index(i_recv_point)) = recv_irr_values(i_recv_point)
         end do
         if (allocated(recv_irr_index)) deallocate (recv_irr_index)
         if (allocated(recv_irr_values)) deallocate (recv_irr_values)
      end do
      if (allocated(send_irr_index)) deallocate (send_irr_index)
      if (allocated(send_irr_values)) deallocate (send_irr_values)

      !end block mapping_mpi

   end subroutine compute_batch_rlsy_full_to_irr_points

   subroutine compute_batch_rlsy_full_to_irr_points_gradient(full_points, irr_points)
      use grids, only: batch_of_points, grid_point, batches
      use rlsy_interface, only: rlsy_h
      use dimensions, only: n_full_points, n_my_batches
      real*8, intent(IN) :: full_points(:, :)
      real*8, intent(OUT) :: irr_points(:, :)
      !mapping: block
      integer :: ipt, ip, ib, i, iop
      integer :: i_my_batch
      integer :: max_index_atom, max_index_radial, max_index_angular

      !mapping: block
      !ipt = 0
      !i_full_points = 0
      !bl6:do i_my_batch = 1, n_my_batches, 1
      !     ! loop over one batch
      !  do i_index = 1, batches(i_my_batch)%size, 1
      !     i_full_points = i_full_points + 1
      !     i_irr_points = batch_irr%map_full_to_irr_ipoint(i_full_points)
      !     iop = batch_irr%map_full_to_irr_op(i_full_points)
      !     !write(*,*) "YY i_full_points, i_irr_points", i_full_points, i_irr_points
      !     if (i_irr_points > 0) then
      !       irr_points(:,i_irr_points) = matmul(rlsy_h%spacegroup%op(iop)%im,&
      !                                    full_points(:,i_full_points))
      !     end if
      !
      !  end do
      !end do bl6
      !end block mapping
      do i = 1, 3
         call compute_batch_rlsy_full_to_irr_points(full_points(i, :), irr_points(i, :))
      end do

   end subroutine compute_batch_rlsy_full_to_irr_points_gradient

   subroutine compute_batch_rlsy_irr_to_full_points(irr_points, full_points)
      use grids, only: batch_of_points, grid_point, batches
      use rlsy_interface, only: rlsy_h
      use dimensions, only: n_full_points, n_my_batches
      use mpi_tasks
      real*8, intent(IN) :: irr_points(:)
      real*8, intent(OUT) :: full_points(:)
      !mapping_mpi: block
      integer :: send_count(n_tasks), send_displ(n_tasks), i_count(n_tasks)
      integer :: recv_count
      integer :: i_my_batch, icpu_irr_points, i_cpu
      integer, allocatable :: send_irr_index(:), recv_irr_index(:), master_full_indexes(:)
      real*8, allocatable :: send_irr_values(:), sendback_irr_values(:)
      real*8, allocatable :: master_full_values(:)

      !mapping: block
      !integer :: ipt,ip,ib,i
      !integer :: i_my_batch
      !integer :: max_index_atom, max_index_radial, max_index_angular
      !ipt = 0
      !i_full_points = 0
      !bl6:do i_my_batch = 1, n_my_batches, 1
      !   ! loop over one batch
      !   do i_index = 1, batches(i_my_batch)%size, 1
      !      i_full_points = i_full_points + 1
      !      i_irr_points = batch_irr%map_full_to_irr_ipoint(i_full_points)
      !      !write(*,*) "YY i_full_points, i_irr_points", i_full_points, i_irr_points
      !      if (i_irr_points > 0) then
      !        full_points(i_full_points) = irr_points(i_irr_points)
      !      end if
      !
      !  end do
      !end do bl6
      !end block mapping

      !mapping_mpi: block

      i_full_points = 0
      send_count = 0
      send_displ = 0
      i_count = 0
      !send_irr_values = 0
      !recv_irr_values = 0

      do i_my_batch = 1, n_my_batches, 1
         do i_index = 1, batches(i_my_batch)%size, 1
            i_full_points = i_full_points + 1
            icpu_irr_points = batch_irr%map_full_to_irr_icpu(i_full_points)
            i_irr_points = batch_irr%map_full_to_irr_ipoint(i_full_points)
            if ((i_irr_points > 0)) then
               send_count(icpu_irr_points + 1) = send_count(icpu_irr_points + 1) + 1
            end if
         end do
      end do

      allocate (send_irr_index(sum(send_count)))
      allocate (master_full_indexes(sum(send_count)))
      allocate (master_full_values(sum(send_count)))
      send_irr_index = 0
      send_displ(1) = 0
      do i_cpu = 2, n_tasks
         send_displ(i_cpu) = send_displ(i_cpu - 1) + send_count(i_cpu - 1)
      end do
      !write(*,*) "YY send_displ myid", send_count, send_displ, myid
      i_full_points = 0
      do i_my_batch = 1, n_my_batches, 1
         do i_index = 1, batches(i_my_batch)%size, 1
            i_full_points = i_full_points + 1
            icpu_irr_points = batch_irr%map_full_to_irr_icpu(i_full_points)
            i_irr_points = batch_irr%map_full_to_irr_ipoint(i_full_points)
            if ((i_irr_points > 0)) then
               i_count(icpu_irr_points + 1) = i_count(icpu_irr_points + 1) + 1
               send_irr_index(send_displ(icpu_irr_points + 1) + i_count(icpu_irr_points + 1)) = &
                  i_irr_points
               master_full_indexes(send_displ(icpu_irr_points + 1) + i_count(icpu_irr_points + 1)) = &
                  i_full_points
            end if
         end do
      end do
      do icpu = 1, n_tasks
         !call MPI_Barrier(mpi_comm_global, mpierr)
         call MPI_scatter(send_count, 1, MPI_INTEGER, &
                          recv_count, 1, MPI_INTEGER, icpu - 1, mpi_comm_global, mpierr)
         if (mpierr /= MPI_SUCCESS) call aims_stop('MPI_scatter error', '4')

         !write(*,*) "YY recv_count", recv_count
         allocate (recv_irr_index(recv_count))
         allocate (sendback_irr_values(recv_count))
         recv_irr_index = 0
         !call MPI_Barrier(mpi_comm_global, mpierr)
         call MPI_scatterv(send_irr_index, send_count, send_displ, MPI_INTEGER, &
                           recv_irr_index, recv_count, MPI_INTEGER, &
                           icpu - 1, mpi_comm_global, mpierr)
         if (mpierr /= MPI_SUCCESS) call aims_stop('MPI_scatterv error', '5')
         sendback_irr_values = 0
         do j_count = 1, recv_count
            sendback_irr_values(j_count) = irr_points(recv_irr_index(j_count))
         enddo
         !call MPI_Barrier(mpi_comm_global, mpierr)
         call MPI_gatherv(sendback_irr_values, recv_count, MPI_DOUBLE_PRECISION, &
                          master_full_values, send_count, send_displ, MPI_DOUBLE_PRECISION, &
                          icpu - 1, mpi_comm_global, mpierr)
         if (mpierr /= MPI_SUCCESS) call aims_stop('MPI_gatherv error', '6')
         do j_count = 1, sum(send_count)
            full_points(master_full_indexes(j_count)) = master_full_values(j_count)
         enddo
         if (allocated(recv_irr_index)) deallocate (recv_irr_index)
         if (allocated(sendback_irr_values)) deallocate (sendback_irr_values)
      end do
      if (allocated(send_irr_index)) deallocate (send_irr_index)
      if (allocated(master_full_indexes)) deallocate (master_full_indexes)
      if (allocated(master_full_values)) deallocate (master_full_values)

      !end block mapping_mpi

   end subroutine compute_batch_rlsy_irr_to_full_points

   subroutine compute_batch_rlsy_irr_to_full_points_gradient(irr_points, full_points)
      use grids, only: batch_of_points, grid_point, batches
      use rlsy_interface, only: rlsy_h
      use dimensions, only: n_full_points, n_my_batches
      real*8, intent(IN) :: irr_points(:, :)
      real*8, intent(OUT) :: full_points(:, :)
      !mapping: block
      real*8 :: full_points_temp(3)
      integer :: ipt, ip, ib, i, iop
      integer :: i_my_batch
      integer :: max_index_atom, max_index_radial, max_index_angular

      !mapping: block
      !ipt = 0
      !i_full_points = 0
      !bl6:do i_my_batch = 1, n_my_batches, 1
      !   ! loop over one batch
      !   do i_index = 1, batches(i_my_batch)%size, 1
      !      i_full_points = i_full_points + 1
      !      i_irr_points = batch_irr%map_full_to_irr_ipoint(i_full_points)
      !      iop = batch_irr%map_full_to_irr_op(i_full_points)
      !      !write(*,*) "YY i_full_points, i_irr_points", i_full_points, i_irr_points
      !      if (i_irr_points > 0) then
      !        full_points(:,i_full_points) = matmul(rlsy_h%spacegroup%op(iop)%m, &
      !                                       irr_points(:,i_irr_points))
      !      end if
      !
      !  end do
      !end do bl6
      !end block mapping
      do i = 1, 3
         call compute_batch_rlsy_irr_to_full_points(irr_points(i, :), full_points(i, :))
      end do
      do i_full_points = 1, n_full_points
         iop = batch_irr%map_full_to_irr_op(i_full_points)
         full_points_temp = full_points(:, i_full_points)
         full_points(:, i_full_points) = matmul(rlsy_h%spacegroup%op(iop)%m, &
                                                full_points_temp(:))
      end do

   end subroutine compute_batch_rlsy_irr_to_full_points_gradient

end module batch_rlsy_interface
