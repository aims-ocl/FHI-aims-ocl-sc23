!****h* FHI-aims/batch_statistics
!*  NAME
!*    batch_statistics
!*  SYNOPSIS
module batch_statistics
!*  PURPOSE
!*    This module stored the statistics calculated for batch-integration-related
!*    quantities.
!*  USES
  use statistics, only: stats_int
  implicit none
!*  AUTHOR
!*    William Huhn (Duke University)
!*  HISTORY
!*    June 2018 - Forked off of json_output.
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is
!*    subject to the terms and conditions of the respective license
!*    agreement.
!*  SOURCE

  private

  public :: synchronize_batch_statistics

  ! batch_size, across all batches
  ! batch_sizes is calculated in partition_grid then discarded, so we intercept
  ! them there
  type(stats_int), public :: batch_sizes_stats

  ! n_max_compute_ham, across MPI tasks
  type(stats_int), public :: n_max_compute_ham_stats

  ! n_hamiltonian_matrix_size/n_local_matrix_size, across MPI tasks
  type(stats_int), public :: ld_ham_integ_stats
  type(stats_int), public :: ld_ham_hpot_stats
  type(stats_int), public :: ld_ham_density_stats

  ! Estimate of number of basis functions touching a batch
  type(stats_int), public :: batch_n_compute_stats

contains

  subroutine synchronize_batch_statistics( )
    use dimensions, only: n_grid_batches, n_max_compute_ham, n_my_batches
    use grids, only: batches
    use load_balancing, only: batch_perm, batch_permutation, n_bp_integ, &
        n_bp_hpot, n_bp_density
    use mpi_tasks, only: myid, n_tasks
    use mpi_utilities, only: batch_task_list
    use statistics, only: calc_stats_array
    use synchronize_mpi_basic, only: sync_int_vector
    implicit none
    integer :: i_grid_batch, i_my_batch, task_for_batch, offset = 0
    integer, allocatable :: n_batches_on_task(:), array(:)

    allocate( array(n_tasks) )

    ! Calculate distribution for various quantities spread across MPI tasks
    array = 0
    array(myid+1) = n_max_compute_ham
    call sync_int_vector(array, n_tasks)
    call calc_stats_array(array, n_max_compute_ham_stats, n_fractiles = 10, &
                          n_bins = 10)

    ! Local matrix size for integration batch permutation
    if (batch_perm(n_bp_integ)%initialized) then
      array = 0
      array(myid+1) = batch_perm(n_bp_integ)%n_local_matrix_size
      call sync_int_vector(array, n_tasks)
      call calc_stats_array(array, ld_ham_integ_stats, n_fractiles = 10, &
                            n_bins = 10)
    end if

    ! Local matrix size for Hartree potential batch permutation
    if (batch_perm(n_bp_hpot)%initialized) then
      array = 0
      array(myid+1) = batch_perm(n_bp_hpot)%n_local_matrix_size
      call sync_int_vector(array, n_tasks)
      call calc_stats_array(array, ld_ham_hpot_stats, n_fractiles = 10, &
                            n_bins = 10)
    end if

    ! Local matrix size for density batch permutation
    if (batch_perm(n_bp_density)%initialized) then
      array = 0
      array(myid+1) = batch_perm(n_bp_density)%n_local_matrix_size
      call sync_int_vector(array, n_tasks)
      call calc_stats_array(array, ld_ham_density_stats, n_fractiles = 10, &
                            n_bins = 10)
    end if

    deallocate( array )

    ! The batch size statistics are computed in partition_grid(), so we don't
    ! need to recompute them here.

    ! Here's where we compute the number of basis functions touching a batch.
    ! This is a bit more difficult than previous statistics, because as far as
    ! I know, there's no good place in the code where we calculate this quantity
    ! and still have the global index for a batch, like we did for batch sizes.
    ! So we have to go through a somewhat awkward construction of extracting the
    ! number of batches per MPI task from batch_task_list, using that knowledge
    ! to determine how many batches all MPI tasks with lower ranks than the
    ! current one have, which then allows us to create a list of batches across
    ! all MPI ranks, ordered by MPI rank.
    allocate( array(n_grid_batches) )

    ! Create a list of how many batches each MPI task has, essentially the
    ! inverse array of batch_task_list
    allocate( n_batches_on_task(0:n_tasks-1) )
    n_batches_on_task = 0
    do i_grid_batch = 1, n_grid_batches
      task_for_batch = batch_task_list(i_grid_batch)
      n_batches_on_task(task_for_batch) = n_batches_on_task(task_for_batch) + 1
    end do
    ! Now use the number of batches all previous MPI tasks have to determine
    ! the current task's offset in the global list of batches
    offset = 0
    if (myid > 0) offset = sum( n_batches_on_task(0:(myid-1)) )
    deallocate( n_batches_on_task )

    ! Now that we know where the batches on this MPI task reside in the global
    ! array of batches, generate statistics on batch_n_compute
    array = 0
    do i_my_batch = 1, n_my_batches
      array(offset + i_my_batch) = batches(i_my_batch)%batch_n_compute
    end do
    call sync_int_vector(array, n_grid_batches)
    call calc_stats_array(array, batch_n_compute_stats, n_fractiles = 10, &
                          n_bins = 10)

    deallocate( array )
  end subroutine synchronize_batch_statistics

  !******
end module
!******
