!------------------------------------------------------------------------------------------------------------

!****s* FHI-aims/integrate_hellam_dynamical_matrix_phonon_reduce_memory
!  NAME integrate_hellam_dynamical_matrix_phonon_reduce_memory
!  SYNOPSIS

subroutine integrate_hellam_dynamical_matrix_phonon_reduce_memory &
     (                                & 
       partition_tab_std,  & 
       rho_std, first_order_rho_std,  &
       i_q_point, j_atom, j_coord, &
       hellman_feynman_dynamical_matrix_delta_part &
     )

!
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use free_atoms
  use spline
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use pbc_lists
  use load_balancing

  implicit none

!  ARGUMENTS


  real*8, target, dimension(n_full_points)        :: partition_tab_std
  real*8, target, dimension(n_spin, n_full_points) :: rho_std
  complex*16, target, dimension(n_full_points) ::  first_order_rho_std

  integer, intent(in) ::  i_q_point
  integer, intent(in) ::  j_atom
  integer, intent(in) ::  j_coord

  complex*16, dimension(3, n_atoms,3,n_atoms),intent(out)   :: hellman_feynman_dynamical_matrix_delta_part
!  INPUTS
!
!  OUTPUT
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


  !  local variables
  real*8, dimension(:,:),  allocatable :: local_rho
  complex*16, dimension(:),  allocatable :: local_first_order_rho
  real*8 :: partition(n_max_batch_size)
  real*8, dimension(:,:), allocatable   :: coords_npoints


  !  counters
  integer :: i_spin
  integer :: i_batch, i_index
  integer :: i_full_points, i_point, n_points

  character*100 :: info_str

  integer :: mpierr
  integer :: info

 ! Load balancing stuff
  integer n_my_batches_work  ! Number of batches actually used
  integer n_full_points_work ! Number of integration points actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used




  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)
  real*8, pointer :: rho(:,:)
  complex*16, pointer :: first_order_rho(:)


  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start
  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all

  call localorb_info("Summing up the delta_hellman_dynamical_matrix", use_unit,'(2X,A)', OL_norm )

  !-----------------------------------------------------------------------------

    n_my_batches_work = n_my_batches
    n_full_points_work = n_full_points
    batches_work => batches
    partition_tab => partition_tab_std
    rho => rho_std
    first_order_rho => first_order_rho_std

  !-----------------------------------------------------------------------------


    allocate( local_rho(n_spin,n_max_batch_size),stat=info)
    call check_allocation(info, 'local_rho                     ')

    allocate( local_first_order_rho(n_max_batch_size),stat=info)
    call check_allocation(info, 'local_first_order_rho                     ')

    allocate (coords_npoints(3,n_max_batch_size),STAT=info)
    call check_allocation(info, 'coords_npoints        ')



   hellman_feynman_dynamical_matrix_delta_part = (0.0d0,0.0d0)


  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()


      ! Now follows the real work: Summing the multipole potential, density, and their
      ! derivatives on the integration grid
      !

      ! Reset grid counter for current Hartree potential center
      i_full_points = 0

      do i_batch = 1, n_my_batches_work

        if(get_batch_weights) time_start = mpi_wtime()

        i_point = 0
        ! loop over one batch
        do i_index = 1, batches_work(i_batch)%size, 1

          ! i_full_points is the index that indicates where we are in the entire grid (for external quanities like rho, potential, ...)
          i_full_points = i_full_points + 1

          ! Only execute if partition_tab is .gt. zero, else
          ! we can run into 1/0 !!!
          if (partition_tab(i_full_points).gt.0.d0) then

            i_point = i_point+1
 
            ! get current integration point coordinate
            coords_npoints(1:3,i_point) = batches_work(i_batch) % points(i_index) % coords(1:3) 

            do i_spin = 1, n_spin, 1
               local_rho(i_spin,i_point) = rho(i_spin,i_full_points)
            enddo


            local_first_order_rho(i_point) = first_order_rho(i_full_points)

            partition(i_point) = partition_tab(i_full_points)

          end if ! end if (partition_tab.gt.0.d0)
        end do  ! end loop over points in a batch

         n_points = i_point

          do i_spin = 1, n_spin, 1

           call  evaluate_hellman_feynman_delta_part_phonon_reduce_memory(   &
                 n_points, partition,      &
                 local_rho(i_spin,1:n_points),local_first_order_rho(1:n_points), &
                 coords_npoints(:,1:n_points), &
                 i_q_point, j_atom,j_coord, &
                 hellman_feynman_dynamical_matrix_delta_part )

          enddo


        if(get_batch_weights) batch_times(i_batch) = batch_times(i_batch) + mpi_wtime() - time_start
      end do ! end loop over batches



  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_world,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for potential: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)





   !-------shanghui begin parallel------  
    if(.not. use_local_index) call sync_vector_complex(hellman_feynman_dynamical_matrix_delta_part,3*n_atoms*3*n_atoms)
   !-------shanghui end parallel------  



  ! Write one last line to bound output
  write(info_str,*) ' '
  call localorb_info(info_str, use_unit,'(A)',OL_norm)



  if(allocated( local_rho            )) deallocate( local_rho            )
  if(allocated( local_first_order_rho      )) deallocate( local_first_order_rho      )
  if(allocated( coords_npoints)) deallocate( coords_npoints)



end subroutine integrate_hellam_dynamical_matrix_phonon_reduce_memory
!******
!------------------------------------------------------------------------------------------------------------
