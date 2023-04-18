!------------------------------------------------------------------------------------------------------------

!****s* FHI-aims/integrate_polarizability
!  NAME
!   integrate_polarizability
!  SYNOPSIS

subroutine integrate_polarizability &
     (partition_tab_std, first_order_rho_std,polarizability)

!  PURPOSE
!    The subroutine get polarizability(1:3)=  int[first_order_rho(1:3,:) r(:) dr ]
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use pbc_lists
  use load_balancing

  implicit none

!  ARGUMENTS


  !-------------------shanghui add here------------------------------------
  real*8, target, dimension(n_full_points)        :: partition_tab_std
  real*8, target, dimension(3,n_spin,n_full_points) :: first_order_rho_std
  real*8, dimension(3,3), intent(inout) :: polarizability
  !------------------shanghui end add------------------------------------- 

!  INPUTS
! o partition_tab_std -- values of partition function
!
!  OUTPUT
!  AUTHOR
!    Honghui Shang @ FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
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


  ! Local variables

!------------------real calc-----------
  integer i_center
  integer i_batch
  integer i_coord, i_coord2
  integer i_index
  integer i_spin
  integer :: i_full_points


  real*8 coord_current(3)

  real*8, pointer :: partition_tab(:)
  real*8, pointer :: first_order_rho(:,:,:)
!-----------------end real calc-----------


!-------------mpi---------------------------
  integer :: info
  character*200 :: info_str

  ! Load balancing stuff
  integer n_my_batches_work  ! Number of batches actually used
  integer n_full_points_work ! Number of integration points actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  integer n_bp

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all
!-------------end mpi---------------------------


  if(use_batch_permutation > 0) then
    call localorb_info("Summing up the Hartree potential with load balancing.", use_unit,'(2X,A)', OL_norm )
  else
    call localorb_info("Summing up the Hartree potential.", use_unit,'(2X,A)', OL_norm )
  endif

  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)

  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    n_full_points_work = batch_perm(n_bp)%n_full_points

    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab

  else
    n_my_batches_work = n_my_batches
    n_full_points_work = n_full_points
    batches_work => batches
    partition_tab => partition_tab_std
    first_order_rho => first_order_rho_std

  endif

  if(get_batch_weights) then
    allocate(batch_times(n_my_batches_work))
    batch_times(:) = 0
  endif

  !-----------------------------------------------------------------------------




  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()

    polarizability = 0.0d0
    i_full_points = 0

      do i_batch = 1, n_my_batches_work

        if(get_batch_weights) time_start = mpi_wtime()

        ! loop over one batch
        do i_index = 1, batches_work(i_batch)%size, 1

          ! i_full_points is the index that indicates where we are in the entire grid (for external quanities like rho, potential, ...)
          i_full_points = i_full_points + 1

          ! Only execute if partition_tab is .gt. zero, else
          ! we can run into 1/0 !!!
          if (partition_tab(i_full_points).gt.0.d0) then

            ! get current integration point coordinate
            coord_current(:) = batches_work(i_batch) % points(i_index) % coords(:)


           do i_coord = 1, 3
            do i_coord2 = 1, 3
            
              do i_spin = 1, n_spin

              ! Nath: The polarizability calculated here should be symmetric, so we could calculate only 6
              ! components instead of 9, but this is a good (and short) check anyway
              polarizability(i_coord,i_coord2)=polarizability(i_coord,i_coord2) & 
              +partition_tab(i_full_points)  & 
              *first_order_rho(i_coord2,i_spin,i_full_points)*coord_current(i_coord)
            
              enddo
             
            enddo
           enddo



          end if ! end if (partition_tab.gt.0.d0)
        end do  ! end loop over points in a batch
        if(get_batch_weights) batch_times(i_batch) = batch_times(i_batch) + mpi_wtime() - time_start
      end do ! end loop over batches



         

  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_world,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for integrate_polarizability: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)


  if(get_batch_weights) then
    call set_batch_weights(n_bp,batch_times)
    deallocate(batch_times)
  endif


  ! Write one last line to bound output
  write(info_str,*) ' '
  call localorb_info(info_str, use_unit,'(A)',OL_norm)


  !     synchronise the polarizability
  !-------shanghui begin parallel------
  if(.not. use_local_index) call sync_integrate_polarizability( polarizability )
  !-------shanghui end parallel------


end subroutine integrate_polarizability
!******
!------------------------------------------------------------------------------------------------------------
