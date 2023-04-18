!****s* FHI-aims/trans_circle_to_cubic
!  NAME
!   trans_circle_to_cubic
!  SYNOPSIS

subroutine trans_circle_to_cubic &
     (partition_tab_std, first_order_potential_circle, first_order_potential_cubic) 

!  PURPOSE
!  
!  Integrates the matrix elements for the first_order overlap(S) matrix,
!  using a fixed basis set. 

!  called by a SCF subroutine

!  shanghui 2012.05.08 : created 
!  shanghui 2013.12.20 : change to phonon_gamma version (p0)
!  shanghui 2014.11.11 : change to phonon only with sparse matrix.
!  
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use basis
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use constants
  use species_data, only: species_name
  use load_balancing
  use pbc_lists
  implicit none

!  ARGUMENTS

  real*8, target, dimension(n_full_points)            :: partition_tab_std
  !shanghui------------------------------------------------------------------------
  real*8, dimension(3, n_centers_in_sc_DFPT,n_full_points),intent(in) :: & 
          first_order_potential_circle       ! could outside of center cubic cell 

  real*8, dimension(3, n_centers_in_sc_DFPT,n_full_points),intent(out) :: & 
          first_order_potential_cubic        ! just inside center cubic cell 

  !shanghui------------------------------------------------------------------------

!  INPUTS
!  o partition_tab_std -- values of partition functions
!  o first_order_potential_circle
!
!  OUTPUT
!  o first_order_potential_cubic 
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

  real*8 coord_current(3), length(3) 
  integer :: cell_index_current(3)
  integer :: i_cell_current, i_cell_in_sc_DFPT, i_cell, i_cell_delta
  integer :: i_center, i_center_trans

  !     optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points

  !     and condensed version of hamiltonian_partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)



  !  counters

  integer :: i_my_batch,i_index
  integer i_point
  integer :: i_full_points

  character*200 :: info_str
  integer ::  info


  ! Load balancing stuff

  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

!  integer ld_hamiltonian  ! leading dimension of hamiltonian in calling routine

  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)

  integer  n_bp
  integer, allocatable :: ins_idx(:)

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all

  ! begin work

  if(use_batch_permutation > 0) then
    write(info_str,'(2X,A)') "trans_circle_to_cubic: batch-based integration with load balancing"
  else
    write(info_str,'(2X,A)') "trans_circle_to_cubic: batch-based integration."
  endif
  call localorb_info(info_str, use_unit,'(A)',OL_norm)


  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)

  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab

    allocate(ins_idx(batch_perm(n_bp)%n_basis_local))

    !ld_hamiltonian = batch_perm(n_bp)%n_local_matrix_size

  else

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std

  endif


  !-----------------------------------------------------------------------------

  ! initialize
  first_order_potential_cubic= 0.0d0

  i_full_points = 0

  ! perform partitioned integration, batch by batch of integration point.
  ! This will be the outermost loop, to save evaluations of the potential.
  ! and the Y_lm functions

  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()

  do i_my_batch = 1, n_my_batches_work, 1

     i_point = 0

     ! loop over one batch
     do i_index = 1, batches_work(i_my_batch)%size, 1

        i_full_points = i_full_points + 1

        if (partition_tab(i_full_points).gt.0.d0) then

           i_point = i_point+1

           ! get current integration point coordinate
           coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)

           length(1:3) = matmul(map_to_center_cell_matrix, coord_current)
           cell_index_current(1:3) =  nint(length(1:3))                                          ! this is R  

 
           if(cell_index_current(1).ne.0.or.cell_index_current(2).ne.0 & 
              .or.cell_index_current(3).ne.0 ) then 
          
             do i_cell_in_sc_DFPT   = 1,n_cells_in_sc_DFPT
             if( cell_index_sc_DFPT(i_cell_in_sc_DFPT,1).eq. &
                 cell_index_current(1) .and. &
                 cell_index_sc_DFPT(i_cell_in_sc_DFPT,2).eq. &
                 cell_index_current(2) .and. &
                 cell_index_sc_DFPT(i_cell_in_sc_DFPT,3).eq. &
                 cell_index_current(3) ) then

                 i_cell_current = i_cell_in_sc_DFPT

            endif
            enddo

           else 
                i_cell_current = 1 
           endif   


           do i_center = 1 , n_centers_in_sc_DFPT 

              i_cell = center_in_sc_DFPT_to_cell_in_sc_DFPT(i_center)
 
              i_cell_delta = cell_add_sc_DFPT(i_cell,i_cell_current)

              i_center_trans = cell_and_atom_to_center_sc_DFPT(i_cell_delta,  &      
                                 center_in_sc_DFPT_to_atom(i_center) )              


              first_order_potential_cubic(1:3, i_center,      i_full_points) = &  ! i_full_points = r
              first_order_potential_circle(1:3,i_center_trans,i_full_points)    ! i_full_points = r + R  
                                                 ! i_center_trans = i_center + R         
           enddo 


        end if
     enddo


  end do ! end loop over batches

  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_world,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for trans_circle_to_cubic: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  if(time_all>time_work*1.3 .and. .not.use_load_balancing) &
    info_str = trim(info_str) // ' => Consider using load balancing!'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)



  if(use_batch_permutation > 0) then
    deallocate(ins_idx)
  endif


end subroutine trans_circle_to_cubic
!******
