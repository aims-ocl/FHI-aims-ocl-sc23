module optics_module

!  PURPOSE
!
!  This file may at some point contain all the needed variable declarations and
!  subroutines for the TDDFT calculation of charge-neutral optical excitations of molecules.
!
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!    This file was written by Jan Kloppenburg
!  SOURCE

  use dimensions
  use runtime_choices
  use physics
  use species_data
  use constants
  use gw_para
  use grids
  use geometry
  use mpi_utilities
  use synchronize_mpi
  use sbt_overlap_aims
  use localorb_io
  use basis
  use prodbas
  use grids
  use hartree_fock
  use timing
  use mpi_tasks
  use mpi_utilities
  use elsi_wrapper, only: aims_elsi_stdevp
  use timing

  implicit none

  real*8, parameter :: sparse_thresh = 1d-12

  integer :: occ_states, unocc_states, occ_limit, unocc_limit, i_proc
  integer :: i_state, a_state, j_state, b_state, low_limit, n_pair_states
  integer :: my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols, num_blocks
  integer :: comm_tag, my_casida_ctxt, num_pcol, num_prow, na_cols, na_rows
  integer :: casida_desc(9), mpi_status(mpi_status_size), i_rot, j_rot
  integer :: iia, jja, iarow, iacol, myrow, mycol, nprow, npcol

  real*8, allocatable, dimension(:)       :: coord_of_center
  real*8, allocatable, dimension(:,:)     :: atom_coord_wrt_center
  real*8, allocatable, dimension(:,:,:,:) :: ovlp_3ks, fxc_3ks, pure_3ks
  real*8, allocatable, dimension(:,:,:,:) :: dipole_moment

  character*128 :: info_str

  integer, external :: numroc
  real*8,  external :: ddot

  intrinsic :: int, nint, dble, floor, min, max, dsqrt, ceiling, mod

  contains

subroutine reset_optics_module()
  implicit none

  if (allocated(coord_of_center)) deallocate(coord_of_center)
  if (allocated(atom_coord_wrt_center)) deallocate(atom_coord_wrt_center)
  if (allocated(ovlp_3ks)) deallocate(ovlp_3ks)
  if (allocated(fxc_3ks)) deallocate(fxc_3ks)
  if (allocated(pure_3ks)) deallocate(pure_3ks)
  if (allocated(dipole_moment)) deallocate(dipole_moment)

end subroutine reset_optics_module

subroutine casida_tddft_2d()
  implicit none

  integer :: i_state_loc, j_state_loc, a_state_loc, b_state_loc, &
             k_cnt, l_cnt, m_cnt, errnum, &
             ia, ja, tmp_int, i_task, sendrecv_dims, &
             len_of_recs, max_na_rows, max_na_cols, i_spin

  integer, allocatable, dimension(:)   :: mpi_bcastproc, tmp_pr, tmp_ps
  integer, allocatable, dimension(:,:) :: tmp1r, tmp2r, tmp1s, tmp2s
  integer, allocatable, dimension(:,:) :: tmp_bcastidx, mpi_bcastidx

  real*8 :: kernel, tmp_real, total_mem, proc_mem, time_comm, clock_time_comm, time_calc, clock_time_calc

  real*8, allocatable, dimension(:) :: par_evals, osc_str, tmp_vecr, tmp_vecs, mpi_vec
  real*8, allocatable, dimension(:,:) :: par_evecs, par_omega, tmom
  real*8, allocatable, dimension(:,:,:,:) :: tmp_o3ks, tmp_f3ks, loc_o3ks
  real*8, allocatable, dimension(:,:,:,:) :: loc_f3ks, metric

  logical :: first_found
  integer :: i_local_start, i_local_end
  integer :: a_local_start, a_local_end
  integer :: ks3_blacs_context
  integer :: ks3_sc_desc(9)
  integer :: blacs_info
  integer :: ks3_nrows
  integer :: ks3_ncols
  integer :: ks3_prows
  integer :: ks3_pcols
  integer :: ks3_my_p_row
  integer :: ks3_my_p_col
  integer :: metric_blacs_context
  integer :: metric_sc_desc(9)
  integer :: metric_nrows
  integer :: metric_ncols
  integer :: metric_prows
  integer :: metric_pcols
  integer :: metric_my_p_row
  integer :: metric_my_p_col

  call get_timestamps(time_neutral_excitation_casida, &
           clock_time_neutral_excitation_casida)

  total_mem = 8d0 * dble(n_pair_states)**2

  ! Some initialization output 
  if(myid==0) then
    if(excited_mode_triplet) then
      write(use_unit,'(2x,a)') '| '
      write(use_unit,'(2x,a)') '| Starting TDDFT excited state calculation for Triplets.'
      write(use_unit,'(2x,a)') '| '
    else
      write(use_unit,'(2x,a)') '| '
      write(use_unit,'(2x,a)') '| Starting TDDFT excited state calculation for Singlets.'
      write(use_unit,'(2x,a)') '| '
    endif
    write(use_unit,'(2x,a,i8.1,a,i8.1,a,i6.1,a,i6.1,a,i5.1,a)') '| Deviding the ', &
      n_pair_states, ' x ', n_pair_states, ' Matrix into ', &
      num_pcol,' colomns and ', num_prow, ' rows on ', n_tasks, ' processors.'
    write(use_unit,'(2x,a,f10.3,a)') '| Memory load of ', total_mem/2d0**30, &
      ' GB devided among all Processors.'
  endif

  proc_mem = 8d0 * dble(na_rows) * dble(na_cols)

  if( (total_mem/2d0**20)/dble(n_tasks).gt.1024d0) then
    write(info_str,'(2x,a,i5.1,a,i8.1,a,i8.1,a,f8.3,a)') '| Processor ', myid, &
      ' gets ', na_rows, ' rows and ', na_cols,' columns  (',proc_mem/2d0**20/1024d0,' GB).'
    call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)
  else
    write(info_str,'(2x,a,i5.1,a,i8.1,a,i8.1,a,f8.3,a)') '| Processor ', myid, &
      ' gets ', na_rows, ' rows and ', na_cols,' columns  (',proc_mem/2d0**20,' MB).'
    call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)
  endif

  if(myid==0) then
     write(use_unit,'(2x,a)') '|'
     write(use_unit,'(2x,a,I5.1)') '| Number of total rotations for matrix building: ', n_tasks
     write(use_unit,'(2x,a)') '|'
  endif 


  ! Ok, start work ...

  ! Get local number of pairstates
  k_cnt = 0
  !  Loop over all occupied states in energy window
  do i_state=low_limit, occ_states
    ! Is this state treated my this process?
    if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
    ! Yes, then loop over all unoccupied states in energy window
    do a_state=occ_states+1, n_states-unocc_limit
      ! Is this state treated on this process?
      if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle
       ! If yes increase the local number of pairstates
       k_cnt = k_cnt + 1
    enddo 
  enddo

  ! Find the maximal number of pairstates by sending all local number of pairstates
  ! to process 0 and compare 
  if(myid==0) then
    do i_proc=1, n_tasks-1
      call mpi_recv(tmp_int, 1, mpi_integer, i_proc, &
        100, mpi_comm_world, mpi_status, errnum)
      if(tmp_int.gt.k_cnt) k_cnt = tmp_int
    enddo
  else
    call mpi_send(k_cnt, 1, mpi_integer, 0, 100, mpi_comm_world, errnum)
  end if
  
  i_state = k_cnt

  ! Do the same again for j/b
  k_cnt = 0
  ! Loop over all occupied states in energy window
  do j_state=low_limit, occ_states
    ! If this state is treated on this process loop over
    ! all unoccupied states in energy window
    if(own_dim1_o3KS(j_state) /= myp1_o3KS) cycle
    do b_state=occ_states+1, n_states-unocc_limit
       ! If this state is treated on this process increase the number 
       ! of pair states by one 
       if(own_dim2_o3KS(b_state) /= myp2_o3KS) cycle
       k_cnt = k_cnt + 1
     enddo 
  enddo

  ! Find the maximal number of pairstates by sending all local number of pairstates
  ! to process 0 and compare 
  if(myid==0) then
    do i_proc=1, n_tasks-1
      call mpi_recv(tmp_int, 1, mpi_integer, i_proc, &
        100, mpi_comm_world, mpi_status, errnum)
      if(tmp_int.gt.k_cnt) k_cnt = tmp_int
    enddo
  else
    call mpi_send(k_cnt, 1, mpi_integer, 0, 100, mpi_comm_world, errnum)
  endif
  
  ! BL: In my opinion this will allways be k_cnt**2
  ! JK: no, this is effectively occ_states * unocc_states
  !     after the maximum of both is found from all cores
  sendrecv_dims = i_state * k_cnt

  call mpi_bcast(sendrecv_dims, 1, mpi_integer, 0, mpi_comm_world, errnum)

  allocate(mpi_bcastidx(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'mpi broadcast index')
  mpi_bcastidx = 0
  allocate(mpi_vec(sendrecv_dims),stat=errnum)
  call check_allocation(errnum, 'mpi broadcast vector')
  mpi_vec = 0d0

  time_omegaresdist = 0d0
  clock_time_omegaresdist = 0d0
  call get_timestamps(time_omegabuild, clock_time_omegabuild)

  ! Now we construct the Casida matrix Omega = (ia,jb) + <ia|f_xc|jb>
  if(myid==0) then
    write(use_unit,'(2x,a)') '|'
    write(use_unit,'(2x,a)') '| Building Omega matrix  < ia | V + f_xc | jb >  now.'
    write(use_unit,'(2x,a)') '|'
    write(use_unit,'(2x,a,en10.1,a)') '| Only final results greater than ', sparse_thresh, ' will be stored.'
    write(use_unit,'(2x,a)') '|'
  endif

  if(.not.excited_mode_triplet) then

    allocate(tmp1r(sendrecv_dims,2),stat=errnum)
    call check_allocation(errnum, ' tmp1r res dist')
    allocate(tmp_vecr(sendrecv_dims),stat=errnum)
    call check_allocation(errnum, ' tmp_vecr res dist')

    allocate(par_omega(na_rows,na_cols),stat=errnum)
    call check_allocation(errnum, 'parallel omega matrix')
    par_omega = 0d0

    allocate(tmp2r(sendrecv_dims,2),stat=errnum)
    call check_allocation(errnum, 'rotation index 2r')
    allocate(tmp1s(sendrecv_dims,2),stat=errnum)
    call check_allocation(errnum, 'rotation index 1s')
    allocate(tmp2s(sendrecv_dims,2),stat=errnum)
    call check_allocation(errnum, 'rotation index 2s')
    allocate(loc_o3ks(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin),stat=errnum)
    call check_allocation(errnum, 'rotation local o3ks')
    allocate(loc_f3ks(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin),stat=errnum)
    call check_allocation(errnum, 'rotation local f3ks')

    tmp1s(:,1) = 0
    tmp2s(:,1) = 0

    k_cnt=1
    do i_state=low_limit, occ_states
      if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
      do a_state=occ_states+1, n_states-unocc_limit
        if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle

        i_state_loc = loc_dim1_o3KS(i_state)
        a_state_loc = loc_dim2_o3KS(a_state)

        tmp1s(k_cnt,1) = i_state
        tmp2s(k_cnt,1) = a_state
        tmp1s(k_cnt,2) = i_state_loc
        tmp2s(k_cnt,2) = a_state_loc

        k_cnt = k_cnt + 1

      enddo 
    enddo

    ! BL: DEBUG OUTPUT for getting the process distribution
    !if(myid==0) then
    !  write(use_unit,'(2x,a)') '|'
    !  write(use_unit,'(2x,a)') '| pure_3ks processor map'
    !  write(use_unit,'(2x,a,I6,a,I6,a,I6,a)') '| ', n_basbas, " x ( ",&
    !    (occ_states - low_limit + 1), " x ",&
    !    (n_states - unocc_limit - (occ_states+1) + 1), " )" 
    !  write(use_unit,'(2x,a,I6,a,I6,a)') '| Process grid: ', num_prow, " x ", num_pcol, " )" 
    !  write(use_unit,'(2x,a)') '|'
    !endif
   !
   !call MPI_BARRIER(mpi_comm_world,errnum)
   !
   !do i_state = low_limit, occ_states
   !  if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
   !  do a_state = occ_states+1, n_states-unocc_limit
   !    if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle
   !
   !    i_state_loc = (i_state - low_limit) &
   !                + (occ_states - low_limit + 1) * (a_state - (occ_states+1))
   !    a_state_loc = mod((a_state - 1),num_prow) &
   !                + num_prow * mod((i_state - 1),num_pcol)
   !
   !    write(use_unit,'(2x,I6,a,I6,a,I6,A,I6,a,I6)') i_state, " ", a_state,&
   !      " ", i_state_loc, " ", myid, " ", a_state_loc 
   !  end do
   !end do
   !
   !call MPI_BARRIER(mpi_comm_world,errnum)
   
  if (casida_setup_scalapack) then
     ! JK: If I understand this correctly, all included in thei if{} block 
     !     is a testing comparision to build the casida matrix 
     !     with scalapack pdgemms?! 
     !     The original routine uses all-to-all communication because
     !     of the distribution of ovlp3KS elements...

     !call aims_stop("casida_setup_scalapack not yet implemented")

     if(myid == 0) then
        write(use_unit,'(2x,a)') '| Entering casida_setup_scalapack part'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)

     ! Determine local start and end of matrix
     ! This shall be used when we reduce the spectrum,
     ! TODO: Do this later, for now we are working on the full matrix
     ! first_found = .false.
     ! do i_state = low_limit, occ_states
     !   if(own_dim1_o3KS(i_state) == myp1_o3KS) then
     !     if (.not. first_found) then
     !       i_local_start = loc_dim1_o3KS(i_state)
     !       first_found = .true.
     !     end if
     !     i_local_end = loc_dim1_o3KS(i_state)
     !   end if
     ! end do
      
     ! first_found = .false.
     ! do a_state = occ_states+1, n_states-unocc_limit
     !   if(own_dim2_o3KS(a_state) == myp2_o3KS) then
     !     if (.not. first_found) then
     !       a_local_start = loc_dim2_o3KS(a_state)
     !       first_found = .true.
     !     end if
     !     a_local_end = loc_dim2_o3KS(a_state)
     !   end if
     ! end do

     !if(myid == 0) then
     !   write(use_unit,'(2x,a)') '| Local indices set'
     !   write(use_unit,'(2x,a)') '| Occupied states:'
     !endif
     !call MPI_BARRIER(mpi_comm_world,errnum)
     !write(use_unit,'(2x,a,i4,a,i4,a,i4)') '| Process ', myid," : ",&
     !  i_local_start, " to ", i_local_end 
     !call MPI_BARRIER(mpi_comm_world,errnum)
     !if(myid == 0) then
     !   write(use_unit,'(2x,a)') '| Unoccupied states:'
     !endif
     !call MPI_BARRIER(mpi_comm_world,errnum)
     !write(use_unit,'(2x,a,i4,a,i4,a,i4)') '| Process ', myid," : ",&
     !  a_local_start, " to ", a_local_end 
     !call MPI_BARRIER(mpi_comm_world,errnum)

     if(myid == 0) then
        write(use_unit,'(2x,a)') '| Setup BLACS descriptors'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)


     ! Setup new blacs grid
     ks3_blacs_context = mpi_comm_world
     call blacs_gridinit(ks3_blacs_context, 'C', 1, n_tasks)
     if(myid == 0) then
        write(use_unit,'(2x,a)') '| ks3 BLACS gridinit done'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)

     call blacs_gridinfo(ks3_blacs_context, ks3_prows, ks3_pcols, &
           ks3_my_p_row, ks3_my_p_col)
    
     if(myid == 0) then
        write(use_unit,'(2x,a)') '| ks3 BLACS info done'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)


     ks3_nrows = numroc(n_basbas, n_basbas, ks3_my_p_row, 0, ks3_prows)
     ks3_ncols = numroc(n_pair_states, 1, ks3_my_p_col, 0, ks3_pcols)

      if(myid == 0) then
        write(use_unit,'(2x,a)') '| KS3 Numroc done'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)


     ! Set up new scalapack descriptor
     call descinit(ks3_sc_desc, n_basbas, n_pair_states, n_basbas, 1, 0, 0,&
           ks3_blacs_context, ks3_nrows, blacs_info)

      if(myid == 0) then
        write(use_unit,'(2x,a)') '| KS3 BLACS Init done'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)


     ! Setup new blacs grid for metric
     metric_blacs_context = mpi_comm_world
     call blacs_gridinit(metric_blacs_context,'C',1,n_tasks)
     if(myid == 0) then
        write(use_unit,'(2x,a)') '| METRIC BLACS grid init done'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)

     call blacs_gridinfo(metric_blacs_context,metric_prows,metric_pcols,&
           metric_my_p_row,metric_my_p_col)

     if(myid == 0) then
        write(use_unit,'(2x,a)') '| METRIC BLACS info done'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)
 
     metric_nrows = numroc(n_states*n_states, n_states*n_states, &
           metric_my_p_row, 0, metric_prows)
     metric_ncols = numroc(n_states*n_states, 1, metric_my_p_col, &
           0, metric_pcols)

     if(myid == 0) then
        write(use_unit,'(2x,a)') '| METRIC Numroc done'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)


     ! Set up new scalapack descriptor
     call descinit(metric_sc_desc, n_pair_states, n_pair_states, &
           n_pair_states, 1, 0, 0,&
           metric_blacs_context, metric_nrows, blacs_info)

     if(myid == 0) then
        write(use_unit,'(2x,a)') '| Metric BLACS Init done'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)


     !call MPI_BARRIER(mpi_comm_world,errnum)
     !write(use_unit,'(2x,a,i4,a,i4,a,i4)') '| Process ', myid," has ",&
     !  (i_local_end - i_local_start + 1) * (a_local_end - a_local_start + 1),&
     !  " cols and may have ", ks3_ncols 
     !call MPI_BARRIER(mpi_comm_world,errnum)

     ! Setup metric to match scalapack distribution
     allocate(metric(n_states, n_states, ndim1_o3KS, ndim2_o3KS), stat=errnum)
     call check_allocation(errnum,' metric for 3KS')

     if(myid == 0) then
        write(use_unit,'(2x,a)') '| Metric initialized'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)


     metric = 0d0
     do i_state = 1, occ_states
       if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
         do a_state = occ_states+1, n_states
           if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle
              i_state_loc = (i_state - 1) &
                          + (occ_states + 1) * (a_state - (occ_states+1))
              a_state_loc = mod((a_state - 1),num_prow) &
                          + num_prow * mod((i_state - 1),num_pcol)
              metric(i_state, a_state, i_state_loc, a_state_loc) = 1d0
         end do
     end do

     if(myid == 0) then
        write(use_unit,'(2x,a)') '| Multiply metric'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)

     do i_state = 0, n_tasks - 1
       if(myid == i_state) then
         write(use_unit,*) n_pair_states, " ", ndim1_o3ks, " ", ndim2_o3ks,&
         " ", n_states 
       endif
       call MPI_BARRIER(mpi_comm_world,errnum)
     end do

       call MPI_BARRIER(mpi_comm_world,errnum)

     ! Multiply metric
     call pdgemm ('N','N', n_basbas, n_states*n_states, n_states*n_states, &
           1d0,&
           ovlp_3ks(1:n_basbas,1:ndim1_o3KS,1:ndim2_o3KS,i_spin), &
           1, 1, ks3_sc_desc,&
           metric(1:n_states,1:n_states,1:ndim1_o3KS,1:ndim2_o3KS), &
           1, 1, metric_sc_desc,&
           0d0,&
           ovlp_3ks, 1, 1, ks3_sc_desc)
     call pdgemm ('N','N', n_basbas, n_states*n_states, n_states*n_states,&
           1d0,&
           loc_o3ks(1:n_basbas,1:ndim1_o3KS,1:ndim2_o3KS,i_spin), &
           1, 1, ks3_sc_desc,&
           metric(1:n_states,1:n_states,1:ndim1_o3KS,1:ndim2_o3KS), &
           1, 1, metric_sc_desc,&
           0d0,&
           loc_o3ks, 1, 1, ks3_sc_desc)
     call pdgemm ('N','N', n_basbas, n_states*n_states, n_states*n_states,&
           1d0,&
           pure_3ks(1:n_basbas,1:ndim1_o3KS,1:ndim2_o3KS,i_spin), &
           1, 1, ks3_sc_desc,&
           metric(1:n_states,1:n_states,1:ndim1_o3KS,1:ndim2_o3KS), &
           1, 1, metric_sc_desc,&
           0d0,&
           pure_3ks, 1, 1, ks3_sc_desc)
     call pdgemm ('N','N', n_basbas, n_states*n_states, n_states*n_states,&
           1d0,&
           loc_f3ks(1:n_basbas,1:ndim1_o3KS,1:ndim2_o3KS,i_spin), &
           1, 1, ks3_sc_desc,&
           metric(1:n_states,1:n_states,1:ndim1_o3KS,1:ndim2_o3KS), &
           1, 1, metric_sc_desc,&
           0d0,&
           loc_f3ks, 1, 1, ks3_sc_desc)
     deallocate(metric)

   if(myid == 0) then
        write(use_unit,'(2x,a)') '| <ia|V|jb>'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)



     ! K = <ia|V|jb>
     call pdgemm ('T','N', n_pair_states, n_pair_states, n_basbas,&
           1d0,&
           ovlp_3ks(1:n_basbas,1:ndim1_o3KS,1:ndim2_o3KS,i_spin), &
           1, 1, ks3_sc_desc,&
           loc_o3ks(1:n_basbas,1:ndim2_o3KS,1:ndim2_o3KS,i_spin), &
           1, 1, ks3_sc_desc,&
           0d0,&
           par_omega, 1, 1, casida_desc)
     
     if(myid == 0) then
        write(use_unit,'(2x,a)') '| <ia|V|jb> calculated'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)


     ! K = K + <ia|f_xc|jb>
     call pdgemm ('T','N', n_pair_states, n_pair_states, n_basbas,&
           1d0,&
           pure_3ks(1:n_basbas,1:ndim1_o3KS,1:ndim2_o3KS,i_spin), 1, 1, ks3_sc_desc,&
           loc_f3ks(1:n_basbas,1:ndim1_o3KS,1:ndim2_o3KS,i_spin), 1, 1, ks3_sc_desc,&
           1d0,&
           par_omega, 1, 1, casida_desc)
     if(myid == 0) then
        write(use_unit,'(2x,a)') '| <ia|f_xc|jb> calculated'
     endif

     call MPI_BARRIER(mpi_comm_world,errnum)


     ! Omega = 2 sqrt(f_i(e_a - e_i)) K_(ia,jb) sqrt(f_j(e_b - e_j)) + delta_(i,j) delta_(a,b)(e_a - e_i)^2
     do i_state = low_limit, occ_states
       do a_state = occ_states + 1, n_states - unocc_limit
         do j_state = low_limit, occ_states
           do b_state = occ_states + 1, n_states - unocc_limit
             ! Get local indices
             call infog2l(mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc, &
                       nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)
          
             ! If these are my elements then scale 2 sqrt(f_i(e_a - e_i)) K_(ia,jb) sqrt(f_j(e_b - e_j))
             if(myrow.eq.iarow.and.mycol.eq.iacol) then
              
               par_omega(iia,jja) = 2d0 * &
                 dsqrt(occ_numbers(i_state,i_spin,n_k_points)) * &
                 dsqrt(ks_eigenvalue(a_state,i_spin,n_k_points) - &
                       ks_eigenvalue(i_state,i_spin,n_k_points)) * &
                 par_omega(iia,jja) * &
                 dsqrt(occ_numbers(j_state,i_spin,n_k_points)) * &
                 dsqrt(ks_eigenvalue(b_state,i_spin,n_k_points) - &
                       ks_eigenvalue(j_state,i_spin,n_k_points))
               ! Is this a main diagonal element? Then add delta_(i,j) delta_(a,b)(e_a - e_i)^2
               if (i_state == j_state .and. a_state == b_state) then
                 par_omega(iia,jja) = par_omega(iia,jja) &
                                    +( ks_eigenvalue(a_state,i_spin,n_k_points) &
                                     - ks_eigenvalue(i_state,i_spin,n_k_points) )**2
               end if
             end if
           end do
         end do
       end do
     end do
     if(myid == 0) then
        write(use_unit,'(2x,a)') '| Omega finalized'
     endif
     call MPI_BARRIER(mpi_comm_world,errnum)
     call blacs_gridexit(ks3_blacs_context)
     call blacs_gridexit(metric_blacs_context)


  else
    ! JK: This is my original implementation, I have not tested any scalapack above

    ! Construct local Casida Matrix OMEGA
    ! Start rotating ovlp_3ks elements
    do i_rot=1, n_tasks
      if(myid==0) write(use_unit,'(2x,a,I5.1)') '| starting rotation: ', i_rot

      if(myid==i_rot-1) then
        tmp1r = tmp1s
        tmp2r = tmp2s
        loc_o3ks = ovlp_3ks
        loc_f3ks = fxc_3ks
      endif

      call get_timestamps(time_comm, clock_time_comm)
      call mpi_bcast(tmp1r, 2*sendrecv_dims, mpi_integer, i_rot-1,&
            mpi_comm_world, errnum)
      call mpi_bcast(tmp2r, 2*sendrecv_dims, mpi_integer, i_rot-1,&
            mpi_comm_world, errnum)
      call mpi_bcast(loc_o3ks, n_basbas*ndim1_o3KS*ndim2_o3KS,&
            mpi_double_precision, i_rot-1, mpi_comm_world, errnum)
      call mpi_bcast(loc_f3ks, n_basbas*ndim1_o3KS*ndim2_o3KS,&
            mpi_double_precision, i_rot-1, mpi_comm_world, errnum)
      call get_times(time_comm, clock_time_comm, time_3ks_comm, &
            clock_time_3ks_comm)

      call get_timestamps(time_calc, clock_time_calc)
      len_of_recs = 1
      do i_spin = 1, n_spin
        do i_state=low_limit, occ_states
          if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
          do a_state=occ_states+1, n_states-unocc_limit
            if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle

            i_state_loc = loc_dim1_o3KS(i_state)
            a_state_loc = loc_dim2_o3KS(a_state)

            do k_cnt=1, sendrecv_dims
 
              if(tmp1r(k_cnt,1).ne.0 .and. tmp2r(k_cnt,1).ne.0) then

                kernel = ddot(n_basbas,&
                              ovlp_3ks(:,i_state_loc,a_state_loc,i_spin),1,&
                              loc_o3ks(:,tmp1r(k_cnt,2),tmp2r(k_cnt,2),i_spin),1)
               
                kernel = kernel &
                       + ddot(n_basbas,&
                              pure_3ks(:,i_state_loc,a_state_loc,i_spin),1,&
                              loc_f3ks(:,tmp1r(k_cnt,2),tmp2r(k_cnt,2),i_spin),1)
                 
                !do l_cnt=1, n_basbas
                !  kernel = kernel + (  ovlp_3ks(l_cnt,i_state_loc,a_state_loc,i_spin)* &
                !     loc_o3ks(l_cnt,tmp1r(k_cnt,2),tmp2r(k_cnt,2),i_spin) + &
                !  pure_3ks(l_cnt,i_state_loc,a_state_loc,i_spin)* &
                !     loc_f3ks(l_cnt,tmp1r(k_cnt,2),tmp2r(k_cnt,2),i_spin) )
                !enddo

                tmp_real = 2d0 * &
                  dsqrt(occ_numbers(i_state,i_spin,n_k_points)) * &
                  dsqrt(ks_eigenvalue(a_state,i_spin,n_k_points) - &
                        ks_eigenvalue(i_state,i_spin,n_k_points)) * &
                  kernel * &
                  dsqrt(occ_numbers(tmp1r(k_cnt,1),i_spin,n_k_points)) * &
                  dsqrt(ks_eigenvalue(tmp2r(k_cnt,1),i_spin,n_k_points) - &
                        ks_eigenvalue(tmp1r(k_cnt,1),i_spin,n_k_points))

                if(mat_idx(i_state,a_state).eq.mat_idx(tmp1r(k_cnt,1),tmp2r(k_cnt,1))) then
                  tmp_real = tmp_real + &
                     ( ks_eigenvalue(a_state,i_spin,n_k_points) &
                     - ks_eigenvalue(i_state,i_spin,n_k_points) )**2
                endif

                if(abs(tmp_real).gt.sparse_thresh) then
                  mpi_bcastidx(len_of_recs,1) = mat_idx(i_state,a_state)
                  mpi_bcastidx(len_of_recs,2) = mat_idx(tmp1r(k_cnt,1),tmp2r(k_cnt,1))
                  mpi_vec(len_of_recs) = tmp_real
                  len_of_recs = len_of_recs + 1
                endif

              endif

            enddo 
          enddo 
        enddo
      enddo
      call get_times(time_calc, clock_time_calc, time_casida_matrix_calculation, &
            clock_time_casida_matrix_calculation)
      
      call get_timestamps(time_comm, clock_time_comm)
      do j_rot=1, n_tasks
        if(myid==j_rot-1) then
          do i_state=1, sendrecv_dims
            tmp1r(i_state,:) = mpi_bcastidx(i_state,:)
            tmp_vecr(i_state) = mpi_vec(i_state)
          enddo
        endif

        call mpi_bcast(tmp1r, 2*sendrecv_dims, mpi_integer, j_rot-1, mpi_comm_world, errnum)
        call mpi_bcast(tmp_vecr, sendrecv_dims, mpi_double_precision, j_rot-1, mpi_comm_world, errnum)

        do i_state = 1, sendrecv_dims
          if(tmp1r(i_state,1).ne.0 .and. tmp1r(i_state,2).ne.0) then
            call infog2l(tmp1r(i_state,1), tmp1r(i_state,2), casida_desc, &
              nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)
        
            if(myrow.eq.iarow.and.mycol.eq.iacol) then
              par_omega(iia,jja) = tmp_vecr(i_state)
            endif

            if(tmp1r(i_state,1).ne.tmp1r(i_state,2)) then
              call infog2l(tmp1r(i_state,2), tmp1r(i_state,1), casida_desc, &
                nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)
              if(myrow.eq.iarow.and.mycol.eq.iacol) then
                par_omega(iia,jja) = tmp_vecr(i_state)
              endif
            endif
          endif
        enddo
      enddo !j_rot
    
      call get_times(time_comm, clock_time_comm, time_omegaresdist, &
            clock_time_omegaresdist)

    enddo  ! I_ROT rotation

    endif

    deallocate(pure_3ks)
    deallocate(fxc_3KS)
    deallocate(ovlp_3KS)
    deallocate(tmp1r)
    deallocate(tmp2r)
    deallocate(tmp1s)
    deallocate(tmp2s)
    deallocate(loc_o3ks)
    deallocate(loc_f3ks)
    deallocate(mpi_vec)
    deallocate(mpi_bcastidx)
    deallocate(tmp_vecr)

    call get_times(time_omegabuild, clock_time_omegabuild)

  else 
    ! it's triplet now, it's is easy, only diagonal energy diffs matter,
    ! all of the above communication is irrelevant for triplets

    call get_timestamps(time_omegabuild, clock_time_omegabuild)

    allocate(par_omega(na_rows,na_cols),stat=errnum)
    call check_allocation(errnum, 'triplet par_omega allocation')
    par_omega = 0d0

    do i_spin = 1, n_spin
      do i_state = low_limit, occ_states
        do a_state = occ_states + 1, n_states - unocc_limit
          call infog2l(mat_idx(i_state,a_state), mat_idx(i_state,a_state), casida_desc, &
                       nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)
          if(myrow.eq.iarow.and.mycol.eq.iacol) then
            par_omega(iia,jja) = ( ks_eigenvalue(a_state,i_spin,n_k_points) - &
                                   ks_eigenvalue(i_state,i_spin,n_k_points) )**2
          endif
        enddo
      enddo
    enddo

    call get_times(time_omegabuild, clock_time_omegabuild)

  endif ! .not.triplet

  allocate(par_evecs(na_rows, na_cols),stat=errnum)
  call check_allocation(errnum,'casida parallel eigenvectors')
  par_evecs = 0d0
  allocate(par_evals(n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida parallel eigenvalues')
  par_evals = 0d0

  ! Solve Casida eigenvalue problem OMEGA A_q = w_q^2 A_q (Equation 1 in paper)
  call get_timestamps(time_omega_solve, clock_time_omega_solve)

  if(myid==0) then
     write(use_unit,'(2x,a)') '|'
     write(use_unit,'(2x,a)') '| Calling ELPA 2stage solver now'
     write(use_unit,'(2x,a)') '|'
  endif

  call aims_elsi_stdevp(n_pair_states, na_rows, na_cols, n_pair_states, &
    par_omega, par_evals, par_evecs, my_casida_ctxt, num_blocks, &
    mpi_comm_global)

  call get_times(time_omega_solve, clock_time_omega_solve)

  deallocate(par_omega)

  ! Get w_q from w_q^2 
  do i_state=1, n_pair_states
    if(par_evals(i_state).lt.0d0) then
      call neg_energy_error(par_evals(i_state))
      par_evals(i_state) = 0d0
    else
      par_evals(i_state) = dsqrt(par_evals(i_state))
    endif
  enddo

  ! Calculate the Oscillator Strengths a_q
  ! a_q = 2/3 (|dip_x S^{1/2} A_q|^2 + |dip_y S^{1/2} A_q|^2 + |dip_z S^{1/2} A_q|^2)  
  allocate(osc_str(n_pair_states),stat=errnum)
  call check_allocation(errnum,'Casida Oscillator Strengths allocation')
  allocate(tmom(n_pair_states,3),stat=errnum)
  call check_allocation(errnum,'moment buffer allocation')
  tmom = 0d0

  do i_spin = 1, n_spin
    ! Loop over all pair states
    do i_state=1, n_pair_states
      ! Loop over all occupied states in energy window
      do j_state=low_limit, occ_states
        ! Loop over all unoccupied states in energy window
        do b_state=occ_states+1, n_states-unocc_limit

          ! Get offset in global matrix
          call infog2l(mat_idx(j_state,b_state), i_state, casida_desc, nprow, &
            npcol, myrow, mycol, iia, jja, iarow, iacol)

          if(myrow.eq.iarow .and. mycol.eq.iacol) then
  
            ! S^{1/2} * A_q = sqrt( 2 * (e_unocc - e_occ) ) * A_q 
            tmp_real = dsqrt(2d0*(ks_eigenvalue(b_state,i_spin,n_k_points) &
              - ks_eigenvalue(j_state,i_spin,n_k_points))) * par_evecs(iia,jja) 
          
            ! dip * S^{1/2} * Aq 
            tmom(i_state,1:3) = tmom(i_state,:) &
              + tmp_real * dipole_moment(j_state,b_state-occ_states,i_spin,1:3)
          
          endif
        enddo
      enddo
    enddo
  enddo

  call sync_matrix(tmom, n_pair_states, 3)

  osc_str(:) = 2d0/3d0 * ( tmom(:,1)**2 + tmom(:,2)**2 + tmom(:,3)**2 )

  if(myid==0) then

    if(excited_mode_triplet) then
  
      write(use_unit,'(2x,a)') '| Writing Triplet excitation energies (TDDFT) '
      write(info_str,'(a)') '. Triplet: '
  
    else

      write(use_unit,'(2x,a)') '| Writing Singlet excitation energies (TDDFT) '
      write(info_str,'(a)') '. Singlet: '
    
    endif

    do i_state=1, num_excited_states
      write(use_unit,'(2x,a,i4.1,a11,f11.4,a,f8.4,a,f12.4)') '| ', i_state, info_str, &
        par_evals(i_state)*hartree, ' eV     [', &
        par_evals(i_state), ' Ha]  -  Oscillator Strength: ', osc_str(i_state)
    
      write(use_unit,'(2x,a,f10.4,a,f10.4,a,f10.4)') &
        '|            Transition Moments    X:', tmom(i_state,1), '   Y:', &
        tmom(i_state,2), '   Z:', tmom(i_state,3)
    
      write(use_unit,'(2x,a)') '|'
    enddo

    write(use_unit,'(2x,a)') '| '
  endif

  if(myid==0) then
    if(excited_mode_triplet) then
      write(info_str,'(a)') 'TDDFT_LR_Spectrum_Triplet.dat'
    else
      write(info_str,'(a)') 'TDDFT_LR_Spectrum_Singlet.dat'
    endif
  open(100,file=info_str)
  do i_state=1, n_pair_states
    write(100,'(f29.14,4x,f29.14)') par_evals(i_state)*hartree, osc_str(i_state)
  enddo
  close(100)
  endif

  ! final cleanup
  deallocate(tmom)
  deallocate(osc_str)
  deallocate(par_evecs)
  deallocate(par_evals)

  if(myid==0) then
    if(excited_mode_triplet) then
      write(use_unit,'(2x,a)') '| End of TDDFT Triplet excitation energy calculation.'
    else
      write(use_unit,'(2x,a)') '| End of TDDFT Singlet excitation energy calculation.'
    endif 
  endif

  call get_times(time_neutral_excitation_casida, &
           clock_time_neutral_excitation_casida)

end subroutine casida_tddft_2d

subroutine casida_tddft()
  implicit none

  integer                             :: k_cnt, errnum, dummy_integer, vecs_found
  integer                             :: print_percent
  integer, allocatable, dimension(:)  :: lapack_ifail, lapack_iwork

  real*8                              :: kernel, lapack_abstol, temp_real
  real*8, allocatable, dimension(:)   :: eig_vals, lapack_space_1d, osc_str
  real*8, allocatable, dimension(:,:) :: omega, eig_vecs, tmom

  if(n_pair_states.lt.10) then
    print_percent = 1
  else
    print_percent = n_pair_states/10
  endif

  allocate(omega(n_pair_states,n_pair_states),stat=errnum)
  call check_allocation(errnum, 'casida omega')
  omega = 0d0

  call get_timestamps(time_omegabuild, clock_time_omegabuild)

  write(use_unit,'(2x,a)') '|'
  write(use_unit,'(2x,a)') '| Building Omega matrix  < ia | V + f_xc | jb >  now.'
  write(use_unit,'(2x,a)') '|'

  if(.not.excited_mode_triplet) then
    do i_state = 1, occ_states
      do a_state = occ_states+1, n_states

        if(mod(mat_idx(i_state,a_state),print_percent)==0) then
          write(use_unit,'(2x,a,i8.1,a,f10.2,a)') '|  < ia |  pair state ', mat_idx(i_state,a_state), '  ', &
            dble(mat_idx(i_state,a_state))/dble(n_pair_states)*100d0, '%'
        endif

        do j_state=1, occ_states
          do b_state=occ_states+1, n_states

            kernel = 0d0
            do k_cnt = 1, n_basbas
      
              kernel = kernel &
                + ovlp_3ks(k_cnt,i_state,a_state,n_spin) * ovlp_3ks(k_cnt,j_state,b_state,n_spin) &
                + pure_3ks(k_cnt,i_state,a_state,n_spin) *  fxc_3ks(k_cnt,j_state,b_state,n_spin) 
        
            enddo

            omega(mat_idx(i_state,a_state),mat_idx(j_state,b_state)) = &
              4d0 &
              * dsqrt(KS_eigenvalue(a_state,n_spin,n_k_points) - KS_eigenvalue(i_state,n_spin,n_k_points)) &
              * kernel &
              * dsqrt(KS_eigenvalue(b_state,n_spin,n_k_points) - KS_eigenvalue(j_state,n_spin,n_k_points))

          enddo
        enddo
      enddo
    enddo
  endif

  do i_state=1, occ_states
    do a_state=occ_states+1, n_states
      omega(mat_idx(i_state,a_state),mat_idx(i_state,a_state)) = &
        omega(mat_idx(i_state,a_state),mat_idx(i_state,a_state)) & 
        + (KS_eigenvalue(a_state,n_spin,n_k_points) - KS_eigenvalue(i_state,n_spin,n_k_points))**2
    enddo
  enddo
  
  write(use_unit,'(2x,a)') '|'
  call get_times(time_omegabuild, clock_time_omegabuild)
  write(use_unit,'(2x,a)') '|'

  !- allocate lapack spaces and call eigensolver
  allocate(lapack_space_1d(8 * n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_workspace')
  lapack_space_1d = 0d0
  allocate(eig_vals(n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_eig_vals')
  eig_vals = 0d0
  allocate(eig_vecs(n_pair_states,n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_eig_vecs')
  eig_vecs = 0d0
  allocate(lapack_iwork(5*n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_iwork')
  lapack_iwork = 0
  allocate(lapack_ifail(n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_ifail')
  lapack_ifail = 0

  lapack_abstol = 2d0 * dlamch()

  call dsyevx( 'V' , 'A' , 'U' , n_pair_states , omega , n_pair_states , &
               dummy_integer, dummy_integer, dummy_integer, dummy_integer, lapack_abstol , &
               vecs_found , eig_vals , eig_vecs , n_pair_states , &
               lapack_space_1d , 8 * n_pair_states , lapack_iwork , lapack_ifail, errnum )

  if((myid==0).and.(errnum.ne.0)) then
    write(info_str,*) 'DSYEVX reportet Error: ',errnum
  endif

  deallocate(lapack_space_1d)
  deallocate(lapack_iwork)
  deallocate(lapack_ifail)
  deallocate(omega)

  !- Omega F = omega^2 F
  do i_state=1, n_pair_states
    if(eig_vals(i_state).lt.0d0) then
      call neg_energy_error(eig_vals(i_state))
      eig_vals(i_state) = 0d0
    else
      eig_vals(i_state) = dsqrt(eig_vals(i_state))
    endif
  enddo

  allocate(tmom(n_pair_states,3),stat=errnum)
  call check_allocation(errnum,'tmom allocation')
  tmom = 0d0

  do i_state=1, n_pair_states
    do j_state=1, occ_states
      do b_state=occ_states+1, n_states
        
        temp_real = dsqrt(2d0*(ks_eigenvalue(b_state,n_spin,n_k_points) &
                             - ks_eigenvalue(j_state,n_spin,n_k_points))) &
                  * eig_vecs(mat_idx(j_state,b_state),i_state)

        tmom(i_state,1:3) = tmom(i_state,1) &
          + temp_real * dipole_moment(j_state,b_state-occ_states,n_spin,1:3)
  
      enddo
    enddo
  enddo

  allocate(osc_str(n_pair_states),stat=errnum)
  call check_allocation(errnum,'osc_str error')
  osc_str = 0d0

  do i_state=1, n_pair_states
    osc_str(i_state) = 2.0d0/3.0d0 * ( tmom(i_state,1)**2 + tmom(i_state,2)**2 + tmom(i_state,3)**2 )
  enddo

  if(excited_mode_triplet) then
    write(use_unit,'(2x,a)') '| Writing Triplet excitation energies (TDDFT) '
    write(info_str,'(a)') '. Triplet: '
  else
    write(use_unit,'(2x,a)') '| Writing Singlet excitation energies (TDDFT) '
    write(info_str,'(a)') '. Singlet: '
  endif
  do i_state=1, num_excited_states
    write(use_unit,'(2x,a,i4.1,a11,f11.4,a,f8.4,a,f12.4,a,f8.4,a,f8.4)') &
       '| ', i_state, info_str, eig_vals(i_state)*hartree, ' eV     [', &
       eig_vals(i_state), ' Ha]  -  Oscillator Strength: ', osc_str(i_state)
    write(use_unit,'(2x,a,f10.4,a,f10.4,a,f10.4)') &
       '|            Transition Moments    X:', tmom(i_state,1), '   Y:', &
       tmom(i_state,2), '   Z:', tmom(i_state,3)
    write(use_unit,'(2x,a)') '|'
  enddo

  if(excited_mode_triplet) then
    write(info_str,'(a)') 'TDDFT_LR_Spectrum_Triplet.dat'
  else
    write(info_str,'(a)') 'TDDFT_LR_Spectrum_Singlet.dat'
  endif
  open(100,file=info_str)
  do i_state=1, n_pair_states
    write(100,'(f29.14,4x,f29.14)') eig_vals(i_state)*hartree, osc_str(i_state)
  enddo
  close(100)

!- cleanup the rest
  deallocate(osc_str)
  deallocate(eig_vals)
  deallocate(eig_vecs)
end subroutine casida_tddft

subroutine casida_rpa_2d()
  implicit none
  integer :: i_state_loc, j_state_loc, a_state_loc, b_state_loc, &
              k_cnt, l_cnt, m_cnt, errnum, &
              ia, ja, tmp_int, i_task, sendrecv_dims, &
              len_of_recs, max_na_rows, max_na_cols

  integer, allocatable, dimension(:) :: mpi_bcastproc, tmp_pr, tmp_ps
  integer, allocatable, dimension(:,:) :: tmp1r, tmp2r, tmp1s, tmp2s, tmp_bcastidx, mpi_bcastidx

  real*8 :: kernel, tmp_real, total_mem, proc_mem, time_comm, clock_time_comm
  real*8, allocatable, dimension(:) :: par_evals, osc_str, tmp_vecr, tmp_vecs, mpi_vec
  real*8, allocatable, dimension(:,:) :: par_evecs, par_omega, tmom
  real*8, allocatable, dimension(:,:,:,:) :: tmp_o3ks, loc_o3ks

  ! ELSI
  real*8, dimension(:,:), allocatable :: dummy

  total_mem = 8d0 * dble(n_pair_states)**2

  if(myid==0) then
    if(excited_mode_triplet) then
      write(use_unit,'(2x,a)') '| '
      write(use_unit,'(2x,a)') '| Starting RPA excited state calculation for Triplets.'
      write(use_unit,'(2x,a)') '| '
    else
      write(use_unit,'(2x,a)') '| '
      write(use_unit,'(2x,a)') '| Starting RPA excited state calculation for Singlets.'
      write(use_unit,'(2x,a)') '| '
    endif
    write(use_unit,'(2x,a,i8.1,a,i8.1,a,i6.1,a,i6.1,a,i5.1,a)') '| Deviding the ', &
      n_pair_states, ' x ', n_pair_states, ' Matrix into ', &
      num_pcol,' colomns and ', num_prow, ' rows on ', n_tasks, ' processors.'
    write(use_unit,'(2x,a,f10.3,a)') '| Memory load of ', total_mem/2d0**30, ' GB devided among all Processors.'
  endif

  proc_mem = 8d0 * dble(na_rows) * dble(na_cols)

  if( (total_mem/2d0**20)/dble(n_tasks).gt.1024d0) then
    write(info_str,'(2x,a,i5.1,a,i8.1,a,i8.1,a,f8.3,a)') '| Processor ', myid, &
      ' gets ', na_rows, ' rows and ', na_cols,' columns  (',proc_mem/2d0**20/1024d0,' GB).'
    call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)
  else
    write(info_str,'(2x,a,i5.1,a,i8.1,a,i8.1,a,f8.3,a)') '| Processor ', myid, &
      ' gets ', na_rows, ' rows and ', na_cols,' columns  (',proc_mem/2d0**20,' MB).'
    call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)
  endif

  if(myid==0) write(use_unit,'(2x,a)') '|'
  if(myid==0) write(use_unit,'(2x,a,I5.1)') '| Number of total rotations for matrix building: ', n_tasks
  if(myid==0) write(use_unit,'(2x,a)') '|'

  k_cnt = 0
  do i_state=low_limit, occ_states
    if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
      do a_state=occ_states+1, n_states-unocc_limit
        if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle

  k_cnt = k_cnt + 1

  enddo; enddo

  if(myid==0) then
  do i_proc=1, n_tasks-1
    call mpi_recv(tmp_int, 1, mpi_integer, i_proc, &
      100, mpi_comm_world, mpi_status, errnum)
    if(tmp_int.gt.k_cnt) k_cnt = tmp_int
  enddo
  else
    call mpi_send(k_cnt, 1, mpi_integer, 0, 100, mpi_comm_world, errnum)
  endif
  i_state = k_cnt

  k_cnt = 0
  do j_state=low_limit, occ_states
    if(own_dim1_o3KS(j_state) /= myp1_o3KS) cycle
      do b_state=occ_states+1, n_states-unocc_limit
        if(own_dim2_o3KS(b_state) /= myp2_o3KS) cycle

    k_cnt = k_cnt + 1

  enddo; enddo

  if(myid==0) then
  do i_proc=1, n_tasks-1
    call mpi_recv(tmp_int, 1, mpi_integer, i_proc, 100, mpi_comm_world, mpi_status, errnum)
    if(tmp_int.gt.k_cnt) k_cnt = tmp_int
  enddo
  else
    call mpi_send(k_cnt, 1, mpi_integer, 0, 100, mpi_comm_world, errnum)
  endif
  sendrecv_dims = i_state * k_cnt

  call mpi_bcast(sendrecv_dims, 1, mpi_integer, 0, mpi_comm_world, errnum)

  allocate(mpi_bcastidx(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'mpi broadcast index')
  mpi_bcastidx = 0
  allocate(mpi_vec(sendrecv_dims),stat=errnum)
  call check_allocation(errnum, 'mpi broadcast vector')
  mpi_vec = 0d0

  time_omegaresdist = 0d0
  clock_time_omegaresdist = 0d0
  call get_timestamps(time_omegabuild, clock_time_omegabuild)

  if(myid==0) then
    write(use_unit,'(2x,a)') '| Building Omega matrix  < ia | V + f_xc | jb >  now.'
    write(use_unit,'(2x,a)') '|'
    write(use_unit,'(2x,a,en10.1,a)') '| Only final results greater than ', sparse_thresh, ' will be stored.'
    write(use_unit,'(2x,a)') '|'
  endif

  if(.not.excited_mode_triplet) then

  allocate(tmp1r(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, ' tmp1r res dist')
  allocate(tmp_vecr(sendrecv_dims),stat=errnum)
  call check_allocation(errnum, ' tmp_vecr res dist')

  allocate(par_omega(na_rows,na_cols),stat=errnum)
  call check_allocation(errnum, 'parallel omega matrix')
  par_omega = 0d0

  allocate(tmp2r(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'rotation index 2r')
  allocate(tmp1s(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'rotation index 1s')
  allocate(tmp2s(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'rotation index 2s')
  allocate(loc_o3ks(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin),stat=errnum)
  call check_allocation(errnum, 'rotation local o3ks')

  tmp1s(:,1) = 0
  tmp2s(:,1) = 0

  k_cnt=1
  do i_state=low_limit, occ_states
    if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
      do a_state=occ_states+1, n_states-unocc_limit
        if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle

    i_state_loc = loc_dim1_o3KS(i_state)
    a_state_loc = loc_dim2_o3KS(a_state)

    tmp1s(k_cnt,1) = i_state
    tmp2s(k_cnt,1) = a_state
    tmp1s(k_cnt,2) = i_state_loc
    tmp2s(k_cnt,2) = a_state_loc

    k_cnt = k_cnt + 1

  enddo; enddo

! Start rotating ovlp_3ks elements
  do i_rot=1, n_tasks
    if(myid==0) write(use_unit,'(2x,a,I5.1)') '| starting rotation: ', i_rot

  if(myid==i_rot-1) then
    tmp1r = tmp1s
    tmp2r = tmp2s
    loc_o3ks = ovlp_3ks
  endif

  call mpi_bcast(tmp1r, 2*sendrecv_dims, mpi_integer, i_rot-1, mpi_comm_world, errnum)
  call mpi_bcast(tmp2r, 2*sendrecv_dims, mpi_integer, i_rot-1, mpi_comm_world, errnum)
  call mpi_bcast(loc_o3ks, n_basbas*ndim1_o3KS*ndim2_o3KS, mpi_double_precision, i_rot-1, mpi_comm_world, errnum)

  len_of_recs = 1
  do i_state=low_limit, occ_states
    if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
      do a_state=occ_states+1, n_states-unocc_limit
        if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle

        i_state_loc = loc_dim1_o3KS(i_state)
        a_state_loc = loc_dim2_o3KS(a_state)

        do k_cnt=1, sendrecv_dims

    if(tmp1r(k_cnt,1).ne.0 .and. tmp2r(k_cnt,1).ne.0) then

      kernel = 0d0
      do l_cnt=1, n_basbas
      kernel = kernel + (  ovlp_3ks(l_cnt,i_state_loc,a_state_loc,n_spin) * &
                           loc_o3ks(l_cnt,tmp1r(k_cnt,2),tmp2r(k_cnt,2),n_spin) )
      enddo

      tmp_real = 2d0 * &
              dsqrt(occ_numbers(i_state,n_spin,n_k_points)) * &
              dsqrt(ks_eigenvalue(a_state,n_spin,n_k_points) - &
                    ks_eigenvalue(i_state,n_spin,n_k_points)) * &
              kernel * &
              dsqrt(occ_numbers(tmp1r(k_cnt,1),n_spin,n_k_points)) * &
              dsqrt(ks_eigenvalue(tmp2r(k_cnt,1),n_spin,n_k_points) - &
                    ks_eigenvalue(tmp1r(k_cnt,1),n_spin,n_k_points))

      if(mat_idx(i_state,a_state).eq.mat_idx(tmp1r(k_cnt,1),tmp2r(k_cnt,1))) then
         tmp_real = tmp_real + &
         ( ks_eigenvalue(a_state,n_spin,n_k_points) &
         - ks_eigenvalue(i_state,n_spin,n_k_points) )**2
      endif

      if(abs(tmp_real).gt.sparse_thresh) then
         mpi_bcastidx(len_of_recs,1) = mat_idx(i_state,a_state)
         mpi_bcastidx(len_of_recs,2) = mat_idx(tmp1r(k_cnt,1),tmp2r(k_cnt,1))
         mpi_vec(len_of_recs) = tmp_real
         len_of_recs = len_of_recs + 1
      endif

    endif

  enddo; enddo; enddo

  call get_timestamps(time_comm, clock_time_comm)
  do j_rot=1, n_tasks
    if(myid==j_rot-1) then
      do i_state=1, sendrecv_dims
         tmp1r(i_state,:) = mpi_bcastidx(i_state,:)
         tmp_vecr(i_state) = mpi_vec(i_state)
      enddo
    endif

    call mpi_bcast(tmp1r, 2*sendrecv_dims, mpi_integer, j_rot-1, mpi_comm_world, errnum)
    call mpi_bcast(tmp_vecr, sendrecv_dims, mpi_double_precision, j_rot-1, mpi_comm_world, errnum)

    do i_state = 1, sendrecv_dims
    if(tmp1r(i_state,1).ne.0 .and. tmp1r(i_state,2).ne.0) then
      call infog2l(tmp1r(i_state,1), tmp1r(i_state,2), casida_desc, &
         nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)
      if(myrow.eq.iarow.and.mycol.eq.iacol) then
         par_omega(iia,jja) = tmp_vecr(i_state)
      endif
      if(tmp1r(i_state,1).ne.tmp1r(i_state,2)) then
         call infog2l(tmp1r(i_state,2), tmp1r(i_state,1), casida_desc, &
                   nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)
         if(myrow.eq.iarow.and.mycol.eq.iacol) then
            par_omega(iia,jja) = tmp_vecr(i_state)
         endif
      endif
    endif 
    enddo
  enddo !j_rot
  call get_times(time_comm, clock_time_comm)
  time_omegaresdist = time_omegaresdist + time_comm
  clock_time_omegaresdist = clock_time_omegaresdist + clock_time_comm

  enddo  ! I_ROT rotation

  deallocate(ovlp_3KS)
  deallocate(tmp1r)
  deallocate(tmp2r)
  deallocate(tmp1s)
  deallocate(tmp2s)
  deallocate(loc_o3ks)
  deallocate(mpi_vec)
  deallocate(mpi_bcastidx)
  deallocate(tmp_vecr)

  call get_times(time_omegabuild, clock_time_omegabuild)

  else ! it's triplet now, it's is easy, only diagonal energy diffs matter,
       ! all of the above communication is irrelevant for triplets

  call get_timestamps(time_omegabuild, clock_time_omegabuild)

  allocate(par_omega(na_rows,na_cols),stat=errnum)
  call check_allocation(errnum, 'triplet par_omega allocation')
  par_omega = 0d0

  do i_state = low_limit, occ_states
    do a_state = occ_states + 1, n_states - unocc_limit
      call infog2l(mat_idx(i_state,a_state), mat_idx(i_state,a_state), casida_desc, &
                   nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)
      if(myrow.eq.iarow.and.mycol.eq.iacol) then
        par_omega(iia,jja) = ( ks_eigenvalue(a_state,n_spin,n_k_points) - &
                               ks_eigenvalue(i_state,n_spin,n_k_points) )**2
      endif
  enddo; enddo

  call get_times(time_omegabuild, clock_time_omegabuild)

  endif ! .not.triplet

  allocate(par_evecs(na_rows, na_cols),stat=errnum)
  call check_allocation(errnum,'casida parallel eigenvectors')
  par_evecs = 0d0
  allocate(par_evals(n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida parallel eigenvalues')
  par_evals = 0d0

  call get_timestamps(time_omega_solve, clock_time_omega_solve)

  if(myid==0) then
     write(use_unit,'(2x,a)') '|'
     write(use_unit,'(2x,a)') '| Calling ELPA 2stage solver now'
     write(use_unit,'(2x,a)') '|'
  endif

  call aims_elsi_stdevp(n_pair_states, na_rows, na_cols, n_pair_states, &
    par_omega, par_evals, par_evecs, my_casida_ctxt, num_blocks, &
    mpi_comm_global)

  call get_times(time_omega_solve, clock_time_omega_solve)

  deallocate(par_omega)

! OF=w^2F
  do i_state=1, n_pair_states
    if(par_evals(i_state).lt.0d0) then
      call neg_energy_error(par_evals(i_state))
      par_evals(i_state) = 0d0
    else
      par_evals(i_state) = dsqrt(par_evals(i_state))
    endif
  enddo

  allocate(osc_str(n_pair_states),stat=errnum)
  call check_allocation(errnum,'Casida Osci Strs')
  allocate(tmom(n_pair_states,3),stat=errnum)
  call check_allocation(errnum,'tmom2 allocation')
  tmom = 0d0

  do i_state=1, n_pair_states
    do j_state=low_limit, occ_states
      do b_state=occ_states+1, n_states-unocc_limit

  call infog2l(mat_idx(j_state,b_state), i_state, casida_desc, nprow, &
               npcol, myrow, mycol, iia, jja, iarow, iacol)

  if(myrow.eq.iarow .and. mycol.eq.iacol) then
    tmom(i_state,1) = tmom(i_state,1) &
       + dsqrt(2d0*(ks_eigenvalue(b_state,n_spin,n_k_points) &
                  - ks_eigenvalue(j_state,n_spin,n_k_points))) * &
      (par_evecs(iia,jja) * dipole_moment(j_state,b_state-occ_states,n_spin,1))
    tmom(i_state,2) = tmom(i_state,2) &
       + dsqrt(2d0*(ks_eigenvalue(b_state,n_spin,n_k_points) &
                  - ks_eigenvalue(j_state,n_spin,n_k_points))) * &
      (par_evecs(iia,jja) * dipole_moment(j_state,b_state-occ_states,n_spin,2))
    tmom(i_state,3) = tmom(i_state,3) &
       + dsqrt(2d0*(ks_eigenvalue(b_state,n_spin,n_k_points) &
                  - ks_eigenvalue(j_state,n_spin,n_k_points))) * &
      (par_evecs(iia,jja) * dipole_moment(j_state,b_state-occ_states,n_spin,3))
  endif
  enddo; enddo; enddo

  call sync_matrix(tmom, n_pair_states, 3)

  osc_str(:) = 2d0/3d0 * ( tmom(:,1)**2 + tmom(:,2)**2 + tmom(:,3)**2 )

  if(myid==0) then
    if(excited_mode_triplet) then
      write(use_unit,'(2x,a)') '| Writing Triplet excitation energies (RPA) '
      write(info_str,'(a)') '. Triplet: '
    else
      write(use_unit,'(2x,a)') '| Writing Singlet excitation energies (RPA) '
      write(info_str,'(a)') '. Singlet: '
    endif
    do i_state=1, num_excited_states
    write(use_unit,'(2x,a,i4.1,a11,f11.4,a,f8.4,a,f12.4)') '| ', &
              i_state, info_str, par_evals(i_state)*hartree, ' eV     [', &
              par_evals(i_state), ' Ha]  -  Oscillator Strength: ', osc_str(i_state)
    write(use_unit,'(2x,a,f10.4,a,f10.4,a,f10.4)') &
       '|            Transition Moments    X:', tmom(i_state,1), '   Y:', &
       tmom(i_state,2), '   Z:', tmom(i_state,3)
    write(use_unit,'(2x,a)') '|'
    enddo
    write(use_unit,'(2x,a)') '| '
  endif

  if(myid==0) then
   if(excited_mode_triplet) then
     write(info_str,'(a)') 'TDDFT_LR_Spectrum_Triplet.dat'
   else
     write(info_str,'(a)') 'TDDFT_LR_Spectrum_Singlet.dat'
   endif
   open(100,file=info_str)
   do i_state=1, n_pair_states
     write(100,'(f29.14,4x,f29.14)') par_evals(i_state)*hartree, osc_str(i_state)
   enddo
   close(100)
  endif

  ! final cleanup
  deallocate(tmom)
  deallocate(osc_str)
  deallocate(par_evecs)
  deallocate(par_evals)

  if(myid==0) then
  if(excited_mode_triplet) then
    write(use_unit,'(2x,a)') '| End of RPA Triplet excitation energy calculation.'
  else
    write(use_unit,'(2x,a)') '| End of RPA Singlet excitation energy calculation.'
  endif; endif
end subroutine casida_rpa_2d

subroutine casida_rpa()
  implicit none
  integer				:: i_state, j_state, a_state, b_state, &
					   n_pair_states, k_cnt, errnum, dummy_integer, vecs_found
  integer, allocatable, dimension(:)	:: lapack_ifail, lapack_iwork

  real*8				:: kernel, lapack_abstol
  real*8, allocatable, dimension(:)	:: eig_vals, lapack_space_1d, osc_str
  real*8, allocatable, dimension(:,:)	:: omega, eig_vecs, tmom

  intrinsic	:: max, int, nint, dsqrt

  allocate(omega(n_pair_states,n_pair_states),stat=errnum)
  call check_allocation(errnum, 'casida omega')
  omega = 0d0

  call get_timestamps(time_omegabuild, clock_time_omegabuild)

  if(.not.excited_mode_triplet) then
  do i_state = 1, occ_states
    do a_state = occ_states+1, n_states
	do j_state=1, occ_states
	  do b_state=occ_states+1, n_states

    kernel = 0d0
    do k_cnt = 1, n_basbas
      kernel = kernel + ( ovlp_3ks(k_cnt,i_state,a_state,n_spin) * ovlp_3ks(k_cnt,j_state,b_state,n_spin) )
    enddo

    omega(mat_idx(i_state,a_state),mat_idx(j_state,b_state)) = &
        4d0 * &
        dsqrt( KS_eigenvalue(a_state,n_spin,n_k_points) - KS_eigenvalue(i_state,n_spin,n_k_points) ) * &
        kernel * &
        dsqrt( KS_eigenvalue(b_state,n_spin,n_k_points) - KS_eigenvalue(j_state,n_spin,n_k_points) )

  enddo; enddo; enddo; enddo
  endif

  do i_state=1, occ_states
  do a_state=occ_states+1, n_states
    omega(mat_idx(i_state,a_state),mat_idx(i_state,a_state)) = omega(mat_idx(i_state,a_state),mat_idx(i_state,a_state)) + &
	  ( KS_eigenvalue(a_state,n_spin,n_k_points) - KS_eigenvalue(i_state,n_spin,n_k_points) )**2
  enddo; enddo

  call get_times(time_omegabuild, clock_time_omegabuild)

!- allocate lapack spaces and call eigensolver
  allocate(lapack_space_1d(8 * n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_workspace')
  lapack_space_1d = 0d0
  allocate(eig_vals(n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_eig_vals')
  eig_vals = 0d0
  allocate(eig_vecs(n_pair_states,n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_eig_vecs')
  eig_vecs = 0d0
  allocate(lapack_iwork(5*n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_iwork')
  lapack_iwork = 0
  allocate(lapack_ifail(n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_ifail')
  lapack_ifail = 0

  lapack_abstol = 2d0 * dlamch()

  call dsyevx( 'V' , 'A' , 'U' , n_pair_states , omega , n_pair_states , &
               dummy_integer, dummy_integer, dummy_integer, dummy_integer, lapack_abstol , &
               vecs_found , eig_vals , eig_vecs , n_pair_states , &
               lapack_space_1d , 8 * n_pair_states , lapack_iwork , lapack_ifail, errnum )

  if((myid==0).and.(errnum.ne.0)) then
    write(info_str,*) 'DSYEVX reportet Error: ',errnum
  endif

  deallocate(lapack_space_1d)
  deallocate(lapack_iwork)
  deallocate(lapack_ifail)
  deallocate(omega)

!- Omega F = omega^2 F
  do i_state=1, n_pair_states
    if(eig_vals(i_state).lt.0d0) then
      call neg_energy_error(eig_vals(i_state))
      eig_vals(i_state) = 0d0
    else
      eig_vals(i_state) = dsqrt(eig_vals(i_state))
    endif
  enddo

  allocate(tmom(n_pair_states,3),stat=errnum)
  call check_allocation(errnum,'tmom allocation')
  tmom = 0d0

  do i_state=1, n_pair_states
    do j_state=1, occ_states
      do b_state=occ_states+1, n_states

  tmom(i_state,1) = tmom(i_state,1) + &
     dsqrt(2d0*(ks_eigenvalue(b_state,n_spin,n_k_points) &
              - ks_eigenvalue(j_state,n_spin,n_k_points))) * &
     ( eig_vecs(mat_idx(j_state,b_state),i_state) * &
       dipole_moment(j_state,b_state-occ_states,n_spin,1) )
  tmom(i_state,2) = tmom(i_state,2) + &
     dsqrt(2d0*(ks_eigenvalue(b_state,n_spin,n_k_points) &
              - ks_eigenvalue(j_state,n_spin,n_k_points))) * &
     ( eig_vecs(mat_idx(j_state,b_state),i_state) * &
       dipole_moment(j_state,b_state-occ_states,n_spin,2) )
  tmom(i_state,3) = tmom(i_state,3) + &
     dsqrt(2d0*(ks_eigenvalue(b_state,n_spin,n_k_points) &
              - ks_eigenvalue(j_state,n_spin,n_k_points))) * &
     ( eig_vecs(mat_idx(j_state,b_state),i_state) * &
       dipole_moment(j_state,b_state-occ_states,n_spin,3) )

  enddo; enddo; enddo

  allocate(osc_str(n_pair_states),stat=errnum)
  call check_allocation(errnum,'osc_str error')
  osc_str = 0d0

  do i_state=1, n_pair_states
    osc_str(i_state) = 2.0d0/3.0d0 * ( tmom(i_state,1)**2 + tmom(i_state,2)**2 + tmom(i_state,3)**2 )
  enddo

  if(excited_mode_triplet) then
    write(use_unit,'(2x,a)') 'Writing Triplet excitation energies (RPA) '
    write(info_str,'(a)') '. Triplet: '
  else
    write(use_unit,'(2x,a)') 'Writing Singlet excitation energies (RPA) '
    write(info_str,'(a)') '. Singlet: '
  endif
  do i_state=1, num_excited_states
    write(use_unit,'(2x,i4.1,a11,f12.4,a,f8.4,a,f12.4,a,f8.4,a,f7.4)') &
       i_state, info_str, eig_vals(i_state), ' Ha      [', &
       eig_vals(i_state)*hartree, ' eV]  - Oscillator Strength: ', &
       osc_str(i_state)
    write(use_unit,'(2x,a,f12.4,a,f12.4,a,f12.4)') &
       '|     Cartesian transition moments:    X: ', tmom(i_state,1), &
       '    Y: ', tmom(i_state,2), '    Z: ',tmom(i_state,3)
    write(use_unit,'(2x,a)') '| '
  enddo
  write(use_unit,*) ' '

  if(excited_mode_triplet) then
    write(info_str,'(a)') 'TDDFT_LR_Spectrum_Triplet.dat'
  else
    write(info_str,'(a)') 'TDDFT_LR_Spectrum_Singlet.dat'
  endif
  open(100,file=info_str)
  do i_state=1, n_pair_states
    write(100,'(f29.14,4x,f29.14)') eig_vals(i_state)*hartree, osc_str(i_state)
  enddo
  close(100)

!- cleanup the rest
  deallocate(osc_str)
  deallocate(eig_vals)
  deallocate(eig_vecs)
  deallocate(tmom)

end subroutine casida_rpa

subroutine casida_tdhf_2d()
  implicit none

  integer				:: i_state, j_state, a_state, b_state, i_proc, &
					   i_state_loc, j_state_loc, a_state_loc, b_state_loc, tmp_int, &
					   n_pair_states, k_cnt, errnum, &
					   my_prow, my_pcol, num_blocks, comm_tag, &
					   my_casida_ctxt, num_pcol, num_prow, na_cols, na_rows, &
					   casida_desc(9), mpi_status(mpi_status_size), &
					   ia, ja, iia, jja, iarow, iacol, myrow, mycol, nprow, npcol

  integer, allocatable, dimension(:)	:: lapack_iwork
  integer, allocatable, dimension(:,:)	:: who_has
  integer, allocatable, dimension(:,:)	:: ovlp_idx, tmp_idx, ovij_idx, ovab_idx

  real*8				:: kernel, tmp_real, triplet_factor
  real*8, allocatable, dimension(:)	:: lapack_workspace, par_evals, o3ks_buffer, osc_str, tau, wr, wi, lapack_scale
  real*8, allocatable, dimension(:,:)	:: par_lmat, par_mmat, par_evecs, par_lpm, par_lmm, par_omega, par_z, tmom

  integer, external :: numroc

  intrinsic	:: nint, floor, min, max, dsqrt, dble, mod

  write(use_unit,*) 'This routine is non-functional.'
  stop

  do num_pcol = nint(dsqrt(dble(n_tasks))),2,-1
    if(mod(n_tasks,num_pcol)==0) exit
  enddo

  num_prow = n_tasks/num_pcol

  my_casida_ctxt = mpi_comm_world
  call blacs_gridinit( my_casida_ctxt, 'C' , num_prow, num_pcol )
  call blacs_gridinfo( my_casida_ctxt, num_prow, num_pcol, my_prow, my_pcol )

  num_blocks = min(n_tasks,floor(dble(n_pair_states)/dble(n_tasks)))

  na_rows = numroc(n_pair_states, num_blocks, my_prow, 0, num_prow)
  na_cols = numroc(n_pair_states, num_blocks, my_pcol, 0, num_pcol)

  if(myid==0) then
    if(excited_mode_triplet) then
      write(use_unit,'(2x,a)') '| '
      write(use_unit,'(2x,a)') '| Starting TDHF excited state calculation for Triplets.'
      write(use_unit,'(2x,a)') '| '
    else
      write(use_unit,'(2x,a)') '| '
      write(use_unit,'(2x,a)') '| Starting TDHF excited state calculation for Singlets.'
      write(use_unit,'(2x,a)') '| '
    endif
    write(use_unit,'(2x,a,i8.1,a,i8.1,a,i6.1,a,i6.1,a,i5.1,a)') '| Deviding the ', &
      n_pair_states, ' x ', n_pair_states, ' Matrix into ', &
      num_pcol,' colomns and ', num_prow, ' rows on ', n_tasks, ' processors.'
      tmp_real = 8 * n_basbas * n_states * n_states
    write(use_unit,'(2x,a,f10.3,a)') '| Memory load of ', tmp_real/2**30, ' GB devided among all Processors.'
  endif
  tmp_real = 8 * n_loc_prodbas * n_states * n_states
  write(info_str,'(2x,a,i5.1,a,i8.1,a,i8.1,a,f12.3,a)') '| Processor ', myid, &
      ' gets ', na_rows, ' rows and ', na_cols,' columns      (',tmp_real/2**20,' MB).'
  call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)
  if(myid==0) write(use_unit,'(2x,a)') '|'
  call descinit( casida_desc, n_pair_states, n_pair_states, &
	num_blocks, num_blocks, 0, 0, my_casida_ctxt, na_rows, errnum )

!map memory layout for communication in the main casida loop
  allocate(who_has(n_pair_states,n_pair_states),stat=errnum)
  call check_allocation(errnum,'who_has matrix alloc')
  who_has = 0

  do i_state = 1, occ_states
    do a_state = occ_states+1, n_states
  do j_state = 1, occ_states
    do b_state = occ_states+1, n_states

  call blacs_gridinfo(casida_desc(2), nprow, npcol, myrow, mycol)

  call infog2l(mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc, nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)

  if(myrow.eq.iarow .and. mycol.eq.iacol) then
    who_has(mat_idx(i_state,a_state), mat_idx(j_state,b_state)) = myid
  endif

  enddo; enddo; enddo; enddo

  call sync_integer_vector(who_has, n_pair_states**2)

!  if(myid==0) write(use_unit,*) 'Building ovlp_idx'

! build array with which proc has which ovlp_3KS element
  allocate(ovlp_idx(occ_states,occ_states+1:n_states),stat=errnum)
  call check_allocation(errnum,'casida ovlp_idx')
  ovlp_idx = -1

  allocate(tmp_idx(occ_states, occ_states+1:n_states),stat=errnum)
  call check_allocation(errnum,'tmp idx ovlp')
  tmp_idx = -1

  do i_state=1, occ_states
    if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
      do a_state=occ_states+1, n_states
	if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle

  ovlp_idx(i_state,a_state) = myid

  enddo; enddo

  tmp_int = size(ovlp_idx(1,:))*size(ovlp_idx(:,1))
  
  if(myid==0) then
  do i_proc=1, n_tasks-1
    call mpi_recv(tmp_idx, tmp_int, mpi_integer, i_proc, 100, mpi_comm_world, mpi_status, errnum)
    do i_state=1, occ_states
      do a_state=occ_states+1, n_states
	  if(tmp_idx(i_state,a_state).ne.-1) then !call aims_stop('double vector entry in ia index detected!')
	    ovlp_idx(i_state,a_state) = tmp_idx(i_state,a_state)
	  endif
    enddo; enddo
  enddo
  else
    call mpi_send(ovlp_idx, tmp_int, mpi_integer, 0, 100, mpi_comm_world, errnum)
  endif

  call mpi_bcast(ovlp_idx, tmp_int, mpi_integer, 0, mpi_comm_world, errnum)

  deallocate(tmp_idx)

  allocate(ovij_idx(occ_states,occ_states),stat=errnum)
  call check_allocation(errnum,'casida ovlp_idx')
  ovij_idx = -1

  do i_state=1, occ_states
    if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
      do j_state=1, occ_states
	if(own_dim2_o3KS(j_state) /= myp2_o3KS) cycle

  ovij_idx(i_state,j_state) = myid

  enddo; enddo

  tmp_int = size(ovij_idx(1,:))*size(ovij_idx(:,1))

  allocate(tmp_idx(occ_states, occ_states),stat=errnum)
  call check_allocation(errnum,'tmp idx ovlp')
  tmp_idx = -1

  if(myid==0) then
  do i_proc=1, n_tasks-1
    call mpi_recv(tmp_idx, tmp_int, mpi_integer, i_proc, 200, mpi_comm_world, mpi_status, errnum)
    do i_state=1, occ_states
      do j_state=1, occ_states
	  if(tmp_idx(i_state,j_state).ne.-1) then
	    ovij_idx(i_state,j_state) = tmp_idx(i_state,j_state)
	  endif
    enddo; enddo
  enddo
  else
    call mpi_send(ovij_idx, tmp_int, mpi_integer, 0, 200, mpi_comm_world, errnum)
  endif

  call mpi_bcast(ovij_idx, tmp_int, mpi_integer, 0, mpi_comm_world, errnum)

  deallocate(tmp_idx)

!  if(myid==0) write(use_unit,*) 'Building ovAB index'

  allocate(ovab_idx(occ_states+1:n_states,occ_states+1:n_states),stat=errnum)
  call check_allocation(errnum,'casida ovlp_idx')
  ovab_idx = -1

  do a_state=occ_states+1, n_states
    if(own_dim1_o3KS(a_state) /= myp1_o3KS) cycle
      do b_state=occ_states+1, n_states
	if(own_dim2_o3KS(b_state) /= myp2_o3KS) cycle

  ovab_idx(a_state,b_state) = myid

  enddo; enddo

  allocate(tmp_idx(occ_states+1:n_states, occ_states+1:n_states),stat=errnum)
  call check_allocation(errnum,'tmp idx ovlp')
  tmp_idx = -1

  tmp_int = size(ovab_idx(:,1))*size(ovab_idx(1,:))

  if(myid==0) then
  do i_proc=1, n_tasks-1
    call mpi_recv(tmp_idx, tmp_int, mpi_integer, i_proc, 300, mpi_comm_world, mpi_status, errnum)
    do a_state=occ_states+1, n_states
      do b_state=occ_states+1, n_states
	  if(tmp_idx(a_state,b_state).ne.-1) then !call aims_stop('double vector entry in ia index detected!')
	    ovab_idx(a_state,b_state) = tmp_idx(a_state,b_state)
	  endif
    enddo; enddo
  enddo
  else
    call mpi_send(ovab_idx, tmp_int, mpi_integer, 0, 300, mpi_comm_world, errnum)
  endif

  call mpi_bcast(ovab_idx, tmp_int, mpi_integer, 0, mpi_comm_world, errnum)

  deallocate(tmp_idx)

!  if(myid==0) write(use_unit,*) 'allocatable matrices for omega build'

  allocate(par_lmat(na_rows, na_cols),stat=errnum)
  call check_allocation(errnum,'casida parallel matrix alloc')
  par_lmat = 0d0

  allocate(par_mmat(na_rows, na_cols),stat=errnum)
  call check_allocation(errnum,'casida parallel matrix alloc')
  par_mmat = 0d0

  allocate(o3ks_buffer(n_basbas), stat=errnum)
  call check_allocation(errnum, 'casida o3ks_buffer')
  o3ks_buffer = 0d0

  allocate(par_lmm(na_rows,na_cols))
  call check_allocation(errnum, 'casida par_lmm')

  allocate(par_lpm(na_rows,na_cols))
  call check_allocation(errnum, 'casida par_lpm')

  if(excited_mode_triplet) then
    triplet_factor = 0d0
  else
    triplet_factor = 2.0d0
  endif

  call get_timestamps(time_omegabuild, clock_time_omegabuild)

  if(myid==0) write(use_unit,'(2x,a)') '| Building omega matrix  < ia | jb >   now.'

  kernel = 0d0

  do i_state=1, occ_states
    do a_state=occ_states+1, n_states

      if(myid==0.and.(mod(mat_idx(i_state,a_state),n_pair_states/10)==0)) then
	write(use_unit,'(2x,a,i8.1,a,f10.2,a)') '|  < ia |  pair state ', mat_idx(i_state,a_state), '  ', &
	dble(mat_idx(i_state,a_state))/dble(n_pair_states)*100d0, '%'
      endif

	do j_state=1, occ_states
	  do b_state=occ_states+1, n_states

  comm_tag = ( a_state + i_state ) * ( b_state + j_state )

!	##################	Calculate (ai|jb) for L 	################
  if( ovlp_idx(i_state,a_state) .eq. myid ) then
    if( ovlp_idx(j_state,b_state) .eq. myid ) then	! proc has all local here to calculate (ia|jb)
      i_state_loc = loc_dim1_o3KS(i_state)
      a_state_loc = loc_dim2_o3KS(a_state)
      j_state_loc = loc_dim1_o3KS(j_state)
      b_state_loc = loc_dim2_o3KS(b_state)
      kernel = 0d0
      do k_cnt = 1, n_basbas
	kernel = kernel + ovlp_3KS(k_cnt,i_state_loc,a_state_loc,n_spin) * ovlp_3KS(k_cnt,j_state_loc,b_state_loc,n_spin)
      enddo
    else							! proc has only (ia| here, needs to get |jb)
      call mpi_recv(o3ks_buffer, n_basbas, mpi_double_precision, &
                    ovlp_idx(j_state,b_state), comm_tag, mpi_comm_world, &
                    mpi_status, errnum)
      i_state_loc = loc_dim1_o3KS(i_state)
      a_state_loc = loc_dim2_o3KS(a_state)
      kernel = 0d0
      do k_cnt = 1, n_basbas
	kernel = kernel + ovlp_3KS(k_cnt,i_state_loc,a_state_loc,n_spin) * o3ks_buffer(k_cnt)
      enddo
    endif
  elseif(myid.eq.ovlp_idx(j_state,b_state)) then
    j_state_loc = loc_dim1_o3KS(j_state)
    b_state_loc = loc_dim2_o3KS(b_state)
    o3ks_buffer(:) = ovlp_3ks(:,j_state_loc,b_state_loc,n_spin)
    call mpi_send(o3ks_buffer, n_basbas, mpi_double_precision, ovlp_idx(i_state,a_state), comm_tag, mpi_comm_world, errnum)
  endif

  if(myid.eq.who_has(mat_idx(i_state,a_state),mat_idx(j_state,b_state))) then
    if(myid.eq.ovlp_idx(i_state,a_state)) then
      call get_local( tmp_real, par_lmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc )
      tmp_real = tmp_real + ( kernel * triplet_factor )
      call set_local( par_lmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc, tmp_real )
    else
      call mpi_recv(kernel, 1, mpi_double_precision, ovlp_idx(i_state,a_state), comm_tag, mpi_comm_world, mpi_status, errnum)
      call get_local( tmp_real, par_lmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc )
      tmp_real = tmp_real + ( kernel * triplet_factor )
      call set_local( par_lmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc, tmp_real )
    endif
  elseif(myid.eq.ovlp_idx(i_state,a_state)) then
    call mpi_send(kernel, 1, mpi_double_precision, &
                  who_has(mat_idx(i_state,a_state),mat_idx(j_state,b_state)), &
                  comm_tag, mpi_comm_world, errnum)
  endif


!	##################	Calculate (ia|jb) for M 	################
  if( ovlp_idx(i_state,a_state) .eq. myid ) then
    if( ovlp_idx(j_state,b_state) .eq. myid ) then	! proc has all local here to calculate (ia|jb)
      i_state_loc = loc_dim1_o3KS(i_state)
      a_state_loc = loc_dim2_o3KS(a_state)
      j_state_loc = loc_dim1_o3KS(j_state)
      b_state_loc = loc_dim2_o3KS(b_state)
      kernel = 0d0
      do k_cnt = 1, n_basbas
	kernel = kernel + ovlp_3KS(k_cnt,i_state_loc,a_state_loc,n_spin) * ovlp_3KS(k_cnt,j_state_loc,b_state_loc,n_spin)
      enddo
    else							! proc has only (ia| here, needs to get |jb)
      call mpi_recv(o3ks_buffer, n_basbas, mpi_double_precision, &
                    ovlp_idx(j_state,b_state), comm_tag, mpi_comm_world, &
                    mpi_status, errnum)
      i_state_loc = loc_dim1_o3KS(i_state)
      a_state_loc = loc_dim2_o3KS(a_state)
      kernel = 0d0
      do k_cnt = 1, n_basbas
	kernel = kernel + ovlp_3KS(k_cnt,i_state_loc,a_state_loc,n_spin) * o3ks_buffer(k_cnt)
      enddo
    endif
  elseif(myid.eq.ovlp_idx(j_state,b_state)) then
    j_state_loc = loc_dim1_o3KS(j_state)
    b_state_loc = loc_dim2_o3KS(b_state)
    o3ks_buffer(:) = ovlp_3ks(:,j_state_loc,b_state_loc,n_spin)
    call mpi_send(o3ks_buffer, n_basbas, mpi_double_precision, ovlp_idx(i_state,a_state), comm_tag, mpi_comm_world, errnum)
  endif

  if(myid.eq.who_has(mat_idx(i_state,a_state),mat_idx(j_state,b_state))) then
    if(myid.eq.ovlp_idx(i_state,a_state)) then
      call get_local( tmp_real, par_mmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc )
      tmp_real = tmp_real + ( kernel * triplet_factor )
      call set_local( par_mmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc, tmp_real )
    else
      call mpi_recv(kernel, 1, mpi_double_precision, ovlp_idx(i_state,a_state), comm_tag, mpi_comm_world, mpi_status, errnum)
      call get_local( tmp_real, par_mmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc )
      tmp_real = tmp_real + ( kernel * triplet_factor )
      call set_local( par_mmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc, tmp_real )
    endif
  elseif(myid.eq.ovlp_idx(i_state,a_state)) then
    call mpi_send(kernel, 1, mpi_double_precision, &
                  who_has(mat_idx(i_state,a_state),mat_idx(j_state,b_state)), &
                  comm_tag, mpi_comm_world, errnum)
  endif


!	##################	Calculate (ib|ja) for B correction	################
  if( ovlp_idx(i_state,b_state) .eq. myid ) then
    if( ovlp_idx(j_state,a_state) .eq. myid ) then	! proc has all local here, calculate (ib|ja) here
      i_state_loc = loc_dim1_o3KS(i_state)
      b_state_loc = loc_dim2_o3KS(b_state)
      j_state_loc = loc_dim1_o3KS(j_state)
      a_state_loc = loc_dim2_o3KS(a_state)
      kernel = 0d0
      do k_cnt = 1, n_basbas
	kernel = kernel + ovlp_3KS(k_cnt,i_state_loc,b_state_loc,n_spin) * ovlp_3KS(k_cnt,j_state_loc,a_state_loc,n_spin)
      enddo
    else							! proc has only (ib| here, needs to get |ja)
      call mpi_recv(o3ks_buffer, n_basbas, mpi_double_precision, &
                    ovlp_idx(j_state,a_state), comm_tag, mpi_comm_world, &
                    mpi_status, errnum)
      i_state_loc = loc_dim1_o3KS(i_state)
      b_state_loc = loc_dim2_o3KS(b_state)
      kernel = 0d0
      do k_cnt = 1, n_basbas
	kernel = kernel + ovlp_3KS(k_cnt,i_state_loc,b_state_loc,n_spin) * o3ks_buffer(k_cnt)
      enddo
    endif
  elseif(myid.eq.ovlp_idx(j_state,a_state)) then
    j_state_loc = loc_dim1_o3KS(j_state)
    a_state_loc = loc_dim2_o3KS(a_state)
    o3ks_buffer(:) = ovlp_3ks(:,j_state_loc,a_state_loc,n_spin)
    call mpi_send(o3ks_buffer, n_basbas, mpi_double_precision, ovlp_idx(i_state,b_state), comm_tag, mpi_comm_world, errnum)
  endif

  if(myid.eq.who_has(mat_idx(i_state,a_state),mat_idx(j_state,b_state))) then
    if(myid.eq.ovlp_idx(i_state,b_state)) then
      call get_local( tmp_real, par_mmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc )
      tmp_real = tmp_real - kernel
      call set_local( par_mmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc, tmp_real )
    else
      call mpi_recv(kernel, 1, mpi_double_precision, ovlp_idx(i_state,b_state), comm_tag, mpi_comm_world, mpi_status, errnum)
      call get_local( tmp_real, par_mmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc )
      tmp_real = tmp_real - kernel
      call set_local(par_mmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc, tmp_real )
    endif
  elseif(myid.eq.ovlp_idx(i_state,b_state)) then
    call mpi_send(kernel, 1, mpi_double_precision, &
                  who_has(mat_idx(i_state,a_state),mat_idx(j_state,b_state)), &
                  comm_tag, mpi_comm_world, errnum)
  endif


!	##################	Calculate (ij|ab) for A correction	################
  if( ovij_idx(i_state,j_state) .eq. myid ) then
    if( ovab_idx(a_state,b_state) .eq. myid ) then	! proc has all local here, calculate (ij|ab) here
      i_state_loc = loc_dim1_o3KS(i_state)
      j_state_loc = loc_dim2_o3KS(j_state)
      a_state_loc = loc_dim1_o3KS(a_state)
      b_state_loc = loc_dim2_o3KS(b_state)
      kernel = 0d0
      do k_cnt = 1, n_basbas
	kernel = kernel + ovlp_3KS(k_cnt,i_state_loc,j_state_loc,n_spin) * ovlp_3KS(k_cnt,a_state_loc,b_state_loc,n_spin)
      enddo
    else							! proc has only (ij| here, needs to get |ab)
      call mpi_recv(o3ks_buffer, n_basbas, mpi_double_precision, &
                    ovab_idx(a_state,b_state), comm_tag, mpi_comm_world, &
                    mpi_status, errnum)
      i_state_loc = loc_dim1_o3KS(i_state)
      j_state_loc = loc_dim2_o3KS(j_state)
      kernel = 0d0
      do k_cnt = 1, n_basbas
	kernel = kernel + ovlp_3KS(k_cnt,i_state_loc,j_state_loc,n_spin) * o3ks_buffer(k_cnt)
      enddo
    endif
  elseif(myid.eq.ovab_idx(a_state,b_state)) then
    a_state_loc = loc_dim1_o3KS(a_state)
    b_state_loc = loc_dim2_o3KS(b_state)
    o3ks_buffer(:) = ovlp_3ks(:,a_state_loc,b_state_loc,n_spin)
    call mpi_send(o3ks_buffer, n_basbas, mpi_double_precision, ovij_idx(i_state,j_state), comm_tag, mpi_comm_world, errnum)
  endif

  if(myid.eq.who_has(mat_idx(i_state,a_state),mat_idx(j_state,b_state))) then
    if(myid.eq.ovij_idx(i_state,j_state)) then
      call get_local( tmp_real, par_lmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc )
      tmp_real = tmp_real - kernel
      call set_local( par_lmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc, tmp_real )
    else
      call mpi_recv(kernel, 1, mpi_double_precision, ovij_idx(i_state,j_state), comm_tag, mpi_comm_world, mpi_status, errnum)
      call get_local( tmp_real, par_lmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc )
      tmp_real = tmp_real - kernel
      call set_local( par_lmat, mat_idx(i_state,a_state), mat_idx(j_state,b_state), casida_desc, tmp_real )
    endif
  elseif(myid.eq.ovij_idx(i_state,j_state)) then
    call mpi_send(kernel, 1, mpi_double_precision, &
                  who_has(mat_idx(i_state,a_state), &
                  mat_idx(j_state,b_state)), comm_tag, mpi_comm_world, errnum)
  endif


  enddo; enddo; enddo; enddo

  deallocate(o3ks_buffer)

!- add diagonal (e_a-e_i) contribution to L
  do i_state=1, occ_states
    do a_state=occ_states+1, n_states

  call pdelget( 'a', 'd', tmp_real, par_lmat, mat_idx(i_state,a_state), mat_idx(i_state,a_state), casida_desc )
  tmp_real = tmp_real + ( KS_eigenvalue(a_state,n_spin,n_k_points) - KS_eigenvalue(i_state,n_spin,n_k_points) )
  call pdelset( par_lmat, mat_idx(i_state,a_state), mat_idx(i_state,a_state), casida_desc, tmp_real )

  enddo; enddo

  call get_times(time_omegabuild, clock_time_omegabuild)

  par_lmm = par_lmat - par_mmat
  par_lpm = par_lmat + par_mmat

  deallocate(par_lmat)
  deallocate(par_mmat)

  allocate(par_omega(na_rows,na_cols),stat=errnum)
  call check_allocation(errnum,'casida parallel matrix alloc')
  par_omega = 0d0

  call pdgemm('n', 'n', n_pair_states, n_pair_states, n_pair_states, 1.0d0, &
	       par_lmm, 1, 1, casida_desc, &
	       par_lpm, 1, 1, casida_desc, &
	       0d0, &
	       par_omega, 1, 1, casida_desc)

!- free mem
  deallocate(par_lpm)
  deallocate(par_lmm)

  allocate(par_evecs(na_rows, na_cols),stat=errnum)
  call check_allocation(errnum,'casida parallel matrix alloc')
  par_evecs = 0d0
  allocate(par_evals(n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida parallel matrix alloc')
  par_evals = 0d0

  allocate(lapack_scale(n_pair_states),stat=errnum)
  call check_allocation(errnum,'scale array')

!  call pdgebal( 'B', n_pair_states, par_omega, casida_desc, lapack_ilo, lapack_ihi, lapack_scale, errnum)

!  if(myid==0) write(use_unit,*) 'bal errnum: ', errnum

  allocate(tau(1+n_pair_states-2),stat=errnum)
  call check_allocation(errnum,'tau solver')

  allocate(lapack_workspace(1))

!  call pdgehrd( n_pair_states, lapack_ilo, lapack_ihi, par_omega, 1, 1, casida_desc, tau, lapack_workspace, -1, errnum)

  tmp_int = int(lapack_workspace(1))

!if(myid==0) write(use_unit,*) 'optimal reduction workspace size: ', lapack_workspace(1), tmp_int

  deallocate(lapack_workspace)
  allocate(lapack_workspace(tmp_int),stat=errnum)
  call check_allocation(errnum,'workspace')

!  call pdgehrd( n_pair_states, 1, n_pair_states, par_omega, 1, 1, casida_desc, tau, lapack_workspace, tmp_int, errnum)

  deallocate(lapack_workspace)

!if(myid==0) write(use_unit,*) 'reduction errnum: ', errnum

  allocate(wr(n_pair_states),stat=errnum)
  call check_allocation(errnum,'real evals')
  allocate(wi(n_pair_states),stat=errnum)
  call check_allocation(errnum,'imag evals')

  allocate(par_z(na_rows,na_cols),stat=errnum)
  call check_allocation(errnum,'par_z')
  par_z = 0d0

  allocate(lapack_workspace(1))
  allocate(lapack_iwork(1))

!  call pdhseqr( 'E', 'N', n_pair_states, lapack_ilo, lapack_ihi, par_omega, casida_desc, wr, wi, par_z, casida_desc, lapack_workspace, -1, lapack_iwork, 1, errnum)

!if(myid==0) write(use_unit,*) 'workspace call errnum: ', errnum
!if(myid==0) write(use_unit,*) 'optimal real space size: ', lapack_workspace(1)
!if(myid==0) write(use_unit,*) 'optimal int space size: ', lapack_iwork(1)


  tmp_int = int(lapack_workspace(1))
  deallocate(lapack_workspace)
  allocate(lapack_workspace(tmp_int),stat=errnum)

  tmp_int = int(lapack_iwork(1))
  deallocate(lapack_iwork)
  allocate(lapack_iwork(tmp_int),stat=errnum)

!  call pdhseqr( 'E', 'N', n_pair_states, lapack_ilo, lapack_ihi, par_omega, casida_desc, wr, wi, &
!		 par_z, casida_desc, lapack_workspace, size(lapack_workspace), lapack_iwork, size(lapack_iwork), errnum)

!if(myid==0) write(use_unit,*) 'run errnum: ', errnum

! OF=w^2F
  do i_state=1, n_pair_states
    if(par_evals(i_state).lt.0d0) then
      call neg_energy_error(par_evals(i_state))
      par_evals(i_state) = 0d0
    else
      par_evals(i_state) = dsqrt(par_evals(i_state))
    endif
  enddo

  allocate(osc_str(n_pair_states),stat=errnum)
  call check_allocation(errnum,'osc_str')

  allocate(tmom(n_pair_states,3),stat=errnum)
  call check_allocation(errnum,'tmom allocation')
  tmom = 0d0

  do i_state=1, n_pair_states
    do j_state=1, occ_states
      do b_state=occ_states+1, n_states

    call pdelget('a','d',tmp_real, par_evecs, mat_idx(j_state,b_state), i_state, casida_desc)

    tmom(i_state,1) = tmom(i_state,1) + (tmp_real * dipole_moment(j_state,b_state-occ_states,n_spin,1))
    tmom(i_state,2) = tmom(i_state,2) + (tmp_real * dipole_moment(j_state,b_state-occ_states,n_spin,2))
    tmom(i_state,3) = tmom(i_state,3) + (tmp_real * dipole_moment(j_state,b_state-occ_states,n_spin,3))
    enddo; enddo;

!   account for spin multiplicity
    tmom(i_state,:) = dsqrt(2d0) * tmom(i_state,:)

    osc_str(i_state) = 2d0/3d0 * par_evals(i_state) * ( tmom(i_state,1)**2 + tmom(i_state,2)**2 + tmom(i_state,3)**2 )

  enddo

  if(myid==0) then
    if(excited_mode_triplet) then
      write(use_unit,'(2x,a)') '| Writing Triplet excitation energies (TDHF) '
      write(info_str,'(a)') '. Triplet: '
    else
      write(use_unit,'(2x,a)') '| Writing Singlet excitation energies (TDHF) '
      write(info_str,'(a)') '. Singlet: '
    endif
    do i_state=1, num_excited_states
      write(use_unit,'(2x,a,i4.1,a11,f12.4,a,f8.4,a,f12.4,a,f8.4,a,f12.4)') &
         '| ', i_state, info_str, par_evals(i_state), ' Ha      [', &
         par_evals(i_state)*hartree, ' eV]', osc_str(i_state)
    enddo
    write(use_unit,'(2x,a)') '| '
  endif

  if(myid==0) then
   if(excited_mode_triplet) then
     write(info_str,'(a)') 'TDDFT_LR_Spectrum_Triplet.dat'
   else
     write(info_str,'(a)') 'TDDFT_LR_Spectrum_Singlet.dat'
   endif
   open(100,file=info_str)
   do i_state=1, n_pair_states
     write(100,'(f29.14,4x,f29.14)') par_evals(i_state)*hartree, osc_str(i_state)
   enddo
   close(100)
  endif

  ! final cleanup
  deallocate(par_evecs)
  deallocate(par_evals)
  deallocate(ovlp_idx)
  deallocate(who_has)
  deallocate(par_omega)

  if(myid==0.and.excited_mode_triplet) write(use_unit,'(2x,a)') '| End of TDHF Triplet excitation energy calculation.'
  if(myid==0.and..not.excited_mode_triplet) write(use_unit,'(2x,a)') '| End of TDHF Singlet excitation energy calculation.'

end subroutine casida_tdhf_2d

subroutine casida_tdhf()
  implicit none

  integer				:: i_state, j_state, a_state, b_state, &
					   k_cnt, l_cnt, errnum

  real*8				:: m_kernel, l_kernel, l_tdhf, m_tdhf, triplet_factor
  real*8, allocatable, dimension(:)	:: eig_vals, lapack_space_1d, osc_str, wr, wi
  real*8, allocatable, dimension(:,:)	:: a_mat, b_mat, apb, amb, omega, vl, vr, tmom

  intrinsic	:: dsqrt

  if(excited_mode_triplet) then
    triplet_factor = 0d0
  else
    triplet_factor = 2d0
  endif

  allocate(a_mat(n_pair_states,n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida A Matrix')
  a_mat = 0d0

  allocate(b_mat(n_pair_states,n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida B Matrix')
  b_mat = 0d0

  call get_timestamps(time_omegabuild, clock_time_omegabuild)

  do i_state = 1, occ_states
    do a_state = occ_states+1, n_states
      do j_state=1, occ_states
	do b_state=occ_states+1, n_states

    l_kernel = 0d0
    m_kernel = 0d0
    l_tdhf = 0d0
    m_tdhf = 0d0
 
    do k_cnt = 1, n_basbas
      l_kernel = l_kernel + ovlp_3KS(k_cnt,i_state,a_state,n_spin) * ovlp_3KS(k_cnt,j_state,b_state,n_spin)
      m_kernel = m_kernel + ovlp_3KS(k_cnt,i_state,a_state,n_spin) * ovlp_3KS(k_cnt,j_state,b_state,n_spin)
      l_tdhf = l_tdhf + ovlp_3KS(k_cnt,i_state,j_state,n_spin) * ovlp_3KS(k_cnt,a_state,b_state,n_spin)
      m_tdhf = m_tdhf + ovlp_3KS(k_cnt,i_state,b_state,n_spin) * ovlp_3KS(k_cnt,j_state,a_state,n_spin)
!    do l_cnt = 1, n_basbas
!      l_kernel = l_kernel + fxc_3ks(k_cnt,i_state,a_state,n_spin) * fxc_matrix(k_cnt,l_cnt) * fxc_3ks(l_cnt,j_state,b_state,n_spin)
!      m_kernel = m_kernel + fxc_3ks(k_cnt,i_state,a_state,n_spin) * fxc_matrix(k_cnt,l_cnt) * fxc_3ks(l_cnt,j_state,b_state,n_spin)
     !enddo
    enddo
    a_mat(mat_idx(i_state,a_state),mat_idx(j_state,b_state)) = triplet_factor * l_kernel - l_tdhf
    b_mat(mat_idx(i_state,a_state),mat_idx(j_state,b_state)) = triplet_factor * m_kernel - m_tdhf

  enddo; enddo; enddo; enddo

  do i_state = 1, occ_states
    do a_state = occ_states+1, n_states

  a_mat(mat_idx(i_state,a_state),mat_idx(i_state,a_state)) = &
    a_mat(mat_idx(i_state,a_state),mat_idx(i_state,a_state)) + &
    ( KS_eigenvalue(a_state,n_spin,n_k_points) - KS_eigenvalue(i_state,n_spin,n_k_points) )

  enddo; enddo

  call get_times(time_omegabuild, clock_time_omegabuild)

  allocate(amb(n_pair_states,n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida L Matrix')
  amb = 0d0

  allocate(apb(n_pair_states,n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida L Matrix')
  apb = 0d0

  apb = a_mat + b_mat
  amb = a_mat - b_mat

  deallocate(a_mat)
  deallocate(b_mat)

  allocate(omega(n_pair_states,n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida L Matrix')
  omega = 0d0

  call dgemm( 'N' , 'N' , n_pair_states , n_pair_states , n_pair_states , &
	     1.0d0 , amb , n_pair_states , apb , n_pair_states , 0d0 , omega , n_pair_states )

  deallocate(apb)
  deallocate(amb)

!- allocate lapack spaces and call eigensolver
  allocate(lapack_space_1d(8 * n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_workspace')
  lapack_space_1d = 0d0
  allocate(wr(n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_eig_vals')
  wr = 0d0
  allocate(wi(n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_eig_vals')
  wi = 0d0
  allocate(vl(n_pair_states,n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_eig_vals')
  vl = 0d0
  allocate(vr(n_pair_states,n_pair_states),stat=errnum)
  call check_allocation(errnum,'casida lapack_eig_vals')
  vr = 0d0

  call dgeev( 'N' , 'V' , n_pair_states, omega, n_pair_states, wr, wi, vl, n_pair_states, &
	      vr, n_pair_states, lapack_space_1d, 8*n_pair_states, errnum )

  if(errnum.ne.0) then
    write(info_str,*) 'DGEEV reportet Error: ',errnum
  endif

  deallocate(lapack_space_1d)
  deallocate(omega)

!- Omega F = omega^2 F
  do i_state=1, n_pair_states
    if(wr(i_state).lt.0d0) then
      call neg_energy_error(wr(i_state))
      wr(i_state) = 0d0
    else
      wr(i_state) = dsqrt(wr(i_state))
    endif
  enddo

  allocate(tmom(n_pair_states,3),stat=errnum)
  call check_allocation(errnum,'tmom allocation')
  tmom = 0d0

  do i_state=1, n_pair_states
    do j_state=1, occ_states
      do b_state=occ_states+1, n_states

	tmom(i_state,1) = tmom(i_state,1) + (vr(mat_idx(j_state,b_state),i_state) * dipole_moment(j_state,b_state-occ_states,n_spin,1))
	tmom(i_state,2) = tmom(i_state,2) + (vr(mat_idx(j_state,b_state),i_state) * dipole_moment(j_state,b_state-occ_states,n_spin,2))
	tmom(i_state,3) = tmom(i_state,3) + (vr(mat_idx(j_state,b_state),i_state) * dipole_moment(j_state,b_state-occ_states,n_spin,3))

  enddo; enddo; enddo

  tmom = sqrt(2.0d0) * tmom

  allocate(osc_str(n_pair_states),stat=errnum)
  call check_allocation(errnum,'osc_str')
  osc_str = 0d0

  do i_state=1, n_pair_states
    osc_str(i_state) = 2.0d0/3.0d0 * wr(i_state) * ( tmom(i_state,1)**2 + tmom(i_state,2)**2 + tmom(i_state,3)**2 )
  enddo

  do j_state = 1, n_pair_states
    do i_state = 1, n_pair_states
    if(wr(j_state).ge.wr(i_state)) then
      l_kernel = wr(i_state)
      wr(i_state) = wr(j_state)
      wr(j_state) = l_kernel
    endif
    enddo
  enddo

  allocate(eig_vals(n_pair_states))

  eig_vals = wr

  wr = 0d0

  do i_state = 1, n_pair_states
  wr(i_state) = eig_vals(n_pair_states-i_state+1)
  enddo

  if(myid==0) then
    if(excited_mode_triplet) then
      write(use_unit,'(2x,a)') 'Writing Triplet excitation energies (TDHF): '
      write(info_str,'(a)') '. Triplet: '
    else
      write(use_unit,'(2x,a)') 'Writing Singlet excitation energies (TDHF): '
      write(info_str,'(a)') '. Singlet: '
    endif

    do i_state=1, num_excited_states
      write(use_unit,'(2x,i4.1,a11,f12.4,a,f8.4,a,f12.4,a,f8.4,a,f12.4)') i_state, info_str, wr(i_state)*hartree, ' eV      [', wr(i_state), ' Ha]', osc_str(i_state)
    enddo
    write(use_unit,*) ' '
  endif

  if(excited_mode_triplet) then
    write(info_str,'(a)') 'TDDFT_LR_Spectrum_Triplet.dat'
  else
    write(info_str,'(a)') 'TDDFT_LR_Spectrum_Singlet.dat'
  endif
  open(100,file=info_str)
  do i_state=1, n_pair_states
    write(100,'(f29.14,4x,f29.14)') wr(i_state)*hartree, osc_str(i_state)
  enddo
  close(100)

!- cleanup the rest
  deallocate(osc_str)
  deallocate(wr,wi)
  deallocate(vr,vl)

  if(excited_mode_triplet) write(use_unit,'(2x,a)') '| End of TDHF Triplet excitation energy calculation.'
  if(.not.excited_mode_triplet) write(use_unit,'(2x,a)') '| End of TDHF Singlet excitation energy calculation.'

end subroutine casida_tdhf

		!---------------------------------------------
		!--------- END main casida routines ----------
		!---------------------------------------------

pure integer function mat_idx(i,a) result(idx)
 implicit none
 integer, intent(in) :: i, a
 idx = (i-low_limit)*(unocc_states-unocc_limit)+(a-occ_states)
end function mat_idx

subroutine get_ks_orbitals_1d()
  implicit none
  integer				:: errnum
  real*8				:: n_bytes_o3KS
  real*8, allocatable, dimension(:,:)	:: ovlp_3fn

  allocate(ovlp_3fn(n_basis_pairs,n_basbas),stat=errnum)
  call check_allocation(errnum,'ovlp_3fxc Casida')
  
  call get_coeff_3fn_v_1d(ovlp_3fn)

  call localorb_info('|', use_unit, '(2x,A)', OL_norm)
  n_bytes_o3KS = 8d0 * n_loc_prodbas * n_states**2 * n_spin
  write(info_str, "(2X,A,I7,A,I7,A,F12.3,A)") &
  & '| The ovlp_3KS matrix takes another', nint(n_bytes_o3KS / 2.d0**20), &
  & ' MB x', n_tasks, ' procs =', n_tasks * n_bytes_o3KS / 2.d0**30, ' GB.'
  call localorb_info(info_str, use_unit, '(A)', OL_norm)
  call localorb_info('|', use_unit, '(2x,A)', OL_norm)

  if (sparse_o3fn) then
  ! coeff_3fn_ten, coulomb_matr_lvl -> ovlp_3KS
    call ovlp3KS_lvl_1d(n_states, KS_eigenvector, &
                        coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
  call aims_stop('no sparse implementation for casida yat!')
  else
  ! coeff_3fn -> ovlp_3KS
    call transform_ovlp3fn(n_states, ks_eigenvector, ovlp_3fn, ovlp_3ks)
  end if

  deallocate(ovlp_3fn)
end subroutine get_ks_orbitals_1d

subroutine get_ks_orbitals_2d()
  implicit none
  integer				:: errnum

  real*8				:: n_bytes_o3KS
  real*8, allocatable, dimension(:,:)	:: ovlp_3fn

  allocate(ovlp_3fn(n_basis_pairs,n_loc_prodbas), stat=errnum)
  call check_allocation(errnum,'ovlp_3fn Casida')

  call get_coeff_3fn_v_2d(ovlp_3fn)

  call get_timestamps(time_ovlp3fn_to_ovlp3KS, clock_time_ovlp3fn_to_ovlp3KS)

  ndim1_o3KS = (n_states-1)/np1_o3KS + 1
  ndim2_o3KS = (n_states-1)/np2_o3KS + 1

  call localorb_info('|', use_unit, '(2x,A)', OL_norm)
  n_bytes_o3KS = 8d0 * n_basbas * ndim1_o3KS * ndim2_o3KS * n_spin
  write(info_str, "(2X,A,I7,A,I7,A,F12.3,A)") &
  & '| The ovlp_3KS matrix takes another', nint(n_bytes_o3KS / 2.d0**20), &
  & ' MB x', n_tasks, ' procs =', n_tasks * n_bytes_o3KS / 2.d0**30, ' GB.'
  call localorb_info(info_str, use_unit, '(A)', OL_norm)
  call localorb_info('|', use_unit, '(2x,A)', OL_norm)

  if (sparse_o3fn) then
  ! coeff_3fn_ten, coulomb_matr_lvl -> ovlp_3KS
    call ovlp3KS_lvl_2d(n_states, KS_eigenvector, &
                      coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
  call aims_stop('No sparse in casida routines yet')
  else
    call transform_ovlp3fn_2(n_states,KS_eigenvector,ovlp_3fn,ovlp_3ks)
  end if

  deallocate(ovlp_3fn)
end subroutine get_ks_orbitals_2d

subroutine get_tddft_orbitals_1d()
  implicit none
  integer				:: k_cnt, l_cnt, errnum

  real*8				:: n_bytes_o3KS
  real*8, allocatable, dimension(:,:)	:: ovlp_3fn, fxc_3fn, fxc_matrix

  call localorb_info('|', use_unit, '(2x,A)', OL_norm)
  n_bytes_o3KS = 2d0 * 8d0 * n_loc_prodbas * n_states**2 * n_spin
  write(info_str, "(2X,A,I7,A,I7,A,F12.3,A)") &
  & '| The ovlp_3KS and fxc_3KS matrices take another', nint(n_bytes_o3KS / 2.d0**20), &
  & ' MB x', n_tasks, ' procs =', n_tasks * n_bytes_o3KS / 2.d0**30, ' GB.'
  call localorb_info(info_str, use_unit, '(A)', OL_norm)
  call localorb_info('|', use_unit, '(2x,A)', OL_norm)

  n_bytes_o3KS = 8d0 * n_basbas**2
  write(info_str, "(2X,A,I7,A,I7,A,F12.3,A)") &
  & '| The f_xc Matrix takes another', nint(n_bytes_o3KS / 2.d0**20), &
  & ' MB x', n_tasks, ' procs =', n_tasks * n_bytes_o3KS / 2.d0**30, ' GB.'
  call localorb_info(info_str, use_unit, '(A)', OL_norm)
  call localorb_info('|', use_unit, '(2x,A)', OL_norm)

  allocate(ovlp_3fn(n_basis_pairs,n_basbas),stat=errnum)
  call check_allocation(errnum,'ovlp_3fxc Casida')
  allocate(fxc_3fn(n_basis_pairs, n_basbas),stat=errnum)
  call check_allocation(errnum,'fxc_3fxc Casida')
  allocate(fxc_matrix(n_basbas, n_basbas),stat=errnum)
  call check_allocation(errnum,'fxc_3fxc Casida')

  call get_timestamps(time_fxc_integr, clock_time_fxc_integr)
  call integrate_fxc_matrix(partition_tab, rho, rho_gradient, l_shell_max, fxc_matrix)
  call get_times(time_fxc_integr, clock_time_fxc_integr)

  call get_coeff_fxc3fn_v_1d(fxc_3fn, ovlp_3fn)

  if (sparse_o3fn) then
  ! coeff_3fn_ten, coulomb_matr_lvl -> ovlp_3KS
    call ovlp3KS_lvl_1d(n_states, KS_eigenvector, &
                        coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
  call aims_stop('no sparse implementation for casida yat!')
  else
  ! coeff_3fn -> ovlp_3KS
    call transform_ovlp3fn(n_states, ks_eigenvector, ovlp_3fn, ovlp_3ks)
    call transform_ovlp3fn(n_states, ks_eigenvector, fxc_3fn, fxc_3ks)
  end if

  call localorb_info('')
  write(info_str,'(2x,a)') '| Contracting intermediate matrix x matrix product:'
  call localorb_info(info_str)
  write(info_str,'(2x,a)') '|     < u | f_xc | v > x < v | jb >'
  call localorb_info(info_str)
  write(info_str,'(2x,a)') '| in auxiliary basis.'
  call localorb_info(info_str)
  call get_timestamps(time_matvec_cont, clock_time_matvec_cont)
  pure_3ks = fxc_3ks
  do k_cnt=1, n_states
    do l_cnt=1, n_states
      call dgemv('n', n_basbas, n_basbas, 1d0, fxc_matrix, n_basbas, &
                 pure_3ks(:,k_cnt,l_cnt,n_spin), 1, 0d0, &
                 fxc_3ks(:,k_cnt,l_cnt,n_spin), 1)
  enddo; enddo
 call get_times(time_matvec_cont, clock_time_matvec_cont)
  write(info_str,'(2x,a)') '|'
  call localorb_info(info_str) 

  deallocate(ovlp_3fn)
  deallocate(fxc_3fn)
  deallocate(fxc_matrix)
end subroutine get_tddft_orbitals_1d

subroutine get_tddft_orbitals_2d()
  implicit none
  integer       :: k_cnt, l_cnt, errnum, i_spin

  character*128 :: info_str

  real*8        :: n_bytes_o3KS

  ! Expansion coefficients for deriving the Coulomb matrix
  real*8, allocatable, dimension(:,:) :: ovlp_3fn
  
  ! Expansion coefficients for deriving the f_xc matrix
  real*8, allocatable, dimension(:,:) :: fxc_3fn
  
  ! fxc_matrix will contain <prod_wave|f_xc|prod_wave> (n_basbas,n_basbas,n_spin)
  real*8, allocatable, dimension(:,:,:) :: fxc_matrix

  call get_timestamps(time_neutral_excitation_orbitals, &
           clock_time_neutral_excitation_orbitals)

  n_bytes_o3KS = occ_states*unocc_states*n_spin*3
  write(info_str, "(2X,A,F8.1,A,I7,A,F12.3,A)") &
  & '| The dipole moments already took ', dble(n_bytes_o3KS / 2.d0**20), &
  & ' MB x', n_tasks, ' procs =', n_tasks * n_bytes_o3KS / 2.d0**30, ' GB.'
  call localorb_info(info_str, use_unit, '(A)', OL_norm)

  call localorb_info('|', use_unit, '(2x,A)', OL_norm)
  n_bytes_o3KS = 2d0 * 8d0 * n_basbas * ndim1_o3KS * ndim2_o3KS * n_spin
  write(info_str, "(2X,A,I7,A,I7,A,F12.3,A)") &
  & '| The ovlp_3KS and fxc_3KS matrices take another', nint(n_bytes_o3KS / 2.d0**20), &
  & ' MB x', n_tasks, ' procs =', n_tasks * n_bytes_o3KS / 2.d0**30, ' GB.'
  call localorb_info(info_str, use_unit, '(A)', OL_norm)
  call localorb_info('|', use_unit, '(2x,A)', OL_norm)

  n_bytes_o3KS = 8d0 * n_basbas**2
  write(info_str, "(2X,A,I7,A,I7,A,F12.3,A)") &
  & '| The f_xc Matrix takes another', nint(n_bytes_o3KS / 2.d0**20), &
  & ' MB x', n_tasks, ' procs =', n_tasks * n_bytes_o3KS / 2.d0**30, ' GB.'
  call localorb_info(info_str, use_unit, '(A)', OL_norm)
  call localorb_info('|', use_unit, '(2x,A)', OL_norm)

  allocate(ovlp_3fn(n_basis_pairs,n_loc_prodbas), stat=errnum)
  call check_allocation(errnum,'get_tddft_orbitals_2d')
  allocate(fxc_3fn(n_basis_pairs,n_loc_prodbas), stat=errnum)
  call check_allocation(errnum,'get_tddft_orbitals_2d')
  allocate(fxc_matrix(n_basbas,n_basbas,n_spin), stat=errnum)
  call check_allocation(errnum,'get_tddft_orbitals_2d')
  
  call get_timestamps(time_fxc_integr, clock_time_fxc_integr)
  ! Calculate <prod_wave|f_xc|prod_wave>
  call integrate_fxc_matrix(partition_tab, rho, rho_gradient, l_shell_max, fxc_matrix)
  call get_times(time_fxc_integr, clock_time_fxc_integr)

  ! Calculate the expansion coefficients for the coulomb and the f_xc matrix
  call get_timestamps(time_coeffs_fxc3fn, clock_time_coeffs_fxc3fn)
  call get_coeff_fxc3fn_v_2d(fxc_3fn, ovlp_3fn)
  call get_times(time_coeffs_fxc3fn, clock_time_coeffs_fxc3fn)

  ! Multiply KS_eigenvectors to transform to KS_eigenstate basis
  call get_timestamps(time_ovlp3fn_to_ovlp3KS, clock_time_ovlp3fn_to_ovlp3KS)
  if (sparse_o3fn) then
  ! coeff_3fn_ten, coulomb_matr_lvl -> ovlp_3KS
    call aims_stop('No sparse in casida routines yet')
    call ovlp3KS_lvl_2d(n_states, KS_eigenvector, &
                      coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
  else
    call transform_ovlp3fn_2(n_states,KS_eigenvector,ovlp_3fn,ovlp_3ks)
    call transform_ovlp3fn_2(n_states,KS_eigenvector,fxc_3fn,pure_3ks)
  end if

  deallocate(fxc_3fn)
  deallocate(ovlp_3fn)
  call get_times(time_ovlp3fn_to_ovlp3KS, clock_time_ovlp3fn_to_ovlp3KS)

  call localorb_info('')
  write(info_str,'(2x,a)') '| Contracting intermediate matrix x matrix product:'
  call localorb_info(info_str)
  write(info_str,'(2x,a)') '|     < u | f_xc | v > x < v | jb >'
  call localorb_info(info_str)
  write(info_str,'(2x,a)') '| in auxiliary basis.'
  call localorb_info(info_str)
  
  call get_timestamps(time_matvec_cont, clock_time_matvec_cont)
  do i_spin = 1, n_spin
  !TODO: BL This seems optimizable (BLAS2 -> BLAS3) maybe dgemm
  ! Multiply f_xc_matrix to 3ks_functions  
  !  do k_cnt=1, ndim1_o3KS
  !    do l_cnt=1, ndim2_o3KS
  !      call dgemv('n', n_basbas, n_basbas, 1d0, fxc_matrix, n_basbas, &
  !               pure_3ks(:,k_cnt,l_cnt,n_spin), 1, 0d0, &
  !               fxc_3ks(:,k_cnt,l_cnt,n_spin), 1)
  !    enddo 
  !   enddo

    call dgemm ('N','N', n_basbas, ndim1_o3KS * ndim2_o3KS, n_basbas, &
        1.0d0, fxc_matrix(1:n_basbas,1:n_basbas,i_spin), n_basbas, &
        pure_3ks(1:n_basbas,1:ndim1_o3KS,1:ndim2_o3KS,i_spin), n_basbas,&
        0.0d0, fxc_3ks(1:n_basbas,1:ndim1_o3KS,1:ndim2_o3KS,i_spin), n_basbas)

  end do

  call get_times(time_matvec_cont, clock_time_matvec_cont)
  write(info_str,'(2x,a)') '|'
  call localorb_info(info_str) 

  deallocate(fxc_matrix)

  call get_times(time_neutral_excitation_orbitals, &
           clock_time_neutral_excitation_orbitals)

end subroutine get_tddft_orbitals_2d

real*8 function dlamch()
  implicit none
  real*8	:: eps, rnd, sfmin, small
  intrinsic	:: epsilon

  rnd = 1.0d0

  if(rnd.eq.1.0d0) then
    eps = epsilon(0d0) * 0.5d0
  else
    eps = epsilon(0d0)
  endif
  sfmin = tiny(0d0)
  small = 1.0d0 / huge(0d0)
  if(small.ge.sfmin) then
    sfmin = small * ( 1.0d0 + eps )
  endif
  dlamch = sfmin
end function dlamch

subroutine get_local( alpha, a, ia, ja, desca )
  implicit none
  integer	:: ia, ja, desca(*)
  real*8	:: alpha, a(*)
  integer	:: iacol, iarow, mycol, myrow, npcol, nprow, iia, jja

  external	:: blacs_gridinfo, infog2l

  call blacs_gridinfo(desca(2), nprow, npcol, myrow, mycol)
  call infog2l(ia, ja, desca, nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)

  alpha = 0d0

  if(myrow.eq.iarow .and. mycol.eq.iacol) then
    alpha = a(iia+(jja-1)*desca(9))
  endif

  return
end subroutine get_local

subroutine set_local( a, ia, ja, desca, alpha )
  implicit none
  integer	:: ia, ja, desca(*)
  real*8	:: alpha, a(*)
  integer	:: iacol, iarow, mycol, myrow, npcol, nprow, iia, jja

  external	:: blacs_gridinfo, infog2l

  call blacs_gridinfo(desca(2), nprow, npcol, myrow, mycol)
  call infog2l(ia, ja, desca, nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)

  if(myrow.eq.iarow .and. mycol.eq.iacol) then
    a(iia+(jja-1)*desca(9)) = alpha
  endif

  return
end subroutine set_local

subroutine set_local_int( a, ia, ja, desca, alpha )
  implicit none
  integer	:: ia, ja, desca(*), a(*), alpha
  integer	:: iacol, iarow, mycol, myrow, npcol, nprow, iia, jja

  external	:: blacs_gridinfo, infog2l

  call blacs_gridinfo(desca(2), nprow, npcol, myrow, mycol)
  call infog2l(ia, ja, desca, nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)

  if(myrow.eq.iarow .and. mycol.eq.iacol) then
    a(iia+(jja-1)*desca(9)) = alpha
  endif

  return
end subroutine set_local_int

subroutine get_coeff_fxc3fn_v_1d(fxc_3fn, ovlp_3fn)
  use timing
  use prodbas
  use sbt_overlap_aims
  implicit none

  real*8, intent(OUT) :: fxc_3fn(n_basis_pairs, n_loc_prodbas)
  real*8, intent(OUT) :: ovlp_3fn(n_basis_pairs, n_loc_prodbas)
  real*8, allocatable :: coulomb_matr(:,:), tmp(:,:)
  integer :: info
  character(*), parameter :: func = 'get_coeff_fxc3fn_v_1d'

  ! Integrate ovlp3fn (T)
  call get_timestamps(time_ovlp3fn, clock_time_ovlp3fn )
  call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, fxc_3fn, ovlp_type_bare_or_hse_coul)
  call get_times(time_ovlp3fn, clock_time_ovlp3fn )
  ovlp_3fn = fxc_3fn

  ! Integrate V
  call get_timestamps(time_coulomb_matr, clock_time_coulomb_matr)
  allocate(coulomb_matr(n_basbas, n_loc_prodbas), stat=info)
  call check_allocation(info, 'coulomb_matr', func)
  allocate(tmp(n_basbas, n_loc_prodbas), stat=info)
  call check_allocation(info, 'coulomb_matr', func)

  if (use_logsbt) then
     call integrate_auxmat_by_atomic_sbt(coulomb_matr, ovlp_type_bare_or_hse_coul, .false.)
  else
     call integrate_coulomb_matr_v0(l_shell_max, coulomb_matr)
  end if
  call get_times(time_coulomb_matr, clock_time_coulomb_matr)
  tmp = coulomb_matr

  call get_timestamps(time_inv_coulomb_matr, clock_time_inv_coulomb_matr)
  if (use_asym_inv_sqrt) then
     ! transposed==.true. because get_v_multi_ovlp3fn() does
     call aims_stop('this procedure is not defined for Casida routine.')
     call asym_inv_sqrt_of_auxmat_scalapack(coulomb_matr, "Coulomb", .true.)
  elseif(use_scalapack) then
     call power_auxmat_scalapack(coulomb_matr, -0.5d0, "Coulomb")
     call power_auxmat_scalapack(tmp, -1d0, "f_xc coefficient")
  else
     call power_auxmat_lapack(coulomb_matr, -0.5d0, "Coulomb")
     call power_auxmat_lapack(tmp, -1d0, "f_xc coefficient")
  endif
  call get_times(time_inv_coulomb_matr, clock_time_inv_coulomb_matr)

  call get_timestamps(time_ovlp_multi, clock_time_ovlp_multi)
  call get_v_multi_fxc3fn(tmp, fxc_3fn)
  call get_v_multi_ovlp3fn(coulomb_matr, ovlp_3fn)
  call get_times(time_ovlp_multi, clock_time_ovlp_multi)

  deallocate(coulomb_matr)
  deallocate(tmp)
end subroutine get_coeff_fxc3fn_v_1d

subroutine get_coeff_fxc3fn_v_2d(fxc_3fn, ovlp_3fn)
  use timing
  use prodbas
  use runtime_choices
  use sbt_overlap_aims
  implicit none

  real*8, intent(out) :: fxc_3fn(n_basis_pairs, n_loc_prodbas)
  real*8, intent(out) :: ovlp_3fn(n_basis_pairs, n_loc_prodbas)
  real*8, allocatable :: coulomb_matr(:,:), coulomb_matr_1d(:,:), fxc_coul_mat(:,:)
  integer :: info
  character(*), parameter :: func = 'get_coeff_fxc3fn_v_2d'

  if(.not.use_scalapack) then
     call aims_stop('Cannot use 2d distribution without scalapack', func)
  endif

  ! Integrate ovlp3fn (T): <wave_i*wave_j|V|prod_wave_mu>
  call get_timestamps(time_ovlp3fn, clock_time_ovlp3fn )
  call integrate_ovlp3fn(l_shell_max, ext_l_shell_max, ovlp_3fn, ovlp_type_bare_or_hse_coul)
  call get_times(time_ovlp3fn, clock_time_ovlp3fn )
  
  ! Store also for fxc
  fxc_3fn = ovlp_3fn
  
  ! Integrate V : <prod_wave|V|prod_wave>
  call get_timestamps(time_coulomb_matr, clock_time_coulomb_matr)
  allocate(coulomb_matr(max_row_2d,max_col_2d), stat=info)
  call check_allocation(info, 'coulomb_matr', func)
  allocate(fxc_coul_mat(max_row_2d,max_col_2d), stat=info)
  call check_allocation(info, 'fxc_coul_mat', func)

  if (use_logsbt) then
     call integrate_auxmat_by_atomic_sbt(coulomb_matr, ovlp_type_bare_or_hse_coul, .true.)
  else
     allocate(coulomb_matr_1d(n_basbas, n_loc_prodbas), stat=info)
     call integrate_coulomb_matr_v0(l_shell_max, coulomb_matr_1d)
     call dist_1d_2d(n_basbas, coulomb_matr_1d, ubound(coulomb_matr_1d, 1), &
     &                         coulomb_matr, ubound(coulomb_matr,1))
     deallocate(coulomb_matr_1d)
  end if
  call get_times(time_coulomb_matr, clock_time_coulomb_matr)
  fxc_coul_mat = coulomb_matr

  ! Evaulate V^-0.5 and V^-1 for deriving the expansion coefficients
  call get_timestamps(time_inv_coulomb_matr, clock_time_inv_coulomb_matr)
  if (use_asym_inv_sqrt) then
     !   C := TV^-0.5.
     call aims_stop('This procedure is not defined or Casida routines')
     call asym_inv_sqrt_of_auxmat_scalapack_2d(coulomb_matr, "Coulomb",.false.)
  else
     call power_auxmat_scalapack_2d(coulomb_matr, -0.5d0, "Coulomb")
     call power_auxmat_scalapack_2d(fxc_coul_mat, -1.0d0, "f_xc coefficient")
  end if
  call get_times(time_inv_coulomb_matr, clock_time_inv_coulomb_matr)

  ! Derive expansion coefficients C := T V^-0.5
  call get_timestamps(time_ovlp_multi, clock_time_ovlp_multi)
  call get_v_multi_ovlp3fn_2(coulomb_matr, ovlp_3fn)
  call get_v_multi_fxc3fn_2(fxc_coul_mat, fxc_3fn)
  call get_times(time_ovlp_multi, clock_time_ovlp_multi)

  deallocate(coulomb_matr)
  deallocate(fxc_coul_mat)
end subroutine get_coeff_fxc3fn_v_2d

subroutine get_v_multi_fxc3fn(inv_coulomb_matr, fxc_3fn)
  use dimensions
  use prodbas
  use mpi_tasks
  use synchronize_mpi
  implicit none

  real*8, dimension(n_basbas, n_loc_prodbas) :: inv_coulomb_matr
  real*8, dimension(n_basis_pairs, n_loc_prodbas) :: fxc_3fn
  real*8, allocatable, dimension(:,:) :: temp_ovlp_matr
  real*8, allocatable, dimension(:,:) :: aux_ovlp_matr

  integer :: i_basis_1
  integer :: i_basis_2
  integer :: i_index_1
  integer :: i_index_2
  integer :: i_prodbas_1
  integer :: n_compute
  integer :: i_task

  call localorb_info('Multiplying fxc_3fn x V^-1 (1d-scalapack)', &
  &                  use_unit, '(2X,A)', OL_norm)

  allocate (temp_ovlp_matr(n_loc_prodbas, n_basis))
  allocate (aux_ovlp_matr(n_basbas, n_basis))

  i_index_1 = 0
  i_index_2 = 0
  do i_basis_1 =1, n_basis, 1

      temp_ovlp_matr (:,:) = 0
      n_compute = n_nghr_basis(i_basis_1)
      do i_basis_2 = 1, n_compute, 1
        i_index_1 = i_index_1 + 1

        temp_ovlp_matr(:,i_basis_2) = &
            fxc_3fn(i_index_1,:)
      enddo

      aux_ovlp_matr (:,:) = 0
      call dgemm('N', 'N', n_basbas, n_compute, &
                  n_loc_prodbas, 1.0d0, &
                  inv_coulomb_matr, n_basbas, &
                  temp_ovlp_matr, n_loc_prodbas, 0.d0, &
                  aux_ovlp_matr, n_basbas &
                  )

      call sync_matrix(aux_ovlp_matr, n_basbas, n_compute)

      temp_ovlp_matr(:,:) =0.d0
      i_prodbas_1 = 0
      do i_task = 1, n_tasks, 1
        do i_basis_2 = 1, n_loc_prodbas, 1

            i_prodbas_1 = map_prodbas(i_basis_2,i_task)
            if(i_prodbas_1.gt.0 .and. myid.eq.i_task-1) then

              temp_ovlp_matr(i_basis_2, 1:n_compute) = &
              aux_ovlp_matr(i_prodbas_1, 1:n_compute)
            endif

        enddo
      enddo

      do i_basis_2 = 1, n_compute, 1

        i_index_2 = i_index_2 + 1
        fxc_3fn(i_index_2,:) = &
          temp_ovlp_matr(:,i_basis_2)
      enddo
!  end of loop over i_basis_1
  enddo

  deallocate (temp_ovlp_matr)
  deallocate (aux_ovlp_matr)

end subroutine get_v_multi_fxc3fn

subroutine get_v_multi_fxc3fn_2(inv_coulomb_matr, fxc_3fn)
  use dimensions
  use prodbas
  use mpi_tasks
  use synchronize_mpi

  implicit none

  real*8, dimension(max_row_2d, max_col_2d) :: inv_coulomb_matr
  real*8, dimension(n_basis_pairs, n_loc_prodbas) :: fxc_3fn

  real*8, allocatable, dimension(:,:) :: aux1, aux2, aux3
  integer, allocatable, dimension(:) :: send_cnt, send_loc, send_dsp, idx
  integer, allocatable, dimension(:) :: recv_cnt, recv_loc, recv_dsp
  integer, parameter :: n_strip = 2048 ! strip size for stripmining multiplication
  integer :: n_off, n_off1, n_len1
  integer :: j, np1, npc, npr, npg, ns, nr, n_col, n_rows, mpierr, desc_fxc_3fn(9)

  call localorb_info('Multiplying fxc_3fn x V^-1 (2d-scalapack)', &
  &                  use_unit, '(2X,A)', OL_norm)

  ! Set up data for mpi_alltoall

  allocate(send_cnt(0:n_tasks-1), send_loc(0:n_tasks-1), send_dsp(0:n_tasks-1))
  allocate(recv_cnt(0:n_tasks-1), recv_loc(0:n_tasks-1), recv_dsp(0:n_tasks-1))
  allocate(idx(0:n_tasks-1))

  ! Get counts how many blocks we send/receive

  recv_cnt(:) = 0
  send_cnt(:) = 0
  do npr = 0, nprow_aux_2d-1
    do j = 0, n_basbas-1
      np1 = MOD(j/nb_aux,n_tasks)         ! Processor in 1D distribution
      npc = MOD(j/nb_aux_2d,npcol_aux_2d) ! Column processor in 2D distribution
      npg = global_id(npr,npc)            ! Global ID in 2D distribution
      if(npc==mypcol_aux_2d .and. npr==myprow_aux_2d) recv_cnt(np1) = recv_cnt(np1) + 1
      if(np1==myid) send_cnt(npg) = send_cnt(npg) + 1
    enddo
  enddo

  ! Set up location where data of remote processor starts

  ns = 0
  nr = 0
  do j=0,n_tasks-1
      send_loc(j) = ns
      recv_loc(j) = nr
      ns = ns + send_cnt(j)
      nr = nr + recv_cnt(j)
  enddo

  ! Scale from rows to variable counts

  send_dsp = send_loc*n_strip
  recv_dsp = recv_loc*n_strip
  send_cnt = send_cnt*n_strip
  recv_cnt = recv_cnt*n_strip

  ! Allocate auxiliary arrays for communictaion/matrix multiply

  allocate(aux1(n_strip, n_loc_prodbas*nprow_aux_2d))
  allocate(aux2(n_strip, max_col_2d))
  allocate(aux3(n_strip, max_col_2d))

  aux1 = 0
  aux2 = 0
  aux3 = 0

  do n_off = 0, n_basis_pairs-1, n_strip*nprow_aux_2d

    ! Redistribute strip of ovlp_3fn from 1D to 2D distribution
    ! =========================================================

    ! Put current strip of ovlp_3fn into send buffers

    idx(:) = send_loc(:)
    do npr = 0, nprow_aux_2d-1

      n_off1 = n_off + npr*n_strip
      if(n_off1>=n_basis_pairs) exit
      n_len1 = MIN(n_strip, n_basis_pairs-n_off1) ! length of current strip

      do j = 1, n_loc_prodbas
        n_col = INDXL2G(j, nb_aux, myid, n_tasks)
        npc = MOD((n_col-1)/nb_aux_2d,npcol_aux_2d) ! Column processor in 2D distribution
        npg = global_id(npr,npc)                    ! Global ID in 2D distribution
        idx(npg) = idx(npg) + 1
        aux1(1:n_len1,idx(npg)) = fxc_3fn(n_off1+1:n_off1+n_len1,j)
      enddo
    enddo

    ! Redistribute strip of ovlp_3fn using mpi_alltoallv

    call mpi_alltoallv(aux1,send_cnt,send_dsp,MPI_REAL8, &
                        aux2,recv_cnt,recv_dsp,MPI_REAL8, &
                        mpi_comm_global,mpierr)

    ! Sort received data to 2D distribution

    idx(:) = recv_loc(:)
    do j = 1, max_col_2d
      n_col = INDXL2G(j, nb_aux_2d, mypcol_aux_2d, npcol_aux_2d)
      np1 = MOD((n_col-1)/nb_aux,n_tasks) ! Processor in 1D distribution
      idx(np1) = idx(np1) + 1
      aux3(:,j) = aux2(:,idx(np1))
    enddo

    ! Calculate ovlp_3fn * inv_coulomb_matr for current strip
    ! =======================================================

    n_rows = min(n_strip*nprow_aux_2d, n_basis_pairs-n_off)

    call descinit(desc_fxc_3fn, n_rows, n_basbas, n_strip, nb_aux_2d, &
                  0, 0, my_blacs_ctxt_aux_2d, n_strip, mpierr)

    call pdgemm('N','N',n_rows,n_basbas,n_basbas,1.d0, &
                aux3, 1, 1, desc_fxc_3fn,                           &
                inv_coulomb_matr, 1, 1, aux_sc_desc_2d, 0.d0,        &
                aux2, 1, 1, desc_fxc_3fn)

    ! Redistribute result to 1D distribution
    ! ======================================

    ! put 2D distributed ovlp_3fn into send buffers

    aux3(:,:) = aux2(:,:)
    idx(:) = recv_loc(:)
    do j = 1, max_col_2d
      n_col = INDXL2G(j, nb_aux_2d, mypcol_aux_2d, npcol_aux_2d)
      np1 = MOD((n_col-1)/nb_aux,n_tasks) ! Processor in 1D distribution
      idx(np1) = idx(np1) + 1
      aux2(:,idx(np1)) = aux3(:,j)
    enddo

    ! Redistribute strip of ovlp_3fn using mpi_alltoallv
    ! recv_cnt/dsp and send_cnt/dsp have the reverse meaning here!

    call mpi_alltoallv(aux2,recv_cnt,recv_dsp,MPI_REAL8, &
                        aux1,send_cnt,send_dsp,MPI_REAL8, &
                        mpi_comm_global,mpierr)

    ! Sort received data to 1D distribution

    idx(:) = send_loc(:)
    do npr = 0, nprow_aux_2d-1

      n_off1 = n_off + npr*n_strip
      if(n_off1>=n_basis_pairs) exit
      n_len1 = MIN(n_strip, n_basis_pairs-n_off1) ! length of current strip

      do j = 1, n_loc_prodbas
        n_col = INDXL2G(j, nb_aux, myid, n_tasks)
        npc = MOD((n_col-1)/nb_aux_2d,npcol_aux_2d) ! Column processor in 2D distribution
        npg = global_id(npr,npc)                    ! Global ID in 2D distribution
        idx(npg) = idx(npg) + 1
        fxc_3fn(n_off1+1:n_off1+n_len1,j) = aux1(1:n_len1,idx(npg))
      enddo

    enddo

  enddo

  deallocate (aux1)
  deallocate (aux2)
  deallocate (aux3)
  deallocate (send_cnt, send_loc, send_dsp, idx)
  deallocate (recv_cnt, recv_loc, recv_dsp)

 contains
  INTEGER FUNCTION INDXL2G( INDXLOC, NB, IPROC, NPROCS )
  ! from Scalapack TOOLS directory with ISRCPROC omitted
  INTEGER INDXLOC, IPROC, NB, NPROCS
  INDXL2G = NPROCS*NB*((INDXLOC-1)/NB) + MOD(INDXLOC-1,NB) + IPROC*NB + 1
  END FUNCTION INDXL2G
end subroutine get_v_multi_fxc3fn_2

subroutine integrate_fxc_pairstates(n_occ, n_unocc, coord_of_center, basis_l_max, KS_eigenvector, dipole_mom)
!	COPY of integrate_dipmom_pairstates.f90
  
  implicit none

  integer n_occ
  integer n_unocc
  integer basis_l_max(n_species)
  real*8 coord_of_center(3)
  real*8 KS_eigenvector(n_basis, n_states, n_spin)
  real*8 dipole_mom(n_occ,n_unocc,n_spin,3)
  integer :: l_ylm_max
  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:,:), allocatable :: ylm_tab
  real*8, dimension(:,:), allocatable :: tmp_dipole_mom
  real*8, dimension(:,:,:), allocatable :: dipole_mom_basis
  real*8, dimension(:,:,:), allocatable :: wave_times_r(:,:,:)
  real*8 coord_current(3, n_max_angular)
  real*8 coord_wrt_center(3, n_max_angular)
  real*8 dist_tab(n_atoms, n_max_angular)
  real*8 i_r(n_atoms, n_max_angular)
  real*8 dir_tab(3,n_atoms, n_max_angular)
  real*8 trigonom_tab(4,n_atoms, n_max_angular)
  real*8 wave(n_basis, n_max_angular)
  real*8 radial_wave(n_basis)
  integer :: n_compute
  integer :: n_compute_onsite
  integer :: n_max_onsite_basis
  integer :: i_basis(n_basis)
  integer :: i_basis_onsite(n_basis)
  integer :: n_points
  real*8 :: partition_tab(n_max_angular)
  real*8 :: partition_tab_2atoms(n_atoms, n_max_angular)
  integer :: partition_type_temp = 1
  integer i_basis_1
  integer i_basis_2
  integer i_atom
  integer i_atom_1
  integer i_radial
  integer i_angular
  integer i_coord
  integer i_spin
  integer i_l
  integer i_m
  integer i_index
  integer i_compute
  integer i_compute_1
  integer i_task

  call get_timestamps(time_neutral_excitation_dipole, &
        clock_time_neutral_excitation_dipole)

  if(myid.eq.0) then
    write(use_unit,'(2x,a)') '|'
    write(use_unit,'(2X,A,A)') &
    "| Integrating the dipole moment assoicated with products of", &
    " KS states ..."
  endif

!     begin with general allocations
    l_ylm_max = l_wave_max

  allocate( ylm_tab( (l_ylm_max+1)**2, n_atoms, &
        n_max_angular) )
  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max) )

!     initialize index_lm

  i_index = 0
  do i_l = 0, l_wave_max, 1
    do i_m = -i_l, i_l
      i_index = i_index+1
      index_lm(i_m,i_l) = i_index
    enddo
  enddo

  n_max_onsite_basis = 0
  do i_atom = 1, n_atoms, 1 
    n_compute = 0 
    do i_basis_1 = 1, n_basis, 1
      if(basis_atom(i_basis_1).eq.i_atom) then
        n_compute = n_compute + 1
      endif
    enddo
    if (n_max_onsite_basis .lt. n_compute) then
        n_max_onsite_basis = n_compute
    endif
  enddo 
!      write(use_unit,*) "n_max_onsite_basis", n_max_onsite_basis

  allocate( tmp_dipole_mom (n_basis,n_max_onsite_basis))
  allocate( dipole_mom_basis (n_basis,n_basis,3))
  allocate( wave_times_r (n_max_onsite_basis,n_max_angular,3))

  dipole_mom_basis(:,:,:) = 0.d0

  i_task = myid + 1
!     perform partitioned integration, atom by atom, and point by point
!     This will be the outermost loop, to save evaluations of the potential.
!     and the Y_lm functions
  do i_atom = 1, n_atoms, 1

    if(myid.eq.0) then
    write(use_unit,*) " | i_atom: ", i_atom
    endif

      do i_radial = 1, n_radial(species(i_atom)), 1

        n_compute = 0
        i_basis = 0

        do i_angular = 1, n_angular(i_radial, species(i_atom)), 1

!     get current integration point coordinate
          do i_coord = 1, 3, 1
              coord_current(i_coord,i_angular) = &
                  coords(i_coord,i_atom ) + &
                  r_angular(i_coord, i_angular, i_radial, &
                  species(i_atom)) * &
                  r_radial(i_radial, species(i_atom))

              coord_wrt_center(i_coord,i_angular) = &
                coord_current(i_coord,i_angular) - &                    
                coord_of_center(i_coord) 
          enddo

!     compute atom-centered coordinates of current integration point,
!     as viewed from all atoms
          call tab_atom_centered_coords &
                ( coord_current(1,i_angular), &
                dist_tab(1,i_angular), i_r(1,i_angular), &
                dir_tab(1,1,i_angular) &
                )

!     evaluate the partition function for two atoms

              call evaluate_partition_tab_2atoms &
                ( i_atom, dist_tab(1:n_atoms,i_angular), &
                  i_r(1:n_atoms,i_angular), &
                  w_radial( i_radial, species (i_atom)), &
                  w_angular( i_angular, i_radial, species (i_atom)), &
                  partition_tab_2atoms(1, i_angular), &
                  partition_type_temp )


          partition_tab(i_angular) = 0.d0
          do i_atom_1 =1, n_atoms, 1
            partition_tab(i_angular) = &
                max( partition_tab(i_angular), &
                    partition_tab_2atoms(i_atom_1, i_angular) )
          enddo

          do i_atom_1 = 1, n_atoms
              if(dist_tab(i_atom_1,i_angular) .lt. 1.e-15) then
                partition_tab(i_angular) = 0.d0
                exit
              endif
          enddo

!     determine which basis functions are relevant at current integration point,
!     and tabulate their indices
          if (partition_tab(i_angular).gt.0.d0) &
                then

              call prune_basis_v1(dist_tab(1,i_angular), n_compute, &
                  i_basis)

          end if
        enddo

        i_basis_onsite(:) = 0
        n_compute_onsite = 0
        do i_compute = 1, n_compute, 1
        if (i_atom .eq. basis_atom(i_basis(i_compute))) then
            n_compute_onsite = n_compute_onsite + 1   
            i_basis_onsite(n_compute_onsite) = i_compute
        endif
        enddo

        if (n_compute.gt.0) then

          n_points = 0
          do i_angular = 1, n_angular(i_radial, species(i_atom)), 1

              if (partition_tab(i_angular).gt.0.d0) &
                  then
!     execute only if partition_tab.gt.0 here, i.e. if the integration point
!     makes sense
                n_points = n_points + 1
!     compute trigonometric functions of spherical coordinate angles
!     of current integration point, viewed from all atoms
                call tab_trigonom &
                      ( dir_tab(1,1,i_angular), &
                      trigonom_tab(1,1,i_angular) &
                      )

!     tabulate distance and Ylm's w.r.t. other atoms
                call tab_wave_ylm &
                      ( trigonom_tab(1,1,i_angular), basis_l_max, &
                      l_ylm_max, &
                      ylm_tab(1,1,i_angular) )

!           tabulate total wave function value for each basis function
                call evaluate_waves_v0 &
                      (i_r(1,i_angular), l_ylm_max, &
                      ylm_tab(1,1,i_angular), &
                      dist_tab(1,i_angular), index_lm, n_compute, &
                      i_basis, radial_wave(1), &
                      wave(1,n_points))

                  do i_coord = 1, 3, 1
                  do i_compute_1 = 1, n_compute_onsite, 1

                      i_compute = i_basis_onsite(i_compute_1)
                      wave_times_r(i_compute_1,n_points,i_coord) = &
                      wave(i_compute,n_points)  * &
                      coord_wrt_center(i_coord, i_angular) 

                  enddo
                  enddo

                  do i_compute = 1, n_compute, 1
                    wave(i_compute,n_points)= wave(i_compute,n_points)  &
                    * partition_tab_2atoms(basis_atom(i_basis(i_compute)), &
                                            i_angular)
                  enddo

!     end if (partition_tab.gt.0)
          end if

!     end angular integration loop
      enddo
!     add the contribution from the current shell

      do i_coord=1, 3, 1
        tmp_dipole_mom(:,:) = 0.d0
        call dgemm('N','T',n_compute,n_compute_onsite,n_points,1.d0, &
                  wave(1,1), n_basis, &
                  wave_times_r(1,1,i_coord), n_max_onsite_basis, 0.d0, &
                  tmp_dipole_mom(1,1), n_basis )
  
        do i_compute = 1, n_compute_onsite,  1
          i_basis_1 = i_basis(i_basis_onsite(i_compute))
          do i_compute_1 = 1, n_compute,  1
            i_basis_2 = i_basis(i_compute_1)
            dipole_mom_basis(i_basis_2,i_basis_1, i_coord) = &
                dipole_mom_basis(i_basis_2, i_basis_1, i_coord) + &
                tmp_dipole_mom(i_compute_1,i_compute)
          enddo
        enddo

      enddo
!     end if (n_compute.gt.0) then
    end if

!       end radial integration loop
  enddo

!     end integration loop over atoms
  enddo

  do i_basis_1 = 1, n_basis, 1
    do i_basis_2 = 1, i_basis_1-1, 1
      if(basis_atom(i_basis_2) .ne. basis_atom(i_basis_1)) then
      dipole_mom_basis(i_basis_2,i_basis_1,:) = &
        dipole_mom_basis(i_basis_2,i_basis_1,:) + &
        dipole_mom_basis(i_basis_1,i_basis_2,:) 

        dipole_mom_basis(i_basis_1,i_basis_2,:) = &
          dipole_mom_basis(i_basis_2,i_basis_1,:) 
      endif
      enddo
    enddo
      
  deallocate(tmp_dipole_mom)
  allocate(tmp_dipole_mom(n_basis,n_occ))

  do i_coord = 1, 3, 1
    do  i_spin = 1, n_spin, 1
      call dgemm('N', 'N', n_basis, n_occ, n_basis, 1.d0, &
                  dipole_mom_basis(1,1,i_coord), n_basis, &
                  KS_eigenvector(1,1,i_spin), n_basis, 0.d0, &      
                  tmp_dipole_mom(1,1), n_basis)
                  
      call dgemm('T', 'N', n_occ, n_unocc, n_basis, 1.d0, &
                  tmp_dipole_mom(1,1), n_basis, &
                  KS_eigenvector(1,n_states-n_unocc+1,i_spin), n_basis, &
                  0.d0, dipole_mom(1,1,i_spin,i_coord), n_occ)

    enddo
  enddo

  if(myid==0) write(use_unit,'(2x,a)') '|'

  if (allocated(ylm_tab)) then
    deallocate(ylm_tab)
  end if
  if (allocated(index_lm)) then
    deallocate(index_lm)
  end if
  if (allocated(tmp_dipole_mom)) then
    deallocate(tmp_dipole_mom)
  end if
  if (allocated(dipole_mom_basis)) then
    deallocate(dipole_mom_basis)
  end if
  if (allocated(wave_times_r)) then
    deallocate(wave_times_r)
  end if

  call get_times(time_neutral_excitation_dipole, clock_time_neutral_excitation_dipole)

  return
end subroutine integrate_fxc_pairstates

subroutine neg_energy_error(value)
  implicit none
  real*8, intent(in) :: value

  if(myid==0) then
    write(use_unit,'(2x,a,f12.4)') '| *****      Negative exitation energy found: ', value
    write(use_unit,'(2x,a)') '| *****  This points to an insufficiently large or badly configured basis set.'
    write(use_unit,'(2x,a)') '| *****  When using Gaussian basis sets, try increasing the radial and angual grid.'
    write(use_unit,'(2x,a)') '| *****  Setting this value to zero now to avoid further errors.'
    write(use_unit,'(2x,a)') '|'
  endif
end subroutine neg_energy_error

end module optics_module
