subroutine output_coulomb_integrals_mo_2d ()
  use dimensions
  use physics
  use prodbas
  use constants
  use runtime_choices
  use synchronize_mpi
  use mpi_tasks
  use timing
  use hartree_fock
  use localorb_io
  implicit none
  real*8, dimension(:, :), allocatable :: ovlp_3fn_coul_int
  real*8, dimension(:, :, :, :), allocatable :: ovlp_3ks_coul_int
  integer :: n_pair_states
  integer :: col_dim_proc, row_dim_proc, start_id_row, start_id_col
  integer :: col_id_2d, row_id_2d, i_id_4d, j_id_4d, k_id_4d, l_id_4d
  integer :: num_pcol, num_prow, my_coul_int_ctxt, my_prow, my_pcol
  integer :: mpi_comm_rows, mpi_comm_cols
  integer :: iia, jja, iarow, iacol
  integer :: num_blocks, na_rows, na_cols
  integer, external :: numroc
  integer :: coul_int_desc(9), nprow, npcol, myrow, mycol, errnum
  real*8, dimension(:, :, :, :), allocatable :: loc_ovlp_3ks_coul_int
  real*8 :: total_mem, proc_mem
  character*128 :: info_str
  integer :: k_cnt, i_spin
  integer :: i_proc, tmp_int, mpi_status(mpi_status_size), sendrecv_dims, i_rot, j_rot
  integer :: i_state, a_state, j_state, b_state
  integer, allocatable, dimension(:,:) :: tmp_bcastidx, mpi_bcastidx
  real*8, allocatable, dimension(:) :: mpi_vec, tmp_vecr
  integer, allocatable, dimension(:,:) :: tmp1r, tmp2r, tmp1s, tmp2s
  integer :: i_state_loc, a_state_loc, len_of_recs, j_state_loc, b_state_loc
  real*8,  external :: ddot
  integer, external :: mat_idx_coul_int
  real*8 :: kernel, tmp_real
  real*8, allocatable, dimension(:,:) ::  par_coul_int
  real*8, allocatable, dimension(:,:) ::  global_coul_int
  real*8, allocatable, dimension(:, :, :, :) :: coul_int_4d
  integer :: i, j, k, l
  integer, dimension(n_tasks) :: na_rows_all_procs, na_cols_all_procs
  integer :: block_size_row, block_size_col
  integer:: global_i_state, global_j_state 
  real*8 :: temp_4ks
  integer :: info
  character(*), parameter :: func = 'output_coulomb_integrals_mo_2d'

!  if (myid == 0) write(use_unit,'(10X,A)') "output_coulomb_integrals_mo_2d"
  allocate(ovlp_3fn_coul_int(n_basis_pairs,n_loc_prodbas), stat=info)
  call check_allocation(info, 'ovlp_3fn_coul_int', func)
  call get_coeff_3fn_v_2d(ovlp_3fn_coul_int)
!  if (myid == 0) write(use_unit,'(10X,A)') "after getting ovlp_3fn_coul_int"
  ndim1_o3ks = (n_states-1)/np1_o3ks + 1
  ndim2_o3ks = (n_states-1)/np2_o3ks + 1

  allocate(ovlp_3ks_coul_int(n_basbas,ndim1_o3ks,ndim2_o3ks,n_spin))
  call transform_ovlp3fn_2(n_states,KS_eigenvector, ovlp_3fn_coul_int, ovlp_3ks_coul_int)
  n_pair_states = n_states * n_states
! Set up "mostly" square processor grid ! 2*1 for 2 cpus, num_prow = 2, num_pcol = 1
  do num_pcol = nint(dsqrt(dble(n_tasks))), 2, -1
    if(mod(n_tasks, num_pcol) == 0) exit
  end do
  num_prow = n_tasks/num_pcol

! Initialize Blacs grid
  my_coul_int_ctxt = mpi_comm_world
  call blacs_gridinit(my_coul_int_ctxt, 'C' , num_prow, num_pcol)
  call blacs_gridinfo(my_coul_int_ctxt, num_prow, num_pcol, my_prow, my_pcol)

! Local dimension of BSE matrix: na_rows, na_cols
  num_blocks = min(n_tasks,floor(dble(n_pair_states)/dble(n_tasks)))
  na_rows = numroc(n_pair_states, num_blocks, my_prow, 0, num_prow)
  na_cols = numroc(n_pair_states, num_blocks, my_pcol, 0, num_pcol)
!  print*, "myid, num_blocks", myid, num_blocks
  call descinit(coul_int_desc, n_pair_states, n_pair_states, num_blocks, num_blocks, 0, 0, my_coul_int_ctxt, na_rows, errnum)
  call blacs_gridinfo(coul_int_desc(2), nprow, npcol, myrow, mycol)
!  do i_state = 0, n_tasks - 1
!    if (myid == i_state) then
!      print*, 'proc, nprow, npcol', myid, nprow, npcol
!      print*, 'proc, myrow, mycol', myid, myrow, mycol
!!      print*, 'proc, na_rows, na_cols', myid, na_rows, na_cols
!!      print*, 'proc, num_prow, num_pcol', myid, num_prow, num_pcol
!    end if
!  end do

  total_mem = 8d0 * dble(n_pair_states)**2
  proc_mem = 8d0 * dble(na_rows) * dble(na_cols)
  write(info_str,'(2x,a,i5.1,a,i8.1,a,i8.1,a,f8.3,a)') '| Processor ', myid, &
    ' gets ', na_rows, ' rows and ', na_cols,' columns  (',proc_mem/2d0**20,' MB).'
  call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)
  if(myid == 0) print*, '| Number of total rotations for matrix building: ', n_tasks

  k_cnt = 0
  !  Loop over all occupied states in energy window
  do i_state = 1, n_states
    ! Is this state treated my this process?
    if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
    ! Yes, then loop over all unoccupied states in energy window
    do a_state = 1, n_states
      ! Is this state treated on this process?
      if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle
         ! If yes increase the local number of pairstates
         k_cnt = k_cnt + 1
    enddo
  enddo
  call MPI_BARRIER(mpi_comm_world,errnum)

! Find the maximal number of pairstates by sending all local number of pairstates to process 0 and compare 
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
  sendrecv_dims = i_state * k_cnt ! why k_cnt**2 ???
  call mpi_bcast(sendrecv_dims, 1, mpi_integer, 0, mpi_comm_world, errnum)

! Matrix allocation
  allocate(mpi_bcastidx(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'mpi broadcast index for coulomb matrix')
  mpi_bcastidx = 0
  allocate(mpi_vec(sendrecv_dims),stat=errnum)
  call check_allocation(errnum, 'mpi broadcast vector for coulomb matrix')
  mpi_vec = 0d0
! Now calculate <ia|V|jb>
  if(myid == 0) print*, 'Building <ia|V|jb>'
  allocate(tmp1s(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'rotation index 1s')
  allocate(tmp2s(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'rotation index 2s')
  tmp1s(:,1) = 0
  tmp2s(:,1) = 0

  k_cnt=1
  do i_state = 1, n_states
    if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
    do a_state = 1, n_states
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

  allocate(tmp1r(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, ' tmp1r res dist')
  allocate(tmp2r(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'rotation index 2r')
  allocate(tmp_vecr(sendrecv_dims),stat=errnum)
  call check_allocation(errnum, ' tmp_vecr res dist')
  allocate(par_coul_int(na_rows,na_cols),stat=errnum)
  call check_allocation(errnum, 'parallel bse matrix')
  par_coul_int = 0d0
  allocate(loc_ovlp_3ks_coul_int(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin),stat=errnum)
  call check_allocation(errnum, 'rotation local o3ks')

! Start rotating ovlp_3ks elements
  do i_rot=1, n_tasks
    if(myid==0) write(use_unit,'(2x,a,I5.1)') '| starting rotation: ', i_rot
    if(myid==i_rot-1) then
      tmp1r = tmp1s
      tmp2r = tmp2s
      loc_ovlp_3ks_coul_int = ovlp_3ks_coul_int
    endif
    call mpi_bcast(tmp1r, 2*sendrecv_dims, mpi_integer, i_rot-1, mpi_comm_world, errnum)
    call mpi_bcast(tmp2r, 2*sendrecv_dims, mpi_integer, i_rot-1, mpi_comm_world, errnum)
    call mpi_bcast(loc_ovlp_3ks_coul_int, n_basbas*ndim1_o3KS*ndim2_o3KS,&
          mpi_double_precision, i_rot-1, mpi_comm_world, errnum)
    len_of_recs = 1
    i_spin = 1
    do i_state = 1, n_states
      if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
      do a_state = 1, n_states
        if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle
        i_state_loc = loc_dim1_o3KS(i_state)
        a_state_loc = loc_dim2_o3KS(a_state)
        do k_cnt=1, sendrecv_dims
          if(tmp1r(k_cnt,1).ne.0 .and. tmp2r(k_cnt,1).ne.0) then
            kernel = ddot(n_basbas,&
                          ovlp_3ks_coul_int(:,i_state_loc,a_state_loc,i_spin),1,&
                          loc_ovlp_3ks_coul_int(:,tmp1r(k_cnt,2),tmp2r(k_cnt,2),i_spin),1)
! give it to par_bse
            tmp_real = kernel
            mpi_bcastidx(len_of_recs,1) = mat_idx_coul_int(i_state,a_state, n_states)
            mpi_bcastidx(len_of_recs,2) = mat_idx_coul_int(tmp1r(k_cnt,1),tmp2r(k_cnt,1),n_states)
            mpi_vec(len_of_recs) = tmp_real
            len_of_recs = len_of_recs + 1
          endif
        enddo
      enddo
    enddo
! j_rot
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
          call infog2l(tmp1r(i_state,1), tmp1r(i_state,2), coul_int_desc, &
            nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)

          if(myrow.eq.iarow.and.mycol.eq.iacol) then
            par_coul_int(iia,jja) = tmp_vecr(i_state)
          endif

          if(tmp1r(i_state,1).ne.tmp1r(i_state,2)) then
            call infog2l(tmp1r(i_state,2), tmp1r(i_state,1), coul_int_desc, &
              nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)
            if(myrow.eq.iarow.and.mycol.eq.iacol) then
              par_coul_int(iia,jja) = tmp_vecr(i_state)
            endif
          endif
        endif
      enddo
    enddo !j_rot
  end do ! i_rot

! get index of the whole matrix
  allocate(global_coul_int(n_pair_states, n_pair_states))
  global_coul_int = 0.0d0
  do i_state = 1, na_rows
    do j_state = 1, na_cols
      global_i_state = (((i_state - 1) / num_blocks) * nprow + myrow) *&
                       num_blocks + mod(i_state - 1, num_blocks) + 1
      global_j_state = (((j_state - 1) / num_blocks) * npcol + mycol) *&
                       num_blocks + mod(j_state - 1, num_blocks) + 1
      global_coul_int(global_i_state, global_j_state) = & 
        par_coul_int(i_state, j_state)
    end do
  end do
  call sync_matrix(global_coul_int, n_pair_states, n_pair_states)
  if (myid == 0) then
    write(use_unit,'(2x, A)') "Writing bielectron Coulomb integrals"
    if (coulomb_integral_format == 'qtmp') then
      open(unit = 99, file = "bielec_mo")
  
      do l = 1, n_states
        do k = 1, n_states
          do j = l, n_states
            do i = max(j, k), n_states
              if (i == j .and. k > l) continue
  !            write(99, *) i, j, k, l, global_coul_int((l - 1) * n_states + k, &
  !                                       (j - 1) * n_states + i)
              temp_4ks = global_coul_int((i - 1) * n_states + k, &
                                         (j - 1) * n_states + l)
              if (abs(temp_4ks) > 1E-16) write(99, *) i, j, k, l, temp_4ks
            end do
          end do
        end do
      end do
    else if (coulomb_integral_format == 'full') then
      open(unit = 99, file = "coulomb_integrals_mo.out")
      do i = 1, n_states
        do j = 1, n_states
          do k = 1, n_states
            do l = 1, n_states
              write(99, *) i, j, k, l, global_coul_int((i - 1) * n_states + j, & 
                                                       (k - 1) * n_states + l)
            end do
          end do
        end do
      end do
    end if
    close(99)
  end if
!  deallocate(coul_int_4d)
  deallocate(global_coul_int)
  deallocate(mpi_bcastidx)
  deallocate(mpi_vec)
  deallocate(tmp_vecr)
  deallocate(tmp1r)
  deallocate(tmp2r)
  deallocate(tmp1s)
  deallocate(tmp2s)
  deallocate(par_coul_int)
  deallocate(ovlp_3fn_coul_int)
end subroutine

pure integer function mat_idx_coul_int(i,a, n_states) result(idx)
 implicit none
 integer, intent(in) :: i, a, n_states
 idx = (i - 1) * n_states + a 
end function mat_idx_coul_int
