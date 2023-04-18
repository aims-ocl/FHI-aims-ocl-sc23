subroutine bse_2d(ndim1_o3ks_bse, ndim2_o3ks_bse, ovlp_3ks_bse, qpe)
  use localorb_io
  use runtime_choices
  use dimensions
  use physics
  use prodbas
  use gw_para
  use mpi_utilities
  use synchronize_mpi
  use elsi_wrapper, only: aims_elsi_stdevp
  use evaluate_polarisability_freq, only: evaluate_polarisability_freq_2
  use constants, only : hartree
!  use scalapack_wrapper
!  use hartree_fock
  use species_data, only : l_shell_max
!  use optics_module, only : integrate_fxc_pairstates, dipole_moment,&
!                            coord_of_center, atom_coord_wrt_center
!  use free_atoms
  use timing

  implicit none
!***********************************************************************************************
!* Arguments:
!* ndim1_o3ks_bse, ndim2_o3ks_bse: local dimensions along the 2 basis direction of ovlp_3ks
!* ovlp_3ks_bse: ovlp_3ks used in this subroutine, to calculate ovlp_4ks =
!*               ovlp_3ks_bse * ovlp_3ks_bse, and w_4ks = ovlp_3ks_bse * w_abf * ovlp_3ks_bse
!* qpe: quasi-particle energies calculated usually by GW, dimension is set to be
!       n_states, which means all empty states are used in constructing BSE matrix by default
!***********************************************************************************************
  integer, intent(IN) :: ndim1_o3ks_bse,ndim2_o3ks_bse
  real*8, dimension(n_basbas,ndim1_o3ks_bse,ndim2_o3ks_bse,n_spin), intent(IN) :: ovlp_3ks_bse
  real*8, dimension(n_states_k(1)), intent(IN) :: qpe
!***********************************************************************************************
!* Variables in other modules and can be used without announcement
!* n_states: integer, number of states, which equals to n_basis
!* n_electrons: double, number of elelectrons
!* n_tasks: number of tasks/processors
!* n_basbas: number of auxilary basis, basbas means two basis function is now
!            described by one auxilary basis
!***********************************************************************************************
!* Local variables:
!*
  integer :: errnum, max_pairs, k_cnt, i_spin, i_k_point
  integer :: occ_states, unocc_states, i_state, a_state, j_state, b_state
  integer :: num_pcol, num_prow, my_bse_ctxt, my_prow, my_pcol
  integer :: n_pair_states, num_blocks, na_rows, na_cols
  integer, external :: numroc
  integer :: bse_desc(9), nprow, npcol, myrow, mycol
  integer :: mpi_comm_rows, mpi_comm_cols
  integer :: iia, jja, iarow, iacol
!  real*8, dimension(:, :), allocatable :: ovlp_3fn_bse
!  real*8, dimension(:, :, :, :), allocatable :: ovlp_3ks_bse, loc_ovlp_3ks_bse, o3ks_for_w
  real*8, dimension(:, :, :, :), allocatable :: loc_ovlp_3ks_bse, o3ks_for_w
  real*8 :: total_mem, proc_mem
  character*128 :: info_str
  integer :: i_proc, tmp_int, mpi_status(mpi_status_size), sendrecv_dims, i_rot, j_rot
  integer, allocatable, dimension(:,:) :: tmp_bcastidx, mpi_bcastidx
  real*8, allocatable, dimension(:) :: mpi_vec, tmp_vecr
  integer, allocatable, dimension(:,:) :: tmp1r, tmp2r, tmp1s, tmp2s
  integer :: i_state_loc, a_state_loc, len_of_recs, j_state_loc, b_state_loc
  real*8,  external :: ddot
  integer, external :: mat_idx
  real*8 :: kernel, tmp_real, bse_alpha
  real*8, allocatable, dimension(:,:) ::  par_bse, par_bse_w, par_evecs
  real*8, allocatable, dimension(:) :: par_evals
! variables for screened coulomb
  real*8, dimension(:,:), allocatable :: screened_coulomb_bse
  integer :: n_blk
  integer, allocatable :: n_homo(:)
  integer :: n_first(n_spin)
!  integer :: i_spin, i_state
  real*8 :: omega_full_i_freq
  real*8  :: occ_numbers_bse(n_states, n_spin)
!  real*8  :: occ_numbers_bse(n_states_k(1),n_spin)
!  integer :: n_high_state
  real*8, dimension(:,:), allocatable :: polar_freq
  integer :: info
  integer :: i_abf, j_abf, i
  real*8, dimension(:), allocatable :: global_w_row
  integer :: abf_cnt, loc_abf_cnt, counter
  real*8, dimension(:, :, :, :), allocatable :: sc_ovlp
  real*8 :: temp
! 2d aux
  integer, allocatable :: basbas2row(:), basbas2col(:)
  character(*), parameter :: deffmt = '2X'
  character(*), parameter :: func = 'bse_2d'
! reduce
  integer :: reduced_occ_states, reduced_unocc_states
! Dummy for ELSI
  real*8, dimension(:,:), allocatable :: dummy
! oscillator
  real*8, allocatable, dimension(:) :: osc_str
  real*8, allocatable, dimension(:,:) :: tmom
  real*8, allocatable, dimension(:,:,:,:) :: dipole_moment
  real*8, allocatable, dimension(:)       :: coord_of_center
  real*8, allocatable, dimension(:,:)     :: atom_coord_wrt_center

  call get_times(time_bse_init, clock_time_bse_init)
  if(myid.eq.0) then
    write(use_unit,*)"-------------------------------------------------"
    write(use_unit,'(2X,A)') 'Start to calculate the coulomb 4ks'
!    print*, 'n_states, n_electrons, n_tasks, n_basbas'
!    print*, n_states, n_electrons, n_tasks, n_basbas
!    print*, 'qpe', qpe
  endif
! singlet or triplet
  if(bse_singlet) then
    bse_alpha = 2
  elseif(bse_triplet) then
    bse_alpha = 0
  else
    call aims_stop('bse_s_t keyword needs to be set either singlet or triplet')
  end if

! bse_reduce_matrix
  if (bse_reduce_matrix) then
    do i = 1, n_states
!    do i = 1, n_states_k(1)
      if (qpe(i) > bse_reduce_occ / hartree) then
        exit
      end if
    end do
    bse_lower_limit = i
    do i = n_states, 1, -1
!    do i = n_states_k(1), 1, -1
      if (qpe(i) < bse_reduce_unocc / hartree) then
        exit
      end if
    end do
    bse_upper_limit = i
!    if(myid == 0) then
!      print*, 'bse_reduce_occ', bse_reduce_occ
!      print*, 'bse_lower_limit', bse_lower_limit
!      print*, 'qpe(bse_lower_limit)', qpe(bse_lower_limit)
!      print*, 'bse_reduce_unocc', bse_reduce_unocc
!      print*, 'bse_upper_limit', bse_upper_limit
!      print*, 'qpe(bse_upper_limit)', qpe(bse_upper_limit)
!    end if
  end if

! bse_lower_limit, bse_upper_limit
  bse_lower_limit = max(1, bse_lower_limit)
  bse_upper_limit = min(bse_upper_limit, n_states)
!  bse_upper_limit = min(bse_upper_limit, n_states_k(1))



! Get number of occ_states, unocc_states, and their product n_pair_states, which is also the
! dimension of mat block A in BSE matrix
  i_spin = 1
  occ_states = 0
  unocc_states = 0
  do i_state = bse_lower_limit, bse_upper_limit
!  do i_state = 1, n_states_k(1)
    if(occ_numbers(i_state, i_spin, n_k_points) < 1e-6) then
      unocc_states = unocc_states + 1
    else
      occ_states = occ_states + 1
    end if
  end do
  reduced_occ_states = occ_states - bse_lower_limit + 1
  reduced_unocc_states = bse_upper_limit - occ_states
!  n_pair_states = occ_states * unocc_states
  n_pair_states = reduced_occ_states * reduced_unocc_states

!  if(myid.eq.0) then
!    print*, 'bse_lower_limit, bse_upper_limit', bse_lower_limit, bse_upper_limit
!    print*, 'occ_states, unocc_states, reduced_occ_states, reduced_unocc_states, n_pair_states'
!!    print*, 'occ_states, unocc_states, n_pair_states'
!    print*, occ_states, unocc_states, reduced_occ_states, reduced_unocc_states, n_pair_states
!!    print*, occ_states, unocc_states, n_pair_states
!  endif
! call initialize_prodbas() ! called in previous subroutines

! Set up "mostly" square processor grid ! 2*1 for 2 cpus, num_prow = 2, num_pcol = 1
  do num_pcol = nint(dsqrt(dble(n_tasks))), 2, -1
    if(mod(n_tasks, num_pcol) == 0) exit
  end do
  num_prow = n_tasks/num_pcol
!  if(myid == 0) print*, 'num_prow, num_pcol', num_prow, num_pcol

! Initialize Blacs grid
  my_bse_ctxt = mpi_comm_world
  call blacs_gridinit(my_bse_ctxt, 'C' , num_prow, num_pcol)
  call blacs_gridinfo(my_bse_ctxt, num_prow, num_pcol, my_prow, my_pcol)

! Local dimension of BSE matrix: na_rows, na_cols
  num_blocks = min(n_tasks,floor(dble(n_pair_states)/dble(n_tasks)))
!  if(myid == 0) print*, 'num_blocks', num_blocks
  na_rows = numroc(n_pair_states, num_blocks, my_prow, 0, num_prow)
  na_cols = numroc(n_pair_states, num_blocks, my_pcol, 0, num_pcol)

  call descinit(bse_desc, n_pair_states, n_pair_states, num_blocks, num_blocks, 0, 0, my_bse_ctxt, na_rows, errnum)
  call blacs_gridinfo(bse_desc(2), nprow, npcol, myrow, mycol)
! nprow, npcol, myrow, mycol correspond to
! num_prow, num_pcol, my_prow, my_pcol
!  if(myid == 0) print*, 'proc0, nprow, npcol, myrow, mycol', nprow, npcol, myrow, mycol
!  if(myid == 1) print*, 'proc1, nprow, npcol, myrow, mycol', nprow, npcol, myrow, mycol

  total_mem = 8d0 * dble(n_pair_states)**2
  proc_mem = 8d0 * dble(na_rows) * dble(na_cols)
  write(info_str,'(2x,a,i5.1,a,i8.1,a,i8.1,a,f8.3,a)') '| Processor ', myid, &
    ' gets ', na_rows, ' rows and ', na_cols,' columns  (',proc_mem/2d0**20,' MB).'
  call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)
  if(myid == 0) print*, '| Number of total rotations for matrix building: ', n_tasks

  call get_times(time_bse_init, clock_time_bse_init)
! Really start to work
! Get local number of pairstates
  call get_times(time_bse_build_V, clock_time_bse_build_V)
  k_cnt = 0
  !  Loop over all occupied states in energy window
  do i_state = bse_lower_limit, occ_states
!  do i_state=1, occ_states
    ! Is this state treated my this process?
    if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
    ! Yes, then loop over all unoccupied states in energy window
    do a_state=occ_states+1, bse_upper_limit
!    do a_state=occ_states+1, n_states
      ! Is this state treated on this process?
      if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle
         ! If yes increase the local number of pairstates
         k_cnt = k_cnt + 1
    enddo
  enddo
  call MPI_BARRIER(mpi_comm_world,errnum)
!  if(myid == 0) print*, 'cpu0, k_cnt, ', k_cnt
!  if(myid == 1) print*, 'cpu1, k_cnt', k_cnt
!  if(myid == 0) print*, 'np1_o3KS, np2_o3KS, myp1_o3KS, myp2_o3KS', np1_o3KS, np2_o3KS, myp1_o3KS, myp2_o3KS
!  if(myid == 1) print*, 'np1_o3KS, np2_o3KS, myp1_o3KS, myp2_o3KS', np1_o3KS, np2_o3KS, myp1_o3KS, myp2_o3KS

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
!  if(myid == 0) print*, 'cpu0, sendrecv_dims', sendrecv_dims
!  if(myid == 1) print*, 'cpu1, sendrecv_dims', sendrecv_dims

! Matrix allocation
  allocate(mpi_bcastidx(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'mpi broadcast index for bse matrix')
  mpi_bcastidx = 0
  allocate(mpi_vec(sendrecv_dims),stat=errnum)
  call check_allocation(errnum, 'mpi broadcast vector for bse matrix')
  mpi_vec = 0d0
! Now calculate <ia|V|jb>
  if(myid == 0) print*, 'Building <ia|V|jb>'
  allocate(tmp1s(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'rotation index 1s')
  allocate(tmp2s(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'rotation index 2s')
  tmp1s(:,1) = 0
  tmp2s(:,1) = 0

!  call aims_stop('bse_develop')
  k_cnt=1
!  do i_state = 1, occ_states
  do i_state = bse_lower_limit, occ_states
    if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
!    do a_state = occ_states+1, n_states
    do a_state = occ_states + 1, bse_upper_limit
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
! check tmp1s, tmp2s
!  if(myid == 1) then
!    print*, 'cpu 1'
!    do k_cnt = 1, sendrecv_dims
!      print*, 'tmp1s(k_cnt, :) , tmp2s(k_cnt, :), k_cnt', k_cnt
!      print*, tmp1s(k_cnt,1), tmp1s(k_cnt,2), tmp2s(k_cnt,1), tmp2s(k_cnt,2)
!    end do
!  end if
  allocate(tmp1r(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, ' tmp1r res dist')
  allocate(tmp2r(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'rotation index 2r')
  allocate(tmp_vecr(sendrecv_dims),stat=errnum)
  call check_allocation(errnum, ' tmp_vecr res dist')
  allocate(par_bse(na_rows,na_cols),stat=errnum)
  call check_allocation(errnum, 'parallel bse matrix')
  par_bse = 0d0
  allocate(loc_ovlp_3ks_bse(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin),stat=errnum)
  call check_allocation(errnum, 'rotation local o3ks')

! Start rotating ovlp_3ks elements
  do i_rot=1, n_tasks
!    call get_times(time_bse_rot, clock_time_bse_rot)

    if(myid==0) write(use_unit,'(2x,a,I5.1)') '| starting rotation: ', i_rot
    if(myid==i_rot-1) then
      tmp1r = tmp1s
      tmp2r = tmp2s
      loc_ovlp_3ks_bse = ovlp_3ks_bse
    endif
    call mpi_bcast(tmp1r, 2*sendrecv_dims, mpi_integer, i_rot-1,&
          mpi_comm_world, errnum)
    call mpi_bcast(tmp2r, 2*sendrecv_dims, mpi_integer, i_rot-1,&
          mpi_comm_world, errnum)
    call mpi_bcast(loc_ovlp_3ks_bse, n_basbas*ndim1_o3KS*ndim2_o3KS,&
          mpi_double_precision, i_rot-1, mpi_comm_world, errnum)
    len_of_recs = 1
    i_spin = 1
!    do i_state=1, occ_states
    do i_state = bse_lower_limit, occ_states
      if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
!      do a_state=occ_states+1, n_states
      do a_state = occ_states + 1, bse_upper_limit
        if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle
        i_state_loc = loc_dim1_o3KS(i_state)
        a_state_loc = loc_dim2_o3KS(a_state)
        do k_cnt=1, sendrecv_dims
          if(tmp1r(k_cnt,1).ne.0 .and. tmp2r(k_cnt,1).ne.0) then
            kernel = ddot(n_basbas,&
                          ovlp_3ks_bse(:,i_state_loc,a_state_loc,i_spin),1,&
                          loc_ovlp_3ks_bse(:,tmp1r(k_cnt,2),tmp2r(k_cnt,2),i_spin),1)
!            print*, 'cpu, kernel, i_state, a_state, tmp1r(k_cnt,1), tmp2r(k_cnt,1)'
!            print*, myid, kernel, i_state, a_state, tmp1r(k_cnt,1), tmp2r(k_cnt,1)
! give it to par_bse
            tmp_real = bse_alpha * kernel
            if(mat_idx(i_state,a_state, occ_states, unocc_states, bse_lower_limit, bse_upper_limit) ==&
mat_idx(tmp1r(k_cnt,1),tmp2r(k_cnt,1), occ_states, unocc_states, bse_lower_limit, bse_upper_limit)) then
              tmp_real = tmp_real + qpe(a_state) - qpe(i_state)
            end if
            mpi_bcastidx(len_of_recs,1) = mat_idx(i_state,a_state, occ_states, unocc_states, &
bse_lower_limit, bse_upper_limit)
            mpi_bcastidx(len_of_recs,2) = mat_idx(tmp1r(k_cnt,1),tmp2r(k_cnt,1), occ_states, unocc_states, &
bse_lower_limit, bse_upper_limit)
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
          call infog2l(tmp1r(i_state,1), tmp1r(i_state,2), bse_desc, &
            nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)

          if(myrow.eq.iarow.and.mycol.eq.iacol) then
            par_bse(iia,jja) = tmp_vecr(i_state)
          endif

          if(tmp1r(i_state,1).ne.tmp1r(i_state,2)) then
            call infog2l(tmp1r(i_state,2), tmp1r(i_state,1), bse_desc, &
              nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)
            if(myrow.eq.iarow.and.mycol.eq.iacol) then
              par_bse(iia,jja) = tmp_vecr(i_state)
            endif
          endif
        endif
      enddo
    enddo !j_rot
!    call get_times(time_bse_rot, clock_time_bse_rot)
!    call output_times(deffmt, 'Time for this rotation', &
!                      time_bse_rot, clock_time_bse_rot)
  end do ! i_rot
!  call output_times(deffmt, 'Time for building <ia|V|jb>', &
!             time_bse_build_V, clock_time_bse_build_V)
! check par_bse
!  if(myid == 0) then
!    print*, 'myid, na_rows, na_cols', myid, na_rows, na_cols
!    open(unit = 99, file = 'par_bse_p0')
!    do i_state = 1, na_rows
!      do j_state = 1, na_cols
!        write(99, *) i_state, j_state, par_bse(i_state, j_state)
!      end do
!    end do
!    close(99)
!  end if
! deallocation
  deallocate(mpi_bcastidx)
  deallocate(mpi_vec)
  deallocate(tmp_vecr)
  deallocate(tmp1r)
  deallocate(tmp2r)
  deallocate(tmp1s)
  deallocate(tmp2s)
!  if (myid == 0) then
!    do abf_cnt = 1, n_basbas
!      do i_state = 1, ndim1_o3KS
!        do j_state = 1, ndim2_o3KS
!           print*, abf_cnt, i_state, j_state, ovlp_3ks_bse(abf_cnt, i_state, j_state, i_spin)
!        end do
!      end do
!    end do
!  end if
  call get_times(time_bse_build_V, clock_time_bse_build_V)
  call MPI_BARRIER(mpi_comm_world,errnum)
!  call aims_stop('bse_develop')
!***************************************************************
!  if(myid.eq.0) then
!    write(use_unit,*)
!    write(use_unit,*)"-------------------------------------------------"
!    write(use_unit,'(2X,A)') "Start to calculate the screened coulomb <ij|W|ab>... "
!  endif
  call get_times(time_bse_w0, clock_time_bse_w0)
  n_blk = (255/nb_aux_2d+1)*nb_aux_2d
!! get n_homo     determine the highest occupied orbital level such complication occurs when some of the orbitals are
!!     either not fully occupied or not fully empty
      occ_numbers_bse(:, :) = occ_numbers(:, :, 1)
      n_first(:) = 1
!      do i_spin = 1, n_spin
       i_spin = 1
       do i_state = 1, n_states
!       do i_state = 1, n_states_k(1)
        if (abs(occ_numbers_bse(i_state,i_spin)-dble(2/n_spin)) &
                        .lt.1.d-8) then
         n_first(i_spin)= i_state + 1
        endif
       enddo
       if(n_first(i_spin) .gt. n_states) then
         n_first(i_spin) = n_states
!       if(n_first(i_spin) .gt. n_states_k(1)) then
!         n_first(i_spin) = n_states_k(1)
       endif
!      enddo
  if(.not.allocated(n_homo)) then
     allocate(n_homo(n_spin))
     n_homo(:) = 0
     i_spin = 1
!     do i_spin = 1, n_spin, 1
       do i_state = 1, n_states, 1
!       do i_state = 1, n_states_k(1), 1
         if(occ_numbers_bse(i_state,i_spin).gt.1.e-6) then
            n_homo(i_spin) = i_state
         endif
       enddo
!      enddo
  endif
      if(myid.eq.0) then
        write(use_unit,'(2X, A,A,4I5)') &
                   "HOMO and first non-fully-occupied", &
                   " orbitals:", n_homo(:), n_first(:)
        write(use_unit,*)
      endif
      allocate(polar_freq(max(1,max_row_2d), max(1,max_col_2d)), stat=info)
      call check_allocation(info, 'polar_freq', func)
      allocate(screened_coulomb_bse(max(1,max_row_2d), max(1,max_col_2d)), stat=info)
      call check_allocation(info, 'screened_coulomb', func)

      omega_full_i_freq = 0
      n_high_state = n_states
!      n_high_state = n_states_k(1)
          call  evaluate_polarisability_freq_2 &
             ( 1, n_homo, n_first, occ_numbers_bse, &
               omega_full_i_freq, &
               KS_eigenvalue, ovlp_3ks_bse, polar_freq &
             )
!    now calculate the screened Coulomb interaction W
          call PDLASET( 'Full', n_basbas, n_basbas, 0.d0, 1.d0, screened_coulomb_bse, 1, 1, aux_sc_desc_2d )
          screened_coulomb_bse(:,:) = screened_coulomb_bse(:,:) - polar_freq(:,:)*2.d0/dble(n_spin)
          call power_auxmat_scalapack_2d(screened_coulomb_bse, -1.d0, '')
!  if (myid == 0) then
!    print*, 'write screened_coulomb_bse, cpu0, max_row_2d, max_col_2d', max_row_2d, max_col_2d
!    open(unit = 99, file = 'w_abf_p0')
!    do i_state = 1, max(1,max_row_2d)
!      do j_state = 1, max(1,max_col_2d)
!        print*, screened_coulomb_bse(i_state, j_state)
!!         print*, i_state, j_state, screened_coulomb_bse(i_state, j_state)
!      end do
!    end do
!    close(99)
!  end if
!  if (myid == 1) then
!    print*, 'write screened_coulomb_bse, cpu0, max_row_2d, max_col_2d', &
!max_row_2d, max_col_2d
!    open(unit = 99, file = 'w_abf_p1')
!    do i_state = 1, max(1,max_row_2d)
!      do j_state = 1, max(1,max_col_2d)
!        print*, screened_coulomb_bse(i_state, j_state)
!      end do
!    end do
!    close(99)
!  end if
!  if (myid == 2) then
!    print*, 'write screened_coulomb_bse, cpu0, max_row_2d, max_col_2d', &
!max_row_2d, max_col_2d
!    open(unit = 99, file = 'w_abf_p2')
!    do i_state = 1, max(1,max_row_2d)
!      do j_state = 1, max(1,max_col_2d)
!        print*, screened_coulomb_bse(i_state, j_state)
!      end do
!    end do
!    close(99)
!  end if
!  if (myid == 3) then
!    print*, 'write screened_coulomb_bse, cpu0, max_row_2d, max_col_2d', &
!max_row_2d, max_col_2d
!    open(unit = 99, file = 'w_abf_p3')
!    do i_state = 1, max(1,max_row_2d)
!      do j_state = 1, max(1,max_col_2d)
!        print*, screened_coulomb_bse(i_state, j_state)
!      end do
!    end do
!    close(99)
!  end if
  call get_times(time_bse_w0, clock_time_bse_w0)
  call MPI_BARRIER(mpi_comm_world,errnum)
! calculate W * ovlp_3ks
! ***********************************************************************
! check ovlp_3ks, global_w_row, and sc_ovlp
!  print*, 'nb_aux_2d', nb_aux_2d
  call get_times(time_bse_w_ovlp_3ks, clock_time_bse_w_ovlp_3ks)
  allocate(sc_ovlp(n_basbas, ndim1_o3KS, ndim2_o3KS, n_spin))
  sc_ovlp = 0
  allocate(global_w_row(n_basbas))
  loc_ovlp_3ks_bse = 0
!  if(myid == 0) print*, 'ovlp_3ks_bse', ovlp_3ks_bse
!      if (myid == 0) print*, 'ovlp_3ks_bse(:, i_state, j_state, i_spin)', ovlp_3ks_bse(:, i_state, j_state, i_spin)
!  if (myid == 0) print*, 'write ovlp_3ks_bse'
!  do abf_cnt = 1, n_basbas
!    do i_state = 1, ndim1_o3KS
!      do j_state = 1, ndim2_o3KS
!!         print*, myid, abf_cnt, i_state, j_state, ovlp_3ks_bse(abf_cnt, i_state, j_state, i_spin)
!      end do
!    end do
!  end do

  i_spin = 1
!  if(myid == 0) open(unit = 99, file = 'w_abf_p0')
!  if(myid == 1) open(un
  allocate(basbas2row(n_basbas), basbas2col(n_basbas), stat=info)
  call check_allocation(info, 'basbas2xxx', func)
  call get_basbas_to_rowcol(.true., basbas2row, basbas2col)
  do abf_cnt = 1, n_basbas
    global_w_row = 0
    do i_state = 1, n_basbas
      if (basbas2row(abf_cnt) * basbas2col(i_state) /= 0) then
        global_w_row(i_state) = screened_coulomb_bse(basbas2row(abf_cnt), basbas2col(i_state))
      end if
    end do
!    global_w_row(1 + myid * nb_aux_2d : (1 + myid) * nb_aux_2d) = screened_coulomb_bse(abf_cnt, 1 : nb_aux_2d)
!    global_w_row(myid * nb_aux_2d + 1:myid * nb_aux_2d + max(1, max_col_2d)) = screened_coulomb_bse(abf_cnt, 1:max(1, max_col_2d))
    call MPI_BARRIER(mpi_comm_world,errnum)
    call sync_vector(global_w_row, n_basbas)
!    if (myid == 0) then
!     do i_state = 1, n_basbas
!!       write(99, *) global_w_row(i_state)
!!       write(99, *) myid, abf_cnt, i_state, global_w_row(i_state)
!       print*, myid, abf_cnt, i_state, global_w_row(i_state)
!     enddo
!    end if
!    if (myid == 1 .and. abf_cnt == 1) then
!     do i_state = 1, n_basbas
!       print*, myid, abf_cnt, i_state, global_w_row(i_state)
!     enddo
!    end if
    do i_state = 1, ndim1_o3KS
      do j_state = 1, ndim2_o3KS
        temp = 0
        do counter = 1, n_basbas
          temp = temp + global_w_row(counter) * ovlp_3ks_bse(counter, i_state, j_state, i_spin)
        end do
        sc_ovlp(abf_cnt, i_state, j_state, i_spin) = temp
!        if (myid == 0) print*, myid, abf_cnt, i_state, j_state, temp
!        if (myid == 1) print*, myid, abf_cnt, i_state, j_state, temp
      end do
    end do
  end do ! abf_cnt
  call get_times(time_bse_w_ovlp_3ks, clock_time_bse_w_ovlp_3ks)
  call get_times(time_bse_build_W, clock_time_bse_build_W)
!  if(myid == 0) close(99)
!***************************************
! Calculating ovlp_sc_olvp
! Get local number of pairstates of (occ, occ)
  k_cnt = 0
!  do i_state = 1, occ_states
  do i_state = bse_lower_limit, occ_states
    if(own_dim1_o3ks(i_state) /= myp1_o3ks) cycle
!    do j_state = 1, occ_states
    do j_state = bse_lower_limit, occ_states
      if(own_dim2_o3ks(j_state) /= myp2_o3ks) cycle
      k_cnt = k_cnt + 1
    end do
  end do
! Find the maximal number of pairstates by sending all local number of
! pairstates to process 0 and compare
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
!*****************************
! Get local number of pairstates of (unocc, unocc)
  k_cnt = 0
!  do a_state = occ_states + 1, n_states
  do a_state = occ_states + 1, bse_upper_limit
    if(own_dim1_o3ks(a_state) /= myp1_o3ks) cycle
    do b_state = occ_states + 1, bse_upper_limit
!    do b_state = occ_states + 1, n_states
      if(own_dim2_o3ks(b_state) /= myp2_o3ks) cycle
      k_cnt = k_cnt + 1
    end do
  end do
! Find the maximal number of pairstates by sending all local number of
! pairstates to process 0 and compare
  if(myid==0) then
    do i_proc=1, n_tasks-1
      call mpi_recv(tmp_int, 1, mpi_integer, i_proc, &
        100, mpi_comm_world, mpi_status, errnum)
      if(tmp_int.gt.k_cnt) k_cnt = tmp_int
    enddo
  else
    call mpi_send(k_cnt, 1, mpi_integer, 0, 100, mpi_comm_world, errnum)
  end if
!  i_state = k_cnt
!  sendrecv_dims = max(i_state * k_cnt, 2)
!  sendrecv_dims = k_cnt * k_cnt
  sendrecv_dims = i_state * k_cnt
  call mpi_bcast(sendrecv_dims, 1, mpi_integer, 0, mpi_comm_world, errnum)
!  if(myid == 0) print*, 'cpu0, k_cnt, ', k_cnt
!  if(myid == 1) print*, 'cpu1, k_cnt', k_cnt
!  if(myid == 0) print*, 'cpu0, sendrecv_dims', sendrecv_dims
!  if(myid == 1) print*, 'cpu1, sendrecv_dims', sendrecv_dims

  allocate(mpi_bcastidx(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'mpi broadcast index for bse matrix W')
  mpi_bcastidx = 0
  allocate(mpi_vec(sendrecv_dims),stat=errnum)
  call check_allocation(errnum, 'mpi broadcast vector for bse matrix W')
  mpi_vec = 0d0
! Now calculate <ij|W|ab>
  if(myid == 0) print*, 'Building <ij|W|ab>'
  allocate(tmp1s(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'rotation index 1s W')
  allocate(tmp2s(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'rotation index 2s W')
  tmp1s(:,1) = 0
  tmp2s(:,1) = 0
  k_cnt = 1
  do i_state = bse_lower_limit, occ_states
!  do i_state = 1, occ_states
    if(own_dim1_o3ks(i_state) /= myp1_o3ks) cycle
    do j_state = bse_lower_limit, occ_states
!    do j_state = 1, occ_states
      if(own_dim2_o3ks(j_state) /= myp2_o3ks) cycle
      i_state_loc = loc_dim1_o3ks(i_state)
      j_state_loc = loc_dim2_o3ks(j_state)
      tmp1s(k_cnt,1) = i_state
      tmp1s(k_cnt,2) = i_state_loc
      tmp2s(k_cnt,1) = j_state
      tmp2s(k_cnt,2) = j_state_loc
      k_cnt = k_cnt + 1
    end do
  end do
! check tmp1s, tmp2s
!  if(myid == 0) then
!    print*, 'cpu 0'
!    do k_cnt = 1, sendrecv_dims
!      print*, 'cpu, tmp1s(k_cnt, :) , tmp2s(k_cnt, :), k_cnt', k_cnt
!      print*, myid, tmp1s(k_cnt,1), tmp1s(k_cnt,2), tmp2s(k_cnt,1), tmp2s(k_cnt,2)
!    end do
!  end if
!  if(myid == 1) then
!    print*, 'cpu 1'
!    do k_cnt = 1, sendrecv_dims
!      print*, 'cpu, tmp1s(k_cnt, :) , tmp2s(k_cnt, :), k_cnt', k_cnt
!      print*, myid, tmp1s(k_cnt,1), tmp1s(k_cnt,2), tmp2s(k_cnt,1), tmp2s(k_cnt,2)
!    end do
!  end if
  if(myid==0) print*, 'Calculate W or ovlp_sc_ovlp from ovlp and sc_ovlp'

  allocate(tmp1r(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, ' tmp1r res dist W')
  allocate(tmp2r(sendrecv_dims,2),stat=errnum)
  call check_allocation(errnum, 'rotation index 2r W')
  allocate(tmp_vecr(sendrecv_dims),stat=errnum)
  call check_allocation(errnum, ' tmp_vecr res dist W')
  allocate(par_bse_w(na_rows,na_cols),stat=errnum)
  call check_allocation(errnum, 'parallel bse matrix W')
  par_bse_w = 0d0

  i_spin = 1
!  if (myid == 0) then
!    print*, 'write ovlp_3ks_bse'
!    open(unit = 99, file = 'ovlp_3ks_bse_p0')
!    do abf_cnt = 1, n_basbas
!      do i_state = 1, ndim1_o3KS
!        do j_state = 1, ndim2_o3KS
!           write(99, *) ovlp_3ks_bse(abf_cnt, i_state, j_state, i_spin)
!!           print*, abf_cnt, i_state, j_state, ovlp_3ks_bse(abf_cnt, i_state, j_state, i_spin)
!        end do
!      end do
!    end do
!    close(99)
!  end if
!  if (myid == 1) then
!    print*, 'write ovlp_3ks_bse'
!    open(unit = 99, file = 'ovlp_3ks_bse_p1')
!    do abf_cnt = 1, n_basbas
!      do i_state = 1, ndim1_o3KS
!        do j_state = 1, ndim2_o3KS
!           write(99, *) ovlp_3ks_bse(abf_cnt, i_state, j_state, i_spin)
!!           print*, abf_cnt, i_state, j_state, ovlp_3ks_bse(abf_cnt, i_state,
!!           j_state, i_spin)
!        end do
!      end do
!    end do
!    close(99)
!  end if
  do i_rot=1, n_tasks
!    call get_times(time_bse_rot, clock_time_bse_rot)
    if(myid==0) write(use_unit,'(2x,a,I5.1)') '| starting rotation: ', i_rot
    if(myid==i_rot-1) then
      tmp1r = tmp1s
      tmp2r = tmp2s
      loc_ovlp_3ks_bse = ovlp_3ks_bse
    endif
    call mpi_bcast(tmp1r, 2*sendrecv_dims, mpi_integer, i_rot-1,&
          mpi_comm_world, errnum)
    call mpi_bcast(tmp2r, 2*sendrecv_dims, mpi_integer, i_rot-1,&
          mpi_comm_world, errnum)
    call mpi_bcast(loc_ovlp_3ks_bse, n_basbas*ndim1_o3KS*ndim2_o3KS,&
          mpi_double_precision, i_rot-1, mpi_comm_world, errnum)
    len_of_recs = 1
    i_spin = 1
  do a_state = occ_states + 1, bse_upper_limit
!    do a_state=occ_states+1, n_states
      if(own_dim1_o3KS(a_state) /= myp1_o3KS) cycle
    do b_state = occ_states + 1, bse_upper_limit
!      do b_state=occ_states+1, n_states
        if(own_dim2_o3KS(b_state) /= myp2_o3KS) cycle
        a_state_loc = loc_dim1_o3KS(a_state)
        b_state_loc = loc_dim2_o3KS(b_state)
        do k_cnt=1, sendrecv_dims
          if(tmp1r(k_cnt,1).ne.0 .and. tmp2r(k_cnt,1).ne.0) then
!            if(myid == 0) print*, 'myid, tmp1r(k_cnt,1), tmp2r(k_cnt,1), a_state, b_state, ovlp, sc_ovlp', &
!                    myid, tmp1r(k_cnt,1), tmp2r(k_cnt,1), a_state, b_state,&
!                    ovlp_3ks_bse(:,tmp1r(k_cnt,2),tmp2r(k_cnt,2),i_spin), &
!                    loc_ovlp_3ks_bse(:,a_state_loc, b_state_loc, i_spin)
            kernel = ddot(n_basbas,loc_ovlp_3ks_bse(:,tmp1r(k_cnt,2),tmp2r(k_cnt,2),i_spin), 1, &
                     sc_ovlp(:,a_state_loc, b_state_loc, i_spin),1)
!            if(myid == 0) print*, 'myid, tmp1r(k_cnt,1), tmp2r(k_cnt,1), a_state, b_state, ovlp, sc_ovlp, kernel', &
!            if(myid == 1) then
!              print*, 'myid, tmp1r(k_cnt,1), tmp2r(k_cnt,1), a_state, b_state, ovlp, sc_ovlp, kernel', &
!                    myid, tmp1r(k_cnt,1), tmp2r(k_cnt,1), a_state, b_state,&
!                    loc_ovlp_3ks_bse(:,tmp1r(k_cnt,2),tmp2r(k_cnt,2),i_spin), &
!                    sc_ovlp(:,a_state_loc, b_state_loc, i_spin), kernel
!              print*, 'myid, tmp1r(k_cnt,2),tmp2r(k_cnt,2), a_state_loc, b_state_loc',&
!                      myid, tmp1r(k_cnt,2),tmp2r(k_cnt,2), a_state_loc, b_state_loc
!              print*, '____________________________________________'
!            end if

!            print*, 'myid, kernel', myid, kernel
!            if(abs(kernel) > 1e-6) print*, 'cpu, kernel, tmp1r(k_cnt,1), tmp2r(k_cnt,1), a_state, b_state', &
!                    myid, kernel, tmp1r(k_cnt,1), tmp2r(k_cnt,1), a_state, b_state
            tmp_real = (-1) * kernel
! get bse_mat index
!            print*, 'myid', myid
!            print*, 'myid, mat_idx', myid, mat_idx(tmp1r(k_cnt,1), a_state), mat_idx(tmp2r(k_cnt,1), b_state)
!            print*, mat_idx(tmp1r(k_cnt,1), a_state), mat_idx(tmp2r(k_cnt,1), b_state)
            mpi_bcastidx(len_of_recs,1) = mat_idx(tmp1r(k_cnt,1), a_state, occ_states, unocc_states, &
bse_lower_limit, bse_upper_limit) ! i, a
            mpi_bcastidx(len_of_recs,2) = mat_idx(tmp2r(k_cnt,1), b_state, occ_states, unocc_states, &
bse_lower_limit, bse_upper_limit) ! j, b
            mpi_vec(len_of_recs) = tmp_real
            len_of_recs = len_of_recs + 1
          endif
        enddo
      enddo
    enddo

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
          call infog2l(tmp1r(i_state,1), tmp1r(i_state,2), bse_desc, &
            nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)

          if(myrow.eq.iarow.and.mycol.eq.iacol) then
            par_bse_w(iia,jja) = tmp_vecr(i_state)
          endif

          if(tmp1r(i_state,1).ne.tmp1r(i_state,2)) then
            call infog2l(tmp1r(i_state,2), tmp1r(i_state,1), bse_desc, &
              nprow, npcol, myrow, mycol, iia, jja, iarow, iacol)
            if(myrow.eq.iarow.and.mycol.eq.iacol) then
              par_bse_w(iia,jja) = tmp_vecr(i_state)
            endif
          endif
        endif
      enddo
    enddo !j_rot
!    call get_times(time_bse_rot, clock_time_bse_rot)
!    call output_times(deffmt, 'Time for this rotation', &
!                      time_bse_rot, clock_time_bse_rot)
  end do ! i_rot
  call get_times(time_bse_build_W, clock_time_bse_build_W)
!  call output_times(deffmt, 'Time for building <ij|W|ab>', &
!             time_bse_build_W, clock_time_bse_build_W)
! check par_bse_w
!  if(myid == 0) then
!  open(unit = 99, file = 'par_bse_w_p0')
!  do i_state = 1, na_rows
!    do j_state = 1, na_cols
!      write(99, *) par_bse_w(i_state, j_state)
!    end do
!  end do
!  close(99)
!  end if
!  if(myid == 1) then
!  open(unit = 99, file = 'par_bse_w_p1')
!  do i_state = 1, na_rows
!    do j_state = 1, na_cols
!      write(99, *) par_bse_w(i_state, j_state)
!    end do
!  end do
!  close(99)
!  end if
!  if(myid == 2) then
!  open(unit = 99, file = 'par_bse_w_p2')
!  do i_state = 1, na_rows
!    do j_state = 1, na_cols
!      write(99, *) par_bse_w(i_state, j_state)
!    end do
!  end do
!  close(99)
!  end if
!  if(myid == 3) then
!  open(unit = 99, file = 'par_bse_w_p3')
!  do i_state = 1, na_rows
!    do j_state = 1, na_cols
!      write(99, *) par_bse_w(i_state, j_state)
!    end do
!  end do
!  close(99)
!  end if

!  do i_state = 1, na_rows
!    do j_state = 1, na_cols
!      if(abs(par_bse_w(i_state, j_state)) > 1e-6) print*, 'myid, par_bse_w', myid, i_state, j_state, par_bse_w(i_state, j_state)
!    end do
!  end do
! add par_bse and par_bse_w
  par_bse = par_bse + par_bse_w
! check par_bse
!  do i_state = 1, na_rows
!    do j_state = 1, na_cols
!      if(abs(par_bse(i_state, j_state)) > 1e-6) print*, 'myid, par_bse',  myid, i_state, j_state, par_bse(i_state, j_state)
!    end do
!  end do
  call MPI_BARRIER(mpi_comm_world,errnum)
!  call aims_stop()
  call get_times(time_bse_construct_mat, clock_time_bse_construct_mat)

! SOLVE BSE by ELPA

  allocate(par_evecs(na_rows, na_cols),stat=errnum)
  call check_allocation(errnum,'bse parallel eigenvectors')
  par_evecs = 0d0
  allocate(par_evals(n_pair_states),stat=errnum)
  call check_allocation(errnum,'bse parallel eigenvalues')
  par_evals = 0d0

  ! Solve Casida eigenvalue problem OMEGA A_q = w_q^2 A_q (Equation 1 in paper)
  call get_timestamps(time_bse_solve_mat, clock_time_bse_solve_mat)
  if(myid==0) then
     write(use_unit,'(2x,a)') '|'
     write(use_unit,'(2x,a)') '| Calling ELPA 2stage solver now'
     write(use_unit,'(2x,a)') '|'
  endif

  call aims_elsi_stdevp(n_pair_states, na_rows, na_cols, n_pair_states, &
          par_bse, par_evals, par_evecs, my_bse_ctxt, num_blocks, &
          mpi_comm_global)
  call get_times(time_bse_solve_mat, clock_time_bse_solve_mat)
!  call output_times(deffmt, 'Time for solving BSE matrix', &
!    time_bse_solve_mat, clock_time_bse_solve_mat)
! check eigenvalues
!  if(myid == 0) print*, 'eigenvalues TDA BSE(lowest 20):', par_evals(:20) * hartree

! Oscillator strength
  call get_timestamps(time_bse_osc, clock_time_bse_osc)
  allocate(coord_of_center(3),stat=errnum)
  call check_allocation(errnum,' coord_of_center allocate')
  allocate(atom_coord_wrt_center(3,n_atoms),stat=errnum)
  call check_allocation(errnum,' atom_coord_of_center allocate')
  call determine_center_of_molecule(coord_of_center, atom_coord_wrt_center)

      allocate(dipole_moment(occ_states,unocc_states,n_spin,3),stat=errnum)
      call check_allocation(errnum,'casida dipole_moment')
!      call integrate_dipmom_pairstates(occ_states, unocc_states, n_first, &
!             coord_of_center, l_shell_max, KS_eigenvector(:,:,:,1), &
!             dipole_moment)
      call integrate_fxc_pairstates(occ_states, unocc_states, coord_of_center, &
              l_shell_max, KS_eigenvector(:,:,:,1), dipole_moment)
!              l_shell_max, KS_eigenvector(:,:,:,i_k_point), dipole_moment)

  allocate(osc_str(n_pair_states),stat=errnum)
  call check_allocation(errnum,'Casida Oscillator Strengths allocation')
  allocate(tmom(n_pair_states,3),stat=errnum)
  call check_allocation(errnum,'moment buffer allocation')

  if(bse_singlet) then
    tmom = 0d0
    do i_spin = 1, n_spin
      ! Loop over all pair states
      do i_state=1, n_pair_states
        ! Loop over all occupied states in energy window
        do j_state=bse_lower_limit, occ_states
  !      do j_state=low_limit, occ_states
          ! Loop over all unoccupied states in energy window
          do b_state=occ_states+1, bse_upper_limit
  !        do b_state=occ_states+1, n_states-unocc_limit

            ! Get offset in global matrix
            call infog2l(mat_idx(j_state,b_state,occ_states, unocc_states, &
              bse_lower_limit, bse_upper_limit), i_state, bse_desc, nprow, &
              npcol, myrow, mycol, iia, jja, iarow, iacol)

            if(myrow.eq.iarow .and. mycol.eq.iacol) then

              ! S^{1/2} * A_q = sqrt( 2 * (e_unocc - e_occ) ) * A_q
              tmp_real = dsqrt(2d0*par_evals(i_state)) * par_evecs(iia,jja)
  !            tmp_real = dsqrt(2d0*par_evals(iia)) * par_evecs(iia,jja)
  !            tmp_real = dsqrt(2d0*(ks_eigenvalue(b_state,i_spin,n_k_points) &
  !              - ks_eigenvalue(j_state,i_spin,n_k_points))) * par_evecs(iia,jja)

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

    if (myid == 0) then
      print*, "sum osc_str", sum(osc_str)
      print*, "n_elec", n_electrons
    end if
  else
    osc_str = 0
  end if
  call get_times(time_bse_osc, clock_time_bse_osc)
! end oscillator strength


      write (info_str,'(2X,A,A)') &
           "Eigenvalues TDA BSE."
      call localorb_info(info_str,use_unit)
  write(info_str,'(2X,A,4X,A,4X,A,4X,A)') &
       "State", "Eigenvalue [Ha]", "Eigenvalue [eV]", "Oscillator Strength"
  call localorb_info(info_str,use_unit)
!               write(info_str,'(2X,I5,6X,F8.5,5X,F14.6,4X,F15.5)') &
!                    i_state,  &
!                    (KS_eigenvalue(i_state, i_spin, i_k_point)), &
!                    (KS_eigenvalue(i_state, i_spin, i_k_point)*hartree)
!               call localorb_info(info_str,use_unit,'(A)',OL_low)

!  if(myid == 0) then
!    print*, 'eigenvalues TDA BSE(lowest 20):'
    do i_state = 1, 20
               write(info_str,'(2X,I5,5X,F14.6,4X,F15.6, 4X, f12.6)') &
                    i_state, par_evals(i_state), par_evals(i_state) * hartree, &
                    osc_str(i_state)
               call localorb_info(info_str,use_unit)

!      print*, par_evals(i_state) * hartree
    end do
write(info_str,'(A)') SEPARATORLINE
          call localorb_info( info_str, use_unit )
!  end if
!  call get_times(time_omega_solve, clock_time_omega_solve)

!  deallocate(par_omega)
!    do i_state=1, occ_states
!      if(own_dim1_o3KS(i_state) /= myp1_o3KS) cycle
!      do a_state=occ_states+1, n_states
!        if(own_dim2_o3KS(a_state) /= myp2_o3KS) cycle
!        i_state_loc = loc_dim1_o3KS(i_state)
!        a_state_loc = loc_dim2_o3KS(a_state)
!        do k_cnt=1, sendrecv_dims
!          if(tmp1r(k_cnt,1).ne.0 .and. tmp2r(k_cnt,1).ne.0) then
!            kernel = ddot(n_basbas,&
!                          ovlp_3ks_bse(:,i_state_loc,a_state_loc,i_spin),1,&
!                          loc_ovlp_3ks_bse(:,tmp1r(k_cnt,2),tmp2r(k_cnt,2),i_spin),1)
!            print*, 'cpu, kernel, i_state, a_state, tmp1r(k_cnt,1), tmp2r(k_cnt,1)'
!            print*, myid, kernel, i_state, a_state, tmp1r(k_cnt,1), tmp2r(k_cnt,1)
!          endif
!        enddo
!      enddo
!    enddo
!  end do
! write sc_ovlp
  i_spin = 1
!  if(myid == 0) then
!    open(unit = 99, file = 'sc_ovlp_proc0')
!      abf_cnt = 1
!!    do abf_cnt = 1, n_basbas
!      do i_state = 1, ndim1_o3KS
!        do j_state = 1, ndim2_o3KS
!          write(99, *) sc_ovlp(abf_cnt, i_state, j_state, i_spin)
!!          write(99, *) abf_cnt, i_state, j_state, sc_ovlp(abf_cnt, i_state, j_state, i_spin)
!        end do
!      end do
!!    end do
!    close(99)
!  end if
!n_basis_pairs,n_loc_prodbas
!  if (myid == 0) then
!    print*, 'write ovlp_3fn_bse'
!    do i_state = 1, n_basis_pairs
!      do j_state = 1, n_loc_prodbas
!        print*, i_state, j_state, ovlp_3fn_bse(i_state, j_state)
!      end do
!    end do
!  end if
!  o3ks_for_w = 0
!  call transform_ovlp3fn_2(n_states,KS_eigenvector, ovlp_3fn_bse, o3ks_for_w)
!  deallocate(global_w_row)
!  deallocate(ovlp_3fn_bse)
!  deallocate(ovlp_3ks_bse)
  deallocate(loc_ovlp_3ks_bse)
  if (allocated(coord_of_center)) deallocate(coord_of_center)
  if (allocated(atom_coord_wrt_center)) deallocate(atom_coord_wrt_center)
  if (allocated(dipole_moment)) deallocate(dipole_moment)
  if (allocated(osc_str)) deallocate(osc_str)
  if (allocated(tmom)) deallocate(tmom)
end subroutine

pure integer function mat_idx(i,a, occ_states, unocc_states, bse_lower_limit, bse_upper_limit) result(idx)
 implicit none
 integer, intent(in) :: i, a, occ_states, unocc_states, bse_lower_limit, bse_upper_limit
! low_limit = 1, unocc_limit = 0
 idx = (i - bse_lower_limit) * (bse_upper_limit - occ_states) + a - occ_states
! idx = (i-low_limit)*(unocc_states-unocc_limit)+(a-occ_states)
! idx = (i - 1) * unocc_states + a - occ_states
end function mat_idx

subroutine integrate_fxc_pairstates(n_occ, n_unocc, coord_of_center, basis_l_max, KS_eigenvector, dipole_mom)
!       COPY of integrate_dipmom_pairstates.f90
  use mpi_utilities
  use localorb_io
  use runtime_choices
  use dimensions
!  use physics
  use prodbas
  use gw_para
  use basis, only : basis_atom
  use grids, only : n_radial, n_angular, r_radial, r_angular, &
                    w_radial, w_angular
  use geometry, only: species, coords
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

      end do ! i_angular

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

!  call get_times(time_neutral_excitation_dipole, clock_time_neutral_excitation_dipole)

!    enddo ! i_radial
!!       end radial integration loop
!
!!     end integration loop over atoms
!  enddo
end subroutine integrate_fxc_pairstates
