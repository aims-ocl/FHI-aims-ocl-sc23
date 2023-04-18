!****h* FHI-aims/cg_scalapack
!  NAME
!    cg
!  SYNOPSIS

module cg_scalapack

!  PURPOSE
!  ??????????????????????????????
!
!  USES

  implicit none

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

  ! Variables for LOPBCG
  real*8, allocatable, private :: ovlp_save(:,:)
  real*8, allocatable, private :: cg_prec(:,:)

  real*8, allocatable, dimension(:,:), private :: another_tmp, tmp

  integer, parameter, private :: dlen_ = 9
  integer, dimension(dlen_), private :: cg_desc
  integer, dimension(dlen_), private :: rz_desc

  integer, private :: mxld, mxcol, npcol, nprow, nb, mb, rsrc, csrc
  integer, private :: my_scalapack_id, my_scalapack_comm_work, my_blacs_ctxt
  integer, private :: mpi_comm_rows, mpi_comm_cols

  integer, dimension(:), allocatable, private :: l_row
  integer, dimension(:), allocatable, private :: l_col

  integer, private :: info

  logical, private :: lopbcg_debug = .false.
  logical, private :: print_times = .false.
  real*8, private :: tt0, tt1, tt2, tt3
  integer, private :: n_cg_cycle = 0

contains

!-----------------------------------------------------------------------------------
!****s* cg_scalapack/initialize_cg_scalapack
!  NAME
!    initialize_cg_scalapack
!  SYNOPSIS
  subroutine initialize_cg_scalapack( ovlp, sc_desc, ext_mxld, ext_mxcol, ext_npcol, ext_nprow, ext_nb, ext_mb, &
       ext_rsrc, ext_csrc, ext_my_blacs_ctxt, &
       ext_my_scalapack_id, ext_my_scalapack_comm_work, ext_mpi_comm_rows, ext_mpi_comm_cols, &
       ext_l_col, ext_l_row, ext_dlen_ )
!  PURPOSE
!    Initializes the LOPCG-method with ScaLAPACK.
!  USES
    use localorb_io, only: localorb_info
    use dimensions, only: n_basis
    use mpi_tasks, only: aims_stop, check_allocation
!  ARGUMENTS
    integer :: ext_mxld, ext_mxcol, ext_npcol, ext_nprow, ext_nb, ext_mb, ext_rsrc, ext_csrc
    integer :: ext_dlen_
    real*8, dimension(ext_mxld,ext_mxcol) :: ovlp
    integer, dimension(ext_dlen_) :: sc_desc
    integer :: ext_my_scalapack_id, ext_my_scalapack_comm_work, ext_my_blacs_ctxt
    integer :: ext_mpi_comm_rows, ext_mpi_comm_cols
    integer, dimension(n_basis) :: ext_l_col, ext_l_row
!  INPUTS
!    o ovlp -- the overlap matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_elem
    real*8 :: shift, shift_all

    if (dlen_ /= ext_dlen_) call aims_stop

    mxld = ext_mxld
    mxcol = ext_mxcol
    npcol = ext_npcol
    nprow = ext_nprow
    nb = ext_nb
    mb = ext_mb
    rsrc = ext_rsrc
    csrc = ext_csrc
    my_blacs_ctxt = ext_my_blacs_ctxt
    my_scalapack_id = ext_my_scalapack_id
    my_scalapack_comm_work = ext_my_scalapack_comm_work
    mpi_comm_rows = ext_mpi_comm_rows
    mpi_comm_cols = ext_mpi_comm_cols

    cg_desc = sc_desc

    if (.not.allocated(l_row)) allocate(l_row(n_basis),stat=info)
    call check_allocation(info, 'l_row                         ')
    if (.not.allocated(l_col)) allocate(l_col(n_basis),stat=info)
    call check_allocation(info, 'l_col                         ')

    l_col = ext_l_col
    l_row = ext_l_row

    if (.not.(allocated(ovlp_save))) then
       allocate(ovlp_save(mxld,mxcol),stat=info)
       call check_allocation(info, 'ovlp_save                        ')    
    end if
    ovlp_save = ovlp
    
    if (.not.allocated(tmp)) then
       allocate(tmp(mxld,mxcol),stat=info)
       call check_allocation(info, 'tmp                   ')    
    end if

    call cg_set_full_matrix_real(ovlp_save)

    
  end subroutine initialize_cg_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* cg_scalapack/initialize_cg_scalapack
!  NAME
!    initialize_cg_scalapack
!  SYNOPSIS
  subroutine initialize_cg_prec_scalapack( ham, KS_eigenvalue )
!  PURPOSE
!    Initializes the LOPCG-method with ScaLAPACK.
!  USES
    use dimensions, only: n_basis, n_states
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks, only: check_allocation
    use runtime_choices, only: prec_inv_ham_ovlp, prec_inv_ovlp, &
        prec_diagonal, cg_preconditioner, lopcg_omega
    use synchronize_mpi_basic, only: get_min_double
!  ARGUMENTS
    implicit none
    real*8, dimension(mxld,mxcol) :: ham
    real*8, dimension(n_states) :: KS_eigenvalue
!  INPUTS
!    o ovlp -- the overlap matrix
!    o ham -- the hamiltonian matrix
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer :: i_elem
    real*8 :: shift, shift_all, alpha
    character*120 :: info_str

    write (info_str,'(2X,A)') "Constructing LOPCG preconidtioner."
    call localorb_info(info_str,use_unit,'(A)')

    if (.not.(allocated(cg_prec))) then
       allocate(cg_prec(mxld,mxcol),stat=info)
       call check_allocation(info, 'cg_prec                        ')    
    end if
    !! cg_prec(:,:) = ham(:,:) + 50.0d0*ovlp_save(:,:)
    select case (cg_preconditioner)
    case (prec_inv_ham_ovlp)
       if (MINVAL(KS_eigenvalue) < 0) then
          shift = - lopcg_omega*MINVAL(KS_eigenvalue)
       else
          shift = 0.0d0
       end if
       cg_prec(:,:) = ham(:,:) + shift*ovlp_save(:,:)
    case (prec_inv_ovlp)
       cg_prec(:,:) = ovlp_save(:,:)
    case (prec_diagonal)
       shift = HUGE(shift)
       do i_elem = 1, n_basis, 1
          if ((l_col(i_elem) > 0) .and. (l_row(i_elem) > 0)) then
             shift = MIN(shift, ham(l_row(i_elem),l_col(i_elem)))
          end if
       end do
       if (shift > 0.0d0) shift = 0.0d0
       call get_min_double(shift_all, shift)

       do i_elem = 1, n_basis, 1
          if ((l_col(i_elem) > 0) .and. (l_row(i_elem) > 0)) then
             cg_prec(l_row(i_elem),l_col(i_elem)) = ham(l_row(i_elem),l_col(i_elem)) - 1.1d0*shift_all
          end if
       end do
!!$       do i_elem = 1, n_basis, 1
!!$          if ((l_col(i_elem) > 0) .and. (l_row(i_elem) > 0)) then
!!$             print *, i_elem, cg_prec(l_row(i_elem),l_col(i_elem))
!!$          end if
!!$       end do

    case default
       do i_elem = 1, n_basis, 1
          if ((l_col(i_elem) > 0) .and. (l_row(i_elem) > 0)) then
             cg_prec(l_row(i_elem),l_col(i_elem)) = 1.0d0
          end if
       end do
    end select

    CALL PDPOTRF('U', n_basis, cg_prec, 1, 1, cg_desc, info )
    if(info /= 0) call cg_scalapack_err_exit(info,"PDPOTRF")
    CALL PDPOTRI('U', n_basis, cg_prec, 1, 1, cg_desc, info )
    if(info /= 0) call cg_scalapack_err_exit(info,"PDPOTRI")
    call cg_set_full_matrix_real(cg_prec(1,1))

!!$             ! For some reason PDPOTRF sometimes fails for small blocksizes
!!$             ! We use therefore own routines
!!$         
!!$             call cholesky_real(n_basis, cg_prec, mxld, nb, mpi_comm_rows, mpi_comm_cols)
!!$             call invert_trm_real(n_basis, cg_prec, mxld, nb, mpi_comm_rows, mpi_comm_cols)

  end subroutine initialize_cg_prec_scalapack
!-----------------------------------------------------------------------------------
!****s* cg_scalapack/cg_solver_scalapack
!  NAME
!    cg_solver_scalapack
!  SYNOPSIS
  subroutine cg_solver_scalapack( n_basis_current, n_states_current, KS_eigenvalue, &
       ham, eigenvec, i_spin, converged)
!  PURPOSE
!    Solves the eigenvalue problem using LOPCG-method with ScaLAPACK.
!  USES
    use dimensions, only: n_basis, n_spin
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks
    use physics, only: occ_numbers, rho_change
    use runtime_choices, only: lopcg_adaptive_tolerance, use_ritz_scalapack, &
        max_cg_iterations, lopcg_use_all_evs, lopcg_tol, lopcg_auto_blocksize, &
        lopcg_block_size, lopcg_skip_rate, lopcg_slide_ev_block
    use synchronize_mpi_basic, only: sync_vector, get_max_double
!  ARGUMENTS
    integer :: n_basis_current, n_states_current
    real*8 KS_eigenvalue (n_states_current)
    integer :: i_spin
    real*8, dimension(mxld,mxcol) :: ham, eigenvec
    logical :: converged
!  INPUTS
!    o n_basis_current -- current number of basis functions
!    o n_states_current -- current number of states
!    o i_spin -- the spin channel
!  OUTPUT
!    o KS_eigenvalue -- Kohn-Sham eigenvalues
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer :: max_ritz_dim, lwork, liwork
    integer :: i_state, i_iter, i_vector, i_row, i_col
    integer :: current_block_size, current_ritz_dim, current_state, shift
    integer :: state_top, state_bottom, n_blocks, total_iter
    integer :: n_occupied_states
    integer :: n_active_states, i_active_state, n_converged_states, n_slide
    integer, dimension(n_states_current) :: occupied_states
    integer :: state_offset

    real*8, allocatable, dimension(:,:) :: ritz_matrix_H, ritz_matrix_S
    real*8, allocatable, dimension(:,:) :: ritz_basis
    real*8, allocatable, dimension(:,:), save :: conj_dir
    real*8, allocatable, dimension(:,:) :: res
    real*8, allocatable, dimension(:,:) :: prec_res
    real*8, allocatable, dimension(:,:) :: Hev, Sev
    real*8, allocatable, dimension(:,:) :: converged_evs
    real*8 :: lopcg_tol_current, temp_norm, temp_dot
    real*8 :: residual_norm(n_states_current)

    logical :: converged_state(n_states_current)
    logical :: slide_block

    real*8, allocatable, dimension(:) :: ritz_values
    real*8, allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork

    integer :: np0, nq0, trilwmin, lwormtr
    integer, external :: numroc

    character*120 :: info_str

    ! This routine is called only from the working set, so no need to check here

    n_cg_cycle = n_cg_cycle + 1
    if (lopcg_skip_rate > 0) then
       if (MOD(n_cg_cycle,lopcg_skip_rate) == 0) then
          write (info_str,'(2X,A,I4,A)') &
               "This is LOPCG cycle ", n_cg_cycle, ". Skipping as requested."
          call localorb_info(info_str,use_unit,'(A)')
          converged = .false.
          return
       end if
    end if

!test
    if ((myid == 0).and.lopbcg_debug) then
       print *, 'ON ENTRY:'
       do i_state = 1, n_states_current, 1
          print *, 'state: ', i_state, ' EV: ', KS_eigenvalue(i_state)
       end do
    end if
!test

    if (lopcg_use_all_evs) then
       max_ritz_dim = n_states_current + 2*lopcg_block_size
    else
       max_ritz_dim = 3*lopcg_block_size
    end if
    
    if (use_ritz_scalapack) then
      liwork = MAX(7*n_basis + 8*NPCOL + 2, n_basis + 2*NB + 2*npcol)

      np0 = NUMROC( MAX(n_basis_current,nb,2), nb, 0, 0, nprow )
      nq0 = NUMROC( MAX(n_basis_current,nb,2), nb, 0, 0, npcol )
      TRILWMIN = 3*n_basis + MAX( NB*( NP0+1 ), 3*NB )
      lwormtr = MAX( (NB*(NB-1))/2, (np0 + nq0)*NB + 2*NB*NB)
      lwork = MAX( 1+6*n_basis_current+2*NP0*NQ0, TRILWMIN, lwormtr ) + 2*n_basis_current
    else
       lwork = 1 + 6*max_ritz_dim + 2*max_ritz_dim*max_ritz_dim
       liwork = 3 + 5*max_ritz_dim
    end if


    if (.not.use_ritz_scalapack) then
       if (.not.allocated(ritz_matrix_H)) then
          allocate(ritz_matrix_H(max_ritz_dim,max_ritz_dim),stat=info)
          call check_allocation(info, 'ritz_matrix_H                 ')    
       end if
       if (.not.allocated(ritz_matrix_S)) then
          allocate(ritz_matrix_S(max_ritz_dim,max_ritz_dim),stat=info)
          call check_allocation(info, 'ritz_matrix_S                 ')    
       end if
    end if

    if (.not.allocated(ritz_basis)) then
       allocate(ritz_basis(mxld,mxcol),stat=info)
       call check_allocation(info, 'ritz_basis                    ')    
    end if
    if (.not.allocated(conj_dir)) then
       allocate(conj_dir(mxld,mxcol),stat=info)
       call check_allocation(info, 'conj_dir                      ')    
    end if
    if (.not.allocated(another_tmp)) then
       allocate(another_tmp(mxld,mxcol),stat=info)
       call check_allocation(info, 'another_tmp                   ')    
    end if
    if (.not.allocated(res)) then
       allocate(res(mxld,mxcol),stat=info)
       call check_allocation(info, 'res                           ')    
    end if
    if (.not.allocated(prec_res)) then
       allocate(prec_res(mxld,mxcol),stat=info)
       call check_allocation(info, 'prec_res                      ')    
    end if
    if (.not.allocated(converged_evs)) then
       allocate(converged_evs(mxld,mxcol),stat=info)
       call check_allocation(info, 'converged_evs                 ')    
    end if
    if (.not.allocated(Hev)) then
       allocate(Hev(mxld,mxcol),stat=info)
       call check_allocation(info, 'Hev                           ')    
    end if
    if (.not.allocated(Sev)) then
       allocate(Sev(mxld,mxcol),stat=info)
       call check_allocation(info, 'Sev                           ')    
    end if
    if (.not.allocated(ritz_values)) then
       allocate(ritz_values(max_ritz_dim),stat=info)
       call check_allocation(info, 'ritz_values                   ')    
    end if
    if (.not.allocated(work)) then
       allocate(work(lwork),stat=info)
       call check_allocation(info, 'work                          ')    
    end if
    if (.not.allocated(iwork)) then
       allocate(iwork(liwork),stat=info)
       call check_allocation(info, 'iwork                          ')    
    end if

    write (info_str,'(2X,A,I4,A,I4)') &
         "Solving the eigenvalue problem with ScaLAPACK+LOPCG method for states ", &
         1, " - ", n_states_current
    call localorb_info(info_str,use_unit,'(A)')
!!$    write (info_str,'(2X,A,I4)') &
!!$         "Maximum Ritz dimension ", max_ritz_dim
!!$    call localorb_info(info_str,use_unit,'(A)')

!!$    min_ev = 1.0d6
!!$    do i_state = 1, n_states_current, 1
!!$       min_ev = MIN(min_ev, KS_eigenvalue(i_state))
!!$    end do
!!$    if (min_ev < 0.0d0) then
!!$       ev_shift = -min_ev + 0.1
!!$    else
!!$       ev_shift = 0.0d0
!!$    end if

    if (lopcg_adaptive_tolerance) then
       lopcg_tol_current = MAX(0.01*(SUM(ABS(rho_change))/dble(n_spin)),lopcg_tol)
       write (info_str,'(2X,A,E10.4)') "Current LOPCG tolerance: ", lopcg_tol_current
       call localorb_info(info_str,use_unit,'(A)')
    else
       lopcg_tol_current = lopcg_tol
    end if

    n_occupied_states = 0       
    do i_state = 1, n_states_current, 1
       if (occ_numbers(i_state,i_spin,1) > 0.0d0) then
          n_occupied_states = n_occupied_states + 1
          occupied_states(n_occupied_states) = i_state
       end if
    end do

    total_iter = 0
    converged_state = .false.

    tt0 = MPI_Wtime()

    ! compute ham*eigenvec as a first part of the gradient
!!$    call PDGEMM('N','N', n_basis_current, n_states_current, n_basis_current, 1.0d0, &
!!$         ham(1,1), 1, 1, cg_desc, eigenvec(1,1), 1, 1, cg_desc, 0.0d0, Hev, 1, 1, cg_desc)
    call PDSYMM('L', 'U', n_basis_current, n_states_current, 1.0d0, ham, 1, 1, cg_desc, &
         eigenvec, 1, 1, cg_desc, 0.0d0, Hev, 1, 1, cg_desc)
!!$    call mult_at_b_real_elpa('N','N', n_basis_current, n_states_current, ham, mxld, eigenvec, mxld, &
!!$         nb, mpi_comm_rows, mpi_comm_cols, Hev, mxld)
    
    ! compute ovlp*eigenvec as a second part of the gradient
!!$    call PDGEMM('N','N', n_basis_current, n_states_current, n_basis_current, 1.0d0, &
!!$         ovlp_save(1,1), 1, 1, cg_desc, eigenvec(1,1), 1, 1, cg_desc, 0.0d0, Sev, 1, 1, cg_desc)
    call PDSYMM('L', 'U', n_basis_current, n_states_current, 1.0d0, ovlp_save, 1, 1, cg_desc, &
         eigenvec, 1, 1, cg_desc, 0.0d0, Sev, 1, 1, cg_desc)
!!$    call mult_at_b_real_elpa('N','N', n_basis_current, n_states_current, ovlp_save, mxld, eigenvec, mxld, &
!!$         nb, mpi_comm_rows, mpi_comm_cols, Sev, mxld)
    
    ! complete the gradient: -(ham*eigenvec - eigenvalue*eigenvec)
    do i_col = 1, n_states_current, 1
       if (l_col(i_col) == 0) cycle
       current_state = i_col
       res(:, l_col(i_col)) = -Hev(:,l_col(i_col)) + KS_eigenvalue(current_state)*Sev(:,l_col(i_col))
    end do

    tt1 = MPI_Wtime()
    if (myid==0 .and. print_times) print *, 'Time initial gradient: ', tt1-tt0
    
    ! check which eigenvectors are already converged
    n_converged_states = 0
    do i_vector = 1, n_states_current, 1
       current_state = i_vector
       temp_dot = 0.0d0
       call PDDOT( n_basis_current, temp_dot, &
            res, 1, i_vector, cg_desc, 1, res, 1, i_vector, cg_desc, 1)
       if (temp_dot == 0.0d0) temp_dot = -1.0e36
       call get_max_double(residual_norm(current_state), temp_dot)
       residual_norm(current_state) = sqrt(residual_norm(current_state))
       
       if (residual_norm(current_state) < lopcg_tol_current) then
          n_converged_states = n_converged_states + 1
          converged_state(current_state) = .true.
          call PDCOPY( n_basis_current, eigenvec(1,1), 1, current_state, cg_desc, 1, &
               converged_evs, 1, n_converged_states, cg_desc, 1)
       end if
    end do

    if (ALL(converged_state(1:n_states_current))) then
       converged = .true.
       write (info_str,'(2X,A)') "LOPCG converged already on entry."
       call localorb_info(info_str,use_unit,'(A)')       
       write (info_str,'(2X,A,E10.4,A,I4)') "Maximum residual: ", &
            MAXVAL(residual_norm(1:n_states_current)), " in state: ", &
            MAXLOC(residual_norm(1:n_states_current))
       call localorb_info(info_str,use_unit,'(A)')
       return
    end if
    
    state_top = 0
    n_blocks = 0
    slide_block = .false.

    ! loop over blocks
    block_loop: do while (state_top < n_states_current)

       tt1 = MPI_Wtime()
       if (myid==0 .and. print_times .and. n_blocks > 0) print *, 'Block: ', n_blocks, ' time: ', tt1-tt0
       tt0 = MPI_Wtime()

       n_blocks = n_blocks + 1

       if (slide_block) then
          state_bottom = state_bottom + n_slide
          state_top = state_top + n_slide
          ! current_block_size = state_top - state_bottom + 1
          slide_block = .false.
          write (info_str,'(2X,A,I6,A)') "Sliding block in LOPCG by ", n_slide, &
               " states."
          call localorb_info(info_str,use_unit,'(A)')

       else
          
          state_bottom = state_top + 1
          
          ! compute the size of the next block
          if (lopcg_auto_blocksize) then
             current_block_size = 1
             if (state_bottom + current_block_size - 1 < n_states_current) then
                do while ((ABS(KS_eigenvalue(state_bottom) - KS_eigenvalue(state_bottom + current_block_size)) < &
                     0.01*ABS(KS_eigenvalue(state_bottom))).and.(current_block_size < lopcg_block_size))
                   current_block_size = current_block_size + 1
                   if (state_bottom + current_block_size - 1 == n_states_current) exit
                end do
             end if
          else
             current_block_size = MIN(lopcg_block_size, n_states_current - state_top)
          end if

          state_top = state_top + current_block_size

       end if

       ! if everything is converged, skip the block
       if (ALL(converged_state(state_bottom:state_top))) then             
          write (info_str,'(2X,A,I6,A)') "All states converged in block ", n_blocks, &
               ". Skipping the block."
          call localorb_info(info_str,use_unit,'(A)')
          cycle block_loop
       end if
       
       i_iter = 0
       ! the main iteration loop
       iter_loop: do while (i_iter <= max_cg_iterations)
          i_iter = i_iter + 1
          ! in the first iteration we don't have the conjugate direction, hence the smaller ritz dimension

          if ((myid == 0).and.lopbcg_debug) then
             print *, 'Iteration: ', i_iter, ' states: ', state_bottom, " - ", &
                  state_top
          end if

          if (i_iter > 1) then

             ! compute ham*eigenvec as a first part of the gradient
!!$             call PDGEMM('N','N', n_basis_current, current_block_size, n_basis_current, 1.0d0, &
!!$                  ham(1,1), 1, 1, cg_desc, eigenvec(1,1), 1, state_bottom, cg_desc, 0.0d0, Hev, 1, 1, cg_desc)
             call PDSYMM('L', 'U', n_basis_current, current_block_size, 1.0d0, ham, 1, 1, cg_desc, &
                  eigenvec, 1, state_bottom, cg_desc, 0.0d0, Hev, 1, 1, cg_desc)

             ! compute ovlp*eigenvec as a second part of the gradient
             call PDGEMM('N','N', n_basis_current, current_block_size, n_basis_current, 1.0d0, &
                  ovlp_save(1,1), 1, 1, cg_desc, eigenvec(1,1), 1, state_bottom, cg_desc, 0.0d0, Sev, 1, 1, cg_desc)
!!$             call PDSYMM('L', 'U', n_basis_current, current_block_size, 1.0d0, ovlp_save, 1, 1, cg_desc, &
!!$                  eigenvec, 1, state_bottom, cg_desc, 0.0d0, Sev, 1, 1, cg_desc)
             
             ! complete the gradient: -(ham*eigenvec - eigenvalue*eigenvec)
             do i_col = 1, current_block_size, 1
                if (l_col(i_col) == 0) cycle
                current_state = state_bottom + i_col - 1
                res(:, l_col(i_col)) = -Hev(:,l_col(i_col)) + KS_eigenvalue(current_state)*Sev(:,l_col(i_col))
             end do

             ! check which eigenvectors are already converged
             do i_vector = 1, current_block_size, 1
                current_state = state_bottom + i_vector - 1
                temp_dot = 0.0d0
                call PDDOT( n_basis_current, temp_dot, &
                     res, 1, i_vector, cg_desc, 1, res, 1, i_vector, cg_desc, 1)
                if (temp_dot == 0.0d0) temp_dot = -1.0e36
                call get_max_double(residual_norm(current_state), temp_dot)
                residual_norm(current_state) = sqrt(residual_norm(current_state))
                
                if ((residual_norm(current_state) < lopcg_tol_current).and. &
                     (.not.(converged_state(current_state)))) then
                   converged_state(current_state) = .true.
                   n_converged_states = n_converged_states + 1
                   call PDCOPY( n_basis_current, eigenvec(1,1), 1, current_state, cg_desc, 1, &
                        converged_evs, 1, n_converged_states, cg_desc, 1)
                end if

             end do

             if (ALL(converged_state(state_bottom:state_top))) cycle block_loop

             if (lopcg_slide_ev_block) then
                n_slide = 0
                do i_state = state_bottom, state_top, 1
                   if (.not.(converged_state(i_state))) exit
                   n_slide = n_slide + 1
                end do
                if ((n_slide > current_block_size/2).and.(n_slide + state_top <= n_states_current)) then
                   if (myid == 0) print *, n_slide, n_states_current, state_bottom
                   slide_block = .true.
                   cycle block_loop
                end if
             end if

             ! precondition
             call PDGEMM('N','N', n_basis_current, current_block_size, n_basis_current, 1.0d0, &
                  cg_prec(1,1), 1, 1, cg_desc, res, 1, 1, cg_desc, 0.0d0, prec_res, 1, 1, cg_desc)

          else

             ! precondition
             call PDGEMM('N','N', n_basis_current, current_block_size, n_basis_current, 1.0d0, &
                  cg_prec(1,1), 1, 1, cg_desc, res, 1, state_bottom, cg_desc, 0.0d0, prec_res, 1, 1, cg_desc)

          end if

          n_active_states = 0
          do i_vector = 1, current_block_size, 1
             current_state = state_bottom + i_vector - 1
             if (converged_state(current_state)) cycle
             n_active_states = n_active_states + 1
          end do

          if ((myid == 0).and.lopbcg_debug) print *, 'n_active_states: ', n_active_states

          if (lopcg_use_all_evs) then
             state_offset = state_bottom - 1
             do i_vector = 1, state_offset + current_block_size, 1
                call PDCOPY( n_basis_current, eigenvec(1,1), 1, i_vector, cg_desc, 1, &
                     ritz_basis, 1, i_vector, cg_desc, 1)
             end do

             i_active_state = 0
             do i_vector = 1, current_block_size, 1
                current_state = state_bottom + i_vector - 1
                !test
                if ((myid == 0).and.lopbcg_debug) then
                   if (.not.(converged_state(current_state))) then
                      print *, 'State: ', current_state, 'resnorm: ', residual_norm(current_state)
                   end if
                end if
                !test
                if (.not.(converged_state(current_state))) then
                   i_active_state = i_active_state + 1
                   
                   ! store the preconditioned gradient to the ritz basis array
                   call PDCOPY( n_basis_current, prec_res(1,1), 1, i_vector, cg_desc, 1, &
                        ritz_basis, 1, state_offset + current_block_size + i_active_state, cg_desc, 1)

                   ! if we have conjugated information, also store this to the ritz basis
                   if (i_iter /= 1) then
                      call PDCOPY( n_basis_current, conj_dir, 1, current_state, cg_desc, 1, &
                           ritz_basis, 1, state_offset + current_block_size + n_active_states + i_active_state, cg_desc, 1)
                   end if
                end if
             end do
          else
             
             state_offset = 0
             i_active_state = 0
             do i_vector = 1, current_block_size, 1
                current_state = state_bottom + i_vector - 1
                !test
                if ((myid == 0).and.lopbcg_debug) then
                   if (.not.(converged_state(current_state))) then
                      print *, 'State: ', current_state, 'resnorm: ', residual_norm(current_state)
                   end if
                end if
                !test
                if (.not.(converged_state(current_state))) then
                   i_active_state = i_active_state + 1
                   
                   ! store the eigenvector to the ritz basis array
                   call PDCOPY( n_basis_current, eigenvec(1,1), 1, current_state, cg_desc, 1, &
                        ritz_basis, 1, i_active_state, cg_desc, 1)

                   ! store the preconditioned gradient to the ritz basis array
                   call PDCOPY( n_basis_current, prec_res(1,1), 1, i_vector, cg_desc, 1, &
                        ritz_basis, 1, n_active_states + i_active_state, cg_desc, 1)

                   ! if we have conjugated information, also store this to the ritz basis
                   if (i_iter /= 1) then
                      call PDCOPY( n_basis_current, conj_dir, 1, current_state, cg_desc, 1, &
                           ritz_basis, 1, 2*n_active_states + i_active_state, cg_desc, 1)
                   end if
                end if
             end do
          end if
          ! now the ritz basis array is ready, there we have:
          ! ritz_basis(:,1:n_active_states): the current eigenvectors
          ! ritz_basis(:,n_active_states:2*n_active_states): the (preconditioned) gradients
          ! ritz_basis(:,2*n_active_states:3*n_active_states): the conjugated directions, if available

          ! this test ensures that if the block is converged, also the block loop is exited
          if (ALL(converged_state(state_bottom:state_top))) then             
             write (info_str,'(2X,A,I6,A)') "All states converged in block ", n_blocks, &
                  ". Skipping the block."
             call localorb_info(info_str,use_unit,'(A)')
             ! KS_eigenvalue(state_bottom:state_top) = KS_eigenvalue(state_bottom:state_top) - ev_shift
             cycle block_loop
          end if

          ! when the vectors in ritz basis span the desired set we can build the corresponding
          ! ritz matrices
          if (lopcg_use_all_evs) then
             if (i_iter == 1) then
                current_ritz_dim = state_offset + current_block_size + n_active_states
             else
                current_ritz_dim = state_offset + current_block_size + 2*n_active_states
             end if
          else
             if (i_iter == 1) then
                current_ritz_dim = 2*n_active_states
             else
                current_ritz_dim = 3*n_active_states
             end if
          end if

          ! lock the Ritz space, i.e., project out the components corresponding to lower eigenvectors
          ! in addition, orthonormalize Ritz vectors

          tt2 = MPI_Wtime()

          if ((.not.lopcg_use_all_evs).and.(n_converged_states > 0)) then
             call project_ritz_basis_scalapack_real( current_ritz_dim, ritz_basis, n_converged_states, converged_evs )
          end if
          call orthonormalize_ritz_basis_scalapack_real( current_ritz_dim, ritz_basis )
          if (lopbcg_debug) call check_rb_orthogonality_real( current_ritz_dim, ritz_basis )

          tt3 = MPI_Wtime()
          if (myid==0 .and. print_times .and. MOD(i_iter,100) == 2) &
               print *, 'Iter: ', i_iter, ' projection time: ', tt3-tt2

          tt2 = MPI_Wtime()

          tmp = 0.0d0
          another_tmp = 0.0d0
          ! construct ritz_matrix_H
          call PDGEMM('N','N', n_basis_current, current_ritz_dim, n_basis_current, 1.0d0, &
               ham(1,1), 1, 1, cg_desc, ritz_basis, 1, 1, cg_desc, 0.0d0, tmp, 1, 1, cg_desc)
!!$          call mult_at_b_real_elpa('N','N', n_basis_current, current_ritz_dim, ham, mxld, ritz_basis, mxld, &
!!$               nb, mpi_comm_rows, mpi_comm_cols, tmp, mxld)

          call PDGEMM('T','N', current_ritz_dim, current_ritz_dim, n_basis_current, 1.0d0, &
               ritz_basis, 1, 1, cg_desc, tmp, 1, 1, cg_desc, 0.0d0, another_tmp, 1, 1, cg_desc)

          if (use_ritz_scalapack) then

             call descinit( rz_desc, current_ritz_dim, current_ritz_dim, mb, nb, rsrc, csrc, &
                  my_blacs_ctxt, MAX(1,mxld), info )

             tmp = another_tmp
             call PDTRAN(current_ritz_dim, current_ritz_dim, 0.5d0, tmp, 1, 1, rz_desc, &
                  0.5d0, another_tmp, 1, 1, rz_desc)

             call PDSYEVD('V','U',current_ritz_dim,another_tmp,1,1,rz_desc,ritz_values, &
                  tmp,1,1,rz_desc, work, lwork, iwork, liwork, info)
             if(info /= 0) then
                if (myid == 0) then
                   print *, 'ERROR in pdsyevd.'
                   print *, 'iter, block, dim, info: ', i_iter, n_blocks, current_ritz_dim, info
                end if
                exit iter_loop
             end if
!!$             if (myid==0) print *, 'Ritz values: '
!!$             if (myid==0) print *, ritz_values

          else

             ritz_matrix_H = 0.0d0
             do i_col = 1, current_ritz_dim, 1
                if (l_col(i_col) == 0) cycle
                do i_row = 1, current_ritz_dim, 1
                   if (l_row(i_row) == 0) cycle
                   ritz_matrix_H(i_row, i_col) = another_tmp(l_row(i_row), l_col(i_col))
                end do
             end do

          ! construct ritz_matrix_S
!!$          ritz_matrix_S = 0.0d0
!!$          call PDGEMM('N','N', n_basis_current, current_ritz_dim, n_basis_current, 1.0d0, &
!!$               ovlp_save(1,1), 1, 1, cg_desc, ritz_basis, 1, 1, cg_desc, 0.0d0, tmp, 1, 1, cg_desc)
!!$          call PDGEMM('T','N', current_ritz_dim, current_ritz_dim, n_basis_current, 1.0d0, &
!!$               ritz_basis, 1, 1, cg_desc, tmp, 1, 1, cg_desc, 0.0d0, another_tmp, 1, 1, cg_desc)
!!$          do i_col = 1, current_ritz_dim, 1
!!$             if (l_col(i_col) == 0) cycle
!!$             do i_row = 1, current_ritz_dim, 1
!!$                if (l_row(i_row) == 0) cycle
!!$                ritz_matrix_S(i_row, i_col) = another_tmp(l_row(i_row), l_col(i_col))
!!$             end do
!!$          end do

             call sync_vector(ritz_matrix_H, max_ritz_dim*max_ritz_dim, my_scalapack_comm_work)
!!$          call sync_vector(ritz_matrix_S, max_ritz_dim*max_ritz_dim, my_scalapack_comm_work)
!!$
             ritz_matrix_H(1:current_ritz_dim,1:current_ritz_dim) = 0.5d0*( &
                  TRANSPOSE(ritz_matrix_H(1:current_ritz_dim,1:current_ritz_dim)) + &
                  ritz_matrix_H(1:current_ritz_dim,1:current_ritz_dim) )

!!$          ritz_matrix_S(1:current_ritz_dim,1:current_ritz_dim) = 0.5d0*( &
!!$               TRANSPOSE(ritz_matrix_S(1:current_ritz_dim,1:current_ritz_dim)) + &
!!$               ritz_matrix_S(1:current_ritz_dim,1:current_ritz_dim) )
          ! solve the ritz problem to find the minimizing vectors

!!$          call dsygvd(1, 'V', 'U', current_ritz_dim, ritz_matrix_H(1:current_ritz_dim,1:current_ritz_dim), &
!!$               current_ritz_dim, ritz_matrix_S(1:current_ritz_dim,1:current_ritz_dim), &
!!$               current_ritz_dim, ritz_values, work, lwork, iwork, liwork, info)

             call dsyevd('V', 'U', current_ritz_dim, ritz_matrix_H(1:current_ritz_dim,1:current_ritz_dim), &
                  current_ritz_dim, ritz_values, work, lwork, iwork, liwork, info)

!!$             if (myid==0) print *, 'Ritz values: '
!!$             if (myid==0) print *, ritz_values

          ! Our worst enemy here is the case when the ritz space becomes linearly dependent.
          ! This usually manifests itself in ritz_matrix_S that becomes indefinite and we get
          ! an error from LAPACK. Using smaller blocksize usually helps here (but this of course
          ! affects the convergence.
             if(info /= 0) then
                if (myid == 0) then
                   print *, 'ERROR in dsygvd.'
                   print *, 'iter, block, dim, info: ', i_iter, n_blocks, current_ritz_dim, info
!!$                do i_row = 1, current_ritz_dim, 1
!!$                   do i_col = 1, current_ritz_dim, 1
!!$                      print *, i_row, i_col, ritz_matrix_S(i_row, i_col)
!!$                   end do
!!$                end do
                end if
!!             cycle block_loop
                exit iter_loop
!!$             call scalapack_err_exit(info,"DSYGV")
             end if

          end if

          tt3 = MPI_Wtime()
          if (myid==0 .and. print_times .and. MOD(i_iter,100)==2) &
               print *, 'Iter: ', i_iter, ' Ritz time: ', tt3-tt2

          tt2 = MPI_Wtime()

          if (lopcg_use_all_evs) then

             if (.not.use_ritz_scalapack) then
                do i_vector = 1, state_offset + current_block_size, 1
                   if (l_col(i_vector) /= 0) then
                      do i_row = 1, current_ritz_dim, 1
                         if (l_row(i_row) /= 0) tmp(l_row(i_row),l_col(i_vector)) = &
                              ritz_matrix_H(i_row, i_vector)
                      end do
                   end if
                end do
             end if

             call PDGEMM('N','N', n_basis_current, state_offset + current_block_size, current_ritz_dim, 1.0d0, &
                  ritz_basis, 1, 1, cg_desc, tmp, 1, 1, cg_desc, 0.0d0, another_tmp, 1, 1, cg_desc)

             do i_vector = 1, state_offset + current_block_size, 1
                call PDCOPY( n_basis_current, another_tmp, 1, i_vector, cg_desc, 1, &
                     eigenvec(1,1), 1, i_vector, cg_desc, 1)
                KS_eigenvalue(i_vector) = ritz_values(i_vector)
             end do

             call PDGEMM('N','N', n_basis_current, n_active_states, current_ritz_dim-current_block_size-state_offset, 1.0d0, &
                  ritz_basis, 1, state_offset + current_block_size + 1, cg_desc, &
                  tmp, state_offset + current_block_size + 1, 1, cg_desc, 0.0d0, &
                  another_tmp, 1, 1, cg_desc)

             i_active_state = 0
             do i_vector = 1, current_block_size, 1
                current_state = state_bottom - 1 + i_vector
                
                if (converged_state(current_state)) cycle
                i_active_state = i_active_state + 1

                call PDCOPY( n_basis_current, another_tmp, 1, i_active_state, cg_desc, 1, &
                     conj_dir, 1, current_state, cg_desc, 1)
                
             end do
             
             ! normalize vectors
             call PDGEMM('N','N', n_basis_current, state_offset + current_block_size, n_basis_current, 1.0d0, &
                  ovlp_save(1,1), 1, 1, cg_desc, eigenvec(1,1), 1, 1, cg_desc, &
                  0.0d0, tmp, 1, 1, cg_desc)
             
             call PDGEMM('N','N', n_basis_current, current_block_size, n_basis_current, 1.0d0, &
                  ovlp_save(1,1), 1, 1, cg_desc, conj_dir, 1, state_bottom, cg_desc, &
                  0.0d0, another_tmp, 1, 1, cg_desc)
             
             do i_vector = 1, state_offset + current_block_size, 1
                current_state = i_vector

                temp_dot = 0.0d0
                call PDDOT( n_basis_current, temp_dot, &
                     eigenvec(1,1), 1, current_state, cg_desc, 1, tmp, 1, i_vector, cg_desc, 1)
                if (temp_dot == 0.0d0) temp_dot = -1.0e36
                call get_max_double(temp_norm, temp_dot)
                
                temp_norm = sqrt(temp_norm)
                if (l_col(current_state) /= 0) then
                   do i_row = 1, n_basis_current, 1
                      if (l_row(i_row) == 0) cycle
                      eigenvec(l_row(i_row), l_col(current_state)) = &
                           eigenvec(l_row(i_row), l_col(current_state)) / temp_norm
                   end do
                end if
                
             end do

             do i_vector = 1, current_block_size, 1
                current_state = state_bottom - 1 + i_vector

                temp_dot = 0.0d0
                call PDDOT( n_basis_current, temp_dot, &
                     conj_dir, 1, current_state, cg_desc, 1, another_tmp, 1, i_vector, cg_desc, 1)
                if (temp_dot == 0.0d0) temp_dot = -1.0e36
                call get_max_double(temp_norm, temp_dot)
                
                temp_norm = sqrt(temp_norm)
                if (l_col(current_state) /= 0) then
                   do i_row = 1, n_basis_current, 1
                      if (l_row(i_row) == 0) cycle
                      conj_dir(l_row(i_row), l_col(current_state)) = &
                           conj_dir(l_row(i_row), l_col(current_state)) / temp_norm
                   end do
                end if

             end do

          else

             ! compute the new conjugated direction
             if (.not.use_ritz_scalapack) then
                do i_vector = 1, n_active_states, 1
                   if (l_col(i_vector) /= 0) then
                      do i_row = 1, current_ritz_dim, 1
                         if (l_row(i_row) /= 0) tmp(l_row(i_row),l_col(i_vector)) = &
                              ritz_matrix_H(i_row, i_vector)
                      end do
                   end if
                end do
             end if

             call PDGEMM('N','N', n_basis_current, n_active_states, current_ritz_dim-n_active_states, 1.0d0, &
                  ritz_basis, 1, n_active_states + 1, cg_desc, &
                  tmp, n_active_states + 1, 1, cg_desc, 0.0d0, &
                  another_tmp, 1, n_active_states + 1, cg_desc)
             
             ! compute first part of the new eigenvector
             call PDGEMM('N','N', n_basis_current, n_active_states, n_active_states, 1.0d0, &
                  ritz_basis, 1, 1, cg_desc, tmp, 1, 1, cg_desc, 0.0d0, &
                  another_tmp, 1, 1, cg_desc)
             
             i_active_state = 0
             do i_vector = 1, current_block_size, 1
                current_state = state_bottom - 1 + i_vector
                
                if (converged_state(current_state)) cycle
                i_active_state = i_active_state + 1
                
                ! complete the new eigenvector
                call PDAXPY(n_basis_current, 1.0d0, another_tmp, 1, n_active_states + i_active_state, cg_desc, 1, &
                     another_tmp, 1, i_active_state, cg_desc, 1)
                
                call PDCOPY( n_basis_current, another_tmp, 1, i_active_state, cg_desc, 1, &
                     eigenvec(1,1), 1, current_state, cg_desc, 1)
                call PDCOPY( n_basis_current, another_tmp, 1, n_active_states + i_active_state, cg_desc, 1, &
                     conj_dir, 1, current_state, cg_desc, 1)
                
                ! update the eigenvalues
                KS_eigenvalue(current_state) = ritz_values(i_active_state)
             end do

             ! normalize vectors
             call PDGEMM('N','N', n_basis_current, current_block_size, n_basis_current, 1.0d0, &
                  ovlp_save(1,1), 1, 1, cg_desc, eigenvec(1,1), 1, state_bottom, cg_desc, &
                  0.0d0, tmp, 1, 1, cg_desc)
             
             call PDGEMM('N','N', n_basis_current, current_block_size, n_basis_current, 1.0d0, &
                  ovlp_save(1,1), 1, 1, cg_desc, conj_dir, 1, state_bottom, cg_desc, &
                  0.0d0, another_tmp, 1, 1, cg_desc)
             
             do i_vector = 1, current_block_size, 1
                current_state = state_bottom - 1 + i_vector
                
                temp_dot = 0.0d0
                call PDDOT( n_basis_current, temp_dot, &
                     eigenvec(1,1), 1, current_state, cg_desc, 1, tmp, 1, i_vector, cg_desc, 1)
                if (temp_dot == 0.0d0) temp_dot = -1.0e36
                call get_max_double(temp_norm, temp_dot)
                
                temp_norm = sqrt(temp_norm)
                if (l_col(current_state) /= 0) then
                   do i_row = 1, n_basis_current, 1
                      if (l_row(i_row) == 0) cycle
                      eigenvec(l_row(i_row), l_col(current_state)) = &
                           eigenvec(l_row(i_row), l_col(current_state)) / temp_norm
                   end do
                end if
                
                temp_dot = 0.0d0
                call PDDOT( n_basis_current, temp_dot, &
                     conj_dir, 1, current_state, cg_desc, 1, another_tmp, 1, i_vector, cg_desc, 1)
                if (temp_dot == 0.0d0) temp_dot = -1.0e36
                call get_max_double(temp_norm, temp_dot)
                
                temp_norm = sqrt(temp_norm)
                if (l_col(current_state) /= 0) then
                   do i_row = 1, n_basis_current, 1
                      if (l_row(i_row) == 0) cycle
                      conj_dir(l_row(i_row), l_col(current_state)) = &
                           conj_dir(l_row(i_row), l_col(current_state)) / temp_norm
                   end do
                end if
                
             end do
             
          end if

          tt3 = MPI_Wtime()
          if (myid==0 .and. print_times .and. MOD(i_iter,100) == 2) &
               print *, 'Iter: ', i_iter, ' vector update time: ', tt3-tt2

!!$          shift = 0
!!$          do i_vector = 1, current_block_size, 1
!!$             current_state = state_bottom - 1 + i_vector
!!$             if (converged_state(current_state)) then
!!$                shift = shift + 1
!!$             else
!!$                exit
!!$             end if
!!$          end do
!!$          state_bottom = state_bottom + shift
!!$          current_block_size = current_block_size - shift
!!$          if (current_block_size <= 0) exit

          total_iter = total_iter + 1

       end do iter_loop

       if (MAXVAL(residual_norm(1:n_states_current)) < lopcg_tol_current) exit block_loop

    end do block_loop

    ! once the lopcg-iteration is over we need to update the eigenvalues since the hamiltonian has
    ! changed from the previous scf-iteration

    tt0 = MPI_Wtime()

    call PDGEMM('N','N', n_basis_current, n_states_current, n_basis_current, 1.0d0, &
         ham(1,1), 1, 1, cg_desc, eigenvec(1,1), 1, 1, cg_desc, 0.0d0, Hev, 1, 1, cg_desc)
    do i_state = 1, n_states_current, 1
       current_state = i_state
       temp_dot = 0.0d0
       call PDDOT( n_basis_current, temp_dot, &
            eigenvec(1,1), 1, current_state, cg_desc, 1, Hev, 1, current_state, cg_desc, 1)
       if (temp_dot == 0.0d0) temp_dot = -1.0e36
       call get_max_double(KS_eigenvalue(current_state), temp_dot)
    end do

    call PDGEMM('N','N', n_basis_current, n_states_current, n_basis_current, 1.0d0, &
         ovlp_save(1,1), 1, 1, cg_desc, eigenvec(1,1), 1, 1, cg_desc, 0.0d0, Sev, 1, 1, cg_desc)
    do i_state = 1, n_states_current, 1
       current_state = i_state
       temp_dot = 0.0d0
       call PDDOT( n_basis_current, temp_dot, &
            eigenvec(1,1), 1, current_state, cg_desc, 1, Sev, 1, current_state, cg_desc, 1)
       if (temp_dot == 0.0d0) temp_dot = -1.0e36
       call get_max_double(temp_norm, temp_dot)
       KS_eigenvalue(current_state) = KS_eigenvalue(current_state) / temp_norm
    end do

    do i_state = 1, n_states_current, 1
       current_state = i_state
       if (l_col(current_state) == 0) cycle
       res(:, l_col(current_state)) = -Hev(:,l_col(current_state)) &
            + KS_eigenvalue(current_state)*Sev(:,l_col(current_state))
    end do
!!$       call PDAXPY( n_basis_current, -KS_eigenvalue(current_state), &
!!$            Sev, 1, current_state, cg_desc, 1, &
!!$            Hev, 1, current_state, cg_desc, 1)
    do i_state = 1, n_states_current, 1
       temp_norm = 0.0d0
!!$       call PDNRM2( n_basis_current, temp_norm, Hev, 1, current_state, cg_desc, 1)
       call PDNRM2( n_basis_current, temp_norm, res, 1, current_state, cg_desc, 1)
       call get_max_double(residual_norm(current_state), temp_norm)
    end do

    tt1 = MPI_Wtime()
    if (myid==0 .and. print_times) print *, 'Time final updates: ', tt1-tt0

    ! store the current eigenvectors for further processing in solve_evp_scalapack
!!    tmp(:,:) = eigenvec_untrafo(:,:,i_spin)
    tmp(:,:) = 0.0d0
!!    eigenvec(:,:,i_spin) = eigenvec_untrafo(:,:,i_spin)

    ! print out some information on how well we did
    if (MAXVAL(residual_norm(1:n_states_current)) >= lopcg_tol_current) then
       converged = .false.
       write (info_str,'(2X,A,I6,A)') "LOPCG didn't converge in ", max_cg_iterations, &
            " iterations per block."
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(2X,A,E10.4,A,I4)') "Maximum remaining residual: ", &
            MAXVAL(residual_norm(1:n_states_current)),  " in state: ", &
            MAXLOC(residual_norm(1:n_states_current))
       call localorb_info(info_str,use_unit,'(A)')
    else
       converged = .true.
       write (info_str,'(2X,A,I6,A,I6,A)') "LOPCG converged in total of ", total_iter, &
            " iterations for ", n_blocks, " blocks."
       call localorb_info(info_str,use_unit,'(A)')       
       write (info_str,'(2X,A,E10.4,A,I4)') "Maximum remaining residual: ", &
            MAXVAL(residual_norm(1:n_states_current)), " in state: ", &
            MAXLOC(residual_norm(1:n_states_current))
       call localorb_info(info_str,use_unit,'(A)')
    end if

!test
    if ((myid == 0).and.lopbcg_debug) then
       print *, 'ON EXIT:'
       do i_state = 1, n_states_current, 1
          print *, 'state: ', i_state, ' EV: ', KS_eigenvalue(i_state), ' res: ', residual_norm(i_state)
       end do
    end if

!    call check_ev_orthogonality_real()

!test

  end subroutine cg_solver_scalapack
!******
!-----------------------------------------------------------------------------------
!****s* cg_scalapack/lock_ritz_space_scalapack
!  NAME
!    lock_ritz_space_scalapack
!  SYNOPSIS
  subroutine lock_ritz_space_scalapack( n_basis_current, n_states_current, ritz_basis, &
       current_state, block_size, active_size, converged, current_iter, cg_desc, i_spin, &
       eigenvec )
!  PURPOSE
!    Locks the Ritz-space for the LOPCG-solver, i.e., it is made orthogonal
!    to the lower eigenspaces.
!  USES
    use synchronize_mpi_basic, only: get_max_double
!  ARGUMENTS
    integer :: n_basis_current, n_states_current
    real*8, dimension(mxld,mxcol) :: ritz_basis, eigenvec
    integer :: current_state, block_size, current_iter, active_size
    integer, dimension(dlen_) :: cg_desc
    integer :: i_spin
    logical :: converged(n_states_current)
!  INPUTS
!    o n_basis_current -- current number of basis functions
!    o ritz_basis -- the Ritz basis vectors
!    o current_state -- the top state of the lower eigenvectors
!    o block_size -- size of the current LOPCG-block
!    o current_iter -- current iteration number
!    o i_spin -- the spin channel
!    o cg_desc -- ScaLAPACK descriptor for the matrices
!  OUTPUT
!    ritz_basis is orthonormalized and otrhogonalized to lower eigenvector
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    integer :: i_state, i_vector, i_vector_2, i_col, i_row
    integer :: ritz_dim
    real*8 :: r, temp_norm, temp_dot
    real*8, dimension(3*active_size) :: r_vec

    ! This routine is called only from the working set, so no need to check here
!!$    if (current_state == 1) return

    if (current_iter == 1) then
       ritz_dim = 2*active_size
    else
       ritz_dim = 3*active_size
    end if

    call PDGEMM('N','N', n_basis_current, ritz_dim, n_basis_current, 1.0d0, &
         ovlp_save(1,1), 1, 1, cg_desc, ritz_basis, 1, 1, cg_desc, 0.0d0, tmp, 1, 1, cg_desc)

    ! stage zero: normalize
    do i_vector = 1, ritz_dim, 1
       temp_dot = 0.0d0
       call PDDOT( n_basis_current, temp_dot, ritz_basis, 1, i_vector, cg_desc, 1, &
            tmp, 1, i_vector, cg_desc, 1)
       if (temp_dot == 0.0d0) temp_dot = -1.0e36
       call get_max_double(temp_norm, temp_dot)
!test
!!       if (myid==0) print *, i_vector, temp_norm
!test
       r_vec(i_vector) = sqrt(temp_norm)
    end do
    do i_col = 1, 2*active_size, 1
       if (l_col(i_col) == 0) cycle
       do i_row = 1, n_basis_current, 1
          if (l_row(i_row) == 0) cycle
          ritz_basis(l_row(i_row), l_col(i_col)) = ritz_basis(l_row(i_row), l_col(i_col))/r_vec(i_col)
       end do
    end do
    if (current_iter /= 1) then
       do i_col = 2*active_size+1, 3*active_size, 1
          if (l_col(i_col) == 0) cycle
          do i_row = 1, n_basis_current, 1
             if (l_row(i_row) == 0) cycle
             ritz_basis(l_row(i_row), l_col(i_col)) = ritz_basis(l_row(i_row), l_col(i_col))/r_vec(i_col)
          end do
       end do    
    end if

!!$    call mult_at_b_real_elpa('U','N', n_basis_current, ritz_dim, ovlp_save(1,1), mxld, &
!!$         ritz_basis, mxld, nb, mpi_comm_rows, mpi_comm_cols, tmp, mxld )
    call PDGEMM('N','N', n_basis_current, ritz_dim, n_basis_current, 1.0d0, &
         ovlp_save(1,1), 1, 1, cg_desc, ritz_basis, 1, 1, cg_desc, 0.0d0, tmp, 1, 1, cg_desc)
    
    ! stage one: project out the space spanned by lower eigenvectors
    do i_state = 1, current_state-1 + block_size, 1
       do i_vector = 1, active_size, 1

          if (.not.converged(i_state)) cycle

          temp_dot = 0.0d0
          call PDDOT( n_basis_current, temp_dot, eigenvec(1,1), 1, i_state, cg_desc, 1, &
               tmp, 1, i_vector, cg_desc, 1)
          if (temp_dot == 0.0d0) temp_dot = -1.0e36
          call get_max_double(r, temp_dot)

          call PDAXPY( n_basis_current, -r, eigenvec(1,1), 1, i_state, cg_desc, 1, &
               ritz_basis, 1, i_vector, cg_desc, 1)

          temp_dot = 0.0d0
          call PDDOT( n_basis_current, temp_dot, eigenvec(1,1), 1, i_state, cg_desc, 1, &
               tmp, 1, active_size + i_vector, cg_desc, 1)
          if (temp_dot == 0.0d0) temp_dot = -1.0e36
          call get_max_double(r, temp_dot)

          call PDAXPY( n_basis_current, -r, eigenvec(1,1), 1, i_state, cg_desc, 1, &
               ritz_basis, 1, active_size + i_vector, cg_desc, 1)

          if (current_iter /= 1) then
             temp_dot = 0.0d0
             call PDDOT( n_basis_current, temp_dot, eigenvec(1,1), 1, i_state, cg_desc, 1, &
                  tmp, 1, 2*active_size + i_vector, cg_desc, 1)
             if (temp_dot == 0.0d0) temp_dot = -1.0e36
             call get_max_double(r, temp_dot)

             call PDAXPY( n_basis_current, -r, eigenvec(1,1), 1, i_state, cg_desc, 1, &
                  ritz_basis, 1, 2*active_size + i_vector, cg_desc, 1)
          end if

       end do
    end do

    call PDGEMM('N','N', n_basis_current, ritz_dim, n_basis_current, 1.0d0, &
         ovlp_save(1,1), 1, 1, cg_desc, ritz_basis, 1, 1, cg_desc, 0.0d0, tmp, 1, 1, cg_desc)

    ! stage two: normalize
    do i_vector = 1, ritz_dim, 1
       temp_dot = 0.0d0
       call PDDOT( n_basis_current, temp_dot, ritz_basis, 1, i_vector, cg_desc, 1, &
            tmp, 1, i_vector, cg_desc, 1)
       if (temp_dot == 0.0d0) temp_dot = -1.0e36
       call get_max_double(temp_norm, temp_dot)
!test
!!       if (myid==0) print *, i_vector, temp_norm
!test
       r_vec(i_vector) = sqrt(temp_norm)
    end do
    do i_col = 1, 2*active_size, 1
       if (l_col(i_col) == 0) cycle
       do i_row = 1, n_basis_current, 1
          if (l_row(i_row) == 0) cycle
          ritz_basis(l_row(i_row), l_col(i_col)) = ritz_basis(l_row(i_row), l_col(i_col))/r_vec(i_col)
       end do
    end do
    if (current_iter /= 1) then
       do i_col = 2*active_size+1, 3*active_size, 1
          if (l_col(i_col) == 0) cycle
          do i_row = 1, n_basis_current, 1
             if (l_row(i_row) == 0) cycle
             ritz_basis(l_row(i_row), l_col(i_col)) = ritz_basis(l_row(i_row), l_col(i_col))/r_vec(i_col)
          end do
       end do    
    end if
    call PDGEMM('N','N', n_basis_current, ritz_dim, n_basis_current, 1.0d0, &
         ovlp_save(1,1), 1, 1, cg_desc, ritz_basis, 1, 1, cg_desc, 0.0d0, tmp, 1, 1, cg_desc)
    
    ! stage three: project out the eigenvecs from the residuals
    do i_vector = 1, active_size, 1
       do i_vector_2 = 1, active_size, 1
          temp_dot = 0.0d0
          call PDDOT( n_basis_current, temp_dot, ritz_basis, 1, i_vector, cg_desc, 1, &
               tmp, 1, active_size + i_vector_2, cg_desc, 1)
          if (temp_dot == 0.0d0) temp_dot = -1.0e36
          call get_max_double(r, temp_dot)

          call PDAXPY( n_basis_current, -r, ritz_basis, 1, i_vector, cg_desc, 1, &
               ritz_basis, 1, active_size + i_vector_2, cg_desc, 1)
          
       end do
    end do

    call PDGEMM('N','N', n_basis_current, ritz_dim, n_basis_current, 1.0d0, &
         ovlp_save(1,1), 1, 1, cg_desc, ritz_basis, 1, 1, cg_desc, 0.0d0, tmp, 1, 1, cg_desc)

    ! stage four: normalize
    do i_vector = 1, ritz_dim, 1
       temp_dot = 0.0d0
       call PDDOT( n_basis_current, temp_dot, ritz_basis, 1, i_vector, cg_desc, 1, &
            tmp, 1, i_vector, cg_desc, 1)
       if (temp_dot == 0.0d0) temp_dot = -1.0e36
       call get_max_double(temp_norm, temp_dot)
!test
!!       if (myid==0) print *, i_vector, temp_norm
!test
       r_vec(i_vector) = sqrt(temp_norm)
    end do
    do i_col = 1, 2*active_size, 1
       if (l_col(i_col) == 0) cycle
       do i_row = 1, n_basis_current, 1
          if (l_row(i_row) == 0) cycle
          ritz_basis(l_row(i_row), l_col(i_col)) = ritz_basis(l_row(i_row), l_col(i_col))/r_vec(i_col)
       end do
    end do
    if (current_iter /= 1) then
       do i_col = 2*active_size+1, 3*active_size, 1
          if (l_col(i_col) == 0) cycle
          do i_row = 1, n_basis_current, 1
             if (l_row(i_row) == 0) cycle
             ritz_basis(l_row(i_row), l_col(i_col)) = ritz_basis(l_row(i_row), l_col(i_col))/r_vec(i_col)
          end do
       end do    
    end if

    call PDGEMM('N','N', n_basis_current, ritz_dim, n_basis_current, 1.0d0, &
         ovlp_save(1,1), 1, 1, cg_desc, ritz_basis, 1, 1, cg_desc, 0.0d0, tmp, 1, 1, cg_desc)

    ! stage five: orthogonalilze within blocks
    do i_vector = 1, active_size, 1
       do i_vector_2 = i_vector+1, active_size, 1
          temp_dot = 0.0d0
          call PDDOT( n_basis_current, temp_dot, ritz_basis, 1, i_vector, cg_desc, 1, &
               tmp, 1, i_vector_2, cg_desc, 1)
          if (temp_dot == 0.0d0) temp_dot = -1.0e36
          call get_max_double(r, temp_dot)

          call PDAXPY( n_basis_current, -r, ritz_basis, 1, i_vector_2, cg_desc, 1, &
               ritz_basis, 1, i_vector, cg_desc, 1)

          temp_dot = 0.0d0
          call PDDOT( n_basis_current, temp_dot, ritz_basis, 1, active_size + i_vector, cg_desc, 1, &
               tmp, 1, active_size + i_vector_2, cg_desc, 1)
          if (temp_dot == 0.0d0) temp_dot = -1.0e36
          call get_max_double(r, temp_dot)

          call PDAXPY( n_basis_current, -r, ritz_basis, 1, active_size + i_vector_2, cg_desc, 1, &
               ritz_basis, 1, active_size + i_vector, cg_desc, 1)

          if (current_iter /= 1) then
             temp_dot = 0.0d0
             call PDDOT( n_basis_current, temp_dot, ritz_basis, 1, 2*active_size + i_vector, cg_desc, 1, &
                  tmp, 1, 2*active_size + i_vector_2, cg_desc, 1)
             if (temp_dot == 0.0d0) temp_dot = -1.0e36
             call get_max_double(r, temp_dot)

             call PDAXPY( n_basis_current, -r, ritz_basis, 1, 2*active_size + i_vector_2, cg_desc, 1, &
                  ritz_basis, 1, 2*active_size + i_vector, cg_desc, 1)
          end if
       end do
    end do

    call PDGEMM('N','N', n_basis_current, ritz_dim, n_basis_current, 1.0d0, &
         ovlp_save(1,1), 1, 1, cg_desc, ritz_basis, 1, 1, cg_desc, 0.0d0, tmp, 1, 1, cg_desc)

    ! stage six: normalize again
    do i_vector = 1, ritz_dim, 1
       temp_dot = 0.0d0
       call PDDOT( n_basis_current, temp_dot, ritz_basis, 1, i_vector, cg_desc, 1, &
            tmp, 1, i_vector, cg_desc, 1)
       if (temp_dot == 0.0d0) temp_dot = -1.0e36
       call get_max_double(temp_norm, temp_dot)
!test
!!       if (myid==0) print *, i_vector, temp_norm
!test
       r_vec(i_vector) = sqrt(temp_norm)
    end do
    do i_col = 1, 2*active_size, 1
       if (l_col(i_col) == 0) cycle
       do i_row = 1, n_basis_current, 1
          if (l_row(i_row) == 0) cycle
          ritz_basis(l_row(i_row), l_col(i_col)) = ritz_basis(l_row(i_row), l_col(i_col))/r_vec(i_col)
       end do
    end do
    if (current_iter /= 1) then
       do i_col = 2*active_size+1, 3*active_size, 1
          if (l_col(i_col) == 0) cycle
          do i_row = 1, n_basis_current, 1
             if (l_row(i_row) == 0) cycle
             ritz_basis(l_row(i_row), l_col(i_col)) = ritz_basis(l_row(i_row), l_col(i_col))/r_vec(i_col)
          end do
       end do    
    end if

  end subroutine lock_ritz_space_scalapack
!******
!-----------------------------------------------------------
!****s* cg_scalapack/project_ritz_basis_scalapack_real
!  NAME
!    project_ritz_basis_scalapack_real
!  SYNOPSIS
subroutine project_ritz_basis_scalapack_real( n_vectors, ritz_basis, n_evs, converged_evs )
!  PURPOSE
!    Orthonormalizes the Ritz basis with ScaLAPACK, real version
!    This routine uses a modified Gram-Schmidt orthogonalization
!    It needs two additional work matrices but is fast
!    because most work is done in matrix-matrix-products
!
!    This routine must be called after the overlap matrix is saved !!
!
!  USES
  use dimensions, only: n_basis
  use localorb_io, only: localorb_info, OL_norm, use_unit
  use mpi_tasks, only: check_allocation
!  ARGUMENTS
  integer :: n_vectors, n_evs
  real*8, dimension(mxld,mxcol) :: ritz_basis, converged_evs
!  INPUTS
!    n_vectors - number of vectors in the Ritz basis
!  OUTPUT
!    -
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  real*8, allocatable, dimension(:,:) :: ovlp_ev, work
  integer :: i_col

  character*200 :: info_str

  if (lopbcg_debug) then
     write(info_str,'(2X,A)') "Projecting out converged states"
     call localorb_info(info_str,use_unit,'(A)',OL_norm)
  end if
     
  allocate(ovlp_ev(mxld,mxcol),stat=info)
  call check_allocation(info, 'ovlp_ev                       ')

  allocate(work(mxld,mxcol),stat=info)
  call check_allocation(info, 'work                          ')

  if(my_scalapack_id<npcol*nprow) then ! The main work is done only on the working set

     ! tmp = ovlp_save * ritz_basis    
     call PDSYMM('L', 'U', n_basis, n_vectors, 1.0d0, ovlp_save, 1, 1, cg_desc, &
          ritz_basis(1,1), 1, 1, cg_desc, &
          0.0d0, tmp, 1, 1, cg_desc)

     ! another_tmp = congerged_evs**T * tmp
     call PDGEMM('T','N',n_evs,n_vectors,n_basis, 1.0d0, converged_evs, 1, 1, cg_desc, &
          tmp, 1, 1, cg_desc, 0.0d0, another_tmp, 1, 1, cg_desc)

     ! ovlp_ev = ovlp_save * converged_evs
     call PDSYMM('L', 'U', n_basis, n_evs, 1.0d0, ovlp_save, 1, 1, cg_desc, &
          converged_evs(1,1), 1, 1, cg_desc, 0.0d0, ovlp_ev, 1, 1, cg_desc)

     ! work = converged_evs**T * ovlp_ev
     call PDGEMM('T','N',n_evs,n_evs,n_basis, 1.0d0, converged_evs, 1, 1, cg_desc, &
          ovlp_ev, 1, 1, cg_desc, 0.0d0, work, 1, 1, cg_desc)

     ! work = work**-1
     CALL PDPOTRF('U', n_evs, work, 1, 1, cg_desc, info )
     if(info /= 0) call cg_scalapack_err_exit(info,"PDPOTRF")
     CALL PDPOTRI('U', n_evs, work, 1, 1, cg_desc, info )
     if(info /= 0) call cg_scalapack_err_exit(info,"PDPOTRI")

     ! M = tmp = work * another_tmp
     call PDGEMM('N','N',n_evs,n_vectors,n_evs, 1.0d0, work, 1, 1, cg_desc, &
          another_tmp, 1, 1, cg_desc, 0.0d0, tmp, 1, 1, cg_desc)

     ! another_tmp = converged_evs * tmp
     call PDGEMM('N','N',n_basis,n_vectors,n_evs, 1.0d0, converged_evs, 1, 1, cg_desc, &
          tmp, 1, 1, cg_desc, 0.0d0, another_tmp, 1, 1, cg_desc)

     do i_col = 1, n_vectors, 1

        if (l_col(i_col) == 0) cycle
        ritz_basis(:,l_col(i_col)) = ritz_basis(:,l_col(i_col)) - another_tmp(:,l_col(i_col))

     end do

     deallocate(ovlp_ev)
     deallocate(work)

  end if

end subroutine project_ritz_basis_scalapack_real
!******
!-----------------------------------------------------------
!****s* cg_scalapack/orthonormalize_ritz_basis_scalapack_real
!  NAME
!    orthonormalize_ritz_basis_scalapack_real
!  SYNOPSIS
subroutine orthonormalize_ritz_basis_scalapack_real( n_vectors, ritz_basis )
!  PURPOSE
!    Orthonormalizes the Ritz basis with ScaLAPACK, real version
!    This routine uses a modified Gram-Schmidt orthogonalization
!    It needs two additional work matrices but is fast
!    because most work is done in matrix-matrix-products
!
!    This routine must be called after the overlap matrix is saved !!
!
!  USES
  use dimensions, only: n_basis
  use localorb_io, only: localorb_info, OL_norm, use_unit
  use mpi_tasks, only: check_allocation
!  ARGUMENTS
  integer :: n_vectors
  real*8, dimension(mxld,mxcol) :: ritz_basis
!  INPUTS
!    n_vectors - number of vectors in the Ritz basis
!  OUTPUT
!    -
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


  real*8, allocatable, dimension(:) :: prod
  real*8, allocatable, dimension(:,:) :: ovlp_ev, work
  real*8 :: dotprod, fact

  integer:: i_spin, i_row, i_col, n_block, i_done

  character*200 :: info_str

  if (lopbcg_debug) then
     write(info_str,'(2X,A)') "Orthonormalizing Ritz basis"
     call localorb_info(info_str,use_unit,'(A)',OL_norm)
  end if

  n_block = npcol*nb

  allocate(prod(mxld),stat=info)
  call check_allocation(info, 'prod                          ')

  allocate(ovlp_ev(mxld,mxcol),stat=info)
  call check_allocation(info, 'ovlp_ev                       ')

  allocate(work(mxld,mxcol),stat=info)
  call check_allocation(info, 'work                          ')


    if(my_scalapack_id<npcol*nprow) then ! The main work is done only on the working set

       ! ovlp_ev = ovlp_save * ritz_basis

       call PDSYMM('L', 'U', n_basis, n_vectors, 1.0d0, ovlp_save, 1, 1, cg_desc, &
            ritz_basis(1,1), 1, 1, cg_desc, &
            0.0d0, ovlp_ev, 1, 1, cg_desc)

      i_done = 0 ! Number of vectors which are orthogonalized against all others

      do i_col = 1, n_vectors

         if(i_col > i_done+1) then

            ! Build the matrix-dot product of ritz_basis(i_col) with all vectors
            ! against which ritz_basis(i_col) is not yet orthogonalized,
            ! these are the ones from i_done+1 to i_col-1

            call PDGEMV('T', n_basis, i_col-1-i_done, 1.0d0, ovlp_ev, 1, i_done+1, cg_desc, &
                        ritz_basis(1,1), 1, i_col, cg_desc, 1, &
                        0.0d0, prod, 1, 1, cg_desc, 1)

            ! Orthogonalize basis against the others, keep ovlp_save * ritz_basis up to date

            ! ritz_basis(i_col) = ritz_basis(i_col) - ritz_basis(i_done+1:i_col-1) * prod

            call PDGEMV('N', n_basis, i_col-1-i_done, -1.0d0, ritz_basis(1,1), 1, i_done+1, cg_desc, &
                        prod, 1, 1, cg_desc, 1, &
                        1.0d0, ritz_basis(1,1), 1, i_col, cg_desc, 1)

            ! The same for ovlp_save * ritz_basis

            call PDGEMV('N', n_basis, i_col-1-i_done, -1.0d0, ovlp_ev, 1, i_done+1, cg_desc, &
                        prod, 1, 1, cg_desc, 1, &
                        1.0d0, ovlp_ev, 1, i_col, cg_desc, 1)

         endif

         ! Now eigenvec(i_col) is orthogonalized against all previous ones, normalize it

         ! dotprod = ovlp_ev(i_col)**T * eigenvec(i_col)

         dotprod = 0
         call PDDOT(n_basis, dotprod, ovlp_ev, 1, i_col, cg_desc, 1, ritz_basis(1,1), 1, i_col, cg_desc, 1)

         if(l_col(i_col) > 0) then

            ! eigenvec(i_col) = eigenvec(i_col)/sqrt(dotprod)

            if(dotprod>0) then
               fact = 1.d0/sqrt(dotprod)
            else
               fact = 1.d0 ! for safety only, should never happen!
            endif

            ritz_basis(:,l_col(i_col)) = ritz_basis(:,l_col(i_col))*fact
            ovlp_ev(:,l_col(i_col)) = ovlp_ev(:,l_col(i_col))*fact

         endif

         ! If i_col-i_done reaches block size, orthogonalize the vectors i_done+1 .. i_col
         ! against all following ones. This can be done with matrix-matrix operations

         if(i_col-i_done == n_block .and. i_col<n_vectors) then

            ! Build the matrix-dot product of ritz_basis(i_done+1..i_col) with ritz_basis(i_col+1 .. n_vectors)

            call PDGEMM('T', 'N', n_block, n_vectors-i_col, n_basis, &
                        1.0d0, ovlp_ev, 1, i_done+1, cg_desc, &
                        ritz_basis(1,1), 1, i_col+1, cg_desc, &
                        0.0d0, work, 1, i_col+1, cg_desc)

            ! Orthogonalize

            call PDGEMM('N', 'N', n_basis, n_vectors-i_col, n_block, &
                 -1.0d0, ritz_basis(1,1), 1, i_done+1, cg_desc, &
                 work, 1, i_col+1, cg_desc, &
                 1.0d0, ritz_basis(1,1), 1, i_col+1, cg_desc)

            i_done = i_done + n_block

         endif

      enddo ! i_col

    endif ! work only on working set

  ! TEST only, may be removed later:
  ! call check_ev_orthogonality_real

  deallocate(prod)
  deallocate(ovlp_ev)
  deallocate(work)

end subroutine orthonormalize_ritz_basis_scalapack_real
!******
!-----------------------------------------------------------
!****s* cg_scalapack/check_rb_orthogonality_real
!  NAME
!    check_rb_orthogonality_real
!  SYNOPSIS
subroutine check_rb_orthogonality_real( n_vectors, ritz_basis )
!  USES
  use dimensions, only: n_basis
  use localorb_io, only: use_unit
  use mpi_tasks, only: myid
  implicit none
!  ARGUMENTS
  integer :: n_vectors
  real*8, dimension(mxld,mxcol) :: ritz_basis
!  PURPOSE
!    Checks the orthogonality of the eigenvectors
!    This is a test routine only!!!!!!!!!!!!

   integer i_col, i_row
   real*8 :: amax, dmax
   real*8, allocatable :: work1(:,:), work2(:,:)

   allocate(work1(mxld,mxcol))
   allocate(work2(mxld,mxcol))

   ! The work in this routine must be done only on the working set
   if(my_scalapack_id>=npcol*nprow) return

   call PDSYMM('L', 'U', n_basis, n_vectors, 1.0d0, ovlp_save, 1, 1, cg_desc, &
        ritz_basis(1,1), 1, 1, cg_desc, &
        0.0d0, work1, 1, 1, cg_desc)

   call PDGEMM('T', 'N', n_vectors, n_vectors, n_basis, &
        1.0d0, ritz_basis(1,1), 1, 1, cg_desc, &
        work1, 1, 1, cg_desc, &
        0.0d0, work2, 1, 1, cg_desc)

   dmax = 0.
   amax = 0.

   do i_col = 1, n_vectors
      do i_row = 1, n_vectors
         if(l_row(i_row) > 0 .and. l_col(i_col)>0) then
            if(i_row==i_col) then
               ! diagonal element
               dmax = max(dmax,abs(1.0-work2(l_row(i_row),l_col(i_col))))
            else
               ! off diagonal element
               amax = max(amax,abs(work2(l_row(i_row),l_col(i_col))))
            endif
         endif
      enddo
   enddo
   
!   print *,'check_ev_orthogonality_real, ID: ',myid,' Errors: ',amax, dmax
   write(use_unit,'(A,I4,A,F10.5,F10.5)') 'check_rb_orthogonality_real, ID: ',myid,' Errors: ',amax, dmax
   deallocate(work1)
   deallocate(work2)

 end subroutine check_rb_orthogonality_real
!******
!-----------------------------------------------------------------------------------
!****s* cg_scalapack/scalapack_err_exit
!  NAME
!    scalapack_err_exit
!  SYNOPSIS
 subroutine cg_scalapack_err_exit(info, name)
!  PURPOSE
!    Exits ScaLAPACK after an error.
!  USES
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks, only: aims_stop
!  ARGUMENTS
    integer :: info
    character*(*) :: name
!  INPUTS
!    o info -- ScaLAPACK error code
!    o name -- name of the routine where error occurred
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


    character*100 :: info_str
!    integer :: mpierr

    write (info_str,'(2X,A,I5,A,A)') 'Error ',info,' in ',name
    call localorb_info(info_str,use_unit,'(A)')

    call aims_stop

!    call MPI_Finalize(mpierr)
!    stop
  end subroutine cg_scalapack_err_exit

!******
!-----------------------------------------------------------------------------------
!****s* cg_solver/set_full_matrix_real
!  NAME
!    set_full_matrix_real
!  SYNOPSIS
  subroutine cg_set_full_matrix_real( mat )
!  PURPOSE
!    Sets the lower half of a distributed matrix from the upper half
!  USES
    use dimensions, only: n_basis
    implicit none
!  ARGUMENTS
    real*8, dimension(mxld, mxcol) :: mat
!  INPUTS
!  OUTPUT
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer i_col, i_row

    ! This routine is called only from the working set, so no need to check here

    call pdtran(n_basis,n_basis,1.d0,mat,1,1,cg_desc,0.d0,tmp,1,1,cg_desc)

    do i_col=1,n_basis-1
       if(l_col(i_col)==0) cycle
       do i_row=i_col+1,n_basis
          if(l_row(i_row)>0) mat(l_row(i_row),l_col(i_col)) = tmp(l_row(i_row),l_col(i_col))
       enddo
    enddo

  end  subroutine cg_set_full_matrix_real
!******
end module cg_scalapack
