!****h* FHI-aims/cg
!  NAME
!    cg
!  SYNOPSIS

module cg

!  PURPOSE
!  ??????????????????????????????
!
!  USES

  use mpi_tasks
  use runtime_choices
  use dimensions
  use separate_core_states
  use timing, only: number_of_loops
  use localorb_io
  use physics, only: rho_change
  use grids, only: batches
  use basis, only: basis_wave_ordered
  use pbc_lists, only: inv_centers_basis_integrals,centers_basis_integrals
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


  real*8, allocatable, dimension(:) :: hartree_overlap_matrix
  real*8, allocatable, dimension(:,:) :: stored_KS_eigenvector
  real*8, allocatable, dimension(:) :: stored_KS_eigenvalue
  real*8, allocatable, dimension(:) :: cholesky_factor

  logical :: factor_cg_overlap = .true.

!******
contains

!-------------------------------------------------------------------------------------------
!****s* cg/cg_solver
!  NAME
!    cg_solver
!  SYNOPSIS

  subroutine cg_solver( n_basis_current, n_states_current, overlap_matrix_current, hamiltonian, &
       KS_eigenvalue, KS_eigenvector, i_spin )

!  PURPOSE
!  ??????????????????????????
!
!  USES

    use physics, only: occ_numbers
    use inner_product, only: S_inner_product

!  ARGUMENTS

    integer :: n_basis_current, n_states_current
    real*8 overlap_matrix_current(n_basis_current*(n_basis_current+1)/2)
    real*8 hamiltonian(n_basis_current*(n_basis_current+1)/2)
    integer :: i_spin
    real*8 KS_eigenvalue (n_states_current)
    real*8 KS_eigenvector (n_basis_current, n_states_current)


! INPUTS
! o  n_basis_current -- ????????????
! o  n_states_current -- ????????????
! o  overlap_matrix_current -- ????????????
! o  hamiltonian -- Hamiltonian matrix
! o  i_spin -- spin index
!
! OUTPUT
! o  KS_eigenvalue -- Kohn-Sham eigenvalues
! o  KS_eigenvector -- Kohn-Sham eigenvectors
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE





    ! locals
    real*8 hamiltonian_store(n_basis_current*(n_basis_current+1)/2)

    real*8, allocatable, dimension(:) :: search_direction
    real*8, allocatable, dimension(:) :: gradient
    real*8, allocatable, dimension(:) :: prec_gradient
    real*8, allocatable, dimension(:) :: temp_vector
    real*8, allocatable, dimension(:,:) :: vectors_y, vectors_x, vectors_w
    real*8, allocatable, dimension(:,:) :: ritz_basis

    real*8, allocatable, dimension(:) :: new_eigenvector

    integer :: i_state, i_state_2, i_index, i_index_2, current_state
    integer :: n_active_states
    integer, dimension(n_states_current) :: active_states
    integer :: n_occupied_states
    integer, dimension(n_states_current) :: occupied_states
    integer, dimension(n_basis_current) :: ipiv

    real*8 :: gamma
    real*8 :: tau1, tau2, tau3, mu1, mu2
    real*8 :: prec_value

    real*8 :: yHx, ySx, yHw, ySw, xHw, xSw, wHw, wSw, xHx, xSx, yHy, ySy
    real*8 :: c0, c1, c2, det

    real*8 :: min_ev, ev_shift, lopcg_tol_current

    real*8 :: alpha = 1.0e-6
    real*8 :: beta = 1.0d0

    real*8, external :: ddot

    integer :: max_ritz_dim, current_ritz_dim
    real*8, allocatable, dimension(:,:) :: ritz_matrix_H
    real*8, allocatable, dimension(:,:) :: ritz_matrix_S
    real*8, allocatable, dimension(:) :: ritz_values
    real*8, allocatable, dimension(:) :: tau
    integer :: lwork
    real*8, allocatable, dimension(:) :: work
    integer :: info

    logical, allocatable, dimension(:) :: converged_state

    real*8 :: residual_norm(n_states_current)

    integer :: n_blocks, remainder_block_size, current_block_size, this_block_size
    integer :: state_top, state_bottom
    integer :: i_iter, i_block, i_vector, i_vector_2
    integer :: total_iter

    logical :: include_state

    character*120 :: info_str

!!$    if (number_of_loops == initial_ev_solutions) then
!!$
!!$       if (.not.allocated(stored_KS_eigenvector)) then
!!$          allocate(stored_KS_eigenvector(n_basis_current, n_states_current))
!!$       end if
!!$       if (.not.allocated(stored_KS_eigenvalue)) then
!!$          allocate(stored_KS_eigenvalue(n_states_current))
!!$       end if
!!$
!!$       stored_KS_eigenvector = KS_eigenvector
!!$       stored_KS_eigenvalue = KS_eigenvalue
!!$
!!$       return
!!$
!!$    end if



    if (.not.allocated(stored_KS_eigenvector)) then
       allocate(stored_KS_eigenvector(n_basis_current, n_states_current),stat=info)
       call check_allocation(info, 'stored_KS_eigenvector         ')
    end if
    stored_KS_eigenvector = 0.0d0

    if (.not.allocated(stored_KS_eigenvalue)) then
       allocate(stored_KS_eigenvalue(n_states_current),stat=info)
       call check_allocation(info, 'stored_KS_eigenvalue          ')
    end if

    if (.not.allocated(cholesky_factor)) then
       allocate(cholesky_factor(n_basis_current*(n_basis_current+1)/2),stat=info)
       call check_allocation(info, 'cholesky_factor               ')
    end if

    max_ritz_dim = 3*lopcg_block_size
    lwork = 3*max_ritz_dim

    ! allocate things
    allocate(gradient(n_basis_current),stat=info)
    call check_allocation(info, 'gradient                      ')

    allocate(prec_gradient(n_basis_current),stat=info)
    call check_allocation(info, 'prec_gradient                 ')

    allocate(search_direction(n_basis_current),stat=info)
    call check_allocation(info, 'search_direction              ')

    allocate(new_eigenvector(n_basis_current),stat=info)
    call check_allocation(info, 'new_eigenvector               ')

    allocate(temp_vector(n_basis_current),stat=info)
    call check_allocation(info, 'temp_vector                   ')

    allocate(vectors_y(n_basis_current,lopcg_block_size),stat=info)
    call check_allocation(info, 'vectors_y                     ')

    allocate(vectors_x(n_basis_current,lopcg_block_size),stat=info)
    call check_allocation(info, 'vectors_x                     ')

    allocate(vectors_w(n_basis_current,lopcg_block_size),stat=info)
    call check_allocation(info, 'vectors_w                     ')

    allocate(ritz_basis(n_basis_current, max_ritz_dim),stat=info)
    call check_allocation(info, 'ritz_basis                    ')

    allocate(ritz_matrix_H(max_ritz_dim,max_ritz_dim),stat=info)
    call check_allocation(info, 'ritz_matrix_H                 ')

    allocate(ritz_matrix_S(max_ritz_dim,max_ritz_dim),stat=info)
    call check_allocation(info, 'ritz_matrix_S                 ')

    allocate(ritz_values(max_ritz_dim),stat=info)
    call check_allocation(info, 'ritz_values                   ')

    allocate(tau(max_ritz_dim),stat=info)
    call check_allocation(info, 'tau                           ')

    allocate(work(lwork),stat=info)
    call check_allocation(info, 'work                          ')

    allocate(converged_state(n_states_current) ,stat=info)
    call check_allocation(info, 'converged_state               ')


!    call integrate_hartree_overlap()

    write (info_str,'(2X,A,I4,A,I4)') &
         "Solving the eigenvalue problem with LOPCG method for states ", &
         n_core_states+1, " - ", n_core_states + n_states_current
    call localorb_info(info_str,use_unit,'(A)')

    hamiltonian_store = hamiltonian

    min_ev = 1.0d6
    do i_state = 1, n_states_current, 1
       min_ev = MIN(min_ev, KS_eigenvalue(i_state))
    end do

    if (min_ev < 0.0d0) then
       ev_shift = -min_ev + 0.1
    else
       ev_shift = 0.0d0
    end if

    if (lopcg_adaptive_tolerance) then
       lopcg_tol_current = MAX(0.01*(SUM(ABS(rho_change))/dble(n_spin)),lopcg_tol)
       write (info_str,'(2X,A,E10.4)') "Current LOPCG tolerance: ", lopcg_tol_current
       call localorb_info(info_str,use_unit,'(A)')
    else
       lopcg_tol_current = lopcg_tol
    end if

    if (factor_cg_overlap) then
       cholesky_factor = overlap_matrix_current
       !    cholesky_factor = hamiltonian - MINVAL(KS_eigenvalue)*overlap_matrix_current
       call dpptrf('U', n_basis_current, cholesky_factor, info)
       factor_cg_overlap = .false.
    end if


    ! compute the gradient

    n_active_states = 0
    n_occupied_states = 0

    do i_state = 1, n_states_current, 1
       n_active_states = n_active_states + 1
       active_states(n_active_states) = i_state
       if (occ_numbers(n_core_states + i_state,i_spin,1) > 0.0d0) then
          n_occupied_states = n_occupied_states + 1
          occupied_states(n_occupied_states) = i_state
       end if
    end do

    total_iter = 0
    i_iter = 0
    converged_state = .false.

    iter_loop: do while (i_iter <= max_cg_iterations)
       i_iter = i_iter + 1

       state_top = 0
       n_blocks = 0

       block_loop: do while (state_top < n_active_states)

          n_blocks = n_blocks + 1

          state_bottom = state_top + 1

          if (lopcg_auto_blocksize) then
             current_block_size = 1
             if (state_bottom + current_block_size - 1 < n_active_states) then
                do while ((ABS(KS_eigenvalue(state_bottom) - KS_eigenvalue(state_bottom + current_block_size)) < &
                     0.01*ABS(KS_eigenvalue(state_bottom))).and. &
                     (state_bottom + current_block_size - 1 < n_active_states).and. &
                     (current_block_size < lopcg_block_size))
                   current_block_size = current_block_size + 1
                end do
             end if
!!$             i_state = state_bottom
!!$             current_block_size = 0
!!$             do while (((state_bottom + current_block_size - 1) < n_active_states) &
!!$                  .and.(current_block_size < lopcg_block_size))
!!$                this_block_size = 0
!!$                do while ((ABS(KS_eigenvalue(i_state) - KS_eigenvalue(i_state + this_block_size)) < &
!!$                     0.01*ABS(KS_eigenvalue(i_state))).and.(i_state + this_block_size < n_active_states))
!!$                   this_block_size = this_block_size + 1
!!$                end do
!!$                if (this_block_size == 0) this_block_size = 1
!!$                if ((current_block_size + this_block_size <= lopcg_block_size).and. &
!!$                     ((state_bottom + current_block_size + this_block_size - 1) <= n_active_states)) then
!!$                   current_block_size = current_block_size + this_block_size
!!$                   i_state = state_bottom + current_block_size
!!$                else
!!$                   exit
!!$                end if
!!$             end do
!!$             if (current_block_size == 0) current_block_size = 1
          else
             current_block_size = MIN(lopcg_block_size, n_active_states - state_top)
          end if
!!$       if (myid == 0) then
!!$          print *, 'In block: ', n_blocks, ' block size: ', current_block_size
!!$       end if

          if (i_iter == 1) then
             current_ritz_dim = 2*current_block_size
          else
             current_ritz_dim = 3*current_block_size
          end if

          state_top = state_top + current_block_size

          if ((myid == 0).and.(i_iter == 1)) then
             print *, 'Iteration 1, in block: ', n_blocks, ' states: ', n_core_states + state_bottom, " - ", &
                  n_core_states + state_top
          end if

          if (ALL(converged_state(state_bottom:state_top))) then
             write (info_str,'(2X,A,I6,A)') "All states converged in block ", n_blocks, &
                  ". Skipping the block."
             call localorb_info(info_str,use_unit,'(A)')
             cycle block_loop
          end if

!!$          min_ev = 1.0d6
!!$          do i_state = state_bottom, state_top, 1
!!$             min_ev = MIN(min_ev, KS_eigenvalue(i_state))
!!$          end do
!!$
!!$          if (min_ev < 0.0d0) then
!!$             ev_shift = -min_ev + 0.1d0
!!$          else
!!$             ev_shift = 0.1d0
!!$          end if

          hamiltonian = hamiltonian_store + ev_shift*overlap_matrix_current
          KS_eigenvalue(state_bottom:state_top) = KS_eigenvalue(state_bottom:state_top) + ev_shift

!!$       iter_loop: do i_iter = 1, max_cg_iterations, 1

          inblock_loop: do i_vector = 1, current_block_size, 1
             i_state = state_bottom + i_vector - 1
             current_state = active_states(i_state)

             ! compute the gradient
             gradient = 0.0d0
             call  DSPMV('U', n_basis_current, 1.0d0, hamiltonian, KS_eigenvector(:,current_state), 1, &
                  0.0d0, gradient(:), 1)
             call  DSPMV('U', n_basis_current, -1.0d0*KS_eigenvalue(current_state), &
                  overlap_matrix_current, KS_eigenvector(:,current_state), 1, &
                  1.0d0, gradient(:), 1)

             ! check if all states in the block are converged
             residual_norm(current_state) = sqrt(ddot(n_basis_current, gradient, 1, gradient, 1))
             if (residual_norm(current_state) < lopcg_tol_current) converged_state(current_state) = .true.
             if (ALL(converged_state(state_bottom:state_top))) exit inblock_loop

             select case(cg_preconditioner)

             case (prec_diagonal)

                i_state_2 = 0
                do i_index = 1, n_basis_current, 1
                   i_state_2 = i_state_2 + i_index
                   prec_value = 2.0d0*occ_numbers(n_core_states+current_state,i_spin,1)* &
                        (hamiltonian(i_state_2) - &
                        KS_eigenvalue(current_state)*overlap_matrix_current(i_state_2))
                   prec_value = MAX(abs(prec_value), 1.0d0)
                   prec_gradient(i_index) = beta*gradient(i_index) / prec_value
                end do

             case (prec_inv_ovlp)

                prec_gradient = gradient
                call DPPTRS( 'U', n_basis_current, 1, cholesky_factor, &
                     prec_gradient, n_basis_current, info )

             case default

                prec_gradient = gradient

             end select

             search_direction = -prec_gradient

             vectors_y(:,i_vector) = stored_KS_eigenvector(:,current_state)
             vectors_x(:,i_vector) = KS_eigenvector(:,current_state)
             vectors_w(:,i_vector) = search_direction
          end do inblock_loop

          if (ALL(converged_state(state_bottom:state_top))) then
             write (info_str,'(2X,A,I6,A)') "All states converged in block ", n_blocks, &
                  ". Skipping the block."
             call localorb_info(info_str,use_unit,'(A)')
             KS_eigenvalue(state_bottom:state_top) = KS_eigenvalue(state_bottom:state_top) - ev_shift
             cycle block_loop
          end if

          call lock_ritz_space( n_basis_current, n_states_current, &
               KS_eigenvector, overlap_matrix_current, vectors_x , &
               vectors_w, vectors_y, active_states, state_bottom, current_block_size, i_iter )

          ! build the Ritz-basis matrix
          i_vector_2 = 0
          do i_vector = 1, current_block_size, 1
             i_vector_2 = i_vector_2 + 1
             ritz_basis(:,i_vector_2) = vectors_x(:,i_vector)
          end do
          do i_vector = 1, current_block_size, 1
             i_vector_2 = i_vector_2 + 1
             ritz_basis(:,i_vector_2) = vectors_w(:,i_vector)
          end do
          do i_vector = 1, current_block_size, 1
             i_vector_2 = i_vector_2 + 1
             ritz_basis(:,i_vector_2) = vectors_y(:,i_vector)
          end do

          ! compute the Ritz matrices (H and S)
          do i_index = 1, current_ritz_dim, 1
             call DSPMV('U', n_basis_current, 1.0d0, hamiltonian, ritz_basis(:,i_index), 1, &
                  0.0d0, temp_vector , 1)
             do i_index_2 = 1, i_index, 1
                ritz_matrix_H(i_index_2, i_index) = &
                     ddot(n_basis_current, ritz_basis(:,i_index_2), 1, temp_vector, 1)
             end do
          end do

          do i_index = 1, current_ritz_dim, 1
             call DSPMV('U', n_basis_current, 1.0d0, overlap_matrix_current, ritz_basis(:,i_index), 1, &
                  0.0d0, temp_vector , 1)
             do i_index_2 = 1, i_index, 1
                ritz_matrix_S(i_index_2, i_index) = &
                     ddot(n_basis_current, ritz_basis(:,i_index_2), 1, temp_vector, 1)
             end do
          end do

          call dsygv(1, 'V', 'U', current_ritz_dim, ritz_matrix_H(1:current_ritz_dim,1:current_ritz_dim), &
               current_ritz_dim, ritz_matrix_S(1:current_ritz_dim,1:current_ritz_dim), &
               current_ritz_dim, ritz_values, work, lwork, info)
          if(info /= 0) then
             print *, myid, 'iter, block, dim, info: ', i_iter, n_blocks, current_ritz_dim, info
          end if

          do i_vector = 1, current_block_size, 1

             current_state = active_states(state_bottom - 1 + i_vector)

             tau(1:current_ritz_dim) = ritz_matrix_H(1:current_ritz_dim,i_vector)

             ! compute and normalize the new eigenvector (= the ritz vector)
             !          new_eigenvector(:) = tau3*vector_y + tau2*vector_x + tau1*vector_w

             call DGEMV('N', n_basis_current, current_ritz_dim, 1.0d0, ritz_basis, n_basis_current, &
                  tau, 1, 0.0d0, new_eigenvector, 1)
             call DSPMV('U', n_basis_current, 1.0d0, overlap_matrix_current, new_eigenvector, 1, &
                  0.0d0, temp_vector, 1)
             new_eigenvector = new_eigenvector / sqrt(ddot(n_basis_current, new_eigenvector, 1, temp_vector, 1))
             KS_eigenvector(:,current_state) = new_eigenvector(:)

             tau(1:current_ritz_dim-current_block_size) = &
                  ritz_matrix_H(current_block_size+1:current_ritz_dim,i_vector)
             call DGEMV('N', n_basis_current, current_ritz_dim-current_block_size, &
                  1.0d0, ritz_basis(1,current_block_size+1), n_basis_current, &
                  tau, 1, 0.0d0, new_eigenvector, 1)
             call DSPMV('U', n_basis_current, 1.0d0, overlap_matrix_current, new_eigenvector, 1, &
                  0.0d0, temp_vector, 1)
             new_eigenvector = new_eigenvector / sqrt(ddot(n_basis_current, new_eigenvector, 1, temp_vector, 1))
             stored_KS_eigenvector(:,current_state) = new_eigenvector(:)

             ! update the eigenvalues
!!$             stored_KS_eigenvalue(current_state) = KS_eigenvalue(current_state)
             call DSPMV('U', n_basis_current, 1.0d0, hamiltonian, KS_eigenvector(:,current_state), 1, &
                  0.0d0, temp_vector, 1)
             KS_eigenvalue(current_state) = &
                  ddot(n_basis_current, KS_eigenvector(:,current_state), 1, temp_vector, 1)
             call DSPMV('U', n_basis_current, 1.0d0, &
                  hamiltonian_store - KS_eigenvalue(current_state)*overlap_matrix_current, &
                  KS_eigenvector(:,current_state), 1, 0.0d0, temp_vector, 1)
             residual_norm(current_state) = sqrt(ddot(n_basis_current, temp_vector, 1, temp_vector, 1))

          end do

          total_iter = total_iter + 1

          KS_eigenvalue(state_bottom:state_top) = KS_eigenvalue(state_bottom:state_top) - ev_shift

       end do block_loop

!       if ((ALL(converged_state(state_bottom:state_top))).and.(block_iter == 0)) then
!!$          if (myid == 0) then
!!$             print *, 'Here: ', i_block
!!$          end if
!!$          do i_vector = 1, current_block_size, 1
!!$             current_state = active_states(state_bottom - 1 + i_vector)
!!$!             stored_KS_eigenvector(:,current_state) = KS_eigenvector(:,current_state)
!!$             stored_KS_eigenvalue(current_state) = KS_eigenvalue(current_state)
!!$             call DSPMV('U', n_basis_current, 1.0d0, hamiltonian, KS_eigenvector(:,current_state), 1, &
!!$                  0.0d0, temp_vector, 1)
!!$             KS_eigenvalue(current_state) = &
!!$                  ddot(n_basis_current, KS_eigenvector(:,current_state), 1, temp_vector, 1) - ev_shift
!!$             call DSPMV('U', n_basis_current, 1.0d0, &
!!$                  hamiltonian - (KS_eigenvalue(current_state) + ev_shift)*overlap_matrix_current, &
!!$                  KS_eigenvector(:,current_state), 1, 0.0d0, temp_vector, 1)
!!$             residual_norm(current_state) = sqrt(ddot(n_basis_current, temp_vector, 1, temp_vector, 1))
!!$          end do
!       end if

!!$       do i_vector = 1, current_block_size, 1
!!$
!!$          current_state = state_bottom - 1 + i_vector
!!$
!!$          stored_KS_eigenvalue(current_state) = KS_eigenvalue(current_state)
!!$          call DSPMV('U', n_basis_current, 1.0d0, hamiltonian, KS_eigenvector(:,current_state), 1, &
!!$               0.0d0, temp_vector, 1)
!!$          KS_eigenvalue(current_state) = &
!!$               ddot(n_basis_current, KS_eigenvector(:,current_state), 1, temp_vector, 1) - ev_shift
!!$          call DSPMV('U', n_basis_current, 1.0d0, hamiltonian - (KS_eigenvalue(current_state) + ev_shift)*overlap_matrix_current, &
!!$               KS_eigenvector(:,i_state), 1, 0.0d0, temp_vector, 1)
!!$          residual_norm(current_state) = sqrt(ddot(n_basis_current, temp_vector, 1, temp_vector, 1))
!!$
!!$       end do
!       call orthonormalize_eigenvectors()

       ! finally, get the eigenvalues
       do i_state = 1, n_states_current, 1
          stored_KS_eigenvalue(i_state) = KS_eigenvalue(i_state)
!!$          call DSPMV('U', n_basis_current, 1.0d0, hamiltonian, KS_eigenvector(:,i_state), 1, &
!!$               0.0d0, temp_vector, 1)
!!$          KS_eigenvalue(i_state) = &
!!$               ddot(n_basis_current, KS_eigenvector(:,i_state), 1, temp_vector, 1) - ev_shift
!!$          call DSPMV('U', n_basis_current, 1.0d0, &
!!$               hamiltonian - (KS_eigenvalue(i_state) + ev_shift)*overlap_matrix_current, &
!!$               KS_eigenvector(:,i_state), 1, 0.0d0, temp_vector, 1)
!!$          residual_norm(i_state) = sqrt(ddot(n_basis_current, temp_vector, 1, temp_vector, 1))
       end do

       if (MAXVAL(residual_norm(1:n_states_current)) < lopcg_tol_current) exit iter_loop
!!$       if (myid == 0) then
!!$          print *,'iteration: ', i_iter
!!$          do i_state = 1, n_states_current, 1
!!$             print *,'i_state: ', i_state, '|| r ||: ', residual_norm(i_state)
!!$          end do
!!$       end if

!!$       if (MAXVAL(residual_norm(occupied_states(1:n_occupied_states))) < lopcg_tol_current) exit

    end do iter_loop ! end loop over ev iterations

    ! finally, get the eigenvalues
    do i_state = 1, n_active_states, 1
       current_state = active_states(i_state)
       stored_KS_eigenvalue(current_state) = KS_eigenvalue(current_state)
       call DSPMV('U', n_basis_current, 1.0d0, hamiltonian_store, KS_eigenvector(:,current_state), 1, &
            0.0d0, temp_vector, 1)
       KS_eigenvalue(current_state) = &
            ddot(n_basis_current, KS_eigenvector(:,current_state), 1, temp_vector, 1)
       call DSPMV('U', n_basis_current, 1.0d0, &
            hamiltonian_store - KS_eigenvalue(current_state)*overlap_matrix_current, &
            KS_eigenvector(:,current_state), 1, 0.0d0, temp_vector, 1)
       residual_norm(current_state) = sqrt(ddot(n_basis_current, temp_vector, 1, temp_vector, 1))
    end do

    if (MAXVAL(residual_norm(active_states(1:n_active_states))) >= lopcg_tol_current) then
       write (info_str,'(2X,A,I6,A)') "LOPCG didn't converge in ", max_cg_iterations, &
            " iterations per block."
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(2X,A,E10.4,A,I4)') "Maximum remaining residual: ", &
            MAXVAL(residual_norm(active_states(1:n_active_states))),  " in state: ", &
            n_core_states + MAXLOC(residual_norm(active_states(1:n_active_states)))
       call localorb_info(info_str,use_unit,'(A)')
    else
!!$       if (remainder_block_size == 0) then
          write (info_str,'(2X,A,I6,A,I6,A)') "LOPCG converged in total of ", total_iter, &
               " iterations for ", n_blocks, " blocks."
          call localorb_info(info_str,use_unit,'(A)')
!!$       else
!!$          write (info_str,'(2X,A,I6,A,I6,A)') "LOPCG converged in total of ", total_iter, &
!!$               " iterations for ", n_blocks+1, " blocks."
!!$          call localorb_info(info_str,use_unit,'(A)')
!!$       end if
       write (info_str,'(2X,A,E10.4,A,I4)') "Maximum remaining residual: ", &
            MAXVAL(residual_norm(active_states(1:n_active_states))), " in state: ", &
            n_core_states + MAXLOC(residual_norm(active_states(1:n_active_states)))
       call localorb_info(info_str,use_unit,'(A)')
    end if

    ! deallocate things

    deallocate(converged_state)
    deallocate(work)
    deallocate(tau)
    deallocate(ritz_values)
    deallocate(ritz_matrix_S)
    deallocate(ritz_matrix_H)
    deallocate(ritz_basis)
    deallocate(vectors_w)
    deallocate(vectors_x)
    deallocate(vectors_y)
    deallocate(temp_vector)
    deallocate(new_eigenvector)
    deallocate(search_direction)
    deallocate(prec_gradient)
    deallocate(gradient)


  end subroutine cg_solver
!******
!------------------------------------------------------------------------------------
!****s* cg/lock_ritz_space
!  NAME
!    lock_ritz_space
!  SYNOPSIS

  subroutine lock_ritz_space( n_basis_current, n_states_current, &
       KS_eigenvector, overlap_matrix_current, vectors1, vectors2, &
       vectors3, active_states, current_state, block_size, current_iter )

!  PURPOSE
!  ???????????????????????????????????
!
!  USES

    use inner_product

!  ARGUMENTS

    integer :: n_basis_current, n_states_current

    real*8 :: KS_eigenvector (n_basis_current, n_states_current)
    real*8 :: overlap_matrix_current(n_basis_current*(n_basis_current+1)/2)

    real*8, dimension(n_basis_current,lopcg_block_size) :: vectors1, vectors2, vectors3
    integer, dimension(n_states_current) :: active_states
    integer :: current_state, block_size, current_iter


! INPUTS
! o n_basis_current -- ????????????
! o n_states_current -- ????????????
! o KS_eigenvector -- ????????????
! o overlap_matrix_current -- ????????????
! o vectors1 -- ????????????
! o vectors2 -- ????????????
! o vectors3 -- ????????????
! o active_states -- ????????????
! o current_state -- ????????????
! o block_size -- ????????????
! o current_iter -- ????????????
!
!  OUTPUT
!    none ??????????????????????????
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




    ! locals
    integer :: i_state, i_vector, i_vector_2
    real*8 :: r

    do i_state = 1, current_state-1, 1
       do i_vector = 1, block_size, 1

          r = S_inner_product( vectors1(:,i_vector), KS_eigenvector(:, active_states(i_state)), &
               overlap_matrix_current, n_basis_current )
          vectors1(:,i_vector) = vectors1(:,i_vector) - r*KS_eigenvector(:, active_states(i_state))

          r = S_inner_product( vectors2(:,i_vector), KS_eigenvector(:, active_states(i_state)), &
               overlap_matrix_current, n_basis_current )
          vectors2(:,i_vector) = vectors2(:,i_vector) - r*KS_eigenvector(:, active_states(i_state))

          if (current_iter /= 1) then
             r = S_inner_product( vectors3(:,i_vector), KS_eigenvector(:, active_states(i_state)), &
                  overlap_matrix_current, n_basis_current )
             vectors3(:,i_vector) = vectors3(:,i_vector) - r*KS_eigenvector(:, active_states(i_state))
          end if
       end do
    end do

    do i_vector = 1, block_size, 1
       r = sqrt(S_inner_product( vectors1(:,i_vector), vectors1(:,i_vector), overlap_matrix_current, n_basis_current))
       vectors1(:,i_vector) = vectors1(:,i_vector) / r
       r = sqrt(S_inner_product( vectors2(:,i_vector), vectors2(:,i_vector), overlap_matrix_current, n_basis_current))
       vectors2(:,i_vector) = vectors2(:,i_vector) / r
       if (current_iter /= 1) then
          r = sqrt(S_inner_product( vectors3(:,i_vector), vectors3(:,i_vector), overlap_matrix_current, n_basis_current))
          vectors3(:,i_vector) = vectors3(:,i_vector) / r
       end if
    end do

    do i_vector = 1, block_size, 1
       do i_vector_2 = i_vector+1, block_size, 1
          r = S_inner_product( vectors1(:,i_vector), vectors1(:,i_vector_2), overlap_matrix_current, n_basis_current)
          vectors1(:,i_vector) = vectors1(:,i_vector) - r*vectors1(:,i_vector_2)
       end do
!!$       do i_vector_2 = 1, block_size, 1
!!$          r = S_inner_product( vectors1(:,i_vector), vectors2(:,i_vector_2), overlap_matrix_current, n_basis_current)
!!$          vectors1(:,i_vector) = vectors1(:,i_vector) - r*vectors2(:,i_vector_2)
!!$       end do
!!$       if (current_iter /= 1) then
!!$          do i_vector_2 = 1, block_size, 1
!!$             r = S_inner_product( vectors1(:,i_vector), vectors3(:,i_vector_2), overlap_matrix_current, n_basis_current)
!!$             vectors1(:,i_vector) = vectors1(:,i_vector) - r*vectors3(:,i_vector_2)
!!$          end do
!!$       end if
    end do

    do i_vector = 1, block_size, 1
       do i_vector_2 = i_vector+1, block_size, 1
          r = S_inner_product( vectors2(:,i_vector), vectors2(:,i_vector_2), overlap_matrix_current, n_basis_current)
          vectors2(:,i_vector) = vectors2(:,i_vector) - r*vectors2(:,i_vector_2)
       end do
!!$       if (current_iter /= 1) then
!!$          do i_vector_2 = 1, block_size, 1
!!$             r = S_inner_product( vectors2(:,i_vector), vectors3(:,i_vector_2), overlap_matrix_current, n_basis_current)
!!$             vectors2(:,i_vector) = vectors2(:,i_vector) - r*vectors3(:,i_vector_2)
!!$          end do
!!$       end if
    end do

    if (current_iter /= 1) then
       do i_vector = 1, block_size, 1
          do i_vector_2 = i_vector+1, block_size, 1
             r = S_inner_product( vectors3(:,i_vector), vectors3(:,i_vector_2), overlap_matrix_current, n_basis_current)
             vectors3(:,i_vector) = vectors3(:,i_vector) - r*vectors3(:,i_vector_2)
          end do
       end do
    end if

    do i_vector = 1, block_size, 1
       r = sqrt(S_inner_product( vectors1(:,i_vector), vectors1(:,i_vector), overlap_matrix_current, n_basis_current))
       vectors1(:,i_vector) = vectors1(:,i_vector) / r
       r = sqrt(S_inner_product( vectors2(:,i_vector), vectors2(:,i_vector), overlap_matrix_current, n_basis_current))
       vectors2(:,i_vector) = vectors2(:,i_vector) / r
       if (current_iter /= 1) then
          r = sqrt(S_inner_product( vectors3(:,i_vector), vectors3(:,i_vector), overlap_matrix_current, n_basis_current))
          vectors3(:,i_vector) = vectors3(:,i_vector) / r
       end if
    end do

  end subroutine lock_ritz_space
!******
!------------------------------------------------------------------------------------
!****s* cg/remove_core_states
!  NAME
!    remove_core_states
!  SYNOPSIS

  subroutine remove_core_states(n_basis_current, n_states_current, KS_eigenvector, &
       stored_KS_eigenvector, overlap_matrix_current )

!  PURPOSE
!  ????????????????????????????????????
!
!  USES

    use physics, only:overlap_matrix
    use separate_core_states
    use inner_product
    implicit none

!  ARGUMENTS

    integer :: n_basis_current, n_states_current
    real*8 :: KS_eigenvector(n_basis_current, n_states_current)
    real*8 :: stored_KS_eigenvector(n_basis_current, n_states_current)
    real*8 :: overlap_matrix_current(n_basis_current*(n_basis_current+1)/2)

!  INPUTS
!  o n_basis_current -- ?????????????
!  o n_states_current -- ?????????????
!  o KS_eigenvector -- ?????????????
!  o stored_KS_eigenvector -- ?????????????
!  o overlap_matrix_current -- ?????????????
!
!  OUTPUT
!    none ????UPDATE HERE?????????
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE






    real*8 :: KS_eigenvector_temp(n_basis), coeff_vector(n_basis_current), load_vector(n_basis_current)
    real*8 :: matrix_temp(n_basis, n_basis)
    real*8 :: matrix_temp2(n_basis, n_basis_current)
    real*8 :: coeff_matrix(n_basis_current, n_basis_current)
    integer :: i_state, i_basis, i_b, i_index, i_basis_1, i_basis_2
    integer :: info
    integer :: ipiv(n_basis_current)
    real*8, allocatable :: ortogonal_basis(:,:)

    real*8 :: r

    i_index = 0
    matrix_temp = 0.d0
    do i_basis_2 = 1,n_basis
       do i_basis_1 = 1,i_basis_2
          i_index = i_index + 1
          matrix_temp(i_basis_1,i_basis_2) =  overlap_matrix(i_index)
          matrix_temp(i_basis_2,i_basis_1) =  overlap_matrix(i_index)
       end do
    end do

    call DGEMM( 'N', 'N',  n_basis ,n_basis_current, n_basis , 1.d0, matrix_temp,n_basis, &
         ortogonal_basis,n_basis,0.d0, matrix_temp2,n_basis)
    call DGEMM( 'T', 'N', n_basis_current, n_basis_current, n_basis, 1.d0, ortogonal_basis, n_basis, &
         matrix_temp2, n_basis, 0.d0, coeff_matrix, n_basis_current)

    do i_state = 1, n_states_current, 1

       KS_eigenvector_temp = 0.0d0
       do i_basis = 1, n_basis, 1
          i_b = 0
          do i_basis_2 = 1, n_basis, 1
             i_b = i_b + 1
             KS_eigenvector_temp(i_basis) = KS_eigenvector_temp(i_basis) + &
                matrix_temp(i_basis, i_basis_2)*KS_eigenvector(i_b,i_state)
          end do
       end do

!!$       CALL DGEMV('N',n_basis,n_basis_current,1.d0, ortogonal_basis, n_basis,  &
!!$            KS_eigenvector(1:n_basis_current,i_state),1 ,0.d0,KS_eigenvector_temp(1:n_basis), 1)

       CALL DGEMV('T',n_basis,n_basis_current,1.d0, ortogonal_basis, n_basis,  &
            KS_eigenvector_temp(1:n_basis),1 ,0.d0, load_vector(1:n_basis_current) ,1)

       call DGESV( n_basis_current, 1, coeff_matrix, n_basis_current, ipiv, load_vector, n_basis_current, info )
       coeff_vector = load_vector

       CALL DGEMV('N',n_basis,n_basis_current,1.d0, ortogonal_basis, n_basis,  &
            coeff_vector(1:n_basis_current),1 ,0.d0, KS_eigenvector_temp(1:n_basis),1)

       i_b = 0
       do i_basis = 1, n_basis, 1
          i_b = i_b + 1
          KS_eigenvector(i_b,i_state) = KS_eigenvector_temp(i_basis)
       end do

       r = sqrt(S_inner_product( KS_eigenvector(:,i_state), KS_eigenvector(:,i_state), &
            overlap_matrix_current, n_basis_current))
       KS_eigenvector(:,i_state) = KS_eigenvector(:,i_state) / r


       KS_eigenvector_temp = 0.0d0
       do i_basis = 1, n_basis, 1
          i_b = 0
          do i_basis_2 = 1, n_basis, 1
             i_b = i_b + 1
             KS_eigenvector_temp(i_basis) = KS_eigenvector_temp(i_basis) + &
                matrix_temp(i_basis, i_basis_2)*stored_KS_eigenvector(i_b, i_state)
          end do
       end do

!!$       CALL DGEMV('N',n_basis,n_basis_current,1.d0, ortogonal_basis, n_basis,  &
!!$            stored_KS_eigenvector(1:n_basis_current,i_state),1 ,0.d0,KS_eigenvector_temp(1:n_basis), 1)

       CALL DGEMV('T',n_basis,n_basis_current,1.d0, ortogonal_basis, n_basis,  &
            KS_eigenvector_temp(1:n_basis),1 ,0.d0, load_vector(1:n_basis_current) ,1)

       call DGESV( n_basis_current, 1, coeff_matrix, n_basis_current, ipiv, load_vector, n_basis_current, info )
       coeff_vector = load_vector

       CALL DGEMV('N',n_basis,n_basis_current,1.d0, ortogonal_basis, n_basis,  &
            coeff_vector(1:n_basis_current),1 ,0.d0, KS_eigenvector_temp(1:n_basis),1)

       i_b = 0
       do i_basis = 1, n_basis, 1
          i_b = i_b + 1
          stored_KS_eigenvector(i_b,i_state) = KS_eigenvector_temp(i_basis)
       end do

       r = sqrt(S_inner_product( stored_KS_eigenvector(:,i_state), stored_KS_eigenvector(:,i_state), &
            overlap_matrix_current, n_basis_current))
       stored_KS_eigenvector(:,i_state) = stored_KS_eigenvector(:,i_state) / r

    end do

  end subroutine remove_core_states
!******
!------------------------------------------------------------------------------------
!****s* cg/cleanup_cg
!  NAME
!    cleanup_cg
!  SYNOPSIS

  subroutine cleanup_cg()

!  PURPOSE
!  Deallocates module variables
!
!  USES
    implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




    if (allocated(stored_KS_eigenvector)) then
       deallocate(stored_KS_eigenvector)
    end if
    if (allocated(stored_KS_eigenvalue)) then
       deallocate(stored_KS_eigenvalue)
    end if
    if (allocated(cholesky_factor)) then
       deallocate(cholesky_factor)
    end if

  end subroutine cleanup_cg
!******
!------------------------------------------------------------------------------------
!****s* cg/integrate_hartree_overlap
!  NAME
!    integrate_hartree_overlap
!  SYNOPSIS

  subroutine integrate_hartree_overlap()

!  PURPOSE
!  ??????????????????????????
!
!  USES

    use physics
    use synchronize_mpi
    use species_data
    implicit none

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE







    !  local variables

    integer basis_l_max (n_species)

    integer :: l_ylm_max
    integer, dimension(:,:), allocatable :: index_lm
    real*8, dimension(:,:), allocatable :: ylm_tab

    real*8 coord_current(3)
    real*8 dist_tab(n_centers_basis_integrals, n_max_angular)
    real*8 dist_tab_sq(n_centers_basis_integrals, n_max_angular)
    real*8 i_r(n_centers_basis_integrals)
    real*8 dir_tab(3, n_centers_basis_integrals, n_max_angular)
    real*8 trigonom_tab(4, n_centers_basis_integrals)

    real*8 radial_wave(n_max_compute_fns_ham)
    real*8 wave(n_max_compute_ham, n_max_angular)
    real*8 V_times_wave(n_max_compute_ham, n_max_angular)

    !    Auxiliary Hamiltonian matrix, to sum up contributions from only a single integration shell
    !     The hope is that such a separate treatment will allow to minimize numerical noise
    !     introduced through ZORA
    real*8, dimension(:,:), allocatable :: matrix_shell

    !     optimal accounting for matrix multiplications: only use points with nonzero components
    integer :: n_points

    !     and condensed version of partition_tabs on angular grids
    real*8 :: partition(n_max_angular)

    !     for pruning of atoms, radial functions, and basis functions, to only the relevant ones ...

    integer :: n_compute_a, n_compute_c
    integer :: i_basis(n_centers_basis_I)

    integer :: n_compute_fns
    integer :: i_basis_fns(n_basis_fns*n_centers_basis_integrals)
    integer :: i_basis_fns_inv(n_basis_fns,n_centers)
    integer :: i_atom_fns(n_basis_fns*n_centers_basis_integrals)

    integer :: n_compute_atoms
    integer :: atom_index(n_centers_basis_integrals)
    integer :: atom_index_inv(n_centers)

    integer :: spline_array_start(n_centers_basis_integrals)
    integer :: spline_array_end(n_centers_basis_integrals)

    !     for splitting of angular shells into "octants"

    integer division_low
    integer division_high

    !  counters

    integer i_basis_1
    integer i_basis_2
    integer i_grid
    integer i_index, i_l, i_m
    integer i_coord
    integer i_center, i_center_L
    integer i_division

    integer i_species, info

    integer i_point
    integer :: i_full_points
    integer :: i_full_points_2

    integer :: i_my_batch

    character*100 :: info_str

    !  begin work

    if (.not.allocated(hartree_overlap_matrix)) then
       allocate(hartree_overlap_matrix(n_hamiltonian_matrix_size),stat=info)
       call check_allocation(info, 'hartree_overlap_matrix        ')
    end if

    write(info_str,'(2X,A,A)')"Integrating hartree overlap matrix."
    call localorb_info(info_str,use_unit,'(A)')

    basis_l_max = l_shell_max

    !     begin with general allocations
    l_ylm_max = l_wave_max

    allocate( ylm_tab( (l_ylm_max+1)**2, n_centers_basis_integrals),stat=info )
    call check_allocation(info, 'ylm_tab                       ')

    allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info )
    call check_allocation(info, 'index_lm                      ')

    allocate ( matrix_shell(n_max_compute_ham,n_max_compute_ham) ,stat=info)
    call check_allocation(info, 'matrix_shell                  ')

    !     initialize

    hartree_overlap_matrix = 0.d0

    !     initialize index_lm

    i_index = 0
    do i_l = 0, l_wave_max, 1
       do i_m = -i_l, i_l
          i_index = i_index+1
          index_lm(i_m,i_l) = i_index
       enddo
    enddo

    i_full_points = 0
    i_full_points_2 = 0

    do i_my_batch = 1, n_my_batches, 1

          n_compute_c = 0
          n_compute_a = 0
          i_basis = 0

          i_point = 0

          ! loop over one batch
          do i_index = 1, batches(i_my_batch)%size, 1

             i_full_points_2 = i_full_points_2 + 1

             if (partition_tab(i_full_points_2).gt.0.d0) then

                i_point = i_point+1

                !     get current integration point coordinate
                coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)

                if(n_periodic > 0)then
                   call map_to_center_cell(coord_current(1:3) )
                end if

                !     compute atom-centered coordinates of current integration point,
                !     as viewed from all atoms
                call tab_atom_centered_coords_p0( coord_current, dist_tab_sq(1,i_point), &
                     dir_tab(1,1,i_point), n_centers_basis_integrals, centers_basis_integrals )

                !    determine which basis functions are relevant at current integration point,
                !     and tabulate their indices

                ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
                call prune_basis_p0( dist_tab_sq(1,i_point), n_compute_a, n_compute_c, &
                     i_basis, n_centers_basis_I, n_centers_basis_integrals, &
                     inv_centers_basis_integrals )

             end if

          enddo ! end loop over one batch

          n_points = i_point

          ! Perform actual integration if more than 0 basis functions
          ! are actually relevant on the present angular shell ...
          if (n_compute_a.gt.0) then

             i_point = 0

             ! loop over one batch of integration points
             do i_index = 1, batches(i_my_batch)%size, 1

                ! Increment the (global) counter for the grid, to access storage arrays
                i_full_points = i_full_points + 1

                if (partition_tab(i_full_points).gt.0.d0) then

                   i_point = i_point+1

                   ! for all integrations
                   partition(i_point) = partition_tab(i_full_points)

                   n_compute_atoms = 0
                   n_compute_fns = 0
                   i_basis_fns_inv = 0

                   ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                   ! Are stored in a compact spline array that can be accessed by spline_vector_waves,
                   ! without any copying and without doing any unnecessary operations.
                   ! The price is that the interface is no longer explicit in terms of physical
                   ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

                   call prune_radial_basis_p0 &
                        ( dist_tab_sq(1,i_point), &
                        dist_tab(1,i_point), &
                        dir_tab(1,1,i_point), &
                        n_compute_atoms, atom_index, atom_index_inv, &
                        n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                        i_atom_fns, spline_array_start, spline_array_end, &
                        n_centers_basis_integrals, centers_basis_integrals)

                   ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                   ! for all atoms which are actually relevant
                   call tab_local_geometry_p0 &
                        ( dist_tab_sq(1, i_point), n_compute_atoms, atom_index, &
                        dir_tab(1,1,i_point), dist_tab(1,i_point),  &
                        i_r )

                   !              compute trigonometric functions of spherical coordinate angles
                   !              of current integration point, viewed from all atoms
                   call tab_trigonom_p0 &
                        ( n_compute_atoms, dir_tab(1,1,i_point),  &
                        trigonom_tab )

                   ! tabulate distance and Ylm's w.r.t. other atoms
                   call tab_wave_ylm_p0 &
                        ( n_compute_atoms, atom_index,  &
                        trigonom_tab, basis_l_max,  &
                        l_ylm_max, &
                        ylm_tab )

                   ! Now evaluate radial functions
                   ! from the previously stored compressed spline arrays
                   call evaluate_radial_functions_p0  &
                        ( spline_array_start, spline_array_end,  &
                        n_compute_atoms, n_compute_fns,   &
                        dist_tab(1,i_point), i_r,  &
                        atom_index, i_basis_fns_inv,  &
                        basis_wave_ordered, radial_wave,  &
                        .false. , n_compute_c, n_max_compute_fns_ham )

                   ! tabulate total wave function value for each basis function
                   call evaluate_waves_p0  &
                        ( l_ylm_max,   &
                        ylm_tab, dist_tab(1,i_point),   &
                        index_lm, n_compute_c,   &
                        i_basis, radial_wave,   &
                        wave(1,i_point), n_compute_atoms,   &
                        atom_index_inv, n_compute_fns,  &
                        i_basis_fns_inv,  n_max_compute_fns_ham )

                   V_times_wave(1:n_compute_c,i_point) = &
                        (hartree_potential(i_full_points) - free_hartree_superpos(i_full_points))* &
                        wave(1:n_compute_c,i_point)

                   !            end if (partition_tab.gt.0)
                end if

                ! end loop over one batch
             enddo

             ! add full non-relativistic contributions and (for relativistic points)
             ! all contributions from the potential to the Hamiltonian matrix elements

             call evaluate_hamiltonian_shell_p1 &
                  ( n_points,  &
                  partition, &
                  n_compute_a, V_times_wave, &
                  n_compute_c, wave(1,1),  &
                  matrix_shell )

             call update_full_matrix_p0 &
                  ( n_compute_c, n_compute_a,  &
                  i_basis, &
                  matrix_shell,  &
                  hartree_overlap_matrix &
                )

          else
             i_full_points = i_full_points + batches(i_my_batch)%size
             !      end if (n_compute.gt.0) then
          end if

          ! end distribution over batches
       !end if

       !     end loop over bathces
    enddo

    !     synchronise the hamiltonian
    call sync_vector( hartree_overlap_matrix, n_hamiltonian_matrix_size )

    if (allocated(ylm_tab)) then
       deallocate(ylm_tab)
    end if
    if (allocated(index_lm)) then
       deallocate(index_lm)
    end if

    if (allocated(matrix_shell)) then
       deallocate( matrix_shell )
    end if

  end subroutine integrate_hartree_overlap
!******

end module cg
