!****h* FHI-aims/wf_extrapolation
!  NAME
!    wf_extrapolation - wave function extrapolation for MD
!  SYNOPSIS
module wf_extrapolation
  !  PURPOSE
  !
  !    Necessary data and subroutines for the wf extrapolation.
  !
  !    What could still be done:
  !      o wf_extrapolation_coefficients() and wf_extrapolation_nikl_coeffs()
  !        could support more-point schemes and least-squares fitting
  !        (solve_LSQ() instead of solve_LEQ()) to stabilize things.  For
  !        ammonia, I did not see much of an improvement in old tests,
  !        but who knows...
  !                
  !                     -- JW
  !  USES

  use dimensions
  use mpi_tasks
  use localorb_io
  use synchronize_mpi
  use runtime_choices
  use general_function
  use scalapack_wrapper, only: dlen_
  implicit none

  !  AUTHOR
  !    FHI-aims team.
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  INPUTS
  !    none
  !  OUTPUT
  !    none
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
  !******

  ! Maximum number of known iterations
  integer :: n_max_wf_steps
  ! Keep track of retrieval in order to allow update_density_densmat()
  ! to unconditionally call wf_extrapolate_density_matrix() without overwriting
  ! subsequent electronic SCF steps:
  integer :: last_saved_wf_step
  integer :: last_extrapolated_wf_step 

  integer :: n_wf_init            ! 1 for 'gen_func', 2 for 'niklasson06'
  integer, parameter, private :: wf_converged = 1
  integer, parameter, private :: wf_initial = 2

  integer, parameter :: n_max_wf_funcs = 10
  integer :: n_wf_funcs
  type(gen_func) :: wf_funcs(n_max_wf_funcs)  ! see general_function.f90

  ! n_dim1, n_dim2 are the first dimensions of the (potentially distributed) matrices,
  ! i.e, n_dim1=n_basis, n_dim2=n_basis if Scalapack is not used,
  ! n_dim1=mxld, n_dim2=mxcol if Scalapack is used (mxld, mxcol from module scalapack_wrapper)
  integer, private :: n_dim1, n_dim2

  ! Scalapack descriptor 
  ! We cannot use sc_desc from scalapack_wrapper since this name clashes with other routines
  integer, dimension(dlen_), private :: sc_desc

  real*8, allocatable :: wf_scalapack_overlap(:,:)
  complex*16, allocatable :: wf_scalapack_overlap_cmplx(:,:)

  real*8, allocatable, private     :: cckdm_save(:,:,:,:,:,:)
  complex*16, allocatable, private :: cckdm_save_cmplx(:,:,:,:,:,:)

  ! size: (n_max_wf_steps, n_wf_init)
  integer, allocatable, private :: pos2step(:,:)
  ! kdm_save(...,i_pos, i_init) is of MD step pos2step(i_pos, i_init).

contains

  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/initialize_wf
  !  NAME
  !    initialize_wf
  !  SYNOPSIS

  subroutine initialize_wf()

    !  PURPOSE
    !    Prepare module data of wf_extrapolation.
    !  USES

    use scalapack_wrapper, only: mxld, mxcol, mb, nb, rsrc, csrc, my_blacs_ctxt

    implicit none

    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character*150 :: info_str
    integer :: info
    character(*), parameter :: func = 'initialize_wf'

    call localorb_info('  Initialize wf_extra')

    last_saved_wf_step = 0
    last_extrapolated_wf_step = 0

    select case (wf_extra_type)
    case(WF_EXTRA_FUNC)
       n_max_wf_steps = n_wf_funcs
       n_wf_init = 1
    case(WF_EXTRA_NIKL)
       n_max_wf_steps = n_wf_funcs + 1
       n_wf_init = 2
    case(WF_EXTRA_NONE)
       n_max_wf_steps = 0
       n_wf_init = 0
    case default
       write(info_str, "('*** initialize_wf: Invalid extrapolation type.')")
       call aims_stop(info_str)
    end select

    ! Set n_dim1, n_dim2

    if (use_scalapack) then
       n_dim1 = mxld
       n_dim2 = mxcol
       call descinit( sc_desc, n_basis, n_basis, mb, nb, rsrc, csrc, &
            my_blacs_ctxt, MAX(1,mxld), info )
       if (real_eigenvectors) then
          allocate(wf_scalapack_overlap(n_dim1, n_dim2))
       else
          allocate(wf_scalapack_overlap_cmplx(n_dim1, n_dim2))
       end if
    else
       n_dim1 = n_basis
       n_dim2 = n_basis
    endif

    if (n_max_wf_steps > 0) then
       if (real_eigenvectors) then
          allocate(cckdm_save(n_dim1, n_dim2, n_spin, n_k_points_task, &
          &                   n_max_wf_steps, n_wf_init), stat=info)
       else
          allocate(cckdm_save_cmplx(n_dim1, n_dim2, n_spin, n_k_points_task, &
          &                         n_max_wf_steps, n_wf_init), stat=info)
       end if
       call check_allocation(info, 'kdm_save')

       allocate(pos2step(n_max_wf_steps, n_wf_init))
       pos2step = 0         ! invalid
    end if

  end subroutine initialize_wf
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extra/cleanup_wf
  !  NAME
  !    cleanup_wf
  !  SYNOPSIS

  subroutine cleanup_wf

    !  PURPOSE
    !    Deallocate module arrays.
    !  USES

    implicit none

    !  ARGUMENTS

    ! none

    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i

    if (allocated(cckdm_save)) deallocate(cckdm_save)
    if (allocated(cckdm_save_cmplx)) deallocate(cckdm_save_cmplx)
    if (allocated(pos2step)) deallocate(pos2step)
    if (allocated(wf_scalapack_overlap)) deallocate(wf_scalapack_overlap)
    if (allocated(wf_scalapack_overlap_cmplx)) deallocate(wf_scalapack_overlap_cmplx)

  end subroutine cleanup_wf
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/wf_save_eigenvectors
  ! NAME
  !   wf_save_eigenvectors
  ! SYNOPSIS

  subroutine wf_save_eigenvectors(occ_numbers, KS_eigenvector, KS_eigenvector_complex)

    ! PURPOSE
    !   Calculate and write out the contra-covariant density matrix and
    !   store it in a module array.
    ! USES

    use physics, only: n_electrons
    use timing
    implicit none

    ! ARGUMENTS

    real*8, intent(IN) :: occ_numbers(n_states, n_spin, n_k_points)
    real*8, intent(IN) :: KS_eigenvector(n_basis, n_states, n_spin, n_k_points_task)
    complex*16, intent(IN) :: KS_eigenvector_complex(n_basis, n_states, n_spin, n_k_points_task)

    ! INPUTS
    !   o occ_numbers -- occupation numbers for density matrix construction
    !   o KS_eigenvector{_complex} -- Eigenvectors to save
    !     These are only used if use_scalapack is not in effect,
    !     otherways eigenvec{_complex} from scalapack_wrapper is used
    ! OUTPUTS
    !   none
    ! AUTHOR
    !   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    ! HISTORY
    !   Release version, FHI-aims (2010).
    ! SOURCE
    real*8, allocatable :: contraco_kdm(:,:,:,:)
    complex*16, allocatable :: contraco_kdm_cmplx(:,:,:,:)
    real*8 :: trace
    real*8 :: time_save, clock_time_save
    character*100 :: info_str
    integer :: this_step

    call localorb_info('  Save eigenvectors for extrapolation')
    call get_timestamps(time_save, clock_time_save)

    ! --- update administrative info

    last_saved_wf_step = last_saved_wf_step + 1
    this_step = last_saved_wf_step

    ! --- allocate

    if (real_eigenvectors) then
       allocate(contraco_kdm(n_dim1, n_dim2, n_spin, n_k_points_task))
       ! complex matrix: dummy only
       allocate(contraco_kdm_cmplx( 1, 1, 1, 1))
    else
       allocate(contraco_kdm_cmplx(n_dim1, n_dim2, n_spin, n_k_points_task))
       ! real matrix: dummy only
       allocate(contraco_kdm(1, 1, 1, 1 ))
    end if

    ! --- calculate & save

    call calculate_contraco_kdm_all(contraco_kdm, contraco_kdm_cmplx, &
    &                               occ_numbers, KS_eigenvector, KS_eigenvector_complex, trace)

    call wf_save_cckdm('converged', this_step, contraco_kdm, contraco_kdm_cmplx)

    ! --- check trace

    if (abs(trace - n_electrons) > 1d-10) then
       write(info_str, "(2X,'Trace of saved ccdm:',ES24.16)") trace
       call localorb_info(info_str, use_unit, '(A)', OL_norm)
    end if

    ! --- deallocate

    if (allocated(contraco_kdm)) deallocate(contraco_kdm)
    if (allocated(contraco_kdm_cmplx)) deallocate(contraco_kdm_cmplx)

    call get_times(time_save, clock_time_save, tot_time_wf_extra_in, tot_clock_time_wf_extra_in)

  end subroutine wf_save_eigenvectors
  !******

  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/wf_extrapolate_eigenvectors
  !  NAME
  !    wf_extrapolate_eigenvectors
  !  SYNOPSIS
  subroutine wf_extrapolate_eigenvectors(occ_numbers, KS_eigenvector, KS_eigenvector_complex)
    !  PURPOSE
    !    Calls wf_extrapolate_cckdm and applies the result to KS_eigenvector.
    !  USES
    use scalapack_wrapper, only: eigenvec, eigenvec_complex, my_k_point, l_col
    use timing
    implicit none
    !  ARGUMENTS
    real*8, intent(IN) :: occ_numbers(n_states, n_spin, n_k_points)
    real*8, intent(INOUT) :: KS_eigenvector(n_basis, n_states, n_spin, n_k_points_task)
    complex*16, intent(INOUT) :: KS_eigenvector_complex(n_basis, n_states, n_spin, n_k_points_task)
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  INPUTS
    !    occ_numbers -- Occupation numbers (for decision between occupied and unoccupied)
    !    KS_eigenvector{_complex} -- Current eigencoefficients
    !      These are only used if use_scalapack is not in effect,
    !      otherways eigenvec{_complex} from scalapack_wrapper is used
    !  OUTPUT
    !    KS_eigenvector{_complex} -- Extrapolated eigencoefficients
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE

    real*8 :: time_extra, clock_time_extra
    real*8 :: half_occ, proj_fac
    real*8, allocatable :: contraco_kdm(:,:,:,:)
    real*8, allocatable :: build_ev(:,:)
    complex*16, allocatable :: contraco_kdm_cmplx(:,:,:,:)
    complex*16, allocatable :: build_ev_cmplx(:,:)
    complex*16 :: proj_fac_cmplx
    integer :: info, mpierr, k_proc
    integer :: i_state, i_spin, i_k_point, i_k
    logical :: is_extrapolated
    character(*), parameter :: func = 'wf_extrapolate_eigenvectors'

    ! --- prepare & allocate

    call get_timestamps(time_extra, clock_time_extra)

    half_occ = 1.d0 / n_spin
    proj_fac = n_spin / 2.d0   ! scale contraco_kdm to be a projector
    proj_fac_cmplx = proj_fac

    ! build_ev needs only to cover (n_basis, n_states) but for now we
    ! allocate it as big as the other matrices
    if (real_eigenvectors) then
       allocate(build_ev(n_dim1, n_dim2), stat=info)
       call check_allocation(info, 'wf_extrapolate_eigenvectors:build_ev')
       allocate(contraco_kdm(n_dim1, n_dim2, n_spin, n_k_points_task), stat=info)
       call check_allocation(info, 'wf_extrapolate_eigenvectors:contraco_kdm')
       ! complex matrix: dummy only
       allocate(contraco_kdm_cmplx(1, 1, 1, 1), stat=info)
       call check_allocation(info, 'wf_extrapolate_eigenvectors:contraco_kdm_cmplx')
    else
       allocate(build_ev_cmplx(n_dim1, n_dim2), stat=info)
       call check_allocation(info, 'wf_extrapolate_eigenvectors:build_ev_cmplx')
       allocate(contraco_kdm_cmplx(n_dim1, n_dim2, n_spin, n_k_points_task), stat=info)
       call check_allocation(info, 'wf_extrapolate_eigenvectors:contraco_kdm_cmplx')
       ! real matrix: dummy only
       allocate(contraco_kdm(1, 1, 1, 1), stat=info)
       call check_allocation(info, 'wf_extrapolate_eigenvectors:contraco_kdm')
    end if

    ! --- get extrapolated contra-covariant density matrix

    call wf_extrapolate_cckdm(contraco_kdm, contraco_kdm_cmplx, is_extrapolated)

    ! --- apply to eigenvectors

    ! Essentially, calculate build_ev := contraco_kdm * KS_eigenvector.
    ! Then, for occupied states, KS_eigenvector := build_ev,
    ! for unoccupied states,     KS_eigenvector := KS_eigenvector - build_ev.

    if (is_extrapolated) then
       i_k = 0
       do i_k_point = 1, n_k_points
          if (use_scalapack) then
             if (my_k_point /= i_k_point) cycle    ! all tasks of this point
          else
             if (myid /= modulo(i_k_point, n_tasks) .or. myid > n_k_points) cycle
          end if
          i_k = i_k + 1

          do i_spin = 1, n_spin
             if (use_scalapack) then
                if (real_eigenvectors) then
                   call pdgemm('N', 'N', n_basis, n_states, n_basis, &
                   &           proj_fac, contraco_kdm(:,:, i_spin, i_k), 1, 1, sc_desc, &
                   &           eigenvec(:,:, i_spin), 1, 1, sc_desc, &
                   &           0.d0, build_ev, 1, 1, sc_desc)
                   do i_state = 1, n_states
                      if(l_col(i_state)==0) cycle
                      if (occ_numbers(i_state, i_spin, i_k_point) > half_occ) then
                         eigenvec(:, l_col(i_state), i_spin) = build_ev(:, l_col(i_state))
                      else
                         eigenvec(:, l_col(i_state), i_spin) &
                         & = eigenvec(:, l_col(i_state), i_spin) - build_ev(:, l_col(i_state))
                      end if
                   end do
                else
                   call pzgemm('N', 'N', n_basis, n_states, n_basis, &
                   &           proj_fac_cmplx, contraco_kdm_cmplx(:,:, i_spin, i_k), 1, 1, sc_desc, &
                   &           eigenvec_complex(:,:, i_spin), 1, 1, sc_desc, &
                   &           (0.d0, 0.d0), build_ev_cmplx, 1, 1, sc_desc)
                   do i_state = 1, n_states
                      if(l_col(i_state)==0) cycle
                      if (occ_numbers(i_state, i_spin, i_k_point) > half_occ) then
                         eigenvec_complex(:, l_col(i_state), i_spin) = build_ev_cmplx(:, l_col(i_state))
                      else
                         eigenvec_complex(:, l_col(i_state), i_spin) &
                         & = eigenvec_complex(:, l_col(i_state), i_spin) &
                         &   - build_ev_cmplx(:, l_col(i_state))
                      end if
                   end do
                end if
             else
                if (real_eigenvectors) then
                   call dgemm('N', 'N', n_basis, n_states, n_basis, &
                   &          proj_fac, contraco_kdm(:,:, i_spin, i_k), n_basis, &
                   &          KS_eigenvector(:,:, i_spin, i_k), n_basis, &
                   &          0.d0, build_ev, n_basis)
                   do i_state = 1, n_states
                      if (occ_numbers(i_state, i_spin, i_k_point) > half_occ) then
                         KS_eigenvector(:, i_state, i_spin, i_k) = build_ev(:, i_state)
                      else
                         KS_eigenvector(:, i_state, i_spin, i_k) &
                         & = KS_eigenvector(:, i_state, i_spin, i_k) - build_ev(:, i_state)
                      end if
                   end do
                else
                   call zgemm('N', 'N', n_basis, n_states, n_basis, &
                   &          proj_fac_cmplx, contraco_kdm_cmplx(:,:, i_spin, i_k), n_basis, &
                   &          KS_eigenvector_complex(:,:, i_spin, i_k), n_basis, &
                   &          (0.d0, 0.d0), build_ev_cmplx, n_basis)
                   do i_state = 1, n_states
                      if (occ_numbers(i_state, i_spin, i_k_point) > half_occ) then
                         KS_eigenvector_complex(:, i_state, i_spin, i_k) = build_ev_cmplx(:, i_state)
                      else
                         KS_eigenvector_complex(:, i_state, i_spin, i_k) &
                         & = KS_eigenvector_complex(:, i_state, i_spin, i_k) &
                         &   - build_ev_cmplx(:, i_state)
                      end if
                   end do
                end if
             endif
          end do
       end do
    end if

    if ( allocated(build_ev) )           deallocate(build_ev)
    if ( allocated(build_ev_cmplx) )     deallocate(build_ev_cmplx)
    if ( allocated(contraco_kdm) )       deallocate(contraco_kdm)
    if ( allocated(contraco_kdm_cmplx) ) deallocate(contraco_kdm_cmplx)

    ! one corner case where the eigenvectors need to be collected explicitly:
    ! normal update + cluster case + LAPACK
    if((.not. use_density_matrix) .and. &
    &  n_periodic==0 .and. &
    &  (.not. use_scalapack) .and. &
    &  n_tasks > 1) then
       k_proc = 1  ! number of processor responsible for the only k-point
       call MPI_Bcast(KS_eigenvector, n_basis*n_states*n_spin, &
       & MPI_DOUBLE_PRECISION, k_proc, mpi_comm_global, mpierr)
       if (mpierr /= MPI_SUCCESS) call aims_stop('MPI_Bcast error', func)
    end if

    call get_times(time_extra, clock_time_extra, &
    &              tot_time_wf_extra_out, tot_clock_time_wf_extra_out)

  end subroutine wf_extrapolate_eigenvectors
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/wf_extrapolate_eigenvectors_diag
  !  NAME
  !    wf_extrapolate_eigenvectors_diag
  !  SYNOPSIS
  subroutine wf_extrapolate_eigenvectors_diag(occ_numbers, KS_eigenvector, KS_eigenvector_complex)
    !  PURPOSE
    !
    !    Gets the extrapolated density matrix using wf_extrapolate_density_matrix()
    !    and diagonalizes the result.
    !
    !    First tests indicate that this is not a good way to proceed.  The
    !    resulting occupation numbers are /very/ smeared, though surprisingly
    !    seem to stay within 0. and max_occ (for alanine dipeptide, the larges
    !    occupation number is about 1.3 instead of 2.).  I'm not completely
    !    sure that these tests were done with the correct overlap matrix,
    !    though.
    !
    !  USES
    use scalapack_wrapper, only: spread_eigenvectors_scalapack, my_k_point
    use physics, only: overlap_matrix, n_electrons
    use timing
    implicit none
    !  ARGUMENTS
    real*8, intent(OUT) :: occ_numbers(n_states, n_spin, n_k_points)
    real*8, intent(OUT) :: KS_eigenvector(n_basis, n_states, n_spin, n_k_points_task)
    complex*16, intent(OUT) :: KS_eigenvector_complex(n_basis, n_states, n_spin, n_k_points_task)
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  INPUTS
    !    none
    !  OUTPUT
    !    occ_numbers -- Extrapolated occupation numbers
    !    KS_eigenvector{_complex} -- Extrapolated eigencoefficients
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE

    real*8 :: time_extra, clock_time_extra
    real*8 :: dummy(1)
    real*8, allocatable :: density_matrix_sparse(:,:)
    real*8, allocatable :: ev_tmp(:,:,:), occ_tmp(:,:)
    complex*16, allocatable :: ev_tmp_cmplx(:,:,:)
    integer :: info, saved_priority, n_occ
    integer :: i_state, j_state, i_spin, i_k_point, i_k
    logical :: is_extrapolated
    character*100 :: info_str
    character(*), parameter :: func = 'wf_extrapolate_eigenvectors_diag'

    ! --- prepare & allocate

    call localorb_info('*** Eigenstates by diagonalizing density matrix.')
    call localorb_info('*** This routine is for testing purposes only and does not perform well.')
    if(use_scalapack) then
       call aims_stop_coll('*** This routine must not be used if use_scalapack is in effect!', func)
    endif

    n_occ = int(n_electrons / 2.d0)

    call get_timestamps(time_extra, clock_time_extra)
    allocate(density_matrix_sparse(n_hamiltonian_matrix_size, n_spin))
    allocate(occ_tmp(n_spin, n_k_points))
    if (real_eigenvectors) then
       allocate(ev_tmp(n_basis, n_spin, n_k_points))
    else
       allocate(ev_tmp_cmplx(n_basis, n_spin, n_k_points))
    end if

    ! --- get density matrix

    call wf_extrapolate_density_matrix(dummy, density_matrix_sparse, .true., is_extrapolated)

    if (is_extrapolated) then
       call solve_KS_eigen(overlap_matrix, density_matrix_sparse, &
       &                   occ_numbers, KS_eigenvector, KS_eigenvector_complex)
       ! occ_numbers are now sorted 0...max_occ; reverse:
       do i_state = 1, n_states/2
          j_state = n_states - i_state + 1
          occ_tmp = occ_numbers(i_state, :,:)
          occ_numbers(i_state, :,:) = occ_numbers(j_state, :,:)
          occ_numbers(j_state, :,:) = occ_tmp
          if (real_eigenvectors) then
             ev_tmp = KS_eigenvector(:, i_state, :,:)
             KS_eigenvector(:, i_state, :,:) = KS_eigenvector(:, j_state, :,:)
             KS_eigenvector(:, j_state, :,:) = ev_tmp
          else
             ev_tmp_cmplx = KS_eigenvector_complex(:, i_state, :,:)
             KS_eigenvector_complex(:, i_state, :,:) = KS_eigenvector_complex(:, j_state, :,:)
             KS_eigenvector_complex(:, j_state, :,:) = ev_tmp_cmplx
          end if
       end do

       if (use_scalapack) then
          ! Put vectors into scalapack-internal, distributed eigenvec arrays.
          call spread_eigenvectors_scalapack(KS_eigenvector, KS_eigenvector_complex)
       end if
    end if
    deallocate(density_matrix_sparse)
    write (info_str,'(2X,A,A)') "=== Writing extrapolated occ_numbers ==="
    call localorb_info(info_str)
    saved_priority = output_priority
    output_priority = 0
    call output_real_eigenfunctions(occ_numbers, KS_eigenvector, occ_numbers)
    output_priority = saved_priority

    occ_numbers(1:n_occ,:,:) = 2.d0 / n_spin
    occ_numbers(n_occ+1:,:,:) = 0.d0

    call get_times(time_extra, clock_time_extra, &
    &              tot_time_wf_extra_out, tot_clock_time_wf_extra_out)

  end subroutine wf_extrapolate_eigenvectors_diag
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/wf_extrapolate_density_matrix
  !  NAME
  !    wf_extrapolate_density_matrix
  !  SYNOPSIS
  subroutine wf_extrapolate_density_matrix(density_matrix, density_matrix_sparse, force_packed, &
  &                                        is_extrapolated)
    !  PURPOSE
    !    Calls wf_extrapolate_cckdm and applies the result to KS_eigenvector.
    !  USES
    use density_matrix_evaluation
    use scalapack_wrapper, only: my_k_point, my_scalapack_id
    use physics, only: overlap_matrix
    use pbc_lists, only: k_weights
    use timing
    implicit none
    !  ARGUMENTS
    real*8, intent(OUT) :: density_matrix(n_centers_basis_T, n_centers_basis_T, n_spin)
    real*8, intent(OUT) :: density_matrix_sparse(n_hamiltonian_matrix_size, n_spin)
    logical, intent(IN) :: force_packed
    logical, intent(OUT) :: is_extrapolated
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  INPUTS
    !    occ_numbers -- Occupation numbers (for constructing the density matrix)
    !    KS_eigenvector{_complex} -- Current eigencoefficients
    !  OUTPUT
    !    KS_eigenvector{_complex} -- Extrapolated eigencoefficients
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE

    logical :: use_packed
    real*8 :: time_extra, clock_time_extra
    real*8, allocatable :: overlap_lpck(:), work_ovl(:,:)
    complex*16, allocatable :: overlap_lpck_cmplx(:)
    real*8, allocatable :: contraco_kdm(:,:,:,:), kdm(:,:)
    complex*16, allocatable :: contraco_kdm_cmplx(:,:,:,:), kdm_cmplx(:,:)
    integer :: info
    integer :: i_k_point, i_k, i_spin, i_state
    integer, allocatable :: ipiv(:)
    character*120 :: info_str

    if(use_scalapack) then
       call localorb_info('*** wf_extrapolate_density_matrix must not be used if use_scalapack is in effect!')
       call aims_stop
    endif

    ! --- prepare & allocate

    if (wf_extrapolation_has_been_done()) then
       ! Has already been done, so this should be a subsequent SCF cycle.
       ! Don't extrapolate in such a case
       write(info_str, "('  * Extrapolation demanded, but already been done')")
       call localorb_info(info_str, use_unit, '(A)', OL_norm)
       is_extrapolated = .false.
       return
    end if

    write(info_str, "('  Extrapolating density matrix')")
    call localorb_info(info_str, use_unit, '(A)', OL_norm)

    call get_timestamps(time_extra, clock_time_extra)

    if (.not. allocated(overlap_matrix) .or. size(overlap_matrix) == 1) then
       call localorb_info('*** Internal error: need overlap_matrix for contraco-dm')
       call aims_stop
    end if

    if (real_eigenvectors) then
       allocate(kdm(n_basis, n_basis), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:kdm')
       allocate(contraco_kdm(n_basis, n_basis, n_spin, n_k_points_task), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:contraco_kdm')
       allocate(overlap_lpck(n_basis*(n_basis+1)/2), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:overlap_lpck')
       ! dummy for complex matrices only
       allocate(kdm_cmplx(1, 1), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:kdm_cmplx')
       allocate(contraco_kdm_cmplx(1, 1, 1, 1), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:contraco_kdm_cmplx')
       allocate(overlap_lpck_cmplx(1))
       call check_allocation(info, 'wf_extrapolate_density_matrix:overlap_lpck_cmplx')
    else
       allocate(ipiv(n_basis), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:ipiv')
       allocate(kdm_cmplx(n_basis, n_basis), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:kdm_cmplx')
       allocate(contraco_kdm_cmplx(n_basis, n_basis, n_spin, n_k_points_task), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:contraco_kdm_cmplx')
       allocate(overlap_lpck_cmplx(n_basis*(n_basis+1)/2))
       call check_allocation(info, 'wf_extrapolate_density_matrix:overlap_lpck_cmplx')
       ! dummy for real matrices only
       allocate(kdm(1, 1), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:kdm')
       allocate(contraco_kdm(1, 1, 1, 1), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:contraco_kdm')
       allocate(overlap_lpck(1), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:overlap_lpck')
    end if
    if (packed_matrix_format == PM_none .and. .not. real_eigenvectors) then
       ! This could be avoided by editing construct_overlap().
       allocate(work_ovl(n_centers_basis_I, n_centers_basis_I), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:work_ovlp')
    else 
       ! dummy allocation if not needed.
       allocate(work_ovl(1, 1), stat=info)
       call check_allocation(info, 'wf_extrapolate_density_matrix:work_ovlp')
    end if

    use_packed = (packed_matrix_format /= PM_none .or. force_packed)
    if (use_packed) then
       density_matrix_sparse = 0.d0
    else
       density_matrix = 0.d0
    end if


    ! --- get extrapolated contra-covariant density matrix

    call wf_extrapolate_cckdm(contraco_kdm, contraco_kdm_cmplx, is_extrapolated)

    ! --- construct overlap matrix, invert, apply to kdm, and accumulate

    ! Essentially, dm := contraco_dm * overlap^-1, density_matrix := dm + dm^T.

    if (is_extrapolated) then
       i_k = 0
       do i_k_point = 1, n_k_points
          if (use_scalapack) then
             if (my_k_point /= i_k_point .or. my_scalapack_id > 0) cycle
          else
             if (myid /= modulo(i_k_point, n_tasks) .or. myid > n_k_points) cycle
          end if
          i_k = i_k + 1

          ! Get Cholesky decomposition of the overlap:
          call construct_overlap(overlap_matrix, overlap_lpck, overlap_lpck_cmplx, i_k_point, work_ovl)
          if (real_eigenvectors) then
             call dpptrf('U', n_basis, overlap_lpck, info)
          else
             ipiv = 0
             call zhptrf('U', n_basis, overlap_lpck_cmplx, ipiv, info)
          end if
          if (info /= 0) then
             write(info_str, "('*** wf_extrapolate_density_matrx: d/zpptrf: info =',I5)") info
             call localorb_info(info_str)
             call aims_stop
          end if

          do i_spin = 1, n_spin
             ! Solve   kdm * overlap = contraco_kdm for kdm -> transpose
             if (real_eigenvectors) then
                kdm = transpose(contraco_kdm(:,:, i_spin, i_k))
                call dpptrs('U', n_basis, n_basis, overlap_lpck, kdm, n_basis, info)
             else
                kdm_cmplx = conjg(transpose(contraco_kdm_cmplx(:,:, i_spin, i_k)))
                call zhptrs('U', n_basis, n_basis, overlap_lpck_cmplx, ipiv, &
                &           kdm_cmplx, n_basis, info)
             end if
             if (info /= 0) then
                write(info_str, "('*** wf_extrapolate_density_matrx: d/zpptrs: info =',I5)") info
                call localorb_info(info_str)
                call aims_stop
             end if

             ! Symmetrize & weight
             if (real_eigenvectors) then
                kdm = (k_weights(i_k_point) * 0.5d0) * (kdm + transpose(kdm))
             else
                kdm_cmplx = (k_weights(i_k_point) * 0.5d0) * &
                &           (kdm_cmplx + conjg(transpose(kdm_cmplx)))
             end if

             ! Accumulate
             call accumulate_k_densmat(density_matrix_sparse(:, i_spin), &
             &                         density_matrix(:,:, i_spin), &
             &                         force_packed, kdm, kdm_cmplx, i_k_point)

          end do
       end do
    end if
    if (use_packed) then
       call sync_vector(density_matrix_sparse, n_hamiltonian_matrix_size)
    else
       call sync_matrix(density_matrix, n_centers_basis_I, n_centers_basis_I)
    end if

    if (allocated(work_ovl)) deallocate(work_ovl)
    ! deallocate including dummy matrices
    deallocate(kdm, kdm_cmplx, contraco_kdm, contraco_kdm_cmplx, overlap_lpck, overlap_lpck_cmplx )
    if (.not.real_eigenvectors) then
       deallocate(ipiv)
    end if

    call get_times(time_extra, clock_time_extra, &
    &              tot_time_wf_extra_out, tot_clock_time_wf_extra_out)

  end subroutine wf_extrapolate_density_matrix
  !******

  !------------------------------------------------------------------------------
  !****s* wf_extra/wf_extrapolation_has_been_done
  !  NAME
  !    wf_extrapolation_has_been_done
  !  SYNOPSIS

  logical function wf_extrapolation_has_been_done()

    !  PURPOSE
    !    Check if extrapolation for this time step has already been done.
    !  USES

    implicit none

    !  ARGUMENTS

    ! none

    !  INPUTS
    !    none
    !  OUTPUTS
    !    o wf_extrapolation_has_been_done -- guess what.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    wf_extrapolation_has_been_done = (last_extrapolated_wf_step > last_saved_wf_step)

  end function wf_extrapolation_has_been_done
  !******


  !------------------------------------------------------------------------------
  !-------------------------------- PRIVATES ------------------------------------
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !****s* wf_extra/wf_i_init
  !  NAME
  !    wf_i_init
  !  SYNOPSIS

  integer function wf_i_init(type_of_guess)

    !  PURPOSE
    !    Return correct offset for last index in kdm_save.
    !  USES

    implicit none

    !  ARGUMENTS

    character*(*), intent(IN) :: type_of_guess

    !  INPUTS
    !    o type_of_guess -- 'initial' or 'converged'
    !  OUTPUTS
    !    o wf_i_init -- 1 for 'converged' and 2 for 'initial'
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'wf_i_init'

    select case (type_of_guess)
    case('converged')
       wf_i_init = wf_converged
    case('initial')
       wf_i_init = wf_initial
    case default
       call aims_stop('invalid type_of_guess', func)
       stop  ! silence compiler
    end select
    if (wf_i_init > n_wf_init) then
       call aims_stop('type_of_guess not allowed with current settings', func)
    end if

  end function wf_i_init
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extra/wf_have_step
  !  NAME
  !    wf_have_step
  !  SYNOPSIS

  logical function wf_have_step(type_of_guess, i_step)

    !  PURPOSE
    !    Check if step i_step can be retrieved
    !  USES

    implicit none

    !  ARGUMENTS

    character*(*), intent(IN) :: type_of_guess
    integer, intent(IN)       :: i_step

    !  INPUTS
    !    o type_of_guess -- 'initial' or 'converged'
    !    o i_step -- MD step number to check
    !  OUTPUTS
    !    o wf_have_step -- .true. if there.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_pos, i_init

    if (i_step > 0) then
       i_init = wf_i_init(type_of_guess)
       i_pos = mod(i_step, n_max_wf_steps) + 1
       wf_have_step = (pos2step(i_pos, i_init) == i_step)
    else
       wf_have_step = .false.
    end if

  end function wf_have_step
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extra/wf_n_steps
  !  NAME
  !    wf_n_steps
  !  SYNOPSIS

  integer function wf_n_steps(type_of_guess, current_step)

    !  PURPOSE
    !    Check how many steps prior to current_step can be retrieved
    !  USES

    implicit none

    !  ARGUMENTS

    character*(*), intent(IN) :: type_of_guess
    integer, intent(IN) :: current_step

    !  INPUTS
    !    o type_of_guess -- 'initial' or 'converged'
    !    o current_step -- current MD step number
    !  OUTPUTS
    !    o wf_n_steps -- number of saved steps
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_step

    do i_step = 1, n_max_wf_steps
       if (.not. wf_have_step(type_of_guess, current_step-i_step)) then
          wf_n_steps = i_step - 1
          return
       end if
    end do
    wf_n_steps = n_max_wf_steps

  end function wf_n_steps
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/wf_save_cckdm
  ! NAME
  !   wf_save_cckdm
  ! SYNOPSIS

  subroutine wf_save_cckdm(type_of_guess, i_step, contraco_kdm, contraco_kdm_cmplx)

    ! PURPOSE
    !   Save the given contra-covariant k-dependent density matrix.
    ! USES

    implicit none

    ! ARGUMENTS

    character(len=*), intent(IN) :: type_of_guess
    integer, intent(IN) :: i_step
    real*8, intent(IN) :: contraco_kdm(n_dim1, n_dim2, n_spin, n_k_points_task)
    complex*16, intent(IN) :: contraco_kdm_cmplx(n_dim1, n_dim2, n_spin, n_k_points_task)

    ! INPUTS
    !   o type_of_guess -- 'initial' or 'converged'
    !   o i_step -- current step number
    !   o contraco_kdm{_cmplx} -- contracovariant k-dependent density matrix to save
    ! OUTPUTS
    !   none
    ! AUTHOR
    !   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    ! HISTORY
    !   Release version, FHI-aims (2010).
    ! SOURCE
    integer :: i_pos, i_init

    ! --- prepare

    i_pos = mod(i_step, n_max_wf_steps) + 1
    i_init = wf_i_init(type_of_guess)

    ! --- actually write

    pos2step(i_pos, i_init) = i_step
    if (real_eigenvectors) then
       cckdm_save(:,:,:,:, i_pos, i_init) = contraco_kdm
    else
       cckdm_save_cmplx(:,:,:,:, i_pos, i_init) = contraco_kdm_cmplx
    end if

  end subroutine wf_save_cckdm
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/wf_retrieve_cckdm
  ! NAME
  !   wf_retrieve_cckdm
  ! SYNOPSIS

  subroutine wf_retrieve_cckdm(type_of_guess, i_step, fac, contraco_kdm, contraco_kdm_cmplx)

    ! PURPOSE
    !   Retrieve the given contra-covariant k-dependent density matrix and add
    !   a multiple of it to contraco_kdm.
    ! USES

    implicit none

    ! ARGUMENTS

    character(len=*), intent(IN) :: type_of_guess
    integer, intent(IN) :: i_step
    real*8, intent(IN) :: fac
    real*8, intent(INOUT) :: contraco_kdm(n_dim1, n_dim2, n_spin, n_k_points_task)
    complex*16, intent(INOUT) :: contraco_kdm_cmplx(n_dim1, n_dim2, n_spin, n_k_points_task)

    ! INPUTS
    !   o type_of_guess -- 'initial' or 'converged'
    !   o i_step -- current step number
    !   o fac -- factor to apply to saved contracovariant density matrix
    !   o contraco_kdm{_cmplx} -- contracovariant density matrix under construction
    ! OUTPUTS
    !   o contraco_kdm{_cmplx} -- contracovariant density matrix under construction
    ! AUTHOR
    !   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    ! HISTORY
    !   Release version, FHI-aims (2010).
    ! SOURCE
    integer :: i_pos, i_init

    ! --- prepare

    i_pos = mod(i_step, n_max_wf_steps) + 1
    i_init = wf_i_init(type_of_guess)
    if (pos2step(i_pos, i_init) /= i_step) then
       call localorb_info('*** wf_retrieve_cckdm: kdm not found')
       call aims_stop
    end if

    ! --- actually get

    if (real_eigenvectors) then
       contraco_kdm = contraco_kdm + fac * cckdm_save(:,:,:,:, i_pos, i_init)
    else
       contraco_kdm_cmplx = contraco_kdm_cmplx + &
       &                    fac * cckdm_save_cmplx(:,:,:,:, i_pos, i_init)
    end if

  end subroutine wf_retrieve_cckdm
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/calculate_contraco_kdm_one_k
  ! NAME
  !   calculate_contraco_kdm_one_k
  ! SYNOPSIS
  subroutine calculate_contraco_kdm_one_k(contraco_kdm, contraco_kdm_cmplx, &
  &                                       occ_numbers, KS_eigenvector, KS_eigenvector_complex, &
  &                                       i_spin, i_k_point, i_k, acc_trace)

    ! PURPOSE
    !   Calculate the contra-covariant density matrix
    !   for one spin channel and one k-point.
    ! USES
    use physics, only: overlap_matrix
    use pbc_lists, only: k_weights
    use density_matrix_evaluation, only: evaluate_k_densmat
    use scalapack_wrapper, only: scalapack_ovlp=>ovlp, scalapack_ovlp_complex=>ovlp_complex, &
           ham, ham_complex, use_ovlp_trafo, n_nonsing_ovlp, &
           construct_dm_scalapack, &
           set_full_matrix_real, set_full_matrix_complex
    implicit none

    ! ARGUMENTS

    real*8, intent(OUT) :: contraco_kdm(n_dim1, n_dim2)
    complex*16, intent(OUT) :: contraco_kdm_cmplx(n_dim1, n_dim2)
    real*8, intent(IN) :: occ_numbers(n_states, n_spin, n_k_points)
    real*8, intent(IN) :: KS_eigenvector(n_basis, n_states, n_spin, n_k_points_task)
    complex*16, intent(IN) :: KS_eigenvector_complex(n_basis, n_states, n_spin, n_k_points_task)
    integer, intent(IN) :: i_spin, i_k_point, i_k
    real*8, intent(INOUT) :: acc_trace

    ! INPUTS
    !   o KS_eigenvector{_complex} -- eigenvectors from which to construct density matrix
    !     These are only used if use_scalapack is not in effect,
    !     otherways eigenvec{_complex} from scalapack_wrapper is used
    !   o i_spin -- spin component
    !   o i_k_point -- which k-point
    !   o i_k -- node-local k-point index
    !   o acc_trace -- real and imaginary part of accumulative trace
    ! OUTPUTS
    !   o contraco_kdm, contraco_kdm_cmplx -- resulting co-contra-variant density matrix
    !   o acc_trace -- real and imaginary part of accumulative trace
    ! AUTHOR
    !   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    ! HISTORY
    !   Release version, FHI-aims (2010).
    ! SOURCE

    real*8, allocatable, dimension(:,:) :: densmat, overlap, work_ovl
    complex*16, allocatable, dimension(:,:) :: densmat_cmplx, overlap_cmplx
    real*8, allocatable, dimension(:) :: overlap_lpck
    complex*16, allocatable, dimension(:) :: overlap_lpck_cmplx

    real*8 :: dummy(1), trace, trace_imag
    complex*16 :: dummy_cmplx(1)
    integer :: i_bas1, i_bas2, i_index, i_col
    character*150 :: info_str

    logical :: use_ovlp_trafo_saved
    integer :: n_nonsing_ovlp_saved


    ! --- set up overlap

    if (.not. use_scalapack) then

       if (.not. allocated(overlap_matrix) .or. size(overlap_matrix) == 1) then
          call localorb_info('*** Internal error: need overlap_matrix for contraco-dm')
          call aims_stop
       end if

       if (real_eigenvectors) then
          allocate(overlap(n_dim1, n_dim2))
       else
          allocate(overlap_cmplx(n_dim1, n_dim2))
       end if

       if (packed_matrix_format == PM_none .and. .not. real_eigenvectors) then
          ! This could be avoided by editing construct_overlap().
          allocate(work_ovl(n_centers_basis_I, n_centers_basis_I))
       else
          ! dummy
          allocate(work_ovl(1,1))
       end if

       if (real_eigenvectors) then
          allocate(overlap_lpck(n_basis*(n_basis+1)/2))
          ! dummy
          allocate(overlap_lpck_cmplx(1))
       else
          allocate(overlap_lpck_cmplx(n_basis*(n_basis+1)/2))
          ! dummy
          allocate(overlap_lpck(1))
       end if
       call construct_overlap(overlap_matrix, overlap_lpck, overlap_lpck_cmplx, i_k_point, work_ovl)
       if (allocated(work_ovl)) deallocate(work_ovl)

       i_index = 0
       do i_bas2 = 1, n_basis
          do i_bas1 = 1, i_bas2
             i_index = i_index + 1
             if (real_eigenvectors) then
                overlap(i_bas1, i_bas2) = overlap_lpck(i_index)
                overlap(i_bas2, i_bas1) = overlap_lpck(i_index)
             else
                overlap_cmplx(i_bas1, i_bas2) = overlap_lpck_cmplx(i_index)
                overlap_cmplx(i_bas2, i_bas1) = conjg(overlap_lpck_cmplx(i_index))
             end if
          end do
       end do
       if (allocated(overlap_lpck)) deallocate(overlap_lpck)
       if (allocated(overlap_lpck_cmplx)) deallocate(overlap_lpck_cmplx)

    endif

    ! --- set up density matrix

    ! In priniciple, it would be more efficient to first calculate
    ! c^HS and cp(c^HS) afterwards, instead of first P=cpc^H and then
    ! PS.

    call check_occs('calculate_contraco_kdm_one_k', occ_numbers, .false.)
    if (real_eigenvectors) then
       allocate(densmat(n_dim1, n_dim2))
       ! dummy for complex matrix
       allocate(densmat_cmplx(1, 1))
    else
       allocate(densmat_cmplx(n_dim1, n_dim2))
       ! dummy for real matrix
       allocate(densmat(1, 1))
    end if

    if (use_scalapack) then

       ! Please note:
       ! - the result of construct_dm_scalapack is stored in ham/ham_complex
       ! - only the upper triangle is stored, thus we have to set the full matrix here
       call construct_dm_scalapack(occ_numbers, i_spin)

       if (real_eigenvectors) then
          densmat(:,:) = ham(:,:,i_spin)
          call set_full_matrix_real(densmat)
       else
          densmat_cmplx(:,:) = ham_complex(:,:,i_spin)
          call set_full_matrix_complex(densmat_cmplx)
       endif
    else
       call evaluate_k_densmat(densmat, densmat_cmplx, occ_numbers, &
       &                       KS_eigenvector, KS_eigenvector_complex, &
       &                       i_spin, i_k_point, i_k)
    end if

    ! --- matrix multiply

    if(use_scalapack) then
       if (real_eigenvectors) then
          call pdgemm('N', 'N', n_basis, n_basis, n_basis, &
          &           1.d0, densmat, 1, 1, sc_desc, &
          &                 wf_scalapack_overlap, 1, 1, sc_desc, &
          &           0.d0, contraco_kdm, 1, 1, sc_desc)
       else
          call pzgemm('N', 'N', n_basis, n_basis, n_basis, &
          &           (1.d0, 0.d0), densmat_cmplx, 1, 1, sc_desc, &
          &                         wf_scalapack_overlap_cmplx, 1, 1, sc_desc, &
          &           (0.d0, 0.d0), contraco_kdm_cmplx, 1, 1, sc_desc)
       end if
    else
       if (real_eigenvectors) then
         call dgemm('N', 'N', n_basis, n_basis, n_basis, &
         &          1.d0, densmat, n_basis, overlap, n_basis, &
         &          0.d0, contraco_kdm, n_basis)
       else
         call zgemm('N', 'N', n_basis, n_basis, n_basis, &
         &          (1.d0, 0.d0), densmat_cmplx, n_basis, overlap_cmplx, n_basis, &
         &          (0.d0, 0.d0), contraco_kdm_cmplx, n_basis)
       end if
    endif
    if (allocated(densmat)) deallocate(densmat)
    if (allocated(densmat_cmplx)) deallocate(densmat_cmplx)

    if (allocated(overlap)) deallocate(overlap)
    if (allocated(overlap_cmplx)) deallocate(overlap_cmplx)

    ! --- check trace

    trace = 0.d0
    trace_imag = 0.d0
    call wf_trace(contraco_kdm, contraco_kdm_cmplx, trace, trace_imag)

    if (abs(trace_imag) > 1d-10 .or. &
    &   abs(trace - sum(occ_numbers(:, i_spin, i_k_point))) > 1d-3) then
       write(info_str, "(2X,A,I5,A,I5,A,I5,2A,ES24.16,A,F16.8)") &
       & '*** Task', myid, ' with i_k_point', i_k_point, ' and i_k', i_k,':', &
       & ' Trace of k-ccdm:', trace, ' instead of', sum(occ_numbers(:, i_spin, i_k_point))
       write(use_unit, "(A)") trim(info_str)  ! No collective write; you never know which nodes...
    end if
    acc_trace = acc_trace + k_weights(i_k_point) * trace

  end subroutine calculate_contraco_kdm_one_k
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/calculate_contraco_kdm_all
  ! NAME
  !   calculate_contraco_kdm_all
  ! SYNOPSIS
  subroutine calculate_contraco_kdm_all(contraco_kdm, contraco_kdm_cmplx, &
  &                                     occ_numbers, KS_eigenvector, KS_eigenvector_complex, &
  &                                     trace)

    ! PURPOSE
    !   Calculate the contra-covariant density matrix.
    ! USES
    use scalapack_wrapper, only: my_k_point, my_scalapack_id
    use density_matrix_evaluation
    implicit none

    ! ARGUMENTS

    real*8, intent(OUT) :: contraco_kdm(n_dim1, n_dim2, n_spin, n_k_points_task)
    complex*16, intent(OUT) :: contraco_kdm_cmplx(n_dim1, n_dim2, n_spin, n_k_points_task)
    real*8, intent(IN) :: occ_numbers(n_states, n_spin, n_k_points)
    real*8, intent(IN) :: KS_eigenvector(n_basis, n_states, n_spin, n_k_points_task)
    complex*16, intent(IN) :: KS_eigenvector_complex(n_basis, n_states, n_spin, n_k_points_task)
    real*8, intent(OUT) :: trace

    ! INPUTS
    !   o KS_eigenvector{_complex} -- eigenvectors from which to construct density matrix
    !     These are only used if use_scalapack is not in effect,
    !     otherways eigenvec{_complex} from scalapack_wrapper is used
    ! OUTPUTS
    !   o contraco_kdm, contraco_kdm_cmplx -- resulting co-contra-variant density matrix
    !   o trace -- Trace of the contra-covariant density matrix (should be n_electrons).
    ! AUTHOR
    !   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    ! HISTORY
    !   Release version, FHI-aims (2010).
    ! SOURCE

    real*8 :: real_dummy(1)
    complex*16 :: complex_dummy(1)
    integer :: i_spin, i_k_point, i_k

    trace = 0d0

    i_k = 0
    do i_k_point = 1, n_k_points
       if (use_scalapack) then
          if (my_k_point /= i_k_point) cycle
       else
          if (myid /= modulo(i_k_point, n_tasks) .or. myid > n_k_points) cycle
       end if
       i_k = i_k + 1
       do i_spin = 1, n_spin
          if (real_eigenvectors) then
             call calculate_contraco_kdm_one_k(contraco_kdm(:,:, i_spin, i_k), complex_dummy, &
             &                                 occ_numbers, KS_eigenvector, KS_eigenvector_complex, &
             &                                 i_spin, i_k_point, i_k, trace)
          else
             call calculate_contraco_kdm_one_k(real_dummy, contraco_kdm_cmplx(:,:, i_spin, i_k), &
             &                                 occ_numbers, KS_eigenvector, KS_eigenvector_complex, &
             &                                 i_spin, i_k_point, i_k, trace)
          end if
       end do
    end do
    if (use_scalapack .and. my_scalapack_id > 0) then
       trace = 0.d0   ! no double counting
    end if
    call sync_real_number(trace)

  end subroutine calculate_contraco_kdm_all
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/wf_trace
  !  NAME
  !    wf_trace
  !  SYNOPSIS

  subroutine wf_trace(kmatrix, kmatrix_cmplx, trace, trace_imag)

    !  PURPOSE
    !    Calculate and output the trace.
    !  USES
    use scalapack_wrapper, only: l_row, l_col, my_scalapack_comm_work
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: kmatrix(n_dim1, n_dim2)
    complex*16, intent(IN) :: kmatrix_cmplx(n_dim1, n_dim2)
    real*8, intent(OUT) :: trace, trace_imag

    !  INPUTS
    !    o kmatrix{_cmplx} -- The matrix of which the trace is to calculate
    !  OUTPUTS
    !    o trace, trace_image -- real and imaginary part of trace
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_bas
    real*8  :: tmp(2)
    complex*16 :: trace_cmplx
    character*150 :: info_str

    trace = 0.d0
    trace_imag = 0.d0
    trace_cmplx = (0.d0, 0.d0)
    if (use_scalapack) then
       do i_bas = 1, n_basis
          if(l_row(i_bas)>0 .and. l_col(i_bas)>0) then
             if (real_eigenvectors) then
                trace = trace + kmatrix(l_row(i_bas), l_col(i_bas))
             else
                trace_cmplx = trace_cmplx + kmatrix_cmplx(l_row(i_bas), l_col(i_bas))
             end if
          endif
       end do

       ! We should use sync_real_number/sync_complex_number for the following, but
       ! - there exists no routine sync_complex_number
       ! - sync_real_number doesn't have a way to specify the communicator
       ! therefore we use sync_vector here

       if (real_eigenvectors) then
          tmp(1) = trace
          call sync_vector(tmp, 1, my_scalapack_comm_work)
          trace = tmp(1)
       else
          tmp(1) = dble(trace_cmplx)
          tmp(2) = aimag(trace_cmplx)
          call sync_vector(tmp, 2, my_scalapack_comm_work)
          trace_cmplx = cmplx(tmp(1), tmp(2))
       endif
    else
       do i_bas = 1, n_basis
          if (real_eigenvectors) then
             trace = trace + kmatrix(i_bas, i_bas)
          else
             trace_cmplx = trace_cmplx + kmatrix_cmplx(i_bas, i_bas)
          end if
       end do
    endif
    if (.not. real_eigenvectors) then
       trace = dble(trace_cmplx)
       trace_imag = aimag(trace_cmplx)
    end if

  end subroutine wf_trace
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/wf_extrapolate_cckdm
  !  NAME
  !    wf_extrapolate_cckdm
  !  SYNOPSIS
  subroutine wf_extrapolate_cckdm(contraco_kdm, contraco_kdm_cmplx, is_extrapolated)
    !  PURPOSE
    !    Performs extrapolation of the eigenvectors,
    !    assuming that the time step is constant.  E.g.:
    !      linear    (order == 1): y_{n+1} = 2*y_{n} -   y_{n-1}
    !      quadratic (order == 2): y_{n+1} = 3*y_{n} - 3*y_{n-1} + y_{n-2}
    !    For higher orders and non-default extrapolation points, see
    !    wf_extrapolation_coefficients().
    !  USES
    implicit none
    !  ARGUMENTS
    real*8, intent(OUT) :: contraco_kdm(n_dim1, n_dim2, n_spin, n_k_points_task)
    complex*16, intent(OUT) :: contraco_kdm_cmplx(n_dim1, n_dim2, n_spin, n_k_points_task)
    logical, intent(OUT) :: is_extrapolated
    !  AUTHOR
    !    FHI-aims team.
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  INPUTS
    !    none
    !  OUTPUT
    !    contraco_kdm{_cmplx} -- extrapolated contra-covariant density matrix (k-dependent)
    !    is_extrapolated -- flag if extrapolation has been successful
    !  SEE ALSO
    !    FHI-aims CPC publication (in copyright notice above)
    !  SOURCE
    logical :: use_niklasson
    integer :: n_initial, n_points, order, info, i_point, n_points_tot, current_step
    character*50 :: fmt, type_of_init_guess
    character*100 :: info_str
    real*8, allocatable :: coeff(:)

    ! --- prepare

    last_extrapolated_wf_step = last_saved_wf_step + 1
    current_step = last_extrapolated_wf_step

    select case (wf_extra_type)
    case(WF_EXTRA_FUNC)
       n_initial = 0
       use_niklasson = .false.
    case(WF_EXTRA_NIKL)
       n_initial = 1
       use_niklasson = .true.
    case default
       ! should have been catched before
       write(info_str, &
       &     "('*** wf_extrapolate_cckdm: Invalid extrapolation type')")
       call aims_stop(info_str)
    end select

    n_points = wf_n_steps('converged', current_step)
    if (use_niklasson) then
       n_points = min(n_points, wf_n_steps('initial', current_step) - 1)
    end if
    n_points_tot = n_points + n_initial
    order = min(n_wf_funcs, n_points-1)

    if (n_points_tot <= 1) then
       ! No use in extrapolation
       is_extrapolated = .false.
       return
    else
       is_extrapolated = .true.
    end if


    if (real_eigenvectors) then
       contraco_kdm = 0d0
    else
       contraco_kdm_cmplx = (0d0, 0d0)
    end if

    ! --- actually extrapolate

    if (n_points_tot > 1) then
       allocate(coeff(n_points))
       if (use_niklasson) then
          call wf_extrapolation_nikl_coeffs(n_points, coeff)
          type_of_init_guess = 'initial' ! The standard
          if (current_step < 2*(n_max_wf_steps+1)) then
             ! Make sure that the first 'initial' has the full order.
             type_of_init_guess = 'converged'
          end if
          call wf_retrieve_cckdm(type_of_init_guess, current_step-(n_points+1), -1d0, &
          &                      contraco_kdm, contraco_kdm_cmplx)
       else
          call wf_extrapolation_gen_coeffs(wf_funcs, n_points, coeff)
          write(fmt, *) "(2X,A,", n_points, "ES10.2)"
          write(info_str, fmt) '| Extrapolation coeffs:', coeff
          call localorb_info(info_str)
       end if

       do i_point = 1, n_points
          call wf_retrieve_cckdm('converged', current_step-i_point, coeff(i_point), &
          &                      contraco_kdm, contraco_kdm_cmplx)
       end do

       ! immediately write extrapolated guess for this iteration
       if (use_niklasson) then
          call wf_save_cckdm('initial', current_step, contraco_kdm, contraco_kdm_cmplx)
       end if

       deallocate(coeff)
    end if

  end subroutine wf_extrapolate_cckdm
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/wf_extrapolation_coefficients
  !  NAME
  !    wf_extrapolation_coefficients
  !  SYNOPSIS

  subroutine wf_extrapolation_coefficients(order, n_points, coeff)

    !  PURPOSE
    !    Obtain extrapolation coefficients of given order for a given
    !    number of previous iterations.
    !    
    !    Given a function f known for f(t-mh) for m = 1, ..., n_points
    !    this procedure returns extrapolation coefficients
    !       f(t) = \sum_{m=1}^M B_m f(t-mh),
    !    where M=n_points and B_m=coeff(m).
    !    The error of the extrapolation is expected to vanish in the
    !    first order orders and as many further odd orders as possible
    !    (in order to approximate time-reversal symmetry).  This procedure
    !    refrains from doing least squares fits to get smaller
    !    coefficients as is done in [1].
    !
    !    The most common examples are linear extrapolation (order=1,
    !    n_points=2), which yields coeff == (/-1, 2/), and quadratic
    !    extrapolation (order=2, n_points=3, coeff=(/1, -3, 3/)).
    ! 
    !    See [1] and [4] for details.
    ! 
    !    [1] J. Kolafa, "Numerical integration of equations of motion ...",
    !        Mol. Simul. 18, 193 (1996)
    !    [2] J. Kolafa, "Time-reversible always stable predictor-corrector ...",
    !        J. Comput. Chem. 25, 335 (2004)
    !    [3] T.D. Kuehne, M. Krack, F.R. Mohamed, and M. Parrinello,
    !        "Efficient and accurate Car-Parrinello-like approach ...",
    !        Phys. Rev. Lett. 98, 066401 (2007)
    !    [4] P. Pulay, G. Fogarasi, "Fock matrix dynamics",
    !        Chem. Phys. Lett. 386, 272 (2004)
    !  USES

    use numerical_utilities
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: order, n_points
    real*8, intent(OUT) :: coeff(n_points)

    !  INPUTS
    !    o order -- order of polynomial expansion
    !    o n_points -- number of known iterations (order <= n_points+1)
    !  OUTPUTS
    !    o coeff -- extrapolation coefficients
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  ALGORITHM
    !    See [1].
    !    The idea is to expand the extrapolation D^p(t) in prior time steps:
    !      D^p(t) := \sum_{m=1}^M B_m D(t - mh)
    !    where M=n_points and B_m=coeff(m).  The k-th order of the
    !    extrapolation error D^p(t) - D(t) is proportional to
    !      C_k = \sum_{m=1}^M B_m (-m)^k - \delta_{k,0}.
    !    Thus, one can choose a set of M indices k\in K which should vanish
    !    by solving the corresponding set of linear equations:
    !      0 = C_k = \sum_m (-m)^k B_m - \delta_{k,0} =: A_{km} B_m - \delta_{k,0}
    !    with
    !      A_{km} = (-m)^k.
    !    Here, K is given as {0, ..., order} + as_much_odd_orders_as_possible,
    !    because odd orders introduce time-irreversibility.
    !
    !    Please note that it is crucial for the definition of "odd" to
    !    Taylor-expand around the time step of the extrapolation in contrast
    !    to what a quick look on [2] might suggest.
    !
    !    For order=1 and n_points=K,
    !      coeff(m)=(-1)**(m+1) * m * binom(2*K, K-m) / binom(2*K-2, K-1),
    !    as noted by [2] and [3].
    !
    !    References are given in the PURPOSE section.
    !
    !  SOURCE

    real*8, allocatable :: A(:,:), rhs(:)
    integer :: i_point, i_ord, i, last_odd
    real*8 :: tt
    character(*), parameter :: func = 'wf_extrapolation_coefficients'

    allocate(A(n_points, n_points), rhs(n_points))

    last_odd = ((order-1) / 2) * 2 + 1

    do i_point = 1, n_points
       tt = - dble(i_point)
       do i_ord = 0, order
          A(i_ord+1, i_point) = tt**(i_ord)
       end do
       do i = 1, n_points - order - 1
          i_ord = last_odd + 2*i
          A(order+1+i, i_point) = tt**(i_ord)
       end do
    end do
    rhs = 0.d0
    rhs(1) = 1.d0
    
    call solve_LEQ(func, n_points, A, rhs, coeff)

    deallocate(A, rhs)

  end subroutine wf_extrapolation_coefficients
  !******
  !------------------------------------------------------------------------------
  !****s* wf_extrapolation/wf_extrapolation_nikl_coeffs
  !  NAME
  !    wf_extrapolation_nikl_coeffs
  !  SYNOPSIS

  subroutine wf_extrapolation_nikl_coeffs(n_points, coeff)

    !  PURPOSE
    !
    !    Obtain extrapolation coefficients for a given number of previous
    !    iterations.  By taking an old extrapolation into account, this scheme
    !    is fully time-reversal and shoul avoid any systematic energy drift.
    ! 
    !    Compute extrapolation coefficients B_m in
    !      D^p(st) := \sum_{m=-(s-1)}^{s-1} B_m D(t + mh) - D^p(-st)
    !    where s=0.5*(n_points+1), and B_m=coeff(s+m).
    !
    !    As quite a thorough (and constant) scf convergence limit is assumed,
    !    stability should be no issue.  Thus, any degrees of freedom can be
    !    used to include higher orders.
    !
    !    Please note that information of n_points+1 iterations is necessary
    !    as the initializing wave functions of the oldest iteration is
    !    needed to actually ensure time reversibility.
    !
    !    [1] A.M.N. Niklasson, C.J. Tymczak, and M. Challacombe,
    !        "Time-reversible Born-Oppenheimer molecular dynamics"
    !        Phys. Rev. Lett. 97, 123001 (2006).
    !
    !  USES

    use numerical_utilities
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_points
    real*8, intent(OUT) :: coeff(n_points)

    !  INPUTS
    !    o n_points -- number of known iterations
    !  OUTPUTS
    !    o coeff -- extrapolation coefficients
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  ALGORITHM
    !    See [1].
    !    The idea is to expand the extrapolation D^p(st) in prior time steps:
    !      D^p(st) := \sum_{m=-(s-1)}^{s-1} B_m D(t - mh) - D^p(-st)
    !    where s = (M+1)/2, M=n_points, B_{-m}=B_m and coeff(i)=B_{s-i}.
    !    The error in D^p(st) then is the error in D^p(-st) plus
    !      \sum_k=0^\infty (sh)^k/k! D^{(k)}(t) C_k
    !    with
    !      C_k = 0                                     | for k odd
    !      C_k = \sum_{m=-(s-1)}^{s-1} (m/s)^k B_m - 2 | for k even.
    !    As many even orders as possible will be removed (C_k == 0).
    !  SOURCE

    real*8, allocatable :: A(:,:), rhs(:), ind_coeff(:)
    integer :: n_indep, i_indep, i, i_ord, i_point
    real*8 :: tt, s
    character(*), parameter :: func = 'wf_extrapolation_coefficients'

    s = dble(n_points+1) / 2.d0
    n_indep = (n_points+1) / 2

    allocate(A(n_indep, n_indep), rhs(n_indep), ind_coeff(n_indep))

    do i_indep = 1, n_indep
       tt = (s - dble(i_indep)) / s
       do i = 1, n_indep
          i_ord = 2*(i-1)
          if (2*i_indep == n_points+1) then ! tt == null
             A(i, i_indep) = tt**(i_ord)
          else
             A(i, i_indep) = 2.d0 * tt**(i_ord)
          end if
       end do
    end do
    rhs = 2.d0
    
    call solve_LEQ(func, n_indep, A, rhs, ind_coeff)
    do i_point = 1, n_points
       if (i_point <= n_indep) then
          coeff(i_point) = ind_coeff(i_point)
       else
          coeff(i_point) = ind_coeff(n_points+1-i_point)
       end if
    end do

    deallocate(A, rhs, ind_coeff)

  end subroutine wf_extrapolation_nikl_coeffs
  !******
  !----------------------------------------------------------------------------
  !****s* wf_extrapolation/wf_extrapolation_gen_coeffs
  !  NAME
  !    wf_extrapolation_gen_coeffs
  !  SYNOPSIS

  subroutine wf_extrapolation_gen_coeffs(funcs, n_points, coeff)

    !  PURPOSE
    !    Obtain extrapolation coefficients of given order for a given
    !    number of previous iterations.  Uses the basis functions given in
    !    funcs to do the extrapolation.
    !  USES

    use numerical_utilities
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_points
    type(gen_func), intent(IN) :: funcs(n_points)
    real*8, intent(OUT) :: coeff(n_points)

    !  INPUTS
    !    o n_points -- number of known iterations (order <= n_points+1)
    !    o funcs -- functions as defined in general_function.f90.
    !  OUTPUTS
    !    o coeff -- resulting extrapolation coefficients.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  ALGORITHM
    !    See [1].
    !    The idea is to expand the extrapolator f~(t) in prior time steps:
    !      f~(t) := \sum_{m=1}^M b_m f(t - mh)
    !    where M=n_points and b_m=coeff(m).  The coefficents are calculated
    !    by a least squares fit of
    !      F = \sum_{m=1}^M | f~(t - mh) - f(t - mh) |^2
    !    where f~(t) is expanded in the basis functions \phi_j
    !      f~(t) = \sum_{j=1}^N c_j \phi_j(t)
    !    given in funcs(:), with N<=M (N<M for any(funcs == 'none'),
    !    Formally, this leads to the least squares problem:
    !      F(c) = \sum_m | \sum_j \phi_j(t_m) c_j - f_m |^2 =: || A c - f ||^2.
    !    The general solution depends linearly on f and can be written as
    !    c = A^+ f, with A^+ being the pseudo-inverse of A.  But as we are
    !    interested in the extrapolated
    !    values of f~(t_0) at the current time step t_0=0, f~ has to be
    !    evaluated:
    !      f~(t_0) = \sum_j \phi_j(t_0) c_j =: phi c = phi A^+ f.
    !    Therefore, the expansion coefficients b_m named above are given by
    !      b_m = (A^+T phi)_m
    !    which is the solution of the least squares problem
    !      A^T b = phi
    !      <=> min \sum_j | \sum_m \phi_j(t_m) b_m - \phi_j(t_0) |^2.
    !    The 'none' basis functions are some kind of 'trick' to get the
    !    desired least squares behaviour.  The result in zero columns
    !    in A = A_{mj} = \phi_j(t_m) which reduce the rank of A.  The least
    !    squares solver should detect this.
    !  SOURCE

    real*8, allocatable :: AT(:,:), phi(:), values(:)
    integer :: m_time, j_func, rank
    real*8 :: tt
    character*150 :: info_str
    character*(*), parameter :: func = 'wf_extrapolation_gen_coeffs'

    allocate(AT(n_points, n_points), phi(n_points))
    allocate(values(0:n_points))

    do j_func = 1, n_points
       call eval_general_function(funcs(j_func), n_points, values, MD_tstep)
       phi(j_func) = values(0)
       AT(j_func,:) = values(1:n_points)
    end do
    call solve_LSQ(func, n_points, n_points, AT, phi, coeff, rank)
    if (rank /= n_points - number_of_nones(funcs)) then
       write(info_str, &
       &     "('** Got rank',I5,' with',I5,' points and',I5,' nones.')") &
       & rank, n_points, number_of_nones(funcs)
    end if
    deallocate(AT, phi, values)

  end subroutine wf_extrapolation_gen_coeffs
  !******

end module wf_extrapolation
!******
