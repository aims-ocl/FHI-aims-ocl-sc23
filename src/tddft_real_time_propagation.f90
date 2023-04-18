module tddft_real_time_propagation
  use physics
  use dimensions
  use mixing
  use runtime_choices
  use grids
  use geometry
  use basis
  use constants
  use localorb_io
  use lapack_wrapper
  use scalapack_wrapper
  use species_data, only: l_shell_max
  use meta_gga_postp

  implicit none

  private 

  public ::                               &
       tddft_real_time_propagation_init,  &
       tddft_real_time_propagation_run,   &
       tddft_real_time_propagation_end

  character(len=100) :: filename, file_prefix
 
  logical :: write_info

  integer ::                                           & 
       i_basis_1, i_basis_2, i_index, i_spin, i_sum,   &
       td_output_every, i_td_step, td_steps, info,     &
       pdsyev_lwork, i_states, i_atom

  integer, dimension (:), allocatable ::               &
       scalapack_ipiv

  real*8 :: dt, cpu_seconds, timing_time_step

  real*8, dimension(:), allocatable ::                 &
       ovlp_eigenvalues, work, pdsyev_work,            & 
       tdks_eigenvalues

  real*8, dimension(:,:), allocatable ::               & 
       overlap_work, overlap_m_sqrt, overlap_p_sqrt,   &
       range_space, null_space, hamiltonian_square,    &
       hamiltonian_extrapolation_1, hamiltonian_step,  &
       hamiltonian_extrapolation_2, hamiltonian_work,  &
       scalapack_work_1, scalapack_work_2,             &
       scalapack_work_3, hamiltonian_ev

  complex*16, dimension(:,:), allocatable ::                & 
       zoverlap_m_sqrt, zoverlap_p_sqrt, zoverlap_matrix,   &
       scalapack_zwork_1, scalapack_zwork_2, propagator_U,  &
       scalapack_zwork_3, tdks_eigenvalues_phases, zwork_1, &
       td_KS_eigenvector, td_KS_gs_eigenvector,             &
       td_KS_initial_eigenvector, td_KS_eigenvector_work 


contains

  ! ---------------------------------------------------------
  subroutine tddft_real_time_propagation_init()
    use aims_memory_tracking, only : aims_allocate, aims_deallocate 

    ! for the moment we have only the parallel version
    if (.not.use_scalapack) then
      if (myid.eq.0) then
        write(use_unit,'(A)') '------------------------------------------------------------'
        write(use_unit,'(A)') 'Warning:  TDDFT real-time propagation only implemented for  '
        write(use_unit,'(A)') '          Scalapack MPI. Stopping FHI-aims.                 '
        write(use_unit,'(A)') '------------------------------------------------------------'
      end if
      call aims_stop
    end if
    ! also we can only treat finite systems at the moment
    if (n_periodic > 0) then
      if (myid.eq.0) then
        write(use_unit,'(A)') '------------------------------------------------------------'
        write(use_unit,'(A)') 'Warning:  TDDFT real-time propagation only implemented for  '
        write(use_unit,'(A)') '            finite systems. Stopping FHI-aims.              '
        write(use_unit,'(A)') '------------------------------------------------------------'
      end if
      call aims_stop
    end if
    
    ! allocate work arrays
    
    allocate(hamiltonian_square(n_basis, n_basis))
    
    allocate(hamiltonian_extrapolation_1(n_basis, n_basis))
    allocate(hamiltonian_extrapolation_2(n_basis, n_basis))
    allocate(hamiltonian_step(n_basis, n_basis))
    
    allocate(scalapack_ipiv(n_basis))
    allocate(scalapack_work_1(mxld, mxcol))
    allocate(scalapack_work_2(mxld, mxcol))
    allocate(scalapack_work_3(mxld, mxcol))
    allocate(scalapack_zwork_1(mxld, mxcol))
    allocate(scalapack_zwork_2(mxld, mxcol))
    allocate(scalapack_zwork_3(mxld, mxcol))
    allocate(hamiltonian_ev(mxld, mxcol))
    
    allocate(hamiltonian_work(n_basis, n_basis))
    allocate(tdks_eigenvalues(n_basis))
    allocate(tdks_eigenvalues_phases(n_basis,n_basis))
    allocate(propagator_U(n_basis, n_basis))
    allocate(td_KS_eigenvector(n_basis, n_basis))
    allocate(td_KS_gs_eigenvector(n_basis, n_basis))
    allocate(td_KS_initial_eigenvector(n_basis, n_basis))
    allocate(td_KS_eigenvector_work(n_basis, n_basis))
    
    allocate(ovlp_eigenvalues(n_basis))
    allocate(zwork_1(n_basis,n_basis))
    allocate(overlap_m_sqrt(n_basis,n_basis))
    allocate(zoverlap_m_sqrt(n_basis,n_basis))
    allocate(overlap_p_sqrt(n_basis,n_basis))
    allocate(zoverlap_p_sqrt(n_basis,n_basis))
    allocate(overlap_work(n_basis,n_basis))
    allocate(zoverlap_matrix(n_basis,n_basis))
    allocate(range_space(n_basis,n_basis))
    allocate(null_space(n_basis,n_basis))

    if(allocated(KS_eigenvector_complex)) call aims_deallocate(KS_eigenvector_complex, "KS_eigenvector_complex")
    call aims_allocate(KS_eigenvector_complex, n_basis, n_states, n_spin, n_k_points_task, "KS_eigenvector_complex")

  end subroutine tddft_real_time_propagation_init

  ! ---------------------------------------------------------
  ! Time-stepping of time-dependent Kohn-Sham equations
  subroutine tddft_real_time_propagation_run()

    use mpi_utilities, only: get_my_task
    implicit none

    integer :: i_coord
    real*8  :: dipole_ele_moment(3)
    real*8  :: dipole_ion_moment(3)
    real*8  :: penalty_energy

    if (myid.eq.0) then
      write(use_unit,'(A)') '------------------------------------------------------------'
      write(use_unit,'(A)') '                 TDDFT real-time propagation                '
      write(use_unit,'(A)') '------------------------------------------------------------'
      write(use_unit,'(A)') ''
    end if

    call get_my_task()
    real_eigenvectors = .false.

    ! ---------------------------------------------------------
    ! Compute square root of overlap
    if (myid.eq.0) then
      write(use_unit,'(A)') ' Computing square root of overlap matrix'
    end if

    call setup_overlap_scalapack(overlap_matrix)

    i_index = 0
    do i_basis_2 = 1, n_basis, 1
      do i_basis_1 = 1, i_basis_2, 1
        i_index = i_index+1
        overlap_work(i_basis_1,i_basis_2) = overlap_matrix(i_index)
        overlap_work(i_basis_2,i_basis_1) = overlap_matrix(i_index)
      enddo
    enddo
    zoverlap_matrix = dcmplx(overlap_work)

    i_index = 0
    do i_basis_2 = 1, n_basis, 1
      do i_basis_1 = 1, i_basis_2, 1
        i_index = i_index+1
        overlap_work(i_basis_1,i_basis_2) = overlap_matrix(i_index)
        overlap_work(i_basis_2,i_basis_1) = overlap_matrix(i_index)
      enddo
    enddo
    call setup_scalapack_rmatrix(overlap_work, scalapack_work_1)
    call scalapack_pdsyev(scalapack_work_1, ovlp_eigenvalues, hamiltonian_ev)
    call get_scalapack_global_rmatrix(hamiltonian_ev, range_space)
    null_space(1:n_basis,1:n_basis) = &
       transpose(range_space(1:n_basis,1:n_basis))    

    ! ---------------------------------------------------------    
    ! compute S^(-1/2) and S^(+1/2) from eigenvalues
    overlap_m_sqrt = M_ZERO
    overlap_p_sqrt = M_ZERO
    do i_basis_1 = 1, n_basis, 1
      do i_basis_2 = 1, n_basis, 1
        do i_sum = 1, n_basis, 1
          overlap_m_sqrt(i_basis_1,i_basis_2) =      &
               overlap_m_sqrt(i_basis_1,i_basis_2) + &
               range_space(i_basis_1, i_sum)/sqrt(ovlp_eigenvalues(i_sum))*null_space(i_sum, i_basis_2) 
          overlap_p_sqrt(i_basis_1,i_basis_2) =      &
               overlap_p_sqrt(i_basis_1,i_basis_2) + &
               range_space(i_basis_1, i_sum)*sqrt(ovlp_eigenvalues(i_sum))*null_space(i_sum, i_basis_2) 
        end do
      end do
    end do
    zoverlap_m_sqrt = dcmplx(overlap_m_sqrt)
    zoverlap_p_sqrt = dcmplx(overlap_p_sqrt)
    
    ! query workspace for pdsyev
    allocate(pdsyev_work(1))
    
    call pdsyev( 'V', 'U', n_basis, scalapack_work_2, &
         1, 1, sc_desc, tdks_eigenvalues, hamiltonian_ev, 1, 1, sc_desc, pdsyev_work, -1, info )
    pdsyev_lwork = pdsyev_work(1)
    
    deallocate(pdsyev_work)
    allocate(pdsyev_work(pdsyev_lwork))
    

    ! ---------------------------------------------------------    
    ! Start time propagation
    ! ---------------------------------------------------------    
    dt              = TDDFT_dt
    td_output_every = TDDFT_output_every
    td_steps = ceiling(TDDFT_time / TDDFT_dt)
    
    if (myid.eq.0) then
      write(use_unit,'(2X, A, f18.12)')  '| Using time step: ', dt
      write(use_unit,'(2X, A, i10)')     '| Number of time steps: ', td_steps
      write(use_unit,'(2X, A, i8, A)')   '| Output will be written every', td_steps , ' steps'
    end if
    
    ! set ground state as initial state
    do i_basis_1 = 1, n_basis, 1
      if (i_basis_1 .le. n_states) then
        td_KS_eigenvector(1:n_basis,i_basis_1) = KS_eigenvector(1:n_basis,i_basis_1,1,1)
      else
        td_KS_eigenvector(1:n_basis,i_basis_1) = M_ZERO  
      end if
    end do
    td_KS_eigenvector_work = td_KS_eigenvector
    td_KS_gs_eigenvector   = td_KS_eigenvector
    
    td_KS_initial_eigenvector = td_KS_eigenvector
    KS_eigenvector_complex(:,:,1,1) = dcmplx(td_KS_eigenvector)
    
    ! used for extrapolation of the hamiltonian
    i_spin = 1
    i_index = 0
    hamiltonian_square = M_ZERO
    do i_basis_2 = 1, n_basis
      do i_basis_1 = 1, i_basis_2
        i_index = i_index + 1
        hamiltonian_square(i_basis_1, i_basis_2) = hamiltonian(i_index, i_spin)
        hamiltonian_square(i_basis_2, i_basis_1) = hamiltonian_square(i_basis_1, i_basis_2)
      end do
    end do
    hamiltonian_extrapolation_1 = hamiltonian_square
    hamiltonian_extrapolation_2 = hamiltonian_square

    info = 0
    if (myid.eq.0) then
      write(use_unit,'(A)') ' Starting time propagation'
    end if
    
    open(55, file='td_eigenvalues.out')
    open(56, file='td_occupations.out')

    ! ---------------------------------------------------------    
    ! main time loop
    do i_td_step = 1, td_steps
    
      !if(i_td_step.eq.200) homogeneous_field = M_ZERO

      if (mod(i_td_step,TDDFT_write_info).eq.0) then
        write_info = .true.
      else
        write_info = .false.
      end if
      
      call cpu_time( cpu_seconds )
      timing_time_step = cpu_seconds
      
      ! ---------------------------------------------------------    
      ! Begin of time step
      if (myid.eq.0.and.write_info) then
        write(use_unit,'(A)') "----------------------------------------------------------------"
        write(use_unit,'(2X, A, i7)')        'Time step:', i_td_step
        write(use_unit,'(2X, A, E20.11, A)') '| Simulation time:', i_td_step*dt, ' a.u.'
      end if
      
      ! predictor step
      hamiltonian_step = &
           - M_THREE/M_FOUR * hamiltonian_extrapolation_1  &
           + M_SEVEN/M_FOUR * hamiltonian_extrapolation_2 
      
      ! perform predictor step
      call tddft_rti(KS_eigenvector_complex, hamiltonian_step)
      
      ! update hamiltonian
      call update_hamiltonian(hamiltonian_step)
      
      ! perform corrector step
      call tddft_rti(KS_eigenvector_complex, hamiltonian_step)

      ! update eigenstates
      td_KS_eigenvector(1:n_basis,1:n_states) = KS_eigenvector_complex(1:n_basis,1:n_states,1,1)
      ! End of time step
      ! ---------------------------------------------------------    
      
      if (write_info) then
        call get_penalty_energy(KS_eigenvector, KS_eigenvector_complex, &
        &                       occ_numbers, condition_penalty, penalty_energy)
        call get_total_energy &
        ( ev_sum, ev_sum_shifted, en_xc, en_pot_xc, en_ion_ion, &
          en_ion_embed, en_density_embed, &
          en_vdw, en_pot_vdw, en_ll_vdw, en_ll_vdw_err, en_lda_c, en_pbe_c, &
          en_elec_free, en_elec_delta, fock_energy, &
          hartree_energy_free, hartree_delta_energy, &
          hartree_multipole_correction, &
          entropy_correction, penalty_energy, total_energy &
        )
      end if

      ! compute dipole moment
      ! call output_dipole_moment()
      call evaluate_moment_p1(M_ONE,rho,partition_tab,dipole_ele_moment,dipole_ion_moment)
      if (myid.eq.0) then
        write(use_unit,'(2X,A,E16.8,A,E30.15,A,3E30.15)') &
             "| Time: ", i_td_step*dt," Energy: ", total_energy," Electronic dipole moment [eAng]:", &
              (dipole_ele_moment(i_coord)*bohr, i_coord=1,3,1)
      end if

          
      ! ---------------------------------------------------------          
      ! compute observables
      call setup_scalapack_full_zmatrix( td_KS_initial_eigenvector, scalapack_zwork_1 )
      call setup_scalapack_full_rmatrix( overlap_work, scalapack_work_2 )
      scalapack_zwork_2 = cmplx(scalapack_work_2)
      call pzgemm ( 'T', 'N', n_basis, n_basis, n_basis, C_ONE, scalapack_zwork_1, 1, 1, sc_desc, &
           scalapack_zwork_2, 1, 1, sc_desc, C_ZERO, scalapack_zwork_3, 1, 1, sc_desc )
      call setup_scalapack_full_zmatrix( td_KS_eigenvector, scalapack_zwork_1 )
      call pzgemm ( 'N', 'N', n_basis, n_basis, n_basis, C_ONE, scalapack_zwork_3, 1, 1, sc_desc, &
           scalapack_zwork_1, 1, 1, sc_desc, C_ZERO, scalapack_zwork_2, 1, 1, sc_desc )
      call get_scalapack_global_zmatrix( scalapack_zwork_2, td_KS_eigenvector_work )
      
      ! eigenvalues      
      call setup_scalapack_full_rmatrix( hamiltonian_step, scalapack_work_2 )
      ! diagonalize transformed predictor hamiltonian
      call pdsyev( 'V', 'U', n_basis, scalapack_work_2, &
           1, 1, sc_desc, tdks_eigenvalues, hamiltonian_ev, 1, 1, sc_desc, pdsyev_work, pdsyev_lwork, info )

      if (td_output_every.ne.0.and.(mod(i_td_step,td_output_every).eq.0.or.i_td_step.eq.1)) then
        write(filename, '(A,I7.7,A)') 'time_step_', i_td_step, '_td_KS_eigenvector.real'
        call scalapack_output_global_matrix(n_basis, real(td_KS_eigenvector), filename)
        write(filename, '(A,I7.7,A)') 'time_step_', i_td_step, '_td_KS_eigenvector.imag'
        call scalapack_output_global_matrix(n_basis, aimag(td_KS_eigenvector), filename)
        write(filename, '(A,I7.7,A)') 'time_step_', i_td_step, '_autocorrelation.out'

        open(45, file=filename)
        do i_basis_1 = 1, n_basis, 1
          write(45, '(5E20.11)') i_td_step*dt, occ_numbers(i_basis_1,1,1), td_KS_eigenvector_work(i_basis_1,i_basis_1)
        end do
        close(45)      

        write(filename, '(A,I7.7,A)') 'time_step_', i_td_step, '_tdks_eigenvalues.out'
        open(45, file=filename)
        do i_basis_1 = 1, n_basis, 1
          write(45, '(5E20.11)') i_td_step*dt, tdks_eigenvalues(i_basis_1), occ_numbers(i_basis_1,1,1)
        end do
        close(45)

        if (use_cube_output) then
          call cpu_time( cpu_seconds )
          timing_time_step = cpu_seconds
          write(file_prefix,'(A,I7.7,A)') 'time_step_', i_td_step, '_'
          ! call output_cube_files_p1(file_prefix) 
          call cpu_time( cpu_seconds )
          if (myid.eq.0.and.write_info) then
            write(use_unit,*) 'Timing for output: ', cpu_seconds - timing_time_step, ' s'
          end if
        end if
      end if

      write(55, '(i8,500E20.11)') i_td_step,i_td_step*dt, tdks_eigenvalues(:)
      write(56, '(i8,500E20.11)') i_td_step,i_td_step*dt, occ_numbers(:,1,1)

      call cpu_time( cpu_seconds )
      if (myid.eq.0.and.write_info) then
        write(use_unit,*) ' Timing for time step: ', cpu_seconds - timing_time_step, ' s'
      end if
 
    end do
    close(55)
    close(56)


    
  contains
    ! ---------------------------------------------------------
    subroutine tddft_rti(KS_eigenvec_cmplx, hamiltonian_rti)
      complex*16, intent(inout) :: KS_eigenvec_cmplx(:, :, :, :)
      real*8,     intent(in)    :: hamiltonian_rti(:, :)

      select case(TDDFT_propagator)
      case('exact_matrix_exponential')
        call tddft_real_time_propagation_exp_midpoint_step(KS_eigenvector_complex, hamiltonian_step)
      case('crank_nicholson')
        call tddft_real_time_propagation_crank_nicholson(KS_eigenvector_complex, hamiltonian_step)
      end select

    end subroutine tddft_rti


    ! ---------------------------------------------------------
    subroutine tddft_real_time_propagation_exp_midpoint_step(KS_eigenvec_cmplx, hamiltonian_rti)
      complex*16, intent(inout) :: KS_eigenvec_cmplx(:, :, :, :)
      real*8,     intent(in)    :: hamiltonian_rti(:, :)

      if (myid.eq.0.and.write_info) then
        write(use_unit,'(2X, A, i7)')        '| Computing S^(-1/2)*H*S^(-1/2)'
      end if
      ! S^(-1/2)*H*S^(-1/2)
      call setup_scalapack_full_rmatrix( overlap_m_sqrt, scalapack_work_1 )
      call setup_scalapack_full_rmatrix( hamiltonian_rti, scalapack_work_2 )
      call pdgemm ( 'N', 'N', n_basis, n_basis, n_basis, M_ONE, scalapack_work_1, 1, 1, sc_desc, &
           scalapack_work_2, 1, 1, sc_desc, M_ZERO, scalapack_work_3, 1, 1, sc_desc )               
      call pdgemm ( 'N', 'N', n_basis, n_basis, n_basis, M_ONE, scalapack_work_3, 1, 1, sc_desc, &
           scalapack_work_1, 1, 1, sc_desc, M_ZERO, scalapack_work_2, 1, 1, sc_desc )                   
      scalapack_zwork_2 = scalapack_work_2
      
      if (myid.eq.0.and.write_info) then
        write(use_unit,'(2X, A, i7)')        '| Diagonalizing S^(-1/2)*H*S^(-1/2)'
      end if
      ! diagonalize transformed predictor hamiltonian
      call pdsyev( 'V', 'U', n_basis, scalapack_work_2, &
           1, 1, sc_desc, tdks_eigenvalues, hamiltonian_ev, 1, 1, sc_desc, pdsyev_work, pdsyev_lwork, info )
      
      if (myid.eq.0.and.write_info) then
        write(use_unit,'(2X, A, i7)')        '| Constructing U(t+delta t/2,t)'
      end if
      ! construct U(t+delta t/2,t)
      tdks_eigenvalues_phases = M_ZERO
      do i_basis_1 = 1, n_basis, 1
        tdks_eigenvalues_phases(i_basis_1,i_basis_1) = exp(-C_IMAG_ONE*tdks_eigenvalues(i_basis_1)*dt/M_TWO)
      end do
      
      call setup_scalapack_full_zmatrix( tdks_eigenvalues_phases, scalapack_zwork_1 )
      scalapack_zwork_2 = dcmplx(hamiltonian_ev)
      call pzgemm ( 'N', 'T', n_basis, n_basis, n_basis, C_ONE, scalapack_zwork_1, 1, 1, sc_desc, &
           scalapack_zwork_2, 1, 1, sc_desc, C_ZERO, scalapack_zwork_3, 1, 1, sc_desc )           
      call pzgemm ( 'N', 'N', n_basis, n_basis, n_basis, C_ONE, scalapack_zwork_2, 1, 1, sc_desc, &
           scalapack_zwork_3, 1, 1, sc_desc, C_ZERO, scalapack_zwork_1, 1, 1, sc_desc )                   
      
      if (myid.eq.0.and.write_info) then
        write(use_unit,'(2X, A, i7)')        '| Constructing S^(-1/2)*U*S^(+1/2)'
      end if
      ! S^(-1/2)*U*S^(+1/2)
      call setup_scalapack_full_zmatrix( zoverlap_m_sqrt, scalapack_zwork_2 )
      call pzgemm ( 'N', 'N', n_basis, n_basis, n_basis, C_ONE, scalapack_zwork_2, 1, 1, sc_desc, &
           scalapack_zwork_1, 1, 1, sc_desc, C_ZERO, scalapack_zwork_3, 1, 1, sc_desc )                
      call setup_scalapack_full_zmatrix( zoverlap_p_sqrt, scalapack_zwork_2 )
      call pzgemm ( 'N', 'N', n_basis, n_basis, n_basis, C_ONE, scalapack_zwork_3, 1, 1, sc_desc, &
           scalapack_zwork_2, 1, 1, sc_desc, C_ZERO, scalapack_zwork_1, 1, 1, sc_desc )                   
      
      if (myid.eq.0.and.write_info) then
        write(use_unit,'(2X, A, i7)')        '| Propagating states'
      end if
      ! propagate states from t to dt/2
      call setup_scalapack_full_zmatrix( td_KS_eigenvector, scalapack_zwork_2 )
      call pzgemm ( 'N', 'N', n_basis, n_basis, n_basis, C_ONE, scalapack_zwork_1, 1, 1, sc_desc, &
           scalapack_zwork_2, 1, 1, sc_desc, C_ZERO, scalapack_zwork_3, 1, 1, sc_desc )
      call get_scalapack_global_zmatrix( scalapack_zwork_3, td_KS_eigenvector_work )
      KS_eigenvec_cmplx(1:n_basis,1:n_states,1,1) = td_KS_eigenvector_work(1:n_basis,1:n_states)

    end subroutine tddft_real_time_propagation_exp_midpoint_step


    ! ---------------------------------------------------------
    subroutine tddft_real_time_propagation_crank_nicholson(KS_eigenvec_cmplx, hamiltonian_rti)
      complex*16, intent(inout) :: KS_eigenvec_cmplx(:, :, :, :)
      real*8,     intent(in)    :: hamiltonian_rti(:, :)

      if (myid.eq.0.and.write_info) then
        write(use_unit,'(2X, A, i7)')        '| Computing |phi> = (S - i H dt/2)|psi>'
      end if

      propagator_U = dcmplx(overlap_work) - C_IMAG_ONE*dt/(M_FOUR*M_TWO)*hamiltonian_rti
      call setup_scalapack_full_zmatrix( propagator_U, scalapack_zwork_1 )
      call setup_scalapack_full_zmatrix( KS_eigenvec_cmplx(1:n_basis,1:n_states,1,1), scalapack_zwork_2 )
      call pzgemm ( 'N', 'N', n_basis, n_basis, n_basis, C_ONE, scalapack_zwork_1, 1, 1, sc_desc, &
           scalapack_zwork_2, 1, 1, sc_desc, C_ZERO, scalapack_zwork_3, 1, 1, sc_desc )

      if (myid.eq.0.and.write_info) then
        write(use_unit,'(2X, A, i7)')        '| Solving (S + i H dt/2)|psi> = |phi>'
      end if

      propagator_U = dcmplx(overlap_work) + C_IMAG_ONE*dt/(M_FOUR*M_TWO)*hamiltonian_rti
      call setup_scalapack_full_zmatrix( propagator_U, scalapack_zwork_1 )
      call pzgesv( n_basis, n_states, scalapack_zwork_1, 1, 1, sc_desc, scalapack_ipiv, scalapack_zwork_3, 1, 1, sc_desc, info )
      call get_scalapack_global_zmatrix( scalapack_zwork_3, zwork_1 )
      KS_eigenvec_cmplx(1:n_basis,1:n_states,1,1) = zwork_1(1:n_basis,1:n_states)

    end subroutine tddft_real_time_propagation_crank_nicholson



    ! ---------------------------------------------------------
    subroutine update_hamiltonian(hamiltonian_rti)
      real*8,     intent(inout)    :: hamiltonian_rti(:, :)
      
      integer :: i_coords

      real_eigenvectors = .false.
      pulay_forces_on = .false.
      gga_forces_on = .false.
      meta_gga_forces_on = .false.
      nlcc_forces_on = .false.

      ! AJL, updated March 2018
       call update_density_and_forces_orbital &
             ( KS_eigenvector, KS_eigenvector_complex, &
             KS_eigenvalue, occ_numbers, &
             partition_tab, hartree_partition_tab, rho, rho_gradient, &
             kinetic_density, hartree_potential, l_shell_max, &
             delta_rho_KS, delta_rho_gradient, delta_kinetic_density, rho_change, &
             hellman_feynman_forces, pulay_forces, &
             gga_forces, nlcc_forces, pulay_forces_on, &
             gga_forces_on, nlcc_forces_on, meta_gga_forces_on)

      if (myid.eq.0.and.write_info) then 
        write(use_unit,'(2X,A,1X,E10.4,1X,E10.4)') &
             "| Change of charge/spin density :", (rho_change(i_spin),i_spin=1,n_spin,1)
      end if
      
      if (n_spin .eq. 1) then
        i_spin = 1
        rho(i_spin,:) = rho(i_spin,:) + delta_rho_KS(:,i_spin)
        if (use_density_gradient) then
          do i_coords = 1, 3
            rho_gradient(i_coords,i_spin,:) = rho_gradient(i_coords,i_spin,:) + delta_rho_gradient(i_coords,:,i_spin)
          end do
          if (use_meta_gga) then
            kinetic_density(i_spin,:) = kinetic_density(i_spin,:) + delta_kinetic_density(:,i_spin)
          end if
        end if
      else
        rho(1,:) = rho(1,:) + M_HALF * ( delta_rho_KS(:,1) +  delta_rho_KS(:,2) )
        rho(2,:) = rho(2,:) + M_HALF * ( delta_rho_KS(:,1) -  delta_rho_KS(:,2) )
        if (use_density_gradient) then
          do i_coords = 1, 3
            rho_gradient(i_coords,1,:) = rho_gradient(i_coords,1,:) + &
                 M_HALF * (delta_rho_gradient(i_coords,:,1) + delta_rho_gradient(i_coords,:,2))
            rho_gradient(i_coords,2,:) = rho_gradient(i_coords,2,:) + &
                 M_HALF * (delta_rho_gradient(i_coords,:,1) - delta_rho_gradient(i_coords,:,2))
          end do
          ! Copied from scf_solver.f90. AJL, March 2018
          if (use_meta_gga) then
            kinetic_density(1,:) = kinetic_density(1,:) + M_HALF * ( delta_kinetic_density(:,1) +  delta_kinetic_density(:,2) )
            kinetic_density(2,:) = kinetic_density(2,:) + M_HALF * ( delta_kinetic_density(:,1) -  delta_kinetic_density(:,2) )
          endif
        end if
      end if
      
     call update_hartree_potential_p1 &
         ( hartree_partition_tab, free_rho_superpos, &
         rho, rho_multipole_old, &
         delta_v_hartree_part_at_zero, &
         delta_v_hartree_deriv_l0_at_zero, &
         multipole_moments, multipole_radius_sq, &
         l_hartree_max_far_distance, rho_radial_multipole_old, &
         outer_potential_radius, AS_stress_on, .true. &
         )

      use_embedding_potential = .false.
      call sum_up_whole_potential_p1 &
           ( delta_v_hartree_part_at_zero, &
                delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
                partition_tab, rho, &
                free_hartree_superpos, free_rho_superpos,  &
                hartree_potential,  &
                hartree_delta_energy, en_elec_delta, hartree_multipole_correction, &
                pot_ion_embed, en_density_embed, &
                multipole_forces, forces_on, multipole_radius_sq, &
                l_hartree_max_far_distance, hellman_feynman_forces, &
                energy_deriv_tress_tensor, rho_multipole_old, &
                outer_potential_radius, AS_stress_on, .true. )

      call integrate_real_hamiltonian_matrix_p2 &
           ( hartree_potential,  rho, rho_gradient, kinetic_density, partition_tab, &
         l_shell_max, en_xc, en_pot_xc, hamiltonian, en_vdw, en_pot_vdw )
      
      i_spin = 1
      i_index = 1
      hamiltonian_square = M_ZERO
      do i_basis_2 = 1, n_basis
        do i_basis_1 = 1, i_basis_2
          hamiltonian_square(i_basis_1, i_basis_2) = hamiltonian(i_index, i_spin)
          hamiltonian_square(i_basis_2, i_basis_1) = hamiltonian_square(i_basis_1, i_basis_2)
          i_index = i_index + 1
        end do
      end do
      
      hamiltonian_extrapolation_1 = hamiltonian_extrapolation_2
      hamiltonian_extrapolation_2 = hamiltonian_square
      hamiltonian_step = hamiltonian_square
      hamiltonian_rti = hamiltonian_square
      
    end subroutine update_hamiltonian
    
  end subroutine tddft_real_time_propagation_run


  ! ---------------------------------------------------------
  subroutine tddft_real_time_propagation_end()
    use aims_memory_tracking, only : aims_deallocate    

    ! cleanup work arrays
    deallocate(hamiltonian_square)
    
    deallocate(hamiltonian_extrapolation_1)
    deallocate(hamiltonian_extrapolation_2)
    deallocate(hamiltonian_step)
    
    deallocate(scalapack_work_1)
    deallocate(scalapack_work_2)
    deallocate(scalapack_work_3)
    deallocate(scalapack_zwork_1)
    deallocate(scalapack_zwork_2)
    deallocate(scalapack_zwork_3)
    deallocate(hamiltonian_ev)
    
    deallocate(hamiltonian_work)
    deallocate(tdks_eigenvalues)
    deallocate(tdks_eigenvalues_phases)
    deallocate(propagator_U)
    deallocate(td_KS_eigenvector)
    deallocate(td_KS_gs_eigenvector)
    deallocate(td_KS_initial_eigenvector)
    deallocate(td_KS_eigenvector_work)
    
    deallocate(zwork_1)
    deallocate(ovlp_eigenvalues)
    deallocate(overlap_m_sqrt)
    deallocate(zoverlap_m_sqrt)
    deallocate(overlap_p_sqrt)
    deallocate(zoverlap_p_sqrt)
    deallocate(overlap_work)
    deallocate(zoverlap_matrix)
    deallocate(range_space)
    deallocate(null_space)

    call aims_deallocate(KS_eigenvector_complex, "KS_eigenvector_complex")
    
  end subroutine tddft_real_time_propagation_end
  
end module tddft_real_time_propagation
