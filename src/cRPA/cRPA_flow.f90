!Manages the control flow
!This is the only module linked explicitly to aims
!Comments always refer to the following line

!Note k_phase_exx have to be allocated before!

!Known bugs:
! if all n_k_points_x,n_k_points_y,n_k_points_z are smaller than 3 then NaN is returned
! for instance 1*1*3 or 3*1*3 works but 2*2*2 does not
MODULE cRPA_flow
    USE runtime_choices
    USE timing

    USE cRPA_view
    USE cRPA_flow_realspace_vector
    USE cRPA_storage
    USE cRPA_parallelism
    USE cRPA_calculation

    USE cRPA_flow_test
    USE cRPA_calculation_integration_test
    USE cRPA_storage_test
    USE cRPA_calculation_test
    USE cRPA_flow_adaptive_grid_test
    USE cRPA_calculation_fouriertrans_test
    USE cRPA_parallelism_test

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: calc_crpa_energy


    INTEGER :: CALC_RESPONSEFUNC_R_SPACE
    PARAMETER (CALC_RESPONSEFUNC_R_SPACE = 1)

    INTEGER :: CALC_RESPONSEFUNC_K_SPACE
    PARAMETER (CALC_RESPONSEFUNC_K_SPACE = 2)

!Indicates wheather to calculate chi_0 in realspace or k-space
    INTEGER :: CALC_RESPONSEFUNC_IN_SPACE
    PARAMETER (CALC_RESPONSEFUNC_IN_SPACE = CALC_RESPONSEFUNC_K_SPACE)

!Tau grid is relevant only for R_SPACE calculation - K_SPACE calc has no tau integration
    INTEGER :: TAU_GRID_FIXED
    PARAMETER (TAU_GRID_FIXED = 1)

    INTEGER :: TAU_GRID_ADAPTIVE
    PARAMETER (TAU_GRID_ADAPTIVE = 2)

!Indicates wheather to use a fixed or an adaptive tau integration grid
    INTEGER :: TAU_GRID_MODE
    PARAMETER (TAU_GRID_MODE = TAU_GRID_ADAPTIVE)

!If fixed integration grid is used: Number of tau points of fixed grid
    INTEGER :: NUM_TAU_MAX_POINTS
    PARAMETER (NUM_TAU_MAX_POINTS = 150)

!If fixed integration grid is used: Upper tau bound
    REAL*8 :: TAU_INTERVAL_BOUND_UPPER
    PARAMETER (TAU_INTERVAL_BOUND_UPPER=10.d0)

!Numerical quadratures order on the subintervals of the summed quadratures
    INTEGER :: INTEGRATION_ORDER
    PARAMETER(INTEGRATION_ORDER = 3)

!Debug settings
    LOGICAL :: DO_UNIT_TESTS
    PARAMETER(DO_UNIT_TESTS = .FALSE.)

!Compare flags have effect only if calculation in R space is used
    LOGICAL :: COMPARE_K_TO_R_SPACE_RESPONSEFUNCTION
    PARAMETER(COMPARE_K_TO_R_SPACE_RESPONSEFUNCTION = .FALSE.)

    LOGICAL :: COMPARE_Romega_TO_R_SPACE_RESPONSEFUNCTION
    PARAMETER(COMPARE_Romega_TO_R_SPACE_RESPONSEFUNCTION = .FALSE.)

    LOGICAL :: PURE_FT_TEST
    PARAMETER(PURE_FT_TEST= .FALSE.)

    LOGICAL :: COMPARE_R_TO_K_SPACE_RESPONSEFUNCTION
    PARAMETER(COMPARE_R_TO_K_SPACE_RESPONSEFUNCTION = .FALSE.)


!see parallelism_storage and parallelism_planing
    type(basbas_distribution),ALLOCATABLE, DIMENSION(:) :: basbas_worklist, &
                                                           basbas_worklist_total
!see parallelism_storage and parallelism_planing
    type(basisfunc_distribution),ALLOCATABLE, DIMENSION(:) :: basis_dist
!see flow_realspace_vectors
    type(realspace_vectors) :: rvecs

!saves the omegas on which chi_0 is calculated
    REAL*8,ALLOCATABLE, DIMENSION(:) :: all_frequencies
    INTEGER :: n_all_frequencies
CONTAINS
!returns rpa correlation energy - this is the function called from aims
    SUBROUTINE calc_crpa_energy(n_omega, omega_points, omega_weights, rpa_c_energy_out)
        INTEGER, INTENT(IN) :: n_omega
        REAL*8, DIMENSION(n_omega), INTENT(IN) :: omega_points
        REAL*8, DIMENSION(n_omega), INTENT(IN) :: omega_weights
        REAL*8, INTENT(OUT) :: rpa_c_energy_out

      PROCEDURE(get_resp_func_for_omega_template),POINTER :: get_resp_func_for_omega
      INTEGER :: i_loop_order

      if(DO_UNIT_TESTS) CALL do_tests()

      CALL prepare_response_func(n_omega, omega_points, omega_weights)

      get_resp_func_for_omega => get_resp_func_in_k_omega

!depending on the method to calculate chi_0 the loop order is interchanged
!if chi_0 is calculated in k-space it is more efficient to loop k outermost
!if chi_0 is calculated in R-space it is more efficient to loop omega outermost
      if(CALC_RESPONSEFUNC_IN_SPACE == CALC_RESPONSEFUNC_R_SPACE) then
         i_loop_order=LOOP_ORDER_FREQ_K
      elseif(CALC_RESPONSEFUNC_IN_SPACE == CALC_RESPONSEFUNC_K_SPACE) then
         i_loop_order=LOOP_ORDER_K_FREQ
      endif

!calculates the RPA correlation energy using method given by get_resp_func_for_omega
      CALL get_crpa_energy(get_resp_func_for_omega, i_loop_order, rpa_c_energy_out)
      CALL finalize_response_function()
    END SUBROUTINE calc_crpa_energy


    SUBROUTINE prepare_response_func(n_omega, omega_points, &
                                                          omega_weights)
        INTEGER, INTENT(IN) :: n_omega
        REAL*8, DIMENSION(n_omega), INTENT(IN) :: omega_points
        REAL*8, DIMENSION(n_omega), INTENT(IN) :: omega_weights

        type(adaptive_grid) :: tau_grid
        type(integration_grid) :: int_grid
        INTEGER :: n_blocks

        CALL write_stdout("Starting calculation response function")

        CALL create_all_realspace_vectors(n_cells_bvk, cell_index_bvk, &
                                          lattice_vector,rvecs)
        CALL write_stdout("Having realspace vectors: " // num2str(rvecs%n_vecs))

        CALL write_stdout("Creating distribution strategy for auxiliary functions")
        n_blocks = n_tasks
        CALL distribute_basbas(n_basbas, n_blocks, n_tasks, &
                               basbas_worklist, basbas_worklist_total)
        CALL print_basbas_distributions(n_tasks, basbas_worklist_total)

        CALL register_communicator_basisfunction(n_blocks, comm_basisfuncs)
        CALL print_basbas_distributions(comm_basisfuncs%n_tasks, basbas_worklist)


        ALLOCATE(basis_dist(comm_basisfuncs%n_tasks))
        CALL distribute_basisfunctions(comm_basisfuncs%n_tasks, n_basis, rvecs, basis_dist)

        CALL prepare_exp_coeffs(rvecs,n_basbas,basis_dist(comm_basisfuncs%myid+1), &
                                basbas_worklist(comm_basisfuncs%myid+1))

        CALL prepare_crpa_energy_calculation(n_omega, omega_points, omega_weights)

        if(CALC_RESPONSEFUNC_IN_SPACE == CALC_RESPONSEFUNC_R_SPACE) then
            CALL prepare_comparison_k_space_polarisabilty()

            CALL prepare_response_function_calculation(rvecs, basis_dist, &
                                                       basbas_worklist(comm_basisfuncs%myid+1), &
                                                       tau_grid, int_grid)

            CALL calculate_response_function_in_R_tau(tau_grid)

            CALL prepare_aux_exp_basis(tau_grid%n_points, tau_grid%points, int_grid)
        elseif(CALC_RESPONSEFUNC_IN_SPACE == CALC_RESPONSEFUNC_K_SPACE) then
            CALL prepare_mo_coeffs()
            CALL calculate_response_function_in_k_omega(n_omega,omega_points)
        endif


    END SUBROUTINE prepare_response_func

    SUBROUTINE get_resp_func_in_k_omega(n_basbas, i_k_point, &
                                        frequency, &
                                        data_out)
        INTEGER, INTENT(IN) :: n_basbas, i_k_point
        REAL*8, INTENT(IN) :: frequency
        DOUBLE COMPLEX, DIMENSION(n_basbas,n_basbas), INTENT(OUT) :: data_out

        if(CALC_RESPONSEFUNC_IN_SPACE == CALC_RESPONSEFUNC_R_SPACE) then
            CALL get_resp_func_in_k_omega_R(n_basbas,i_k_point,frequency, data_out)
        elseif(CALC_RESPONSEFUNC_IN_SPACE == CALC_RESPONSEFUNC_K_SPACE) then
            CALL get_resp_func_in_k_omega_k(basbas_worklist, &
                                            maxval(basbas_worklist_total(:)%n_rows), &
                                            n_all_frequencies, &
                                            all_frequencies, &
                                            i_k_point, &
                                            n_basbas, &
                                            frequency, &
                                            data_out)
        endif

    END SUBROUTINE get_resp_func_in_k_omega

    SUBROUTINE prepare_exp_coeffs(rvecs, n_basbas, my_basis_dist, my_basbas_dist)
        type(realspace_vectors), INTENT(IN) :: rvecs
        INTEGER, INTENT(IN) :: n_basbas
        type(basisfunc_distribution),INTENT(IN) :: my_basis_dist
        type(basbas_distribution),INTENT(IN) :: my_basbas_dist


        CALL write_stdout("Preparing expansion coefficient integration")

        CALL register_communicator_exp_coeffs(n_tasks, comm_exp_coeffs)
        CALL init_exp_coeff_storage(exp_coeffs_local_offsite,n_basbas,n_cells,my_basis_dist)
        CALL calc_lvl_tricoeff(exp_coeffs_local_offsite,my_basbas_dist)
        CALL condense_exp_coeffs(exp_coeffs_local_offsite, &
                                 exp_coeffs_local)
        CALL finalize_exp_coeff_storage(exp_coeffs_local_offsite)

        CALL communicate_exp_coeffs_headers(exp_coeffs_local)

        CALL set_norms(exp_coeffs_local)

    END SUBROUTINE prepare_exp_coeffs

    SUBROUTINE prepare_mo_coeffs()
        CALL register_communicator_mo_coeffs(1, comm_mo_coeffs)
        CALL init_mo_coefficient_storage(mo_coeffs_local, n_k_points, &
                                         n_basis, n_states,n_spin)
        CALL setup_communication_mo_coeffs(mo_coeffs_local, &
                                           n_k_points_task, &
                                           KS_eigenvector_complex)
    END SUBROUTINE prepare_mo_coeffs

    SUBROUTINE calculate_response_function_in_k_omega(n_omega, omega_points)
        INTEGER, INTENT(IN) :: n_omega
        REAL*8,DIMENSION(n_omega), INTENT(IN) :: omega_points


        ALLOCATE(all_frequencies(n_omega))
        all_frequencies(:) = omega_points(:)
        n_all_frequencies = n_omega

    END SUBROUTINE calculate_response_function_in_k_omega

    SUBROUTINE register_communicator_basisfunction(n_blocks, communicator_basisfuncs)
        INTEGER, INTENT(IN) :: n_blocks
        type(custom_communicator), INTENT(INOUT) :: communicator_basisfuncs

        INTEGER :: n_rank_block
        INTEGER, ALLOCATABLE, DIMENSION(:) :: rank_block

        CALL get_cpus_basbas_row_block(n_tasks, myid, &
                                       n_blocks, &
                                       n_rank_block, rank_block)

        CALL create_custom_communicator(n_rank_block, rank_block, &
                                        communicator_basisfuncs)
        DEALLOCATE(rank_block)

        CALL register_custom_communicator(communicator_basisfuncs)
!        CALL print_custom_communicator(communicator_basisfuncs)

    END SUBROUTINE register_communicator_basisfunction


    SUBROUTINE register_communicator_exp_coeffs(n_tasks, communicator_expcoeffs)
        INTEGER, INTENT(IN) :: n_tasks
        type(custom_communicator), INTENT(INOUT) :: communicator_expcoeffs

        INTEGER :: i_rank
        INTEGER, DIMENSION(n_tasks) :: ranks

        do i_rank = 0,n_tasks-1
            ranks(i_rank+1) = i_rank
        enddo

        CALL create_custom_communicator(n_tasks, ranks, &
                                        communicator_expcoeffs)

        CALL register_custom_communicator(communicator_expcoeffs)
!        CALL print_custom_communicator(communicator_expcoeffs)

    END SUBROUTINE register_communicator_exp_coeffs


    SUBROUTINE register_communicator_mo_coeffs(n_blocks, communicator_mocoeffs)
        INTEGER, INTENT(IN) :: n_blocks
        type(custom_communicator), INTENT(INOUT) :: communicator_mocoeffs

        INTEGER :: n_rank_block
        INTEGER, ALLOCATABLE, DIMENSION(:) :: rank_block

        CALL get_cpus_basbas_row_block(n_tasks, myid, &
                                       n_blocks, &
                                       n_rank_block, rank_block)

        CALL create_custom_communicator(n_rank_block, rank_block, &
                                        communicator_mocoeffs)
        DEALLOCATE(rank_block)

        CALL register_custom_communicator(communicator_mocoeffs)
!        CALL print_custom_communicator(communicator_basisfuncs)

    END SUBROUTINE register_communicator_mo_coeffs

    SUBROUTINE prepare_response_function_calculation(rvecs, &
                                                     basis_dist, worklist, &
                                                     a_grid, &
                                                     int_grid)
        type(realspace_vectors) :: rvecs
        type(basisfunc_distribution), DIMENSION(n_tasks), INTENT(IN) :: basis_dist
        type(basbas_distribution), INTENT(IN) :: worklist
        type(adaptive_grid), INTENT(OUT) :: a_grid
        type(integration_grid),INTENT(OUT) :: int_grid


        INTEGER :: n_partitions, n_rest_points
        type(integration_intervals) :: int_intervals

        CALL write_stdout("Preparing response function calculation")

        n_partitions = NUM_TAU_MAX_POINTS
        CALL factorize_integer_product(INTEGRATION_ORDER,NUM_TAU_MAX_POINTS, &
                                       n_partitions, n_rest_points)
!        CALL create_integration_intervals_logarithmic(0.d0, &
!                                                      TAU_INTERVAL_BOUND_UPPER, &
!                                                      n_partitions, &
!                                                      int_intervals)

        CALL create_integration_intervals_uniform(0.d0, &
                                                  TAU_INTERVAL_BOUND_UPPER, &
                                                  n_partitions, &
                                                  int_intervals)

!        CALL create_integration_intervals_uniform_unsymmetric(minval(a_grid%points), &
!                                                              maxval(a_grid%points), &
!                                                              n_partitions, &
!                                                              int_intervals)
!        CALL create_grid_old(minval(a_grid%points), &
!                             maxval(a_grid%points), &
!                             NUM_TAU_MAX_POINTS/2, tau_grid)


        CALL create_integration_grid_gauss_legendre(int_intervals,INTEGRATION_ORDER, &
                                                    int_grid)
!        CALL create_integration_grid_trapezodial(int_intervals, int_grid)

        CALL free_integration_intervals(int_intervals)

        CALL create_adaptive_grid(int_grid%n_points, int_grid%points,int_grid%points, &
                                  a_grid)


        CALL init_green_func_storage(green_funcs_local, rvecs,n_spin)
        CALL init_green_func_storage(green_funcs_local_minus_tau, rvecs,n_spin)

        CALL init_response_function_storage(rvecs,worklist)

        CALL write_stdout("Finished preparation response function calculation")

    END SUBROUTINE prepare_response_function_calculation

    SUBROUTINE calculate_response_function_in_R_tau(tau_grid)
        type(adaptive_grid),INTENT(IN) :: tau_grid

        REAL*8 :: cur_tau, cur_crpa_energy, resp_func_time, resp_func_time_wall

        INTEGER :: i_tau

        PROCEDURE(get_signal_template),POINTER :: calc_resp_func
        type(adaptive_grid) :: a_grid

        CALL get_timestamps(resp_func_time, resp_func_time_wall)
        CALL write_stdout("Starting response function calculation in imaginary time and realspace")

        if(TAU_GRID_MODE == TAU_GRID_ADAPTIVE) then
            calc_resp_func => calculate_and_score_tau_point

            CALL create_adaptive_grid_adaptivly(1.d-5, &
                                                1.d-2, &
                                                calc_resp_func, &
                                                a_grid)
        elseif(TAU_GRID_MODE == TAU_GRID_FIXED) then
stop
            do i_tau = 1, tau_grid%n_points
                cur_tau = tau_grid%points(i_tau)
                CALL calculate_and_score_tau_point(cur_tau,cur_crpa_energy)
            end do
        endif

        CALL write_stdout("Finished response function calculation in imaginary time and realspace")
        CALL get_times(resp_func_time, resp_func_time_wall)
        CALL output_times("2X","Overall time for resp func",resp_func_time, resp_func_time_wall)

        CALL write_stdout("Overall resp chksum:"//num2str(chksum_response_function()))


!        CALL write_response_function()

        CALL print_response_function()
    END SUBROUTINE calculate_response_function_in_R_tau

    SUBROUTINE calculate_and_score_tau_point(cur_tau,cur_crpa_energy)
       REAL*8,INTENT(IN) :: cur_tau
       REAL*8,INTENT(OUT):: cur_crpa_energy
       DOUBLE COMPLEX, DIMENSION(rvecs%n_vecs, &
                         basbas_worklist(1)%n_columns, &
                         basbas_worklist(1)%n_rows) :: result_this_tau
       INTEGER :: i,j, i_tau

       CALL write_stdout("Calculating tau point" &
                              //num2str(response_func%tau_grid%n_points))

       if(.NOT. PURE_FT_TEST) then
           CALL calculate_tau_point(rvecs, basbas_worklist(1), &
                                    cur_tau,result_this_tau(:,:,:))
       endif

       CALL add_response_function_tau_point(dble(result_this_tau(:,:,:)), cur_tau, i_tau)


       if(COMPARE_K_TO_R_SPACE_RESPONSEFUNCTION) then
            CALL compare_methods()
       endif

       CALL symmetry_greensfunc_test(green_funcs_local)
       CALL symmetry_greensfunc_test(green_funcs_local_minus_tau)

       cur_crpa_energy = chksum_response_function_for_point(i_tau) / &
                         chksum_response_function_for_point(1)

       CALL write_stdout("   Current signal:"//num2str(cur_crpa_energy))
    CONTAINS
       SUBROUTINE compare_methods()
          INTEGER :: i,j,i_row_total, i_tau, i_vec
          REAL*8 :: value_k_space, value_R_space,max_diff, &
                    value_k_space_transpose, value_R_space_transpose
          COMPLEX*16, ALLOCATABLE,DIMENSION(:,:,:) :: responsefunc_from_k_space

          ALLOCATE(responsefunc_from_k_space(n_basbas,n_basbas, rvecs%n_vecs))
          CALL get_kspace_polarisability_for_tau(cur_tau,rvecs, &
                                                 responsefunc_from_k_space)

          CALL find_point(response_func%tau_grid, cur_tau,i_tau)

          max_diff = 0.d0

          if(PURE_FT_TEST) then
              do i = 1,basbas_worklist(comm_basisfuncs%myid+1)%n_rows
                i_row_total = i + basbas_worklist(comm_basisfuncs%myid+1)%n_offset_rows
                do j = i,basbas_worklist(comm_basisfuncs%myid+1)%n_rows
                    do i_vec = 1, rvecs%n_vecs
                    value_k_space = dble(responsefunc_from_k_space(i_row_total,j,i_vec)) * 2.d0
                        !if(isnan(value_k_space)) then
                        !  value_k_space = 0.d0
                        !else 
                    if(abs(value_k_space) < 1.d-20) then
                       value_k_space = 0.d0
                    endif

                    response_func%element(i_vec,j,i,i_tau) = value_k_space
                 enddo
              enddo
           enddo
        endif

          write(use_unit,*) "No rvecs",rvecs%n_vecs

          do i_vec = 1, rvecs%n_vecs, 8
              do i = 1,basbas_worklist(comm_basisfuncs%myid+1)%n_rows
                 i_row_total = i + basbas_worklist(comm_basisfuncs%myid+1)%n_offset_rows
                 do j = i,basbas_worklist(comm_basisfuncs%myid+1)%n_rows
                      value_k_space = responsefunc_from_k_space(i_row_total,j,i_vec)
                      value_k_space_transpose = responsefunc_from_k_space(j,i_row_total,i_vec)
                      value_R_space = response_func%element(i_vec,j,i,i_tau) / 2.d0 !cf. add resp ele

                      if(i>=basbas_worklist(comm_basisfuncs%myid+1)%n_offset_rows .AND. &
                         i<basbas_worklist(comm_basisfuncs%myid+1)%n_offset_rows + &
                            basbas_worklist(comm_basisfuncs%myid+1)%n_rows) then
                        value_R_space_transpose = response_func%element(i_vec,i,j,i_tau) / 2.d0 !cf. add resp ele
                      else
                        value_R_space_transpose = 0.d0
                      endif

                      max_diff = max(max_diff, dble(value_k_space) - value_R_space)

                      if(i/=j .AND. abs(value_R_space)> 1.d-20 ) then
                        write(use_unit,*) i_vec, i,j, &
                                   responsefunc_from_k_space(i_row_total,j,i_vec), &
                                   !value_k_space, value_k_space_transpose, &
                                   value_R_space, value_R_space_transpose
!                                   dble(value_k_space/value_R_space)
                      endif
                 enddo
              enddo
          enddo


         CALL write_debug("---Max diff:"//num2str(max_diff))

         DEALLOCATE(responsefunc_from_k_space)

       END SUBROUTINE compare_methods

    END SUBROUTINE calculate_and_score_tau_point


    SUBROUTINE compare_rpa_correlation_energy_realspace(n_basbas,n_kpoints, k_space_polar,omega)
        INTEGER, INTENT(IN) :: n_basbas, n_kpoints
        COMPLEX*16, DIMENSION(n_basbas,n_basbas,n_kpoints), INTENT(IN) :: k_space_polar
        REAL*8, INTENT(IN) :: omega

        INTEGER :: i_k, i_row, i_column
        COMPLEX*16, DIMENSION(n_basbas,n_basbas) :: diff
        REAL*8 :: max_diff, max_diff_row(n_basbas), max_val_row(n_basbas)
        COMPLEX*16, DIMENSION(n_basbas,n_basbas,n_kpoints) :: data_copy
        REAL*8 :: max_rel_error, min_rel_error


        if(.NOT. COMPARE_R_TO_K_SPACE_RESPONSEFUNCTION) then
           return
        endif

        CALL get_resp_func_in_k_omega_R(n_basbas,n_kpoints,omega, data_copy)


        do i_row = 1,n_basbas
           do i_column = 1,n_basbas

            max_rel_error=0.d0
            min_rel_error=100000.d0
            do i_k = 1, n_kpoints

                !not my k point
                if(ALL(abs(data_copy(:,:,i_k))<1.d-20)) cycle

                diff(:,:) = k_space_polar(:,:,i_k) - data_copy(:,:,i_k)


    !                    write(use_unit,*) i_row,i_column, &
    !                               k_space_polar(i_row,i_column,i_k), &
    !                               data_copy(i_row,i_column,i_k), &
    !                               dble(k_space_polar(i_row,i_column,i_k)) / &
    !                               dble(data_copy(i_row,i_column,i_k))
                max_rel_error = max(max_rel_error, &
                                    abs(dble(diff(i_row,i_column)) / &
                                        dble(k_space_polar(i_row,i_column,i_k))))
                min_rel_error = min(min_rel_error, &
                                    abs(dble(diff(i_row,i_column)) / &
                                        dble(k_space_polar(i_row,i_column,i_k))))

                if(i_k==1) then
                  if(abs(k_space_polar(i_row,i_column,i_k))>=1.d-8) then
                  write(use_unit,*) i_row,i_column, &
                             (k_space_polar(i_row,i_column,i_k)), &
                             dble(data_copy(i_row,i_column,i_k)) , &
                             dble(data_copy(i_column,i_row,i_k))
                  endif
                endif

            enddo

!            if(max_rel_error > 1.d-3) then
!                write(use_unit,*) i_row,i_column, max_rel_error, min_rel_error
!            endif

!                max_val_row(i_row) = maxval(abs(k_space_polar(i_row,:,i_k)))
!                max_diff_row(i_row) = maxval(abs(diff(i_row,:)))
!                write(use_unit,*) "Max diff row ",i_row, &
!                           max_val_row(i_row), max_diff_row(i_row)
            enddo
!            write(use_unit,*) "Max diff ", maxval(max_val_row(:)), &
!                         maxval(max_diff_row(:))

        enddo

    END SUBROUTINE compare_rpa_correlation_energy_realspace

    SUBROUTINE prepare_crpa_energy_calculation(n_omega, omega_points,omega_weights)
        INTEGER,INTENT(IN) :: n_omega
        REAL*8, DIMENSION(n_omega),INTENT(IN) :: omega_points,omega_weights

        type(integration_grid) :: omega_grid
        type(integration_intervals) :: int_intervals
        CALL create_integration_intervals_uniform(0.d0, &
                                                  200.d0, &
                                                  100/INTEGRATION_ORDER, &
                                                  int_intervals)

!        CALL create_integration_grid_gauss_legendre(int_intervals,INTEGRATION_ORDER, &
!                                                    omega_grid)

        CALL create_integration_grid(n_omega,omega_points,omega_weights, &
                                     omega_grid)
        CALL free_integration_intervals(int_intervals)

        CALL set_energy_omega_grid(omega_grid)
        CALL free_integration_grid(omega_grid)

    END SUBROUTINE prepare_crpa_energy_calculation

    SUBROUTINE prepare_comparison_k_space_polarisabilty()

        if(.NOT. COMPARE_K_TO_R_SPACE_RESPONSEFUNCTION .AND. &
           .NOT. COMPARE_Romega_TO_R_SPACE_RESPONSEFUNCTION) then
            return
        endif

        CALL calculate_kspace_polarisability &
             (n_low_state, occ_numbers, &
              KS_eigenvalue, KS_eigenvector, &
              KS_eigenvector_complex, rvecs)
    END SUBROUTINE prepare_comparison_k_space_polarisabilty

    SUBROUTINE prepare_aux_exp_basis(n_points, points, int_grid)
        INTEGER, INTENT(IN) :: n_points
        REAL*8, DIMENSION(n_points), INTENT(IN) :: points
        type(integration_grid), INTENT(IN) :: int_grid

        REAL*8, ALLOCATABLE, DIMENSION(:) :: mink_diff, mean_energies
        INTEGER :: n_mink_diff, n_mean_energies,k


        CALL write_stdout("Preparing auxiliary exp basis for "//num2str(n_points))

        CALL get_possible_energy_differences(occ_numbers,KS_eigenvalue, &
                                             -1.d-1, &
                                             n_mink_diff,mink_diff)
        k = min(n_points,80)

        CALL get_at_most_n_distinct_energies(n_mink_diff,mink_diff, &
                                             -1000.d0, -1.d-8, &
                                             k, n_mean_energies, mean_energies)

        CALL create_exp_basis(n_mean_energies, mean_energies, &
                              aux_exp_basis)

        DEALLOCATE(mean_energies)

        CALL write_stdout("Created auxiliary exp basis")


!        CALL get_at_most_n_distinct_energies(n_mink_diff,mink_diff, &
!                                             -1.d-2, 1000.d0, &
!                                             k, n_mean_energies, mean_energies)
!
!        CALL create_exp_basis(n_mean_energies, mean_energies, &
!                              aux_exp_basis2)
!
!        DEALLOCATE(mean_energies)


        DEALLOCATE(mink_diff)

        CALL prepare_fouriertransform_by_exp_basis(n_points,points, aux_exp_basis, &
                                                   int_grid%weights)
        CALL write_stdout("Prepared ft by auxiliary exp basis")


        CALL print_exp_basis(aux_exp_basis)


!        CALL prepare_fouriertransform_by_exp_basis(n_points,points, aux_exp_basis2, &
!                                                   int_grid%weights)
!
!        CALL print_exp_basis(aux_exp_basis2)

    END SUBROUTINE prepare_aux_exp_basis

    SUBROUTINE finalize_response_function()
        CALL finalize_exp_coeff_storage(exp_coeffs_local)

        if(CALC_RESPONSEFUNC_IN_SPACE == CALC_RESPONSEFUNC_R_SPACE) then
            CALL finalize_green_func_storage(green_funcs_local)
        endif

        if(CALC_RESPONSEFUNC_IN_SPACE == CALC_RESPONSEFUNC_K_SPACE) then
            CALL finalize_mo_coeffs_storage(mo_coeffs_local)
            CALL finalize_calculation_responsefunc_k()
        endif

        DEALLOCATE (basbas_worklist)
        DEALLOCATE (basbas_worklist_total)
    END SUBROUTINE finalize_response_function

    SUBROUTINE do_tests()
        CALL get_realspace_vectors_index_inversion_test()
        CALL get_realspace_vectors_index_all_possible_differences_test()
        CALL discrete_fouriertransform_test()
        CALL adaptive_grid_test()
        CALL distribute_basbas_test()
        CALL fouriertransform_test()
        CALL get_basbas_row_blocks_test()
        CALL distribute_basbas_row_blocks_test()
        CALL get_n_good_mean_energies_test()
        CALL exp_basis_test()
    END SUBROUTINE do_tests

END MODULE cRPA_flow
