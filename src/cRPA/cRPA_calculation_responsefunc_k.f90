!Calculates response function chi_0 in k-space
module cRPA_calculation_responsefunc_k
    USE cRPA_view
    USE cRPA_storage
    USE cRPA_parallelism
    USE cRPA_calculation_energy
    use evaluate_polarisability_kspace_mod

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: calculate_k_point, &
              get_resp_func_in_k_omega_k, &
              finalize_calculation_responsefunc_k

    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:) :: M_kq_my, &
                                                       M_kq_my_cont, &
                                                       M_kq_requested, &
                                                       M_weight, &
                                                       M_work

    DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: M_column, &
                                                     M_row, &
                                                     M_my, &
                                                     mo_coeff_k, &
                                                     mo_coeff_q, &
                                                     ec_column_k, &
                                                     ec_column_q, &
                                                     ft_buffer_row, &
                                                     ft_buffer_column,&
                                                     moc_ec, &
                                                     moc_buffer, &
                                                     moc_buffer_prod

!buffers for basbas pair calculation
    DOUBLE COMPLEX,ALLOCATABLE, DIMENSION(:,:)  :: ft_buffer_prod, &
                                                   ft_buffer, &
                                                   ft_buffer_offsite, &
                                                   res_vec


    type(result_last_calculation) :: result_last_calc
CONTAINS
    SUBROUTINE calculate_k_point(i_k_point_start, n_points, work_dist, &
                                 max_n_rows_over_all_cpus, &
                                 n_frequencies, frequencies)
       INTEGER, INTENT(IN) :: i_k_point_start, n_points
       type(basbas_distribution), DIMENSION(comm_basisfuncs%n_tasks), INTENT(IN) :: work_dist
       INTEGER, INTENT(IN) :: max_n_rows_over_all_cpus
       INTEGER, INTENT(IN) :: n_frequencies
       REAL*8,DIMENSION(n_frequencies), INTENT(IN):: frequencies

       REAL*8 :: tau_point_time, tau_point_time_clock
       REAL*8 :: chksum_my_exp_coeffs,chksum_my_green_func

       INTEGER :: i, &
                  mpi_result

       CALL write_stdout("----------------------------------------------------")
       CALL write_stdout("Calculating k_point="//num2str(i_k_point_start) &
                                               //num2str(i_k_point_start+n_points-1))

       CALL get_timestamps(tau_point_time, tau_point_time_clock)

       CALL init_result_last_calculation(result_last_calc, &
                                         i_k_point_start, n_points , &
                                         n_frequencies, frequencies, &
                                         work_dist(comm_basisfuncs%myid+1))

       CALL allocate_response_function_buffers(mo_coeffs_local%n_basis, &
                                               exp_coeffs_local%n_max_all_cells_filled_rows, &
                                               exp_coeffs_local%n_basbas_local, &
                                               exp_coeffs_local%n_basbas_local_max, &
                                               mo_coeffs_local%n_k_points)

       CALL calculate_my_basbas_distribution(i_k_point_start,n_points, &
                                             n_frequencies, &
                                             frequencies, &
                                             max_n_rows_over_all_cpus, &
                                             work_dist(1), &
                                             result_last_calc%data)

       CALL deallocate_response_function_buffers()

       CALL MPI_BARRIER(mpi_comm_global, mpi_result)

       CALL get_times(tau_point_time, tau_point_time_clock)
       CALL output_times("2X","Time for this k point",tau_point_time, tau_point_time_clock)

!       chksum_my_exp_coeffs = chksum_exp_coeffs(exp_coeffs_local)
!       chksum_my_green_func = chksum_mo_coeffs(mo_coeffs_local)

!       CALL print_all_chksums(chksum_3d_matrix(abs(result_this_k_point(:,:,:))), &
!                              chksum_my_exp_coeffs, &
!                              chksum_my_green_func)


    END SUBROUTINE calculate_k_point

    SUBROUTINE allocate_response_function_buffers(n_basis, &
                                                  n_max_extend_exp_coeff_row,&
                                                  n_basbas_local, &
                                                  n_basbas_local_max, &
                                                  n_k_points)
        INTEGER, INTENT(IN) :: n_basis,n_max_extend_exp_coeff_row, &
                               n_basbas_local, n_basbas_local_max, n_k_points
        CHARACTER(*), PARAMETER :: funcname = 'allocate_response_function_buffers'
        INTEGER :: allocation_info

        ALLOCATE(mo_coeff_k(n_basis,n_states,n_spin), stat=allocation_info)
        CALL check_allocation(allocation_info, 'mo_coeff_k', funcname)

        ALLOCATE(mo_coeff_q(n_basis,n_states,n_spin), stat=allocation_info)
        CALL check_allocation(allocation_info, 'mo_coeff_q', funcname)

        ALLOCATE(moc_buffer(n_basis,n_states,n_spin), stat=allocation_info)
        CALL check_allocation(allocation_info, 'moc_buffer', funcname)

        ALLOCATE(M_kq_my(n_states, n_states,n_spin,n_basbas_local), stat=allocation_info)
        CALL check_allocation(allocation_info, 'M_kq_my', funcname)

        ALLOCATE(M_kq_my_cont(n_basbas_local*2,n_states, n_states,n_spin), stat=allocation_info)
        CALL check_allocation(allocation_info, 'M_kq_my_cont', funcname)


        ALLOCATE(M_kq_requested(n_states, n_states,n_spin,n_basbas_local_max), stat=allocation_info)
        CALL check_allocation(allocation_info, 'M_kq_requested', funcname)

        ALLOCATE(M_weight(n_frequencies,n_states, n_states,n_spin), stat=allocation_info)
        CALL check_allocation(allocation_info, 'M_weight', funcname)


        ALLOCATE(M_work(n_states, n_states,n_spin,n_frequencies), stat=allocation_info)
        CALL check_allocation(allocation_info, 'M_work', funcname)

        ALLOCATE(M_row(n_states, n_states,n_spin), stat=allocation_info)
        CALL check_allocation(allocation_info, 'M_row', funcname)

        ALLOCATE(M_column(n_states, n_states,n_spin), stat=allocation_info)
        CALL check_allocation(allocation_info, 'M_column', funcname)

        ALLOCATE(M_my(n_states, n_states,n_spin), stat=allocation_info)
        CALL check_allocation(allocation_info, 'M_my', funcname)

        ALLOCATE(ec_column_k(2*n_max_extend_exp_coeff_row, n_states,n_spin), stat=allocation_info)
        CALL check_allocation(allocation_info, 'ec_column_k', funcname)

        ALLOCATE(ec_column_q(2*n_max_extend_exp_coeff_row, n_states,n_spin), stat=allocation_info)
        CALL check_allocation(allocation_info, 'ec_column_q', funcname)


        ALLOCATE(moc_ec(n_basis,n_states,n_spin), stat=allocation_info)
        CALL check_allocation(allocation_info, 'moc_ec', funcname)


        ALLOCATE(ft_buffer_row(n_max_extend_exp_coeff_row, n_basis,n_k_points+1), stat=allocation_info)
        CALL check_allocation(allocation_info, 'ft_buffer_row', funcname)

        ALLOCATE(ft_buffer_column(n_max_extend_exp_coeff_row, n_basis,n_k_points+1), stat=allocation_info)
        CALL check_allocation(allocation_info, 'ft_buffer_column', funcname)

        ALLOCATE(ft_buffer(n_max_extend_exp_coeff_row, n_basis), stat=allocation_info)
        CALL check_allocation(allocation_info, 'ft_buffer', funcname)

        ALLOCATE(ft_buffer_offsite(n_max_extend_exp_coeff_row, n_basis), stat=allocation_info)
        CALL check_allocation(allocation_info, 'ft_buffer_offsite', funcname)


        ALLOCATE(ft_buffer_prod(2*n_max_extend_exp_coeff_row, n_basis), stat=allocation_info)
        CALL check_allocation(allocation_info, 'ft_buffer_prod', funcname)

        ALLOCATE(moc_buffer_prod(2*n_max_extend_exp_coeff_row, n_states,n_spin), stat=allocation_info)
        CALL check_allocation(allocation_info, 'moc_buffer_prod', funcname)

        ALLOCATE(res_vec(n_basbas_local,n_frequencies), stat=allocation_info)
        CALL check_allocation(allocation_info, 'res_vec', funcname)

    END SUBROUTINE allocate_response_function_buffers

    INTEGER FUNCTION k_point_minus(i_k_point)
        INTEGER, INTENT(IN) :: i_k_point
        INTEGER :: i_k_run
        DOUBLE COMPLEX :: scalprod_conjg

        do i_k_run=1,n_k_points
          scalprod_conjg= sum( K_phase_exx(:,i_k_point)*conjg(conjg(K_phase_exx(:,i_k_run))) )
          if(abs(scalprod_conjg-DCMPLX(n_k_points,0.d0))<1.d-5) then
             k_point_minus = i_k_run
             return
          endif
        enddo

        stop "error in k_point_minus"
    END FUNCTION

    SUBROUTINE get_kq_pair_weight(i_k_point_start, n_points, &
                                  i_k,i_q, kq_pair_weight)
       INTEGER,INTENT(IN) :: i_k_point_start, n_points,i_k, i_q
       DOUBLE COMPLEX, INTENT(OUT) :: kq_pair_weight

       INTEGER :: i_k_minus_q, i_k_run, i_q_run,i_k_run_end, i_k_minus,i_q_minus
       LOGICAL :: swapped_pair_before
       DOUBLE COMPLEX :: prod

       i_k_minus_q = kq_point_list(i_k,i_q)

       if(i_k_minus_q < i_k_point_start &
          .OR. &
          i_k_point_start + n_points -1< i_k_minus_q ) then
           kq_pair_weight=DCMPLX(0.d0,0.d0)
           return
       endif

!       kq_pair_weight=DCMPLX(1.d0,0.d0)
!       return


       i_k_minus = k_point_minus(i_k)
       i_q_minus = k_point_minus(i_q)
       swapped_pair_before = .FALSE.

       do i_q_run=1,i_q
          if(i_q_run < i_q) then
            i_k_run_end = n_k_points
          else
            i_k_run_end = i_k -1
          endif

          do i_k_run=1,i_k_run_end
             if(i_q_run == i_k_minus .AND. i_k_run == i_q_minus) then
                swapped_pair_before=.TRUE.
             endif
          enddo
       enddo

       if(swapped_pair_before) then
         kq_pair_weight=DCMPLX(0.d0,0.d0)
       else
         if(i_k_minus == i_q .AND. i_q_minus == i_k) then
            kq_pair_weight=DCMPLX(0.5d0,0.d0)
         else
             kq_pair_weight=DCMPLX(1.d0,0.d0)
         endif
       endif

    END SUBROUTINE get_kq_pair_weight

    SUBROUTINE calculate_my_basbas_distribution(i_k_point_start, n_points, &
                                                n_frequencies, &
                                                frequencies, &
                                                max_n_rows_over_all_cpus, &
                                                work_dist, &
                                                result_this_k_point)
        INTEGER, INTENT(IN) :: i_k_point_start, n_points, &
                               n_frequencies
        REAL*8, DIMENSION(n_frequencies) :: frequencies
        INTEGER, INTENT(IN) :: max_n_rows_over_all_cpus
        type(basbas_distribution),INTENT(IN) :: work_dist
        DOUBLE COMPLEX, DIMENSION(n_frequencies,n_points, n_basbas, &
                                  work_dist%n_rows), &
                                             INTENT(OUT) :: result_this_k_point


        INTEGER :: i_freq, mpi_result, i_k,i_q, &
                   i_basbas_row, i_basbas_row_my

        DOUBLE COMPLEX, DIMENSION(n_frequencies) :: frequencies_cmplx
        DOUBLE COMPLEX :: kq_pair_weight
        DOUBLE COMPLEX :: basbas_tmp(n_basbas,n_basbas)
        DOUBLE COMPLEX :: res_before(2),res_after(2)
        DOUBLE PRECISION :: time_communication,time_contraction, &
                            time_multiplication, time_weight, &
                            time_m_kq_all

        FORALL(i_freq=1:n_frequencies) frequencies_cmplx(i_freq) = DCMPLX(0.d0,frequencies(i_freq))

        result_this_k_point(:,:,:,:) = DCMPLX(0.d0,0.d0)

        time_contraction = 0.d0
        time_multiplication = 0.d0
        time_communication = 0.d0
        time_weight = 0.d0
        time_m_kq_all = 0.d0

        do i_q=1,n_k_points
            CALL communicate_mo_coeffs(mo_coeffs_local,i_q,mo_coeff_q)

            do i_k=1,n_k_points

                CALL get_kq_pair_weight(i_k_point_start, n_points, i_k,i_q,&
                                        kq_pair_weight)

                if(abs(kq_pair_weight) < 1.d-8) CYCLE

                CALL communicate_mo_coeffs(mo_coeffs_local,i_k,mo_coeff_k)

                CALL triggered_timer(time_weight)
                CALL build_M_weights(KS_eigenvalue, &
                                     occ_numbers, &
                                     i_k,i_q, &
                                     n_frequencies, &
                                     frequencies_cmplx(:), &
                                     M_weight(:,:,:,:))
                CALL triggered_timer(time_weight)

                CALL triggered_timer(time_m_kq_all)
                CALL calculate_M_kq_all_basbas_sparse(exp_coeffs_local, &
                                                      i_k,i_q,&
                                                      mo_coeff_k, &
                                                      conjg(mo_coeff_q), &
                                                      M_kq_my)
                CALL triggered_timer(time_m_kq_all)

                CALL my_evaluate_polarisability_kspace( KS_eigenvalue, &
                                                       occ_numbers, &
                                                       i_k_point_start, n_points, &
                                                       n_frequencies, frequencies_cmplx, &
                                                       work_dist, &
                                                       i_k,i_q, &
                                                       kq_pair_weight, &
                                                       M_weight, &
                                                       M_kq_my, &
                                                       time_communication, &
                                                       time_multiplication, &
                                                       time_contraction, &
                                                       result_this_k_point(:,:,:,:))
            enddo
        enddo

        do i_k=1,n_points
            do i_freq = 1,n_frequencies
                basbas_tmp(:,:) = DCMPLX(0.d0,0.d0)
                do i_basbas_row_my=1,work_dist%n_rows
                    i_basbas_row = i_basbas_row_my+work_dist%n_offset_rows
                    basbas_tmp(i_basbas_row,:) = conjg(result_this_k_point(i_freq,i_k,:,i_basbas_row_my))
                enddo

                CALL sync_vector_complex(basbas_tmp,size(basbas_tmp),mpi_comm_global)

                do i_basbas_row_my=1,work_dist%n_rows
                    i_basbas_row = i_basbas_row_my+work_dist%n_offset_rows

                    result_this_k_point(i_freq,i_k,:,i_basbas_row_my) = &
                        result_this_k_point(i_freq,i_k,:,i_basbas_row_my) + &
                        basbas_tmp(:,i_basbas_row)
                enddo
            enddo
        enddo


      CALL write_stdout("Constructing weights:"//num2str(time_weight))
      CALL write_stdout("Constructing M all:"//num2str(time_m_kq_all))
      CALL write_stdout("Communication M:"//num2str(time_communication))
      CALL write_stdout("Multiplication:"//num2str(time_multiplication))
      CALL write_stdout("Contraction:"//num2str(time_contraction))
    CONTAINS

  SUBROUTINE build_ec_column_k_q(this, &
                                 i_basbas, &
                                 i_k,i_q, &
                                 moc_k, &
                                 moc_q, &
                                 i_total_start,i_total_end, &
                                 i_total_start_off,i_total_end_off, &
                                 ec_column_k, &
                                 ec_column_q)
       type(expansion_coefficients) :: this
       INTEGER, INTENT(IN) ::   i_basbas, i_k,i_q

       DOUBLE COMPLEX,DIMENSION(n_basis,n_states,n_spin), INTENT(IN) :: moc_k, &
                                                                        moc_q

       INTEGER, INTENT(IN) ::   i_total_start,i_total_end, &
                                i_total_start_off,i_total_end_off


       DOUBLE COMPLEX,DIMENSION(2*this%n_max_all_cells_filled_rows, &
                                 n_states,n_spin), INTENT(OUT) :: ec_column_k, &
                                                                  ec_column_q

       INTEGER :: i_prod_buffer_rows, &
                  n_prod_buffer_row_dim, i_k_point

       n_prod_buffer_row_dim = 2*this%n_max_all_cells_filled_rows
       CALL decompress_offsite_sparse(this, i_basbas, ft_buffer_offsite(:,:))

       CALL ft_sparse(this, i_basbas,i_k, ft_buffer(:,:))
       CALL set_ft_buffer_prod(this, &
                               i_total_start, i_total_end,conjg(ft_buffer(:,:)), &
                               i_total_start_off, i_total_end_off,ft_buffer_offsite(:,:), &
                               i_prod_buffer_rows)

       call zgemm('N', 'N', i_prod_buffer_rows, n_states*n_spin, n_basis, DCMPLX(1.d0,0.d0), &
                  ft_buffer_prod(:,:), n_prod_buffer_row_dim, &
                  conjg(moc_k(:,:,:)), n_basis, DCMPLX(0.d0,0.d0), &
                  ec_column_k, n_prod_buffer_row_dim)


       CALL ft_sparse(this, i_basbas,i_q, ft_buffer(:,:))
       CALL set_ft_buffer_prod(this, &
                               i_total_start, i_total_end,conjg(ft_buffer(:,:)), &
                               i_total_start_off, i_total_end_off,ft_buffer_offsite(:,:), &
                               i_prod_buffer_rows)

       call zgemm('N', 'N', i_prod_buffer_rows, n_states*n_spin, n_basis, DCMPLX(1.d0,0.d0), &
                  ft_buffer_prod(:,:), n_prod_buffer_row_dim, &
                  conjg(moc_q(:,:,:)), n_basis, DCMPLX(0.d0,0.d0), &
                  ec_column_q, n_prod_buffer_row_dim)

       END SUBROUTINE

    SUBROUTINE calculate_M_kq_all_basbas_sparse(this,i_k,i_q,moc_k,moc_q_conjg, M_kq)
       type(expansion_coefficients) :: this
       INTEGER, INTENT(IN) :: i_k,i_q
       DOUBLE COMPLEX,DIMENSION(n_basis,n_states,n_spin), INTENT(IN) :: moc_k,moc_q_conjg

       DOUBLE COMPLEX, DIMENSION(n_states,n_states,n_spin,this%n_basbas_local), INTENT(OUT) :: M_kq
       INTEGER ::i_spin, &
                 i_total_start,i_total_end, &
                 i_total_start_off,i_total_end_off
       INTEGER :: i_basbas, &
                  n_prod_buffer_row_dim, &
                  i_moc_buffer_prod_rows, &
                  i_id_iterations


      do i_id_iterations = 1, exp_coeffs_local%n_basbas_local

        i_basbas = this%basbas_local_by_id(comm_exp_coeffs%myid, &
                                           i_id_iterations)
        if(i_basbas == EXP_COEFFICIENT_UNASSIGNED) CYCLE

        i_total_start = this%coeff(i_basbas)%all_cells_first_filled_row_onsite
        i_total_end = this%coeff(i_basbas)%all_cells_last_filled_row_onsite

        i_total_start_off = this%coeff(i_basbas)%all_cells_first_filled_column_offsite
        i_total_end_off = this%coeff(i_basbas)%all_cells_last_filled_column_offsite


        n_prod_buffer_row_dim = 2*this%n_max_all_cells_filled_rows

        CALL build_ec_column_k_q(this, &
                                 i_basbas, &
                                 i_k,i_q, &
                                 moc_k,conjg(moc_q_conjg),&
                                 i_total_start,i_total_end, &
                                 i_total_start_off,i_total_end_off, &
                                 ec_column_k, &
                                 ec_column_q)

        CALL set_moc_buffer_prod(moc_k, &
                                 i_total_start, i_total_end, &
                                 i_total_start_off, i_total_end_off, &
                                 i_moc_buffer_prod_rows)

        do i_spin = 1, n_spin
                call zgemm('T', 'N', n_states, n_states, i_moc_buffer_prod_rows, DCMPLX(1.d0,0.d0), &
                           moc_buffer_prod(:,:,i_spin), n_prod_buffer_row_dim, &
                           ec_column_q(:,:,i_spin), n_prod_buffer_row_dim, &
                           DCMPLX(0.d0,0.d0), &
                           M_kq(:,:,i_spin,i_id_iterations), n_states)
        enddo

        CALL set_moc_buffer_prod(moc_q_conjg, &
                                 i_total_start, i_total_end, &
                                 i_total_start_off, i_total_end_off, &
                                 i_moc_buffer_prod_rows)

        do i_spin = 1, n_spin
                call zgemm('C', 'N', n_states, n_states, i_moc_buffer_prod_rows, DCMPLX(1.d0,0.d0), &
                           ec_column_k(:,:,i_spin), n_prod_buffer_row_dim, &
                           moc_buffer_prod(:,:,i_spin), n_prod_buffer_row_dim, &
                           DCMPLX(1.d0,0.d0), &
                           M_kq(:,:,i_spin,i_id_iterations), n_states)
        enddo
     enddo

    END SUBROUTINE calculate_M_kq_all_basbas_sparse

    subroutine my_evaluate_polarisability_kspace &
    ( KS_eigenvalue, &
      occ_numbers, &
      i_k_point_start, n_points, &
      n_frequencies, frequencies_cmplx, &
      work_dist, &
      i_k,i_q,kq_pair_weight, &
      M_weight, &
      M_kq_my, &
      time_communication, &
      time_multiplication, &
      time_contraction, &
      polar_kspace)


      REAL*8, dimension(n_states,n_spin,n_k_points),INTENT(IN) :: KS_eigenvalue
      REAL*8, dimension(n_states,n_spin,n_k_points),INTENT(IN) :: occ_numbers
      INTEGER, INTENT(IN) :: i_k_point_start, n_points, n_frequencies
      DOUBLE COMPLEX, DIMENSION(n_frequencies), INTENT(IN) :: frequencies_cmplx
      type(basbas_distribution),INTENT(IN) :: work_dist
      INTEGER,INTENT(IN) :: i_k,i_q
      DOUBLE COMPLEX,INTENT(IN) :: kq_pair_weight
      DOUBLE COMPLEX,DIMENSION(n_frequencies,n_states, n_states,n_spin),INTENT(IN) :: M_weight
      DOUBLE COMPLEX,DIMENSION(n_states, n_states,n_spin,exp_coeffs_local%n_basbas_local),INTENT(IN) :: M_kq_my
      DOUBLE PRECISION, INTENT(INOUT) :: time_communication,time_multiplication, time_contraction
      DOUBLE COMPLEX,INTENT(INOUT) :: polar_kspace(n_frequencies,n_points,n_basbas,work_dist%n_rows)



      integer :: info, mpierr
      character(*), parameter :: func = 'evaluate_polarisability_kspace.f90'
      character*150 :: info_str

      integer i_spin
      INTEGER :: i_k_minus_q
      INTEGER :: i_basbas_column,i_basbas_column_next, &
                 i_basbas_row, i_basbas_row_my, i_basbas_my, &
                 i_freq, i_store, i_id, i_state1, i_state2,  &
                 i_id_iterations,i_id_iterations_row, &
                 h_transaction, i, &
                 n_M_kq_requested
!      DOUBLE COMPLEX :: basbas_tmp(n_basbas,n_basbas,n_frequencies), result_tmp(n_frequencies)
!      basbas_tmp(:,:,:) = (0.d0,0.d0)

      i_k_minus_q = kq_point_list(i_k,i_q)
      i_store = i_k_minus_q - i_k_point_start + 1


      do i_spin = 1,n_spin
         do i_state2=1,n_states
            do i_state1=1,n_states
               do i_basbas_my = exp_coeffs_local%n_basbas_local+1, &
                                2*exp_coeffs_local%n_basbas_local
                  M_kq_my_cont(i_basbas_my,i_state1,i_state2,i_spin)= &
                                  CONJG(M_kq_my(i_state1,i_state2,i_spin,i_basbas_my-exp_coeffs_local%n_basbas_local))
               enddo
            enddo
          enddo
      enddo


      do i_id = 0, comm_exp_coeffs%n_tasks-1
         n_M_kq_requested = COUNT(exp_coeffs_local%basbas_to_id == i_id)

         if(i_id == comm_exp_coeffs%myid) then
             M_kq_requested(:,:,:,1:n_M_kq_requested) = M_kq_my(:,:,:,:)
         endif

         CALL triggered_timer(time_communication)

         CALL communicate_M_basbas_multiple(exp_coeffs_local, &
                                            i_id, &
                                            n_M_kq_requested, &
                                            M_kq_requested)

         CALL triggered_timer(time_communication)

         do i_id_iterations=1,n_M_kq_requested

             i_basbas_column = exp_coeffs_local%basbas_local_by_id(i_id, &
                                                                   i_id_iterations)
!             i_basbas_column_next = get_i_basbas_column_next(exp_coeffs_local, &
!                                                             i_id,i_id_iterations)

             if(i_basbas_column == EXP_COEFFICIENT_UNASSIGNED) CYCLE




             CALL triggered_timer(time_multiplication)

!             do i_freq = 1,n_frequencies
!                do i_spin = 1,n_spin
!                    do i_state2=1,n_states
!                        do i_state1=1,n_states
!                            if(abs(M_tmp(i_freq,i_state1,i_state2,i_spin) - &
!                                     M_kq_requested(i_state1,i_state2,i_spin,i_id_iterations)* &
!                                     M_weight(i_state1,i_state2,i_spin,i_freq)) > 1.d-5 ) &
!                                     write(use_unit,*) i_freq,i_state1,i_state2,i_spin, &
!                                                M_tmp(i_freq,i_state1,i_state2,i_spin), &
!                                                M_kq_requested(i_state1,i_state2,i_spin,i_id_iterations)* &
!                                                M_weight(i_state1,i_state2,i_spin,i_freq)
!                      enddo
!                    enddo
!                enddo
!             enddo

!             do i_freq = 1,n_frequencies
!                    M_tmp(i_freq,:,:,:) = M_tmp(i_freq,:,:,:)-&
!                             M_kq_requested(:,:,:,i_id_iterations)*M_weight(:,:,:,i_freq)
!             enddo


             do i_spin = 1,n_spin
                do i_state2=1,n_states
                    do i_state1=1,n_states
                      do i_basbas_my = 1,exp_coeffs_local%n_basbas_local
                        M_kq_my_cont(i_basbas_my,i_state1,i_state2,i_spin) = &
                           M_kq_requested(i_state1, i_state2, &
                                          i_spin,i_id_iterations) * &
                           M_kq_my_cont( &
                              exp_coeffs_local%n_basbas_local+i_basbas_my, &
                              i_state1, i_state2, i_spin)
                      enddo
                    enddo
                enddo
             enddo


             CALL triggered_timer(time_multiplication)

             CALL triggered_timer(time_contraction)

             CALL zgemm('N','T', &
                        exp_coeffs_local%n_basbas_local, &
                        n_frequencies, &
                        n_states**2 * n_spin, &
                        DCMPLX(1.d0,0.d0), &
                        M_kq_my_cont, exp_coeffs_local%n_basbas_local*2, &
                        M_weight, n_frequencies,&!n_states**2 * n_spin, &
                        DCMPLX(0.d0,0.d0), &
                        res_vec,exp_coeffs_local%n_basbas_local)

!             res_vec(:,:) = DCMPLX(0.d0,0.d0)
!             do i_spin = 1,n_spin
!                 do i_state2=1,n_states
!                    do i_state1=1,n_states
!                        do i_freq = 1,n_frequencies
!                          do i_basbas_my = 1,exp_coeffs_local%n_basbas_local
!                            res_vec(i_basbas_my, i_freq) = res_vec(i_basbas_my,i_freq) + &
!                                              M_kq_my_cont(i_basbas_my,i_state1,i_state2,i_spin)*&
!                                              M_kq_requested(i_state1,i_state2,i_spin,i_id_iterations)* &
!                                              M_weight(i_freq,i_state1,i_state2,i_spin)
!                          enddo
!                        enddo
!                    enddo
!                 enddo
!             enddo

!res_vec_real(:,:) = dble(res_vec(:,:))
!if( abs(dble(res_vec(1,1)) - res_vec_real(1,1))>1.d-4)  write(use_unit,*) res_vec(1,1),res_vec_real(1,1)

             CALL triggered_timer(time_contraction)


!              do i=1,exp_coeffs_local%n_basbas_local
!                 M_row(:,:,:)=M_kq_my(:,:,:,i)
!                 M_my(:,:,:)=M_column(:,:,:)
!                 CALL contract_M_old(KS_eigenvalue, &
!                                     occ_numbers, &
!                                     n_frequencies, &
!                                     frequencies_cmplx, &
!                                     i_k,i_q, &
!                                     M_my, M_row, &
!                                     result_tmp)
!                   do i_freq = 1,n_frequencies
!                     res_vec(i,i_freq) = result_tmp(i_freq)
!                   enddo
!              enddo

           do i_freq = 1,n_frequencies

                do i_id_iterations_row = 1, exp_coeffs_local%n_basbas_local_max

                     i_basbas_row = exp_coeffs_local%basbas_local_by_id(comm_exp_coeffs%myid, &
                                                                        i_id_iterations_row)

                     if(i_basbas_row == EXP_COEFFICIENT_UNASSIGNED) CYCLE

                     i_basbas_row_my = i_basbas_row - work_dist%n_offset_rows
                     polar_kspace(i_freq,i_store,i_basbas_column,i_basbas_row_my) = &
                     polar_kspace(i_freq,i_store,i_basbas_column,i_basbas_row_my) + &
                     kq_pair_weight*res_vec(i_id_iterations_row,i_freq)
                     !(DCMPLX(res_vec_real(i_id_iterations_row,i_freq*2-1),0.d0) + &
                     ! DCMPLX(res_vec_real(i_id_iterations_row,i_freq*2),0.d0))


!                     if(abs(kq_pair_weight-DCMPLX(2.d0,0.d0))<1.d-5) then
!                         basbas_tmp(i_basbas_row,i_basbas_column,i_freq)= &
!                                basbas_tmp(i_basbas_row,i_basbas_column,i_freq) &
!                                +conjg(res_vec(i_id_iterations_row,i_freq))
!                     endif

                enddo
           enddo
       enddo
    enddo

!    if(abs(kq_pair_weight-DCMPLX(2.d0,0.d0))<1.d-5) then
!        CALL sync_vector_complex(basbas_tmp,size(basbas_tmp),mpi_comm_global)
!
!        do i_freq = 1,n_frequencies
!            do i_basbas_row_my=1,work_dist%n_rows
!              polar_kspace(i_freq,i_store,:,i_basbas_row_my) = &
!                     polar_kspace(i_freq,i_store,:,i_basbas_row_my) + &
!                     basbas_tmp(:,i_basbas_row_my+work_dist%n_offset_rows,i_freq)
!            enddo
!        enddo
!    endif

    end subroutine my_evaluate_polarisability_kspace
    END SUBROUTINE

    SUBROUTINE build_M_weights(KS_eigenvalue, &
                              occ_numbers, &
                              i_k,i_q, &
                              n_frequencies, &
                              frequency_cmplx, &
                              M_weights)
      REAL*8, dimension(n_states,n_spin,n_k_points),INTENT(IN) :: KS_eigenvalue
      REAL*8, dimension(n_states,n_spin,n_k_points),INTENT(IN) :: occ_numbers
      INTEGER, INTENT(IN) :: i_k,i_q,n_frequencies
      DOUBLE COMPLEX, DIMENSION(n_frequencies), INTENT(IN) :: frequency_cmplx
      DOUBLE COMPLEX, DIMENSION(n_frequencies,n_states,n_states,n_spin), INTENT(OUT) :: M_weights

      REAL*8 :: e_diff1, e_diff2,&
                occ_diff1, occ_diff2,&
                cache_energy, cache_occ_number

      INTEGER :: i_spin, i_state,i_state_1, i_freq, mpi_result
            do i_spin = 1, n_spin, 1
              do i_state_1 = 1,n_states
                cache_energy = KS_eigenvalue(i_state_1,i_spin,i_q)
                cache_occ_number = occ_numbers(i_state_1,i_spin,i_q)
                do i_state = 1,n_states

                    e_diff1 = KS_eigenvalue(i_state,i_spin,i_k) &
                              - cache_energy
                    occ_diff1 = occ_numbers(i_state,i_spin,i_k) &
                                - cache_occ_number

                  do i_freq = 1,n_frequencies
                     M_weights(i_freq,i_state,i_state_1,i_spin) = &
                           occ_diff1 * DCMPLX(k_weights(i_k),0.d0) &
                            / &
                          (e_diff1 - frequency_cmplx(i_freq))
                  enddo
                enddo
               enddo
        enddo
    END SUBROUTINE build_M_weights

subroutine contract_M_old &
( KS_eigenvalue, &
  occ_numbers, &
  n_frequencies, &
  frequencies, &
  i_k,i_q, &
  M_column, M_row, &
  result_contract)


  real*8, dimension(n_states,n_spin,n_k_points) :: KS_eigenvalue
  real*8, dimension(n_states,n_spin,n_k_points) :: occ_numbers
  INTEGER :: n_frequencies
  DOUBLE COMPLEX, DIMENSION(n_frequencies) :: frequencies
  DOUBLE COMPLEX, DIMENSION(n_states,n_states,n_spin) :: M_column,M_row
  DOUBLE COMPLEX :: result_contract(n_frequencies)

  integer, allocatable :: n_homo_k(:,:) ! the HOMO level at a given k-point q
  REAL*8 :: e_diff1, e_diff2,&
            occ_diff1, occ_diff2,&
            cache_energy, cache_occ_number


  integer :: max_n_homo
  integer i_state
  integer i_state_1,i_k,i_q, a_k_point
  integer i_spin,info
  INTEGER :: i_basbas_column,i_basbas_my, &
             i_freq, i_store, i_id, i_id_iterations

  allocate(n_homo_k(n_k_points,n_spin),stat=info)
  call check_allocation(info, 'n_homo_k', 'contract old')

          n_homo_k(:,:) = 0
          do i_spin = 1, n_spin, 1
            do a_k_point = 1, n_k_points, 1
              do i_state = 1, n_states
                if(occ_numbers(i_state,i_spin,a_k_point) .gt. 1.e-12) then
                 n_homo_k(a_k_point,i_spin) = i_state
                endif
              enddo
            enddo
          enddo
          n_homo(1)=maxval(n_homo_k(:,1), n_k_points)
          n_homo(n_spin)=maxval(n_homo_k(:,n_spin), n_k_points)
          max_n_homo = max(n_homo(1), n_homo(n_spin))

          result_contract=DCMPLX(0.d0,0.d0)

          M_column(:,:,:) = CONJG(M_row(:,:,:)) * &
                                  M_column(:,:,:) * &
                                  DCMPLX(k_weights(1),0.d0)

          do i_spin = 1, n_spin, 1
            do i_state_1 = n_homo(i_spin)+1, n_states, 1
            cache_energy = KS_eigenvalue(i_state_1,i_spin,i_q)
            cache_occ_number = occ_numbers(i_state_1,i_spin,i_q)
              do i_state = 1, n_homo(i_spin), 1

                e_diff1 = KS_eigenvalue(i_state,i_spin,i_k) &
                          - cache_energy
                occ_diff1 = occ_numbers(i_state,i_spin,i_k) &
                            - cache_occ_number


                     result_contract(:) = &
                       result_contract(:) + &
                       occ_diff1  *  &
                      ( &
                        M_column(i_state,i_state_1,i_spin))  &
                        / &! * M_local(i_state,i_state_1,i_spin) / &
                       (e_diff1 - frequencies(:))!* k_weights(i_k)
              enddo
            enddo

            do i_state = 1, n_homo(i_spin), 1
              cache_energy = KS_eigenvalue(i_state,i_spin,i_q)
              cache_occ_number = occ_numbers(i_state,i_spin,i_q)
              do i_state_1 = n_homo(i_spin)+1, n_states, 1
                e_diff2 = KS_eigenvalue(i_state_1,i_spin,i_k) &
                          - cache_energy


                 occ_diff2 = occ_numbers(i_state_1,i_spin,i_k) &
                             - cache_occ_number

                 result_contract(:) = &
                          result_contract(:) + &
                           occ_diff2 * &
                             ( &
                               M_column(i_state_1,i_state,i_spin))  &
                               / &! *M_local(i_state_1,i_state,i_spin) / &
                            (e_diff2 - frequencies(:))!* k_weights(i_k)
              enddo
            enddo
          enddo

  deallocate(n_homo_k)

end subroutine contract_M_old

    SUBROUTINE map_total_pos_to_rel_pos(i_total_start, i_total_end, &
                                        i_total_start_off, i_total_end_off, &
                                        i_rel_start, i_rel_end, &
                                        i_rel_start_off, i_rel_end_off)
        INTEGER, INTENT(IN) :: i_total_start, i_total_end, &
                               i_total_start_off, i_total_end_off
        INTEGER, INTENT(OUT) :: i_rel_start, i_rel_end, &
                                i_rel_start_off, i_rel_end_off

        INTEGER :: abs_start, abs_end

        if(i_total_start_off/=0) then
            abs_start = min(i_total_start, i_total_start_off)
            abs_end = max(i_total_end, i_total_end_off)

            i_rel_start = i_total_start - abs_start + 1
            i_rel_start_off = i_total_start_off - abs_start + 1

            i_rel_end = i_rel_start + (i_total_end - i_total_start)
            i_rel_end_off = i_rel_start_off + (i_total_end_off - i_total_start_off)
        else
           i_rel_start = 1
           i_rel_end = i_rel_start + (i_total_end - i_total_start)

           i_rel_start_off = 0
           i_rel_end_off = 0
        endif
    END SUBROUTINE map_total_pos_to_rel_pos


    SUBROUTINE set_ft_buffer_prod(this, &
                                  i_total_start, i_total_end,ec_ft, &
                                  i_total_start_off, i_total_end_off,ec_ft_offsite, &
                                  i_prod_buffer_rows)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: i_total_start, i_total_end
        DOUBLE COMPLEX,DIMENSION(this%n_max_all_cells_filled_rows,n_basis), INTENT(IN) :: ec_ft
        INTEGER, INTENT(IN) :: i_total_start_off, i_total_end_off
        DOUBLE COMPLEX,DIMENSION(this%n_max_all_cells_filled_rows,n_basis), INTENT(IN) :: ec_ft_offsite
        INTEGER, INTENT(OUT) :: i_prod_buffer_rows

        INTEGER :: n_row_filled, &
                   n_row_filled_off, &
                   i_rel_start, i_rel_end, &
                   i_rel_start_off, i_rel_end_off


        ft_buffer_prod(:,:) = DCMPLX(0.d0,0.d0)
        CALL map_total_pos_to_rel_pos(i_total_start, i_total_end, &
                                      i_total_start_off, i_total_end_off, &
                                      i_rel_start, i_rel_end, &
                                      i_rel_start_off, i_rel_end_off)

        n_row_filled = i_total_end - i_total_start + 1
        ft_buffer_prod(i_rel_start:i_rel_end,:)= ec_ft(1:n_row_filled,:)

        if(i_rel_start_off>0) then
           n_row_filled_off = i_total_end_off - i_total_start_off + 1
           ft_buffer_prod(i_rel_start_off:i_rel_end_off,:)= ft_buffer_prod(i_rel_start_off:i_rel_end_off,:) &
                                                    +ec_ft_offsite(1:n_row_filled_off,:)
        endif

        i_prod_buffer_rows = max(i_rel_end, i_rel_end_off)

!write(use_unit,*) "map ft", i_total_start, i_total_end, &
!           i_total_start_off, i_total_end_off, &
!           i_rel_start, i_rel_end, &
!           i_rel_start_off, i_rel_end_off, i_prod_buffer_rows

    END SUBROUTINE set_ft_buffer_prod

    SUBROUTINE set_moc_buffer_prod(moc, &
                                   i_total_start1, i_total_end1, &
                                   i_total_start2, i_total_end2, &
                                   i_prod_buffer_rows)
       DOUBLE COMPLEX,DIMENSION(n_basis,n_states,n_spin), INTENT(IN) :: moc
       INTEGER, INTENT(IN) :: i_total_start1, i_total_end1, &
                              i_total_start2, i_total_end2

       INTEGER, INTENT(OUT) :: i_prod_buffer_rows

       INTEGER :: i_rel_start, i_rel_end, &
                  i_rel_start_off, i_rel_end_off, &
                  n_diff


!       moc_buffer_prod(:,:,:) = DCMPLX(0.d0,0.d0)
       CALL map_total_pos_to_rel_pos(i_total_start1, i_total_end1, &
                                     i_total_start2, i_total_end2, &
                                     i_rel_start, i_rel_end, &
                                     i_rel_start_off, i_rel_end_off)

       moc_buffer_prod(i_rel_start:i_rel_end,:,:) = moc(i_total_start1:i_total_end1,:,:)

!write(use_unit,*) "map moc", myid,i_total_start1, i_total_end1, &
!           i_total_start2, i_total_end2, &
!           i_rel_start, i_rel_end, &
!           i_rel_start_off, i_rel_end_off

       if(i_rel_start_off>0) then
           if(i_rel_start_off<i_rel_start) then
write(use_unit,*) "map moc less", myid,i_total_start1, i_total_end1, &
           i_total_start2, i_total_end2, &
           i_rel_start, i_rel_end, &
           i_rel_start_off, i_rel_end_off
               n_diff = i_rel_start_off - i_rel_start
               moc_buffer_prod(i_rel_start_off:i_rel_start-1,:,:) = &
                                       moc(i_total_start2:i_total_start2+n_diff-1,:,:)
           endif

           if(i_rel_end_off>i_rel_end) then
write(use_unit,*) "map moc greater", myid,i_total_start1, i_total_end1, &
           i_total_start2, i_total_end2, &
           i_rel_start, i_rel_end, &
           i_rel_start_off, i_rel_end_off
               n_diff = i_rel_end_off - i_rel_end
               moc_buffer_prod(i_rel_end+1:i_rel_end_off,:,:) = &
                                       moc(i_total_end2-n_diff+1:i_total_end2,:,:)

           endif
       endif

       i_prod_buffer_rows = max(i_rel_end,i_rel_end_off)

    END SUBROUTINE set_moc_buffer_prod


    SUBROUTINE ft_sparse(this,i_basbas,i_k_point, ft_out)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: i_basbas,i_k_point
        DOUBLE COMPLEX,INTENT(OUT):: ft_out(this%n_max_all_cells_filled_rows,n_basis)

        INTEGER :: i,j, &
                   i_accept, &
                   i_drop, &
                   i_pos, &
                   i_total_pos, &
                   i_dic_entry, &
                   i_run, i_cell,i_sparse_pos, &
                   i_first_filled_row

        DOUBLE COMPLEX :: weight_k

        ft_out(:,:) = (0.d0,0.d0)

        i_first_filled_row = this%coeff(i_basbas)%all_cells_first_filled_row_onsite
            do i_cell = 1, this%n_cells
                weight_k = k_phase_exx(i_cell,i_k_point)
                i_pos = this%coeff(i_basbas)%cell_start(i_cell)

                i = 0
                j = 1
                do i_dic_entry = this%coeff(i_basbas)%dictionary_start(i_cell), &
                                 this%coeff(i_basbas)%dictionary_end(i_cell), 2
                    i_accept = this%coeff(i_basbas)%dictionary(i_dic_entry)
                    i_drop = this%coeff(i_basbas)%dictionary(i_dic_entry+1)

                    do i_run = 1,i_accept
                        i = i+1
                        if(i>n_basis) then
                             i = 1
                             j=j+1
                        endif
    !c_compressed is real!

                        i_sparse_pos = i-i_first_filled_row+1

                        ft_out(i_sparse_pos,j) = ft_out(i_sparse_pos,j) + &
                                                    this%coeff(i_basbas)%c_compressed(i_pos) &
                                                    * weight_k


                        i_pos = i_pos + 1
                    enddo

                    i_total_pos = (j-1)*n_basis + i
                    i_total_pos = i_total_pos + i_drop

                    j = 1 + (i_total_pos-1) / n_basis
                    i = 1 + mod(i_total_pos-1,n_basis)
               enddo
            enddo


    END SUBROUTINE ft_sparse

    SUBROUTINE decompress_offsite_sparse(this,i_basbas, offsite)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: i_basbas
        DOUBLE COMPLEX,INTENT(OUT):: offsite(this%n_max_all_cells_filled_rows,n_basis)

        INTEGER :: i,j, &
                   i_accept, &
                   i_drop, &
                   i_pos, &
                   i_total_pos, &
                   i_dic_entry, &
                   i_run, i_cell,i_sparse_pos, &
                   i_first_filled_column

        offsite(:,:) = DCMPLX(0.d0,0.d0)
        i_first_filled_column = this%coeff(i_basbas)%all_cells_first_filled_column_offsite
        i_cell = this%n_cells+1

        i_pos = this%coeff(i_basbas)%cell_start(i_cell)

        i = 0
        j = 1
        do i_dic_entry = this%coeff(i_basbas)%dictionary_start(i_cell), &
                             this%coeff(i_basbas)%dictionary_end(i_cell), 2
                i_accept = this%coeff(i_basbas)%dictionary(i_dic_entry)
                i_drop = this%coeff(i_basbas)%dictionary(i_dic_entry+1)

                do i_run = 1,i_accept
                    i = i+1
                    if(i>n_basis) then
                         i = 1
                         j=j+1
                    endif
!c_compressed is real!

                    i_sparse_pos = j-i_first_filled_column+1

                    offsite(i_sparse_pos,i) = this%coeff(i_basbas)%c_compressed(i_pos)

                    i_pos = i_pos + 1
                enddo

                i_total_pos = (j-1)*n_basis + i
                i_total_pos = i_total_pos + i_drop

                j = 1 + (i_total_pos-1) / n_basis
                i = 1 + mod(i_total_pos-1,n_basis)
           enddo

    END SUBROUTINE decompress_offsite_sparse

     SUBROUTINE ft_and_store_sparse(this,i_basbas, ft_out)
        type(expansion_coefficients) :: this
        INTEGER, INTENT(IN) :: i_basbas
        DOUBLE COMPLEX,INTENT(OUT):: ft_out(this%n_max_all_cells_filled_rows,n_basis, n_k_points+1)

        INTEGER :: i_k_point
        INTEGER :: i,j, &
                   i_accept, &
                   i_drop, &
                   i_pos, &
                   i_total_pos, &
                   i_dic_entry, &
                   i_run, i_cell,i_sparse_pos, &
                   i_first_filled_row, &
                   i_first_filled_column

        DOUBLE COMPLEX :: weight_k

 !       i_basbas_row_save = i_basbas
        ft_out(:,:,:) = (0.d0,0.d0)

        i_first_filled_row = this%coeff(i_basbas)%all_cells_first_filled_row_onsite
        do i_k_point=1,n_k_points
            do i_cell = 1, this%n_cells
                weight_k = k_phase_exx(i_cell,i_k_point)
                i_pos = this%coeff(i_basbas)%cell_start(i_cell)

                i = 0
                j = 1
                do i_dic_entry = this%coeff(i_basbas)%dictionary_start(i_cell), &
                                 this%coeff(i_basbas)%dictionary_end(i_cell), 2
                    i_accept = this%coeff(i_basbas)%dictionary(i_dic_entry)
                    i_drop = this%coeff(i_basbas)%dictionary(i_dic_entry+1)

                    do i_run = 1,i_accept
                        i = i+1
                        if(i>n_basis) then
                             i = 1
                             j=j+1
                        endif
    !c_compressed is real!

                        i_sparse_pos = i-i_first_filled_row+1

                        ft_out(i_sparse_pos,j,i_k_point) = ft_out(i_sparse_pos,j,i_k_point) + &
                                                    this%coeff(i_basbas)%c_compressed(i_pos) &
                                                    * weight_k


                        i_pos = i_pos + 1
                    enddo

                    i_total_pos = (j-1)*n_basis + i
                    i_total_pos = i_total_pos + i_drop

                    j = 1 + (i_total_pos-1) / n_basis
                    i = 1 + mod(i_total_pos-1,n_basis)
               enddo
            enddo
        enddo


        i_first_filled_column = this%coeff(i_basbas)%all_cells_first_filled_column_offsite

        do i_cell = this%n_cells+1, this%n_cells+1

            i_pos = this%coeff(i_basbas)%cell_start(i_cell)

            i = 0
            j = 1
            do i_dic_entry = this%coeff(i_basbas)%dictionary_start(i_cell), &
                             this%coeff(i_basbas)%dictionary_end(i_cell), 2
                i_accept = this%coeff(i_basbas)%dictionary(i_dic_entry)
                i_drop = this%coeff(i_basbas)%dictionary(i_dic_entry+1)

                do i_run = 1,i_accept
                    i = i+1
                    if(i>n_basis) then
                         i = 1
                         j=j+1
                    endif
!c_compressed is real!

                    i_sparse_pos = j-i_first_filled_column+1

                    ft_out(i_sparse_pos,i,n_k_points+1) = this%coeff(i_basbas)%c_compressed(i_pos)

                    i_pos = i_pos + 1
                enddo

                i_total_pos = (j-1)*n_basis + i
                i_total_pos = i_total_pos + i_drop

                j = 1 + (i_total_pos-1) / n_basis
                i = 1 + mod(i_total_pos-1,n_basis)
           enddo
        enddo

    END SUBROUTINE ft_and_store_sparse

    SUBROUTINE deallocate_response_function_buffers()

        DEALLOCATE(mo_coeff_k)
        DEALLOCATE(mo_coeff_q)
        DEALLOCATE(moc_buffer)

        DEALLOCATE(M_kq_my)
        DEALLOCATE(M_kq_my_cont)

        DEALLOCATE(M_kq_requested)
        DEALLOCATE(M_weight)
        DEALLOCATE(M_work)
        DEALLOCATE(M_row)
        DEALLOCATE(M_column)
        DEALLOCATE(M_my)

        DEALLOCATE(ec_column_k)
        DEALLOCATE(ec_column_q)

        DEALLOCATE(moc_ec)

        DEALLOCATE(ft_buffer_row)
        DEALLOCATE(ft_buffer_column)
        DEALLOCATE(ft_buffer)
        DEALLOCATE(ft_buffer_offsite)

        DEALLOCATE(ft_buffer_prod)
        DEALLOCATE(moc_buffer_prod)

        DEALLOCATE(res_vec)
    END SUBROUTINE deallocate_response_function_buffers


    SUBROUTINE print_all_chksums(chksum_my_response_func, &
                                 chksum_my_exp_coeffs, &
                                 chksum_my_green_func)

        REAL*8,INTENT(IN) :: chksum_my_response_func, &
                             chksum_my_exp_coeffs, &
                             chksum_my_green_func

        REAL*8 :: chksum_response_func_sum, chksum_exp_coeffs_sum,chksum_green_func_sum
        INTEGER ::mpi_result

!        CALL write_debug("Chksum: (resp,exp,green)"//num2str(chksum_my_response_func)&
!                         //num2str(chksum_my_exp_coeffs)//num2str(chksum_my_green_func))


        CALL MPI_reduce(chksum_my_response_func,chksum_response_func_sum,1, MPI_REAL8, MPI_SUM, &
                        0, mpi_comm_global, mpi_result)


        CALL MPI_reduce(chksum_my_exp_coeffs,chksum_exp_coeffs_sum,1, MPI_REAL8, MPI_SUM, &
                        0, comm_basisfuncs%mpi_comm, mpi_result)


         CALL MPI_reduce(chksum_my_green_func,chksum_green_func_sum,1, MPI_REAL8, MPI_SUM, &
                         0, comm_basisfuncs%mpi_comm, mpi_result)


         CALL write_stdout("Total chksum (resp)"//num2str(chksum_response_func_sum))
         CALL write_stdout("Total chksum (exp)"//num2str(chksum_exp_coeffs_sum))
         CALL write_stdout("Total chksum (green)"//num2str(chksum_green_func_sum))

     END SUBROUTINE print_all_chksums

!this routines calculates k_points starting from i_k_point up to i_k_point+n_points
!the result is temporarily saved. The purpose is to reduce amount of discrete fourier transforms
    SUBROUTINE get_resp_func_in_k_omega_k(basbas_worklist, &
                                          max_n_rows_over_all_cpus, &
                                          n_frequencies_all, &
                                          frequencies_all, &
                                          i_k_point, &
                                          n_basbas, &
                                          frequency, &
                                          data_out)
        type(basbas_distribution), DIMENSION(comm_basisfuncs%n_tasks) :: basbas_worklist
        INTEGER, INTENT(IN) :: max_n_rows_over_all_cpus, n_frequencies_all
        REAL*8, DIMENSION(n_frequencies_all), INTENT(IN) :: frequencies_all


        INTEGER, INTENT(IN) :: i_k_point, n_basbas
        REAL*8,INTENT(IN) :: frequency
        DOUBLE COMPLEX,DIMENSION(n_basbas,n_basbas), INTENT(OUT) :: data_out

        INTEGER :: n_points, &
                   i_k_point_start
        LOGICAL, SAVE :: is_first_call = .TRUE.
        REAL*8, SAVE :: debug_compare_freq_threshold = 1.d-10

        if(is_first_call .OR. &
           (.NOT. has_result_last_calculation(result_last_calc, &
                                             i_k_point, frequency))) then

            i_k_point_start = i_k_point
            n_points = min(n_k_points - i_k_point + 1,5)


            if(.NOT. is_first_call) then
                CALL free_result_last_calculation(result_last_calc)
            endif
!CALL write_debug("Iam here"//num2str(i_k_point_start)&
!//num2str(n_points)&
!//num2str(max_n_rows_over_all_cpus)&
!//num2str(n_frequencies_all))
            CALL calculate_k_point(i_k_point_start, n_points, basbas_worklist, &
                                   max_n_rows_over_all_cpus, &
                                   n_frequencies_all, &
                                   frequencies_all)
            is_first_call = .FALSE.
        endif

        CALL get_result(result_last_calc, i_k_point,frequency, n_basbas, data_out)

        if(i_k_point == -2) then
            if(frequency>debug_compare_freq_threshold) then
                CALL write_stdout("DEBUG COMPARE"//num2str(frequency))
                CALL DEBUG_COMPARE(frequency,i_k_point)
                debug_compare_freq_threshold = frequency * 10.d0
            endif
        endif
CONTAINS
    SUBROUTINE DEBUG_COMPARE(omega, i_k_point)
            REAL*8, INTENT(IN) :: omega
            INTEGER, INTENT(IN) :: i_k_point
            INTEGER :: n_rows, n_columns, i_freq
            complex*16, dimension(n_basbas,n_basbas,n_k_points) :: polar_kspace

            CALL evaluate_polarisability_kspace( 1, omega, KS_eigenvalue, KS_eigenvector, &
! from Xinguo
!            CALL evaluate_polarisability_kspace(n_low_state, omega, KS_eigenvalue, KS_eigenvector, &
                                                KS_eigenvector_complex, &
                                                occ_numbers, polar_kspace)


            do n_rows=1,n_basbas
                do n_columns=1,n_basbas
if(abs(polar_kspace(n_rows,n_columns, i_k_point) &
                            -data_out(n_rows,n_columns))<1.d-4) CYCLE
if(myid==0) write(use_unit,*) i_k_point, n_rows,n_columns, &
                        polar_kspace(n_rows,n_columns, i_k_point), &
                        data_out(n_rows,n_columns), &
                        abs(polar_kspace(n_rows,n_columns, i_k_point) &
                            -data_out(n_rows,n_columns))
                enddo
            enddo
        end subroutine

    END SUBROUTINE get_resp_func_in_k_omega_k

    SUBROUTINE finalize_calculation_responsefunc_k()
        CALL free_result_last_calculation(result_last_calc)
    END SUBROUTINE

    SUBROUTINE triggered_timer(time_stop)
        DOUBLE PRECISION,INTENT(INOUT) :: time_stop
        DOUBLE PRECISION,SAVE :: start_time,start_time_clock
        DOUBLE PRECISION,SAVE :: stop_time,stop_time_clock
        LOGICAL,SAVE :: is_start_call = .TRUE.

        if(is_start_call) then
            call get_timestamps(start_time, start_time_clock)
            is_start_call = .FALSE.
        else
            call get_timestamps(stop_time, stop_time_clock)

            time_stop = time_stop + stop_time-start_time
            is_start_call = .TRUE.
        endif
    END SUBROUTINE triggered_timer

end module cRPA_calculation_responsefunc_k
