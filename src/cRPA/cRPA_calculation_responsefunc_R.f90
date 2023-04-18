!Calculates response function chi_0 in realspace
module cRPA_calculation_responsefunc_R
    USE cRPA_view
    USE cRPA_storage
    USE cRPA_parallelism
    USE cRPA_calculation_greensfunc
    USE cRPA_calculation_integration
    USE cRPA_calculation_energy

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: calculate_tau_point, &
              get_resp_func_in_k_omega_R

! I removed most of the screening routines to reduce complexity
! Routines are partially only commented
!      - check wheater they are up to date before you rely on them!
!    DOUBLE PRECISION :: SCREENING_TOLERANCE
!    PARAMETER (SCREENING_TOLERANCE = 1.d-20)

    DOUBLE COMPLEX, EXTERNAL :: zdotu

    DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:) :: max_row_reduced_column_rl_vecs

!buffers for basbas pair calculation
    DOUBLE COMPLEX,ALLOCATABLE, DIMENSION(:,:) :: res_column_bl, &
                                                  res_column_br, &
                                                  res_reduced_column_lr, &
                                                  res_reduced_column_rl, &
                                                  res_row_lb, &
                                                  res_row_rb

    DOUBLE COMPLEX,ALLOCATABLE, DIMENSION(:,:,:) :: res_reduced_column_lr_vecs_part1, &
                                                    res_reduced_column_rl_vecs_part_1, &
                                                    res_reduced_column_rl_vecs_part_2, &
                                                    res_reduced_column_rl_vecs_part3, &
                                                    res_reduced_column_rl_vecs_part4, &
                                                    ec_local_onsite, &
                                                    ec_local_offsite, &
                                                    ec_requested_onsite, &
                                                    ec_requested_offsite

    INTEGER :: i_basbas_row_save_part1, &
               i_basbas_row_save_part2, &
               i_basbas_row_save_part3, &
               i_basbas_row_save_part4

CONTAINS
    SUBROUTINE calculate_tau_point(rvecs, work_dist,cur_tau, &
                                   result_this_tau)
       type(realspace_vectors), INTENT(IN) :: rvecs
       type(basbas_distribution), INTENT(IN) :: work_dist
       REAL*8, INTENT(IN) :: cur_tau

       DOUBLE COMPLEX, DIMENSION(rvecs%n_vecs, &
                         work_dist%n_columns, &
                         work_dist%n_rows),INTENT(OUT) :: result_this_tau

       REAL*8 :: tau_point_time, tau_point_time_clock
       REAL*8 :: chksum_my_exp_coeffs,chksum_my_green_func

       INTEGER :: i, mpi_result

       INTEGER, DIMENSION(rvecs%n_vecs):: i_vec_inversion
       INTEGER, DIMENSION(rvecs%n_vecs,rvecs%n_vecs)::  i_vec_diff_1_minus_2, &
                                                        i_vec_diff_2_minus_1

       CALL write_stdout("----------------------------------------------------")
       CALL write_stdout("Calculating tau="//num2str(cur_tau))

       CALL get_realspace_vectors_index_inversion(rvecs, i_vec_inversion)
       CALL get_realspace_vectors_index_all_possible_differences(rvecs, &
                                                                 i_vec_diff_1_minus_2)
       i_vec_diff_2_minus_1(:,:) = transpose(i_vec_diff_1_minus_2)

       CALL get_timestamps(tau_point_time, tau_point_time_clock)

       result_this_tau(:,:,:) = 0.d0

       CALL calculate_local_greensfunc(rvecs, cur_tau)


       CALL allocate_response_function_for_basbas_pair_buffers(n_basis, &
                                                               n_basis, &
                                                               n_basis, &
                                                               rvecs)

       CALL calculate_basbas_pair_of_distribution(rvecs, &
                                                  i_vec_inversion, &
                                                  i_vec_diff_1_minus_2, &
                                                  i_vec_diff_2_minus_1, &
                                                  work_dist, &
                                                  result_this_tau(:,:,:))


       CALL deallocate_response_function_for_basbas_pair_buffers()

       CALL MPI_BARRIER(mpi_comm_global, mpi_result)

       CALL get_times(tau_point_time, tau_point_time_clock)
       CALL output_times("2X","Time for this tau point",tau_point_time, tau_point_time_clock)


       chksum_my_green_func = chksum_green_func(green_funcs_local)
       chksum_my_exp_coeffs = chksum_exp_coeffs(exp_coeffs_local)

       CALL print_all_chksums(chksum_3d_matrix(abs(result_this_tau(:,:,:))), &
                              chksum_my_exp_coeffs, &
                              chksum_my_green_func)

    END SUBROUTINE calculate_tau_point

    SUBROUTINE allocate_response_function_for_basbas_pair_buffers(n_basis, &
                                                                  n_basis_local, &
                                                                  n_basis_requested, &
                                                                  rvecs)
        INTEGER, INTENT(IN) :: n_basis, n_basis_local,n_basis_requested
        type(realspace_vectors), INTENT(IN) :: rvecs

!TODO check alloc
        ALLOCATE(ec_local_onsite(n_basis, n_basis_local,rvecs%n_vecs))
        ALLOCATE(ec_local_offsite(n_basis, n_basis_local,rvecs%n_vecs))
        ALLOCATE(ec_requested_onsite(n_basis, n_basis_requested,rvecs%n_vecs))
        ALLOCATE(ec_requested_offsite(n_basis, n_basis_requested,rvecs%n_vecs))


        ALLOCATE(res_column_bl(n_basis, n_basis_local))
        ALLOCATE(res_column_br(n_basis, n_basis_requested))

        ALLOCATE(res_row_lb(n_basis_local, n_basis))
        ALLOCATE(res_row_rb(n_basis_requested, n_basis))


        ALLOCATE(res_reduced_column_rl(n_basis_requested, n_basis_local))
        ALLOCATE(res_reduced_column_lr(n_basis_local, n_basis_requested))
        ALLOCATE(res_reduced_column_lr_vecs_part1(n_basis_local,n_basis_requested, rvecs%n_vecs))
        ALLOCATE(res_reduced_column_rl_vecs_part_1(n_basis_requested, n_basis_local, rvecs%n_vecs))
        ALLOCATE(res_reduced_column_rl_vecs_part_2(n_basis_requested, n_basis_local, rvecs%n_vecs))
        ALLOCATE(res_reduced_column_rl_vecs_part3(n_basis_requested, n_basis_local, rvecs%n_vecs))
        ALLOCATE(res_reduced_column_rl_vecs_part4(n_basis_requested, n_basis_local, rvecs%n_vecs))
        ALLOCATE(max_row_reduced_column_rl_vecs(rvecs%n_vecs))

        res_reduced_column_lr = 0
        res_reduced_column_rl = 0
        res_reduced_column_rl_vecs_part3 = 0
        max_row_reduced_column_rl_vecs(:) = 0.d0


    END SUBROUTINE allocate_response_function_for_basbas_pair_buffers

    SUBROUTINE calculate_basbas_pair_of_distribution(rvecs, &
                                                     i_vec_inversion, &
                                                     i_vec_diff_1_minus_2, &
                                                     i_vec_diff_2_minus_1, &
                                                     work_dist, &
                                                     result_basbas)
        type(realspace_vectors), INTENT(IN) :: rvecs
        INTEGER, DIMENSION(rvecs%n_vecs),INTENT(IN):: i_vec_inversion
        INTEGER, DIMENSION(rvecs%n_vecs,rvecs%n_vecs),INTENT(IN):: &
                                                         i_vec_diff_1_minus_2, &
                                                         i_vec_diff_2_minus_1

        type(basbas_distribution),INTENT(IN) :: work_dist
        DOUBLE COMPLEX, DIMENSION(rvecs%n_vecs,work_dist%n_columns,work_dist%n_rows),INTENT(INOUT) :: result_basbas


        INTEGER :: i_basbas_row, i_basbas_column, i_my_basbas_row, i_spin, &
                   exp_coeff_vecs_local, exp_coeff_vecs_requested, i_part, &
                   i_local_or_requested


        do i_basbas_row = work_dist%n_offset_rows + 1, &
                          work_dist%n_offset_rows + work_dist%n_rows

           i_my_basbas_row = i_basbas_row - work_dist%n_offset_rows

           result_basbas(:,:,i_my_basbas_row) = 0.d0

           CALL load_exp_column(exp_coeffs_local,i_basbas_row,&
                                i_vec_inversion, &
                                EXP_COEFFICIENT_LOCAL, &
                                exp_coeff_vecs_local, &
                                ec_local_onsite, &
                                ec_local_offsite)

            do i_spin = 1, n_spin

                i_basbas_row_save_part1 = -1
                i_basbas_row_save_part2 = -1
                i_basbas_row_save_part3 = -1
                i_basbas_row_save_part4 = -1

                do i_basbas_column = 1, work_dist%n_columns

                    CALL communicate_exp_coeffs_compressed(exp_coeffs_local,i_basbas_column)

                    if(exp_coeffs_local%coeff(i_basbas_column)%is_stored_locally) then
                        i_local_or_requested = EXP_COEFFICIENT_LOCAL
                    else
                        i_local_or_requested = EXP_COEFFICIENT_REQUESTED
                    endif

                    CALL load_exp_column(exp_coeffs_local,i_basbas_column, &
                                         i_vec_inversion, &
                                         i_local_or_requested, &
                                         exp_coeff_vecs_requested, &
                                         ec_requested_onsite, &
                                         ec_requested_offsite)
                    do i_part = 1,4
                        CALL calculate_response_function_for_basbas_pair(i_basbas_row, &
                                                                         i_spin, &
                                                                         i_basbas_column, &
                                                                         rvecs, &
                                                                         i_part, &
                                                                         i_vec_inversion, &
                                                                         i_vec_diff_1_minus_2, &
                                                                         i_vec_diff_2_minus_1, &
                                                                         green_funcs_local%g, &
                                                                         green_funcs_local_minus_tau%g, &!green_funcs_requested%g, &
                                                                         ec_local_onsite, &
                                                                         ec_local_offsite, &
                                                                         ec_requested_onsite, &
                                                                         ec_requested_offsite, &
                                                                         result_basbas(:,:,i_my_basbas_row))
                    enddo
                enddo
            enddo
        enddo

    END SUBROUTINE

        SUBROUTINE calculate_response_function_for_basbas_pair(i_basbas_row, &
                                                               i_spin, &
                                                               i_basbas_column, &
                                                               rvecs, &
                                                               i_part, &
                                                               i_vec_inversion, &
                                                               i_vec_diff_1_minus_2, &
                                                               i_vec_diff_2_minus_1, &
                                                               g_local, &
                                                               g_requested, &
                                                               exp_column_local, &
                                                               exp_column_local_offsite, &
                                                               exp_column_requested, &
                                                               exp_column_requested_offsite, &
                                                               result_basbas)

            INTEGER, INTENT(IN) :: i_basbas_row,i_spin, i_basbas_column
            type(realspace_vectors), INTENT(IN) :: rvecs
            INTEGER, INTENT(IN) :: i_part
            INTEGER, DIMENSION(rvecs%n_vecs),INTENT(IN):: i_vec_inversion
            INTEGER, DIMENSION(rvecs%n_vecs,rvecs%n_vecs),INTENT(IN):: &
                                                            i_vec_diff_1_minus_2, &
                                                            i_vec_diff_2_minus_1


            DOUBLE COMPLEX,DIMENSION(n_basis, &
                             n_basis,rvecs%n_vecs,n_spin), INTENT(IN) :: g_local, &
                                                                         g_requested


            DOUBLE COMPLEX,DIMENSION(n_basis, &
                                     n_basis, &
                             rvecs%n_vecs), INTENT(IN) :: exp_column_local, &
                                                          exp_column_local_offsite, &
                                                          exp_column_requested, &
                                                          exp_column_requested_offsite

            DOUBLE COMPLEX, DIMENSION(rvecs%n_vecs,n_basbas), INTENT(INOUT) :: result_basbas
            INTEGER :: i_vec_R, i_vec_S, i_vec_T, i_vec_minusSplusR, i_vec_SminusR

            type(realspace_vector) :: vec_R,vec_S, vec_T, vec_sum,vec_tmp, &
                                      vec_minusR, vec_minusS, vec_minusT, &
                                      vec_minusSplusR

            REAL*8 :: prefactor, max_norm_green_requested_exp_requested

            INTEGER ::  i,i_green_func

!            max_norm_green_requested_exp_requested = &
!                        maxval(green_funcs_local_minus_tau%g_norm_column(:,i_spin)) * &
!                        maxval(loaded_row_norm_max_column(:))

            if(i_part == 1) then
              CALL calculate_part1()
            endif

            if(i_part == 2) then
              CALL calculate_part2()
            endif

            if(i_part == 3) then
              CALL calculate_part3()
            endif


            if(i_part == 4) then
              CALL calculate_part4()
            endif


        CONTAINS

        SUBROUTINE calculate_part1()
          DOUBLE COMPLEX :: tmp

          CALL calculate_part1_tmp(i_basbas_row)

          do i_vec_R = 1,rvecs%n_vecs


            CALL calculate_part_1_lr()

            tmp = trace_part_reduced_columns_transposed(n_basis,n_basis, &
                                                        res_reduced_column_lr(:,:), &
                                                        res_reduced_column_lr_vecs_part1(:,:,i_vec_R))


            result_basbas(i_vec_R,i_basbas_column) = result_basbas(i_vec_R,i_basbas_column) &
                                                     + 1.d0 * tmp
!write(use_unit,*) i_spin,i_vec_R,i_basbas_column,result_basbas(i_vec_R,i_basbas_column), tmp
          enddo

        END SUBROUTINE calculate_part1

        SUBROUTINE calculate_part2_tmp(i_basbas_row_cur)
            INTEGER, INTENT(IN) :: i_basbas_row_cur

            INTEGER :: i_vec_R_tmp,i_vec_Rminus_T

            if(i_basbas_row_save_part2 == i_basbas_row_cur) return

            res_reduced_column_rl_vecs_part_2(:,:,:) = 0.d0
            do i_vec_R_tmp=1,rvecs%n_vecs
                do i_vec_T = 1, rvecs%n_vecs
                    i_vec_Rminus_T = i_vec_diff_1_minus_2(i_vec_T,i_vec_R_tmp)

                    CALL matrix_product_row_column(n_basis, &
                                                   exp_column_local(:,:, i_vec_inversion(i_vec_T)), &
                                                   n_basis, &
                                                   g_local(:,:,i_vec_Rminus_T,i_spin), &
                                                   res_reduced_column_rl_vecs_part_2(:,:,i_vec_R_tmp))
                enddo
            enddo
            i_basbas_row_save_part2 = i_basbas_row_cur
        END SUBROUTINE calculate_part2_tmp

        SUBROUTINE calculate_part2()
          DOUBLE COMPLEX :: tmp

          CALL calculate_part2_tmp(i_basbas_row)

          do i_vec_R = 1,rvecs%n_vecs
            res_column_bl(:,:) = 0
            do i_vec_S = 1, rvecs%n_vecs
                i_vec_SminusR = i_vec_diff_2_minus_1(i_vec_S,i_vec_R)
                CALL matrix_product_row_column(n_basis, &
                                               exp_column_requested(:,:, i_vec_SminusR), &
                                               n_basis, &
                                               g_requested(:,:,i_vec_S,i_spin), &
                                               res_column_bl(:,:))

            enddo

            tmp = trace_part_reduced_columns(n_basis,n_basis, &
                                             res_reduced_column_rl_vecs_part_2(:,:,i_vec_R),&
                                             res_column_bl(:,:))


            result_basbas(i_vec_R,i_basbas_column) = result_basbas(i_vec_R,i_basbas_column) &
                                                     + tmp
          enddo
        END SUBROUTINE calculate_part2

        SUBROUTINE calculate_part1_tmp(i_basbas_row_cur)
            INTEGER, INTENT(IN) :: i_basbas_row_cur
            INTEGER :: i_vec_R_tmp


            if(i_basbas_row_save_part1 == i_basbas_row_cur) return

!            res_reduced_column_rl_vecs_part_1(:,:,:) = 0.d0
            res_reduced_column_lr_vecs_part1(:,:,:) = 0.d0

            do i_vec_R_tmp = 1, rvecs%n_vecs
                do i_vec_T = 1, rvecs%n_vecs

                    i_green_func = i_vec_diff_2_minus_1(i_vec_T,i_vec_R_tmp)

!                    if(is_part_1_rl_screened_level2()) cycle

!                    CALL matrix_product_row_column(n_basisfunc_extend_requested, &
!                                                   g_requested(:,:,i_green_func,i_spin), &
!                                                   n_basisfunc_extend_local, &
!                                                   exp_column_local(:,:, i_vec_T), &
!                                                   res_reduced_column_rl_vecs_part_1(:,:,i_vec_R_tmp))

                    CALL matrix_product_row_column_transposed &
                                                  (n_basis, &
                                                   g_requested(:,:,i_green_func,i_spin), &
                                                   n_basis, &
                                                   exp_column_local_offsite(:,:, i_vec_T), &
                                                   res_reduced_column_lr_vecs_part1(:,:,i_vec_R_tmp))
                enddo

!write(use_unit,*) "max abs diff", &
!                maxval(abs(res_reduced_column_rl_vecs_part_1(:,:,i_vec_R_tmp) - &
!                      transpose(res_reduced_column_lr_vecs(:,:,i_vec_R_tmp))))

            enddo

            i_basbas_row_save_part1 = i_basbas_row_cur
        END SUBROUTINE calculate_part1_tmp

        LOGICAL FUNCTION is_part_1_rl_screened_level2()
!            if(green_funcs_requested%g_norm_row(i_green_func,i_spin) *&
!               loaded_row_norm_max_column(i_vec_T) &
!               < SCREENING_TOLERANCE) then
!               is_part_1_rl_screened_level2 = .TRUE.
!
!               return
!            endif


            is_part_1_rl_screened_level2 = .FALSE.
        END FUNCTION
        SUBROUTINE calculate_part_1_lr()
            res_reduced_column_lr(:,:) = 0
            do i_vec_S = 1, rvecs%n_vecs
                i_vec_SminusR = i_vec_diff_1_minus_2(i_vec_S,i_vec_R)

                i_green_func =  i_vec_inversion(i_vec_S)

!                if(is_part_1_lr_screened_level2()) cycle


                CALL matrix_product_row_column(n_basis, &
                                               g_local(:,:,i_green_func,i_spin), &
                                               n_basis, &
                                               exp_column_requested_offsite(:,:,(i_vec_SminusR)), &
                                               res_reduced_column_lr)
            enddo
        END SUBROUTINE calculate_part_1_lr


        LOGICAL FUNCTION is_part_1_lr_screened_level2()
!            if(green_funcs_local%g_norm_row(i_green_func,i_spin) * &
!               loaded_column_norm_max_column(i_vec_SminusR) &
!               < SCREENING_TOLERANCE) then
!               is_part_1_lr_screened_level2 = .TRUE.
!
!               return
!            endif

            is_part_1_lr_screened_level2 = .FALSE.
        END FUNCTION

        SUBROUTINE calculate_part3_tmp(i_basbas_row_cur)
            INTEGER, INTENT(IN) :: i_basbas_row_cur

            if(i_basbas_row_save_part3 == i_basbas_row_cur) return

            res_reduced_column_rl_vecs_part3(:,:,:) = 0

            do i_vec_S = 1, rvecs%n_vecs
               do i_vec_T = 1, rvecs%n_vecs
                 i_green_func = i_vec_diff_2_minus_1(i_vec_T,i_vec_S)

                 if(is_part3_screened_level3()) cycle

                 CALL matrix_product_row_column(n_basis, &
                                                g_requested(:,:,i_green_func,i_spin), &
                                                n_basis, &
                                                exp_column_local_offsite(:,:, i_vec_T), &
                                                res_reduced_column_rl_vecs_part3(:,:,i_vec_S))
              enddo

              CALL get_norm_max_row(n_basis ,&
                                    n_basis ,&
                                    abs(res_reduced_column_rl_vecs_part3(:,:,i_vec_S)), &
                                    max_row_reduced_column_rl_vecs(i_vec_S))

            enddo

            i_basbas_row_save_part3 = i_basbas_row_cur

        END SUBROUTINE calculate_part3_tmp

        SUBROUTINE calculate_part3()
          DOUBLE COMPLEX:: tmp,tmp2

          CALL calculate_part3_tmp(i_basbas_row)

          do i_vec_R = 1, rvecs%n_vecs

!            res_column_bl(:,:) = 0
            res_row_lb(:,:) = 0.d0

            do i_vec_S = 1, rvecs%n_vecs
                i_vec_minusSplusR = i_vec_diff_2_minus_1(i_vec_S,i_vec_R)

                if(is_part3_screened_level2()) cycle

!                CALL matrix_product_column_reduced_row(n_basisfunc_extend_requested,&
!                                               exp_column_requested(:,:,i_vec_minusSplusR), &
!                                               n_basisfunc_extend_local, &
!                                               res_reduced_column_rl_vecs(:,:,i_vec_S),&
!                                               res_column_bl)

                CALL matrix_product_column_reduced_row_transposed(n_basis,&
                                               exp_column_requested(:,:,i_vec_minusSplusR), &
                                               n_basis, &
                                               res_reduced_column_rl_vecs_part3(:,:,i_vec_S),&
                                               res_row_lb)

            enddo

            i_green_func = i_vec_inversion(i_vec_R)
!            i_green_func = i_vec_R

            tmp2 = trace_part_row_column_transposed(n_basis, &
                                                    g_local(:,:,i_green_func,i_spin),res_row_lb)


            result_basbas(i_vec_R,i_basbas_column) &
                                   = result_basbas(i_vec_R,i_basbas_column) + &
                                     1.d0 * tmp2
          enddo
        END SUBROUTINE calculate_part3

        SUBROUTINE calculate_part4_tmp(i_basbas_row_cur)
            INTEGER, INTENT(IN) :: i_basbas_row_cur

            INTEGER :: i_vec_S_tmp,i_vec_TminusS

            if(i_basbas_row_save_part4 == i_basbas_row_cur) return

            res_reduced_column_rl_vecs_part4(:,:,:) = 0.d0

            do i_vec_S_tmp = 1,rvecs%n_vecs
             do i_vec_T = 1, rvecs%n_vecs

                i_vec_TminusS = i_vec_diff_1_minus_2(i_vec_T,i_vec_S_tmp)
                CALL matrix_product_row_column(n_basis, &
                                               exp_column_local(:,:,i_vec_inversion(i_vec_T)), &
                                               n_basis, &
                                               g_local(:,:,i_vec_TminusS,i_spin), &
                                               res_reduced_column_rl_vecs_part4(:,:,i_vec_S_tmp))
             enddo
            enddo

            i_basbas_row_save_part4 = i_basbas_row_cur
        END SUBROUTINE calculate_part4_tmp

        SUBROUTINE calculate_part4()
           INTEGER::i_vec_minusSplusR
           DOUBLE COMPLEX :: tmp

           CALL calculate_part4_tmp(i_basbas_row)

           do i_vec_R = 1, rvecs%n_vecs
             res_column_bl(:,:) = 0.d0

             do i_vec_S = 1,rvecs%n_vecs
                    i_vec_minusSplusR = i_vec_diff_2_minus_1(i_vec_S,i_vec_R)

                    CALL matrix_product_row_column(n_basis, &
                                                   res_reduced_column_rl_vecs_part4(:,:,i_vec_S), &
                                                   n_basis, &
                                                   exp_column_requested_offsite(:,:, i_vec_inversion(i_vec_minusSplusR)), &
                                                   res_column_bl(:,:))

             enddo

             tmp = trace_part_reduced_columns(n_basis,&
                                              n_basis, &
                                              res_column_bl(:,:), &
                                              g_requested(:,:,i_vec_R,i_spin))

!            write(use_unit,*) "ohoh-- ",i_vec_R,i_basbas_row,i_basbas_column,tmp2,tmp
!
             result_basbas(i_vec_R,i_basbas_column) &
                                    = result_basbas(i_vec_R,i_basbas_column) + &
                                      1.d0 * tmp

          enddo
        END SUBROUTINE calculate_part4

       SUBROUTINE calculate_part_3()
           INTEGER::i_vec_TminusR
           DOUBLE COMPLEX :: tmp

            do i_vec_R = 1, rvecs%n_vecs
            res_reduced_column_rl_vecs_part3(:,:,1:2) = 0.d0
            do i_vec_T = 1, rvecs%n_vecs
                res_reduced_column_rl_vecs_part3(:,:,1) = 0.d0

                do i_vec_S = 1, rvecs%n_vecs
                    i_vec_TminusR = i_vec_diff_1_minus_2(i_vec_T,i_vec_S)
                    i_vec_minusSplusR = i_vec_diff_2_minus_1(i_vec_S,i_vec_R)

                    CALL matrix_product_row_column(n_basis, &
                                                   g_local(:,:,i_vec_TminusR,i_spin), &
                                                   n_basis, &
                                                   transpose(exp_column_requested(:,:, i_vec_minusSplusR)), &
                                                   res_reduced_column_rl_vecs_part3(:,:,1))
                enddo

                CALL matrix_product_row_column(n_basis, &
                                               (exp_column_local(:,:,i_vec_inversion(i_vec_T))), &
                                               n_basis, &
                                               res_reduced_column_rl_vecs_part3(:,:,1), &
                                               res_reduced_column_rl_vecs_part3(:,:,2))
            enddo


            tmp = trace_part_reduced_columns(n_basis,&
                                             n_basis, &
                                             res_reduced_column_rl_vecs_part3(:,:,2), &
                                             g_requested(:,:,i_vec_R,i_spin))

!            write(use_unit,*) "ohoh-- ",i_vec_R,i_basbas_row,i_basbas_column,tmp2,tmp
!
            result_basbas(i_vec_R,i_basbas_column) &
                                   = result_basbas(i_vec_R,i_basbas_column) + &
                                     1.d0 * tmp



            enddo
        END SUBROUTINE calculate_part_3

        LOGICAL FUNCTION is_part3_screened_level1()
!i_vec_inversion(i_vec_R)
!            if(green_funcs_local%g_norm_column(i_vec_R,i_spin) * &
!               maxval(loaded_column_norm_max_row(:)) * &
!               max_norm_green_requested_exp_requested * &
!               exp_coeff_vecs_requested &
!               < SCREENING_TOLERANCE) then
!               is_part3_screened_level1 = .TRUE.
!               return
!            endif

            is_part3_screened_level1 = .FALSE.
        END FUNCTION

        LOGICAL FUNCTION is_part3_screened_level2()

!            if(loaded_column_norm_max_row(i_vec_minusSplusR) * &
!               max_row_reduced_column_rl_vecs(i_vec_S) &
!               < SCREENING_TOLERANCE) then
!               is_part3_screened_level2 = .TRUE.
!
!               return
!            endif

            is_part3_screened_level2 = .FALSE.
        END FUNCTION


        LOGICAL FUNCTION is_part3_screened_level3()
!            if(green_funcs_requested%g_norm_column(i_green_func,i_spin) *&
!               loaded_row_norm_max_column(i_vec_T) &
!               < SCREENING_TOLERANCE) then
!               is_part3_screened_level3 = .TRUE.
!
!               return
!            endif

            is_part3_screened_level3 = .FALSE.
        END FUNCTION


        SUBROUTINE matrix_product_row_column(n_dim_row,block_row, &
                                             n_dim_column, block_column, &
                                             res)
            INTEGER, INTENT(IN) :: n_dim_row, n_dim_column
            DOUBLE COMPLEX, DIMENSION(n_dim_row,n_basis),INTENT(IN) :: block_row
            DOUBLE COMPLEX, DIMENSION(n_basis,n_dim_column),INTENT(IN) :: block_column
            DOUBLE COMPLEX, DIMENSION(n_dim_row,n_dim_column),INTENT(INOUT) :: res

!            INTEGER :: i_row, i_column, i_basis
!
!            do i_column = 1, n_dim_column
!                do i_row = 1, n_dim_row
!                    do i_basis = 1,n_basis
!                        res(i_row, i_column) = res(i_row, i_column) + &
!                                               block_row(i_row,i_basis) * &
!                                               block_column(i_basis,i_column)
!                    enddo
!                enddo
!            enddo

            CALL zgemm('N','N',&
                       n_dim_row,n_dim_column,n_basis, &
                       DCMPLX(1.0d0,0.d0), &
                       block_row,n_dim_row, &
                       block_column,n_basis, &
                       DCMPLX(1.0d0,0.d0), &
                       res, n_dim_row)
        END SUBROUTINE matrix_product_row_column

       SUBROUTINE matrix_product_row_column_transposed(n_dim_row,block_row, &
                                             n_dim_column, block_column, &
                                             res)
            INTEGER, INTENT(IN) :: n_dim_row, n_dim_column
            DOUBLE COMPLEX, DIMENSION(n_dim_row,n_basis),INTENT(IN) :: block_row
            DOUBLE COMPLEX, DIMENSION(n_basis,n_dim_column),INTENT(IN) :: block_column
            DOUBLE COMPLEX, DIMENSION(n_dim_column,n_dim_row),INTENT(INOUT) :: res

            CALL zgemm('T','T',&
                       n_dim_row,n_dim_column,n_basis, &
                       DCMPLX(1.0d0,0.d0), &
                       block_column,n_basis, &
                       block_row,n_dim_row, &
                       DCMPLX(1.0d0,0.d0), &
                       res, n_dim_column)
        END SUBROUTINE matrix_product_row_column_transposed

        SUBROUTINE matrix_product_column_reduced_row(n_dim_column,block_column, &
                                                     n_dim_row,reduced_row, res_row)
            INTEGER, INTENT(IN) :: n_dim_column,n_dim_row
            DOUBLE COMPLEX, DIMENSION(n_basis,n_dim_column),INTENT(IN) :: block_column
            DOUBLE COMPLEX, DIMENSION(n_dim_column,n_dim_row),INTENT(IN) :: reduced_row
            DOUBLE COMPLEX, DIMENSION(n_basis,n_dim_row),INTENT(INOUT) :: res_row

            INTEGER :: i_row, i_column, i_basis

!            do i_row = 1, n_dim_row
!                do i_basis = 1,n_basis
!                    do i_column = 1, n_dim_column
!                        res_row(i_basis, i_row) = res_row(i_basis, i_row) + &
!                                               block_column(i_basis,i_column) * &
!                                               reduced_row(i_column,i_row)
!                    enddo
!                enddo
!            enddo

            CALL zgemm('N','N',&
                       n_basis,n_dim_row,n_dim_column, &
                       DCMPLX(1.0d0,0.d0), &
                       block_column,n_basis, &
                       reduced_row,n_dim_column, &
                       DCMPLX(1.0d0,0.d0), &
                       res_row, n_basis)
        END SUBROUTINE matrix_product_column_reduced_row


        SUBROUTINE matrix_product_column_reduced_row_transposed(n_dim_column,block_column, &
                                                                n_dim_row,reduced_row, res_column)
            INTEGER, INTENT(IN) :: n_dim_column,n_dim_row
            DOUBLE COMPLEX, DIMENSION(n_basis,n_dim_column),INTENT(IN) :: block_column
            DOUBLE COMPLEX, DIMENSION(n_dim_column,n_dim_row),INTENT(IN) :: reduced_row
            DOUBLE COMPLEX, DIMENSION(n_dim_row,n_basis),INTENT(INOUT) :: res_column

            INTEGER :: i_row, i_column, i_basis

            CALL zgemm('T','T',&
                       n_basis,n_dim_row,n_dim_column, &
                       DCMPLX(1.0d0,0.d0), &
                       reduced_row,n_dim_column, &
                       block_column,n_basis, &
                       DCMPLX(1.0d0,0.d0), &
                       res_column, n_dim_row)
        END SUBROUTINE matrix_product_column_reduced_row_transposed

        SUBROUTINE matrix_product_column_square(n_dim_column,block_column, &
                                                square, res_column)
            INTEGER :: n_dim_column
            DOUBLE COMPLEX, DIMENSION(n_basis,n_dim_column),INTENT(IN) :: block_column
            DOUBLE COMPLEX, DIMENSION(n_dim_column,n_dim_column),INTENT(IN) :: square
            DOUBLE COMPLEX, DIMENSION(n_basis,n_dim_column),INTENT(INOUT) :: res_column

            CALL zgemm('N','N',&
                       n_basis,n_dim_column,n_dim_column, &
                       DCMPLX(1.0d0,0.d0), &
                       block_column,n_basis, &
                       square,n_dim_column, &
                       DCMPLX(1.0d0,0.d0), &
                       res_column, n_basis)
        END SUBROUTINE matrix_product_column_square

        DOUBLE COMPLEX FUNCTION trace_part_row_column(n_dim,block_row,block_column)
            INTEGER, INTENT(IN) :: n_dim
            DOUBLE COMPLEX, DIMENSION(n_dim,n_basis),INTENT(IN) :: block_row
            DOUBLE COMPLEX, DIMENSION(n_basis,n_dim),INTENT(IN) :: block_column

!            trace_part_row_column = 0.d0
!            do i = 1, n_dim
!                trace_part_row_column = trace_part_row_column + &
!                                        ddot(n_basis, &
!                                             block_row(i,:),1, &
!                                             block_column(:,i),1)
!            enddo

            trace_part_row_column = zdotu(n_basis*n_dim, &
                                          block_row(:,:),1, &
                                          transpose(block_column(:,:)),1)


        END FUNCTION trace_part_row_column

       DOUBLE COMPLEX FUNCTION trace_part_row_column_transposed(n_dim,block_row,block_column_transposed)
            INTEGER, INTENT(IN) :: n_dim
            DOUBLE COMPLEX, DIMENSION(n_dim,n_basis),INTENT(IN) :: block_row, &
                                                                   block_column_transposed

            trace_part_row_column_transposed = zdotu(n_basis*n_dim, &
                                                     block_row(:,:),1, &
                                                     block_column_transposed(:,:),1)

!            trace_part_row_column_transposed = scalprod_stable(n_basis*n_dim, &
!                                                               block_row(:,:), &
!                                                               block_column_transposed(:,:))

        END FUNCTION trace_part_row_column_transposed


        DOUBLE COMPLEX FUNCTION trace_part_reduced_columns(n_dim_local,n_dim_requested, &
                                                   column_lr, &
                                                   column_rl)
            INTEGER, INTENT(IN) :: n_dim_local,n_dim_requested
            DOUBLE COMPLEX, DIMENSION(n_dim_local,n_dim_requested),INTENT(IN) :: column_lr
            DOUBLE COMPLEX, DIMENSION(n_dim_requested,n_dim_local),INTENT(IN) :: column_rl

!            trace_part_reduced_columns = 0
!            do i = 1, n_dim_local
!                trace_part_reduced_columns = trace_part_reduced_columns + &
!                                             ddot(n_dim_requested, &
!                                                  column_lr(i,:),1, &
!                                                  column_rl(:,i),1)
!            enddo

            trace_part_reduced_columns = zdotu(n_dim_requested * n_dim_local, &
                                               column_lr(:,:),1, &
                                               transpose(column_rl(:,:)),1)

        END FUNCTION trace_part_reduced_columns

        DOUBLE COMPLEX FUNCTION trace_part_reduced_columns_transposed(n_dim_local,n_dim_requested, &
                                                   column_lr, &
                                                   column_rl_transposed)
            INTEGER, INTENT(IN) :: n_dim_local,n_dim_requested
            DOUBLE COMPLEX, DIMENSION(n_dim_local,n_dim_requested),INTENT(IN) :: column_lr
            DOUBLE COMPLEX, DIMENSION(n_dim_local,n_dim_requested),INTENT(IN) :: column_rl_transposed

            trace_part_reduced_columns_transposed = zdotu(n_dim_requested * n_dim_local, &
                                                          column_lr(:,:),1, &
                                                          column_rl_transposed(:,:),1)

!            trace_part_reduced_columns_transposed = scalprod_stable(n_dim_requested * n_dim_local, &
!                                                                    column_lr(:,:), &
!                                                                    column_rl_transposed(:,:))


        END FUNCTION trace_part_reduced_columns_transposed

END SUBROUTINE

    SUBROUTINE deallocate_response_function_for_basbas_pair_buffers()

        DEALLOCATE(ec_local_onsite)
        DEALLOCATE(ec_local_offsite)
        DEALLOCATE(ec_requested_onsite)
        DEALLOCATE(ec_requested_offsite)

        DEALLOCATE(res_column_bl)
        DEALLOCATE(res_column_br)
        DEALLOCATE(res_row_lb)
        DEALLOCATE(res_row_rb)

        DEALLOCATE(res_reduced_column_rl)
        DEALLOCATE(res_reduced_column_lr)
        DEALLOCATE(res_reduced_column_lr_vecs_part1)
        DEALLOCATE(res_reduced_column_rl_vecs_part_1)
        DEALLOCATE(res_reduced_column_rl_vecs_part_2)
        DEALLOCATE(res_reduced_column_rl_vecs_part3)
        DEALLOCATE(res_reduced_column_rl_vecs_part4)
        DEALLOCATE(max_row_reduced_column_rl_vecs)

    END SUBROUTINE deallocate_response_function_for_basbas_pair_buffers


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

    SUBROUTINE get_resp_func_in_k_omega_R(n_basbas,i_k_point,omega, data_out)
        INTEGER, INTENT(IN) :: n_basbas, i_k_point
        REAL*8, INTENT(IN) :: omega
        DOUBLE COMPLEX, DIMENSION(n_basbas,n_basbas), INTENT(INOUT) :: data_out

        CALL calculate_response_function_in_k_omega(omega,n_tasks,i_k_point,data_out)
    END SUBROUTINE get_resp_func_in_k_omega_R

    SUBROUTINE calculate_response_function_in_k_omega(omega,n_tasks, i_k_point, data_out)
        REAL*8,INTENT(IN) :: omega
        INTEGER, INTENT(IN) :: n_tasks, i_k_point
        DOUBLE COMPLEX, DIMENSION(n_basbas,n_basbas), INTENT(INOUT) :: data_out


        type(basbas_distribution), ALLOCATABLE,DIMENSION(:) :: basbas_worklist, &
                                                               basbas_worklist_total
        DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:) :: data_tmp
        INTEGER :: i_basbas,i_row,mpierr
        REAL*8, SAVE :: LAST_OMEGA = -1.d0

        CALL distribute_basbas(n_basbas, n_tasks, n_tasks, &
                               basbas_worklist,basbas_worklist_total)

        if(abs(omega-LAST_OMEGA) > 1.d-7) then
            LAST_OMEGA= omega
!        CALL write_stdout("Starting response function calculation in k and omega")
            CALL calculate_resp_func_for_one_frequency_in_k(omega)
        endif


        ALLOCATE(data_tmp(n_basbas,n_basbas))
        data_tmp(:,:) = DCMPLX(0.d0,0.d0)

        do i_row = 1,basbas_worklist_total(myid+1)%n_rows
            data_tmp(basbas_worklist_total(myid+1)%n_offset_rows + i_row,:) = &
                                                response_func%one_omega_in_k(:,i_row,i_k_point)
        enddo

        CALL MPI_ALLREDUCE(data_tmp, &
                           data_out, size(data_tmp), &
                           MPI_COMPLEX16, MPI_SUM, mpi_comm_global , mpierr)

        DEALLOCATE(data_tmp)
        DEALLOCATE(basbas_worklist)
        DEALLOCATE(basbas_worklist_total)

    END SUBROUTINE calculate_response_function_in_k_omega

end module cRPA_calculation_responsefunc_R
