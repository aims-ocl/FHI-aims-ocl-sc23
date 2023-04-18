!Storage and processing of response function chi_0 if calculated in realspace
MODULE cRPA_storage_responsefunc_R
    USE dimensions
    USE numerical_utilities


    USE cRPA_view
    USE cRPA_parallelism_storage
    USE cRPA_flow_realspace_vector
    USE cRPA_calculation_integration
    USE cRPA_calculation_fouriertrans
    USE cRPA_calculation_expbasis
    USE cRPA_flow_adaptive_grid
    USE cRPA_storage_norms
    USE cRPA_calculation_energy
    IMPLICIT NONE

    type response_function
        ! rvecs, basbas_column, basbas_row, tau
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: element
        DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: buffer_ft

        DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:) :: one_omega_my_k_points, &
                                                     one_omega_in_k

        INTEGER :: n_rvecs
        type(realspace_vectors) :: rvecs
        type(adaptive_grid) :: tau_grid
        type(basbas_distribution) :: worklist
    end type

    type(response_function) :: response_func
    type(exp_basis) :: aux_exp_basis
    type(exp_basis) :: aux_exp_basis2


    LOGICAL :: is_response_func_storage_inited = .FALSE.

CONTAINS

!response function storage
    SUBROUTINE init_response_function_storage(rvecs, worklist)
        type(realspace_vectors) :: rvecs
        type(basbas_distribution), INTENT(IN) :: worklist

        !response_func
        if(is_response_func_storage_inited) then
            write(use_unit,*) "Response function storage already inited!"
            stop
        end if

        CALL write_stdout("Initializing response function storage")

        ALLOCATE(response_func%buffer_ft(rvecs%n_vecs, &
                                         worklist%n_columns, &
                                         worklist%n_rows))

        ALLOCATE(response_func%one_omega_in_k(worklist%n_columns, &
                                              worklist%n_rows, &
                                              rvecs%n_vecs))

        response_func%n_rvecs = rvecs%n_vecs
        response_func%rvecs = rvecs
        CALL copy_matrixrow(worklist, response_func%worklist)

        is_response_func_storage_inited = .TRUE.
    END SUBROUTINE init_response_function_storage


    SUBROUTINE set_response_function_tau_grid(tau_grid)
        type(adaptive_grid),INTENT(IN) :: tau_grid

        if(ALLOCATED(response_func%element)) then
            DEALLOCATE(response_func%element)
        endif

 !TODO check alloc
        ALLOCATE(response_func%element(response_func%rvecs%n_vecs, &
                                       response_func%worklist%n_columns, &
                                       response_func%worklist%n_rows, &
                                       tau_grid%n_points))
        response_func%element(:,:,:,:) = 0.d0

        CALL clone_adaptive_grid(tau_grid, response_func%tau_grid)
    END SUBROUTINE set_response_function_tau_grid

    SUBROUTINE add_response_function_tau_point(data,tau,i_pos)
        REAL*8, DIMENSION(response_func%rvecs%n_vecs, &
                          response_func%worklist%n_columns, &
                          response_func%worklist%n_rows),INTENT(IN) :: data
        REAL*8, INTENT(IN) :: tau
        INTEGER,INTENT(OUT) :: i_pos

        INTEGER :: i_pos_added
        REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: buffer
        type(adaptive_grid) :: a_grid
        REAL*8 :: prefactor

        if(.NOT. ALLOCATED(response_func%element)) then
            CALL create_adaptive_grid(1,(/tau/),(/0.d0/), &
                                      a_grid)

            CALL set_response_function_tau_grid(a_grid)
        endif

        if(.NOT. has_point(response_func%tau_grid,tau)) then
!TODO check alloc
            ALLOCATE(buffer(response_func%rvecs%n_vecs, &
                          response_func%worklist%n_columns, &
                          response_func%worklist%n_rows, &
                          response_func%tau_grid%n_points))

            buffer(:,:,:,:) = response_func%element(:,:,:,:)
            DEALLOCATE(response_func%element)


            CALL add_point_to_adaptive_grid(tau,0.d0, &
                                            response_func%tau_grid, &
                                            i_pos_added)
            i_pos = i_pos_added

!TODO check alloc
            ALLOCATE(response_func%element(response_func%rvecs%n_vecs, &
                                           response_func%worklist%n_columns, &
                                           response_func%worklist%n_rows, &
                                           response_func%tau_grid%n_points))

            response_func%element(:,:,:,1:i_pos-1) = buffer(:,:,:,1:i_pos-1)

            response_func%element(:,:,:,i_pos+1:response_func%tau_grid%n_points) &
                                 = buffer(:,:,:,i_pos:response_func%tau_grid%n_points-1)

            DEALLOCATE(buffer)
        else
            CALL find_point(response_func%tau_grid, tau, i_pos)
        endif

        prefactor = -2.d0

        !CALL check_for_nan()

        response_func%element(:,:,:,i_pos) = prefactor * data(:,:,:)
    !CONTAINS
        !SUBROUTINE check_for_nan()
        !    INTEGER :: i_vec, i_row, i_column

        !    do i_row = 1, response_func%worklist%n_rows
        !        do i_column = 1,response_func%worklist%n_columns
        !            do i_vec = 1,response_func%rvecs%n_vecs
        !                if(isnan(data(i_vec,i_column,i_row))) then
                            !response_func%element(i_vec,i_column,i_row,i_pos) = 0.d0
        !                    stop "nan on storage response func"
        !                endif
        !            enddo
        !        enddo
        !    enddo
        !END SUBROUTINE check_for_nan

    END SUBROUTINE add_response_function_tau_point

    SUBROUTINE calculate_resp_func_for_one_frequency_in_k(omega)
        REAL*8, INTENT(IN) :: omega

        CALL fouriertransform_resp_func_for_one_frequency(omega)
        CALL print_chksum_buffer_ft()

        CALL fouriertransform_resp_func_one_freq_R_to_k()
        CALL print_chksum_one_omega_in_k()

    END SUBROUTINE calculate_resp_func_for_one_frequency_in_k

    SUBROUTINE fouriertransform_resp_func_for_one_frequency(omega)
        REAL*8, INTENT(IN) :: omega

        INTEGER :: i_vector,i,j

        CALL write_stdout("Fourier transforming for omega"//num2str(omega))

        do i_vector = 1, response_func%rvecs%n_vecs
           CALL fouriertransform_resp_func_for_one_frequency_one_vector(omega,i_vector)
        enddo
    END SUBROUTINE fouriertransform_resp_func_for_one_frequency

    SUBROUTINE fouriertransform_resp_func_for_one_frequency_one_vector(omega,i_vector)
        REAL*8, INTENT(IN) :: omega
        INTEGER, INTENT(IN) :: i_vector

        INTEGER :: i_row, n_splined_data
        type(integration_grid) :: grid_splined, grid_lincont, int_grid
        REAL*8 :: max_dist_between_points,pi, &
                  max_interval_start, &
                  min_interval_end
        REAL*8, ALLOCATABLE, DIMENSION(:,:) :: data_splined, &
                                               data_lincont
        REAL*8, ALLOCATABLE, DIMENSION(:) :: tmp_buffer, tmp_buffer2

!i_omega and old_omega only for debug/compare purpose
        INTEGER :: i_omega = 0
        REAL*8,SAVE :: old_omega = 0.d0

        if( abs(old_omega-omega) > 1.d-10) then
           old_omega = omega
           i_omega = i_omega +1
           if(myid==0) write(use_unit,*) "i_omega",i_omega
        endif

        pi = acos(-1.d0)
        max_dist_between_points = 1.d0/omega!(2.d0*pi) / (omega)

        CALL create_integration_grid_from_adaptive_grid(response_func%tau_grid, &
                                                        int_grid)

        max_interval_start = 0.1d0!(-pi) / (omega)
        min_interval_end = pi / omega

!TODO refactor!
        ALLOCATE(tmp_buffer(response_func%worklist%n_columns))
        ALLOCATE(tmp_buffer2(response_func%worklist%n_columns))


        do i_row = 1,response_func%worklist%n_rows

if(omega<0.1d10) then
             CALL spline_extrapolate(int_grid, &
                                     response_func%worklist%n_columns, &
                                     response_func%element(i_vector,:,i_row,:), &
                                     max_dist_between_points, &
                                     grid_splined, &
                                     data_splined)
             CALL linear_extrapolate_boundaries(max_interval_start, &
                                                min_interval_end, &
                                                max_dist_between_points, &
                                                response_func%worklist%n_columns, &
                                                grid_splined, &
                                                data_splined, &
                                                grid_lincont, &
                                                data_lincont)
            DEALLOCATE(data_splined)
            CALL free_integration_grid(grid_splined)

!            if(i_vector == 1 .OR. i_vector == 10) then
!                CALL write_response_function_splined(i_vector, omega, &
!                                                     i_row, &
!                                                     grid_lincont, data_lincont)
!            endif
!
            CALL fouriertransform_x_to_t(grid_lincont, &
                                         response_func%worklist%n_columns, &
                                         data_lincont, &
                                         omega, &
                                         response_func%buffer_ft(i_vector,:,i_row))
            DEALLOCATE(data_lincont)
!!!            CALL free_integration_grid(grid_lincont)
tmp_buffer(:) =  response_func%buffer_ft(i_vector,:,i_row)
tmp_buffer2(:) = 0.d0
elseif(omega < -0.01d10) then
!            CALL fouriertransform_by_integration(int_grid%n_points, int_grid%points, &
!                                               response_func%worklist%n_columns, &
!                                               response_func%element(i_vector,:,i_row,:), &
!                                               omega, &
!                                               aux_exp_basis, &
!                                               tmp_buffer)


            CALL spline_cos_integration(int_grid, &
                                        response_func%worklist%n_columns, &
                                        response_func%element(i_vector,:,i_row,:), &
                                        min_interval_end, &
                                        omega, &
                                        tmp_buffer)

!            do n_splined_data=1,response_func%worklist%n_columns
!            response_func%buffer_ft(i_vector,n_splined_data,i_row)=DCMPLX(tmp_buffer(n_splined_data),0.d0)
!            enddo
    tmp_buffer2(:) = 0.d0
endif
!else
!            CALL fouriertransform_by_exp_basis(int_grid%n_points, int_grid%points, &
!                                               response_func%worklist%n_columns, &
!                                               response_func%element(i_vector,:,i_row,:), &
!                                               omega, &
!                                               aux_exp_basis, &
!                                               tmp_buffer)
!
!            tmp_buffer2(:) = 0.d0
!endif
            do n_splined_data=1,response_func%worklist%n_columns
               response_func%buffer_ft(i_vector,n_splined_data,i_row)= &
                   DCMPLX(tmp_buffer(n_splined_data)+tmp_buffer2(n_splined_data),0.d0)

!             write(use_unit,*) dble(response_func%buffer_ft(i_vector,n_splined_data,i_row)),tmp_buffer(n_splined_data)
!debug
!            if(1==1) then
!               if(abs( polar_Rspace(i_row+response_func%worklist%n_offset_rows,n_splined_data, i_vector,i_omega))>1.d-3) then
!               if(myid==0) write(use_unit,*)  i_vector, i_row, n_splined_data,&
!                   polar_Rspace(i_row+response_func%worklist%n_offset_rows,n_splined_data, i_vector,i_omega), &
!                   dble(response_func%buffer_ft(i_vector,n_splined_data,i_row)), &
!                   dble(polar_Rspace(i_row+response_func%worklist%n_offset_rows,n_splined_data, i_vector,i_omega))/ &
!                   dble(response_func%buffer_ft(i_vector,n_splined_data,i_row))
!               endif
!            endif

            enddo

!            CALL free_integration_grid(grid_lincont)

        enddo

        DEALLOCATE(tmp_buffer)
        DEALLOCATE(tmp_buffer2)

    END SUBROUTINE fouriertransform_resp_func_for_one_frequency_one_vector

    SUBROUTINE fouriertransform_resp_func_one_freq_R_to_k()

        INTEGER :: i_row

        do i_row = 1, response_func%worklist%n_rows
            CALL discrete_fouriertransform_R_to_k(response_func%rvecs, &
                                                  response_func%worklist%n_columns, &
                                                  transpose(response_func%buffer_ft(:,:,i_row)), &
                                                  response_func%one_omega_in_k(:,i_row,:))
        enddo

    END SUBROUTINE fouriertransform_resp_func_one_freq_R_to_k


    REAL*8 FUNCTION chksum_response_function_for_point(i_tau_point)

       INTEGER :: i_tau_point
       chksum_response_function_for_point = &
                                         chksum_3d_matrix(response_func%element(:,:,:,i_tau_point))! &


       CALL sync_real_number(chksum_response_function_for_point)
    END FUNCTION chksum_response_function_for_point


    REAL*8 FUNCTION chksum_response_function()

       INTEGER :: i_tau
       chksum_response_function = 0.d0
       do i_tau = 1, response_func%tau_grid%n_points
           chksum_response_function = chksum_response_function + &
                                      chksum_response_function_for_point(i_tau)
       enddo

    END FUNCTION chksum_response_function

    SUBROUTINE init_response_func_my_k_points()
!TODO check alloc
LOGICAL, save::bfirst = .TRUE.
if(bfirst) then
        bfirst =.FALSE.
        ALLOCATE(response_func%one_omega_my_k_points(n_basbas,n_basbas, n_k_points_task))
endif
    END SUBROUTINE init_response_func_my_k_points

    SUBROUTINE print_response_function()

        REAL*8 :: my_min_norm_column, &
                  my_min_norm_row, &
                  my_max_norm_column, &
                  my_max_norm_row, &
                  min_norm_column, &
                  min_norm_row, &
                  max_norm_column, &
                  max_norm_row, &
                  percent_small_elements, &
                  norm, cur_tau
        INTEGER :: i_vec,i_tau, &
                   i_elements_total, &
                   i_elements_positive, i_elements_negative, i_elements_small

        CALL write_stdout("Response function info")

        do i_tau = 1, response_func%tau_grid%n_points
            cur_tau = response_func%tau_grid%points(i_tau)
            my_max_norm_column = 0.d0
            my_min_norm_column = 10.d10
            my_max_norm_row = 0.d0
            my_min_norm_row = 10.d10

            do i_vec = 1, response_func%rvecs%n_vecs
                CALL get_norm_max_column(response_func%worklist%n_rows, &
                                         response_func%worklist%n_rows, &
                                         transpose(response_func%element(i_vec,:,:,i_tau)),&
                                         norm)

                my_max_norm_column = max(my_max_norm_column, &
                                         norm)

                my_min_norm_column = min(my_min_norm_column, &
                                         norm)

                CALL get_norm_max_row(response_func%worklist%n_rows, &
                                      response_func%worklist%n_rows, &
                                      transpose(response_func%element(i_vec,:,:,i_tau)),&
                                      norm)
                my_max_norm_row = max(my_max_norm_row, &
                                      norm)

                my_min_norm_row = min(my_min_norm_row, &
                                      norm)


            enddo

            CALL get_max_over_all_cpus(my_max_norm_column,max_norm_column)
            CALL get_min_over_all_cpus(my_min_norm_column,min_norm_column)

            CALL get_max_over_all_cpus(my_max_norm_row,max_norm_row)
            CALL get_min_over_all_cpus(my_min_norm_row,min_norm_row)

            CALL write_stdout("For tau"//num2str(cur_tau))
            CALL write_stdout("    norm max column row" //num2str(max_norm_column) &
                                                         //num2str(max_norm_row))

            CALL write_stdout("    norm min column row" //num2str(min_norm_column) &
                                                        //num2str(min_norm_row))

        enddo

        !i_elements_total = count(.NOT. isnan(response_func%element(:,:,:,:)))
        i_elements_total = size(response_func%element(:,:,:,:))
        i_elements_positive = count(response_func%element(:,:,:,:)>=0.d0)
        i_elements_negative = count(response_func%element(:,:,:,:)<0.d0)
        i_elements_small = count(abs(response_func%element(:,:,:,:))<1.d-10)

        CALL sync_integer(i_elements_total)
        CALL sync_integer(i_elements_positive)
        CALL sync_integer(i_elements_negative)
        CALL sync_integer(i_elements_small)

        percent_small_elements = 100.d0*dble(i_elements_small)/dble(i_elements_total)

        CALL write_stdout("Elements")
        CALL write_stdout("--- total"//num2str(i_elements_total))
        CALL write_stdout("--- positive"//num2str(i_elements_positive))
        CALL write_stdout("--- negative"//num2str(i_elements_negative))
        CALL write_stdout("--- small"//num2str(i_elements_small))
        CALL write_stdout("--- small in percent"//num2str(percent_small_elements))
        CALL write_stdout("-----------------------------------")

    END SUBROUTINE print_response_function

    SUBROUTINE print_chksum_buffer_ft()
        INTEGER :: mpi_result
        REAL*8 :: chksum_my_part, chksum_total

       chksum_my_part = chksum_my_part + &
                        chksum_3d_matrix(dble(response_func%buffer_ft(:,:,:))) + &
                        chksum_3d_matrix(aimag(response_func%buffer_ft(:,:,:)))

        CALL MPI_reduce(chksum_my_part,chksum_total,1, MPI_REAL8, MPI_SUM, &
                        0, mpi_comm_global, mpi_result)

        CALL write_stdout("Chksum (buffer_ft):"//num2str(chksum_total))

    END SUBROUTINE print_chksum_buffer_ft

    SUBROUTINE print_chksum_one_omega_in_k()
        INTEGER :: mpi_result
        REAL*8 :: chksum_my_part, chksum_total

       chksum_my_part = chksum_my_part + &
                        chksum_3d_matrix(dble(response_func%one_omega_in_k(:,:,:))) + &
                        chksum_3d_matrix(aimag(response_func%one_omega_in_k(:,:,:)))

        CALL MPI_reduce(chksum_my_part,chksum_total,1, MPI_REAL8, MPI_SUM, &
                        0, mpi_comm_global, mpi_result)

        CALL write_stdout("Chksum (resp_for_one_omega):"//num2str(chksum_total))

    END SUBROUTINE print_chksum_one_omega_in_k

    SUBROUTINE print_chksum_resp_func_my_k_points()
       INTEGER :: mpi_result, i_row
       REAL*8 :: chksum_my_part, chksum_total

       chksum_my_part = 0.d0

       chksum_my_part = chksum_my_part + &
                        chksum_3d_matrix(dble(response_func%one_omega_my_k_points(:,:,:))) + &
                        chksum_3d_matrix(aimag(response_func%one_omega_my_k_points(:,:,:)))

!CALL write_debug(num2str(chksum_my_part))
        CALL MPI_reduce(chksum_my_part,chksum_total,1, MPI_REAL8, MPI_SUM, &
                        0, mpi_comm_global, mpi_result)


        CALL write_stdout("Chksum (resp_for_all_k):"//num2str(chksum_total))
    END SUBROUTINE print_chksum_resp_func_my_k_points

    SUBROUTINE write_response_function()
       INTEGER :: i_vec
       CHARACTER(len=256) :: filename

       do i_vec = 1, response_func%rvecs%n_vecs

            if(i_vec /= 1 .AND. &
               i_vec /= response_func%rvecs%n_vecs/2 .AND. &
               i_vec /= response_func%rvecs%n_vecs) then
                cycle
            endif

            filename = trim("respele")//"_v"//trim(num2str(i_vec))

            CALL write_matrix_array_to_file(response_func%worklist%n_rows, &
                                            response_func%worklist%n_columns, &
                                            response_func%tau_grid%n_points, &
                                            response_func%element(i_vec,:,:,:), &
                                            response_func%tau_grid%points(:), &
                                            filename)
       enddo
    END SUBROUTINE write_response_function

    SUBROUTINE write_response_function_splined(i_vec, omega, &
                                               i_row, &
                                               grid_splined, splined_data)
       INTEGER, INTENT(IN) :: i_vec, i_row
       REAL*8, INTENT(IN) :: omega
       type(integration_grid),INTENT(IN) :: grid_splined
       REAL*8, INTENT(IN),DIMENSION(response_func%worklist%n_columns, &
                                    grid_splined%n_points) :: splined_data
       INTEGER :: i_column
       CHARACTER(len=256) :: filename, filename_trimed

       do i_column = 1, response_func%worklist%n_columns

           filename = trim("splined")//"_v"//trim(num2str(i_vec)) &
                                     //"_"//trim(num2str(i_row)) &
                                     //"_"//trim(num2str(i_column)) &
                                     //"_o"//trim(num2str(omega))

           CALL remove_all_spaces(filename,filename_trimed)


           CALL write_two_arrays_to_file(grid_splined%n_points, &
                                         grid_splined%points(:), &
                                         splined_data(i_column,:), &
                                         filename_trimed)
       enddo
    END SUBROUTINE write_response_function_splined


    SUBROUTINE finalize_response_func_my_k_points()
!        DEALLOCATE(response_func%one_omega_my_k_points)
    END SUBROUTINE finalize_response_func_my_k_points

    SUBROUTINE finalize_response_function_storage()
        if(is_response_func_storage_inited) then
            DEALLOCATE(response_func%element)
            DEALLOCATE(response_func%buffer_ft)
            DEALLOCATE(response_func%one_omega_in_k)

            if(ALLOCATED(response_func%one_omega_my_k_points)) then
                DEALLOCATE(response_func%one_omega_my_k_points)
            endif

            is_response_func_storage_inited = .FALSE.
        endif
    END SUBROUTINE finalize_response_function_storage

END MODULE cRPA_storage_responsefunc_R
