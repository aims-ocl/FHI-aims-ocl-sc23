MODULE cRPA_calculation_fouriertrans_test
    USE cRPA_calculation_fouriertrans
    implicit none

CONTAINS
    SUBROUTINE discrete_fouriertransform_test()
        COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: data_R, &
                                                   data_k,&
                                                   forward_backward
        INTEGER :: i,n, n_elements
        type(realspace_vectors) :: rvecs

        CALL write_stdout("Testing: discrete_fouriertransform (R->k->R)")

        CALL create_all_realspace_vectors(n_cells_bvk, cell_index_bvk, &
                                          lattice_vector, rvecs)

        if(n_k_points /= rvecs%n_vecs) then
            stop "discrete_fouriertransform test failed! Unequal number of k and R vecs"
        endif

        n_elements = 100

        ALLOCATE(data_R(n_elements,rvecs%n_vecs))
        ALLOCATE(data_k(n_elements,n_k_points))
        ALLOCATE(forward_backward(n_elements,rvecs%n_vecs))
        forward_backward(:,:) = DCMPLX(0.d0, 0.d0)

        CALL fill_test_data(n_elements,rvecs%n_vecs,data_R)

        CALL discrete_fouriertransform_R_to_k(rvecs,n_elements, data_R, data_k)
        CALL discrete_fouriertransform_k_to_R(rvecs,n_elements, data_k, forward_backward)

        if(.NOT. is_data_equal(n_elements, rvecs%n_vecs, data_R, forward_backward)) then
            write(use_unit,*) "discrete_fouriertransform test failed! forward backward unequal"
            do n=1,n_elements
                do i = 1, rvecs%n_vecs
                    write(use_unit,*) n,i,data_R(n,i), data_R(n,i)-forward_backward(n,i)
                enddo
            enddo
            stop "discrete_fouriertransform test failed! forward backward unequal"
        endif

        CALL write_stdout("Testing: discrete_fouriertransform (k->R->k)")
        data_k(:,:) = DCMPLX(0.d0,0.d0)
        forward_backward(:,:) = DCMPLX(0.d0, 0.d0)
        CALL fill_test_data(n_elements, n_k_points,data_k)

        CALL discrete_fouriertransform_k_to_R(rvecs,n_elements, data_k, data_R)
        CALL discrete_fouriertransform_R_to_k(rvecs,n_elements, data_R, forward_backward)


        if(.NOT. is_data_equal(n_elements, n_k_points, data_k, forward_backward)) then
            write(use_unit,*) "discrete_fouriertransform test failed! backward forward unequal"
            do n=1,n_elements
                do i = 1, rvecs%n_vecs
                    write(use_unit,*) n,i,data_k(n,i), forward_backward(n,i)
                enddo
            enddo
            stop "discrete_fouriertransform test failed! backward forward unequal"
        endif


        DEALLOCATE(data_R)
        DEALLOCATE(data_k)
        DEALLOCATE(forward_backward)

        do n=1,rvecs%n_vecs
            do i=1,n_k_points
                if(myid==0) write(use_unit,*) "k_phase_exx", n,i, k_phase_exx(n,i)
            enddo
        enddo


        do i=1,n_k_points
            CALL write_stdout("k_weights"//num2str(i)//num2str(k_weights(i)))
        enddo

        CALL write_stdout("Test discrete_fouriertransform: ok!")
    CONTAINS
        SUBROUTINE fill_test_data(n_elements,n_vecs,data_R)
            INTEGER, INTENT(IN) :: n_elements
            INTEGER, INTENT(IN) :: n_vecs
            COMPLEX*16, DIMENSION(n_elements,n_vecs) :: data_R
            INTEGER :: i,n

            do n=1, n_elements
                do i = 1,n_vecs
                    data_R(n,i) = DCMPLX(real(n)*real(i),real(n/2)*real(-i))
                enddo
            enddo
        END SUBROUTINE fill_test_data

        LOGICAL FUNCTION is_data_equal(n_elements,n_vecs,data_R, forward_backward)
            INTEGER, INTENT(IN) :: n_elements,n_vecs
            COMPLEX*16, DIMENSION(n_elements,n_vecs) :: data_R, forward_backward
            REAL*8 :: max_diff
            INTEGER :: n

            is_data_equal = .FALSE.

            max_diff = 0.d0
            do n=1,n_elements
                max_diff =max(max_diff, &
                              maxval(abs(data_R(n,:)-forward_backward(n,:))))
            enddo

            write(use_unit,*) "Max diff is ",max_diff

            if(max_diff < 1.d-8) then
                is_data_equal = .TRUE.
            endif
        END FUNCTION is_data_equal

    END SUBROUTINE discrete_fouriertransform_test

!notice this is not a "stop on failure" test
    SUBROUTINE fouriertransform_test()
        REAL*8 :: interval_start, interval_end
        INTEGER :: n_partitions, n_order
        type(integration_intervals) :: int_intervals
        type(integration_grid) :: int_grid_forward, int_grid_backward
        INTEGER :: i
        REAL*8 :: error, &
                  ft_diff

        REAL*8, ALLOCATABLE, DIMENSION(:,:) :: f_of_x
        COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: f_of_t, f_of_tx

        interval_start = 0.d0
        interval_end = 2.d0
        n_order = 2
        n_partitions = 500 !int(ceiling(abs(interval_end - interval_start) / 20000.d0))

        CALL create_integration_intervals_uniform(interval_start, &
                                                  interval_end, n_partitions, &
                                                  int_intervals)

        CALL create_integration_grid_gauss_legendre(int_intervals,n_order, &
                                                    int_grid_forward)

        CALL free_integration_intervals(int_intervals)

        CALL create_integration_intervals_uniform(interval_start*200, &
                                                  interval_end*200, n_partitions*100, &
                                                  int_intervals)


        CALL create_integration_grid_gauss_legendre(int_intervals,n_order, &
                                                    int_grid_backward)

        CALL free_integration_intervals(int_intervals)

        ALLOCATE(f_of_x(int_grid_forward%n_points,1))
        ALLOCATE(f_of_tx(int_grid_forward%n_points,1))
        ALLOCATE(f_of_t(int_grid_backward%n_points,1))

        do i=1, int_grid_forward%n_points
            f_of_x(i,1) = exp(-abs(int_grid_forward%points(i)))
        enddo

       !     CALL print_integration_grid(int_grid)
        write(use_unit,*) "FT forward"!, f_of_x(1:6)
        do i=1, int_grid_backward%n_points
            CALL fouriertransform_x_to_t(int_grid_forward,1, f_of_x, &
                                         int_grid_backward%points(i), f_of_t(i,1))
        enddo

        write(use_unit,*) "FT backward"!,f_of_t(1:5)
        do i=1, int_grid_forward%n_points
           CALL fouriertransform_t_to_x(int_grid_backward, &
                                        1, f_of_t, &
                                        int_grid_forward%points(i), f_of_tx(i,1))
        enddo

        write(use_unit,*) "FT diff"!, f_of_x(1:6)
        ft_diff = 0.d0
        do i=1, int_grid_forward%n_points
             ft_diff = max(ft_diff, &
                           abs(f_of_x(i,1) - dble(f_of_tx(i,1))/abs(f_of_x(i,1))))
             write(use_unit,*) f_of_x(i,1) , dble(f_of_tx(i,1))
        enddo

        write(use_unit,*) "   FT forward, backward diff:", ft_diff
        write(use_unit,*) "Is it small?!"


        DEALLOCATE(f_of_x)
        DEALLOCATE(f_of_t)
        DEALLOCATE(f_of_tx)
        CALL free_integration_grid(int_grid_forward)
        CALL free_integration_grid(int_grid_backward)
    END SUBROUTINE fouriertransform_test

end MODULE cRPA_calculation_fouriertrans_test
