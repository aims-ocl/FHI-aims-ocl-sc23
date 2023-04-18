MODULE cRPA_calculation_integration_test
    USE cRPA_calculation_integration
    USE cRPA_calculation_fouriertrans
    implicit none

CONTAINS
    SUBROUTINE integration_gauss_legendre_test()
        REAL*8 :: interval_start, interval_end
        INTEGER :: n_partitions, n_order
        type(integration_intervals) :: int_intervals
        type(integration_grid) :: int_grid
        INTEGER :: i,i_mode
        REAL*8 :: integral_num, integral_analy, error, &
                  integral_num_old, error_old,ft_diff
        REAL*8, ALLOCATABLE, DIMENSION(:) :: x,w

        REAL*8, ALLOCATABLE, DIMENSION(:,:) :: f_of_x
        COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: f_of_t

        interval_start = -3.d0
        interval_end = 3.d0
        write(use_unit,*) "interval_start/end", interval_start, interval_end

        integral_analy = (1.d0-exp(-interval_end))+(1.d0-exp(-abs(interval_start)))
                            !exp(interval_end)- exp(interval_start)


        do n_order = 2,4
             do i_mode = 1,2
             n_partitions = 1
             error = 10.d0
             error_old = 10.d0
             do while( (abs(error) > 1.d-8 &
                        .OR. n_partitions*n_order < 100) &
                        .AND. n_partitions*n_order < 1000)
                if(i_mode == 1) then
                    CALL create_integration_intervals_uniform(interval_start, &
                                                              interval_end, n_partitions, &
                                                              int_intervals)
                endif

                if(i_mode == 2) then
                    CALL create_integration_intervals_logarithmic(interval_start,interval_end, &
                                                              n_partitions, int_intervals)
                endif
!                CALL print_inte<gration_intervals(int_intervals)


                CALL create_integration_grid_gauss_legendre(int_intervals,n_order, &
                                                            int_grid)

!                CALL print_integration_intervals(int_intervals)


                CALL free_integration_intervals(int_intervals)

!                CALL print_integration_grid(int_grid)


                integral_num = 0.d0
                do i=1, int_grid%n_points
                    integral_num = integral_num + &
                                   exp(-abs(int_grid%points(i))) * int_grid%weights(i)
                enddo

                error = integral_analy-integral_num

                CALL free_integration_grid(int_grid)

!old method
                ALLOCATE(x(-n_order*n_partitions/2:n_order*n_partitions/2))
                ALLOCATE(w(-n_order*n_partitions/2:n_order*n_partitions/2))

                CALL get_tau_grid_and_weights_old(interval_start,interval_end, &
                                                  n_order*n_partitions/2,x,w)

                integral_num_old = 0.d0
                do i=-n_order*n_partitions/2, n_order*n_partitions/2
                    integral_num_old = integral_num_old + &
                                        exp(-abs(x(i))) * w(i)
                enddo

                error_old = integral_analy - integral_num_old

                DEALLOCATE(x)
                DEALLOCATE(w)


                n_partitions = n_partitions +1
            enddo

            if(i_mode == 1) then
                write(use_unit,*) "uniform"
            endif

            if(i_mode == 2) then
                write(use_unit,*) "logarithmic"
            endif


            write(use_unit,*) "     order, n_partitions", n_order, n_partitions
            write(use_unit,*) "         Numerically ", integral_num, "analytically ", integral_analy
            write(use_unit,*) "         Difference ", error, error_old
            write(use_unit,*) ""

            CALL fouriertest()

            enddo
        enddo

    CONTAINS
        SUBROUTINE fouriertest()

            if(i_mode == 1) then
                CALL create_integration_intervals_uniform(interval_start, &
                    interval_end, n_partitions, &
                    int_intervals)
            endif

            if(i_mode == 2) then
                CALL create_integration_intervals_logarithmic(interval_start,interval_end, &
                    n_partitions, int_intervals)
            endif

            CALL create_integration_grid_gauss_legendre(int_intervals,n_order, &
                                                            int_grid)


            ALLOCATE(f_of_t(int_grid%n_points,1))
            ALLOCATE(f_of_x(int_grid%n_points,1))

            do i=1, int_grid%n_points
                f_of_x(i,1) = exp(-abs(int_grid%points(i)))!exp(-int_grid%points(i)**2/5.d0)
            enddo

       !     CALL print_integration_grid(int_grid)
            write(use_unit,*) "FT forward"!, f_of_x(1:6)
            do i=1, int_grid%n_points
                CALL fouriertransform_x_to_t(int_grid,1, f_of_x, &
                                             int_grid%points(i), f_of_t(i,1))
            enddo
            f_of_x(:,:) = 0.d0

            write(use_unit,*) "FT backward"!,f_of_t(1:5)
            do i=1, int_grid%n_points
!change to complex res value
!                CALL fouriertransform_t_to_x(int_grid, &
!                                             1, f_of_t, &
!                                             int_grid%points(i), f_of_x(i,1))
            enddo

            write(use_unit,*) "FT diff"!, f_of_x(1:6)
            ft_diff = 0.d0
            do i=1, int_grid%n_points
                ft_diff = ft_diff + &
                          abs(f_of_x(i,1) - exp(-abs(int_grid%points(i))))!exp(-int_grid%points(i)**2/5.d0))
            enddo


            write(use_unit,*) "   FT forward, backward diff:", ft_diff / sum(abs(exp(-abs(int_grid%points(:)))))
            write(use_unit,*) ""
            write(use_unit,*) ""
            write(use_unit,*) ""
            write(use_unit,*) ""


            DEALLOCATE(f_of_x)
            DEALLOCATE(f_of_t)
            CALL free_integration_grid(int_grid)
        END SUBROUTINE fouriertest

    END SUBROUTINE integration_gauss_legendre_test


    SUBROUTINE get_tau_grid_and_weights_old(x1,x2,n_points,x,weight)
        REAL*8, intent(in)  :: x1,x2       ! ranges
        INTEGER,          intent(in)  :: n_points    ! number of points
        REAL*8, intent(out) :: x(-n_points:n_points), &
                               weight(-n_points:n_points) ! abcsissa and weights
        ! internal vars
        real*8 h, a, t_0
        !real*8 t_0
        integer :: i,j,m

          t_0 = 0.01
          h = 1./real(n_points)*log ((x2-x1)/t_0)
          a = 1.05

          do i=1, n_points, 1
            x(i) = t_0 * (exp(i*h)- 1.)
            x(-i) = -x(i)
            weight(i) = h * t_0 * exp(i*h)
            weight(-i) = weight(i)
          enddo

          x(0)=0
          weight(0)=t_0 * h

          return
    END SUBROUTINE get_tau_grid_and_weights_old

END MODULE cRPA_calculation_integration_test
