!Calculates integrals numerically - mainly with quadratues
MODULE cRPA_calculation_integration
    USE cRPA_view
    USE cRPA_flow_adaptive_grid
    USE spline
    implicit none

!a connected interval
    type integration_interval
        REAL*8 :: start
        REAL*8 :: end
        INTEGER :: n_points
        REAL*8, ALLOCATABLE, DIMENSION(:) :: points
        REAL*8, ALLOCATABLE, DIMENSION(:) :: weights
    end type integration_interval

!composition of intervals
    type integration_intervals
        REAL*8 :: total_interval_start
        REAL*8 :: total_interval_end
        INTEGER :: n_partitions

        type(integration_interval), ALLOCATABLE, DIMENSION(:) :: interval
    end type

    type integration_grid
        INTEGER :: n_points
        REAL*8, ALLOCATABLE, DIMENSION(:) :: points
        REAL*8, ALLOCATABLE, DIMENSION(:) :: weights
        COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: weights_complex
    end type integration_grid

    !order value
    REAL*8,DIMENSION(4,4) :: GAUSS_LEGENDRE_ABSCISSAS, &
                             GAUSS_LEGENDRE_WEIGHTS

    PARAMETER(GAUSS_LEGENDRE_ABSCISSAS = reshape( (/ &
              0.d0, -sqrt(1.d0/3.d0), -sqrt(3.d0/5.d0), -sqrt(3.d0/7.d0 + (2.d0/7.d0)*sqrt(6.d0/5.d0)), &
              0.d0,  sqrt(1.d0/3.d0),  0.d0           , -sqrt(3.d0/7.d0 - (2.d0/7.d0)*sqrt(6.d0/5.d0)), &
              0.d0,  0.d0           ,  sqrt(3.d0/5.d0), sqrt(3.d0/7.d0 - (2.d0/7.d0)*sqrt(6.d0/5.d0)), &
              0.d0,  0.d0           ,  0.d0          ,  sqrt(3.d0/7.d0 + (2.d0/7.d0)*sqrt(6.d0/5.d0))/), &
              (/4,4/)))

    PARAMETER(GAUSS_LEGENDRE_WEIGHTS = reshape( (/ &
              2.d0,  1.d0, 5.d0/9.d0, (18.d0-sqrt(30.d0))/36.d0, &
              0.d0,  1.d0, 8.d0/9.d0, (18.d0+sqrt(30.d0))/36.d0, &
              0.d0,  0.d0, 5.d0/9.d0, (18.d0+sqrt(30.d0))/36.d0, &
              0.d0,  0.d0, 0.d0     , (18.d0-sqrt(30.d0))/36.d0/), &
              (/4,4/)))

!GAUSS HERMITE
!abscissa, weight, total_weight
!1   0   1.77245385091   1.77245385091
!
!1   -0.707106781187 0.886226925453  1.46114118266
!2   0.707106781187  0.886226925453  1.46114118266
!
!1 -1.22474487139  0.295408975151  1.32393117521
!2   0   1.1816359006    1.1816359006
!3   1.22474487139   0.295408975151  1.32393117521
!
!1   -1.65068012389  0.0813128354472 1.2402258177
!2   -0.524647623275 0.804914090006  1.05996448289
!3   0.524647623275  0.804914090006  1.05996448289
!4   1.65068012389   0.0813128354472 1.2402258177

CONTAINS

    SUBROUTINE create_integration_intervals_uniform_unsymmetric(interval_start, interval_end, &
                                                                n_partitions, int_intervals)
        REAL*8, INTENT(IN) :: interval_start, &
                              interval_end
        INTEGER, INTENT(IN) :: n_partitions
        type(integration_intervals), INTENT(OUT) :: int_intervals


        INTEGER :: n_partitions_positiv, n_partitions_negativ
        type(integration_intervals) :: log_positiv, log_negativ

        n_partitions_positiv = n_partitions / 2
        n_partitions_positiv = max(n_partitions_positiv,1)
        n_partitions_negativ = n_partitions_positiv

        CALL create_integration_intervals_uniform(0.d0,interval_end, &
                                                  n_partitions_positiv, &
                                                  log_positiv)

        CALL create_integration_intervals_uniform(interval_start,0.d0, &
                                                  n_partitions_negativ, &
                                                  log_negativ)

        CALL merge_integration_intervals(log_negativ,log_positiv, int_intervals)


    END SUBROUTINE create_integration_intervals_uniform_unsymmetric

    SUBROUTINE create_integration_intervals_uniform(interval_start, interval_end, &
                                                    n_partitions, int_intervals)
        REAL*8, INTENT(IN) :: interval_start, &
                              interval_end
        INTEGER, INTENT(IN) :: n_partitions
        type(integration_intervals), INTENT(OUT) :: int_intervals

        INTEGER :: i
        REAL*8 :: int_width
        REAL*8, DIMENSION(n_partitions) :: points_start, &
                                           points_end

        int_width = (interval_end - interval_start)/real(n_partitions)
        points_start(1) = interval_start

        do i=1, n_partitions-1
            points_end(i) = interval_start + real(i) * int_width
            points_start(i+1) = points_end(i)
        enddo

        points_end(n_partitions) = interval_end

        CALL create_integration_intervals(interval_start, interval_end, &
                                          n_partitions, &
                                          points_start, points_end, &
                                          int_intervals)


    END SUBROUTINE create_integration_intervals_uniform

    SUBROUTINE create_integration_intervals_logarithmic(interval_start, &
                                                        interval_end, &
                                                        n_partitions, &
                                                        int_intervals)
        REAL*8, INTENT(IN) :: interval_start, &
                              interval_end
        INTEGER, INTENT(IN) :: n_partitions
        type(integration_intervals), INTENT(OUT) :: int_intervals

        INTEGER :: n_partitions_positiv, n_partitions_negativ
        type(integration_intervals) :: log_positiv, log_negativ

        n_partitions_positiv = n_partitions ! / 2
!        n_partitions_positiv = max(n_partitions_positiv,1)
!        n_partitions_negativ = n_partitions_positiv

        CALL create_integration_intervals_logarithmic_positiv(interval_end, &
                                                              n_partitions_positiv, &
                                                              int_intervals)!log_positiv)

!        CALL create_integration_intervals_logarithmic_positiv(abs(interval_start), &
!                                                              n_partitions_negativ, &
!                                                              log_negativ)
!        CALL mirror_integration_intervals(log_negativ)

!        CALL merge_integration_intervals(log_negativ,log_positiv, int_intervals)


    END SUBROUTINE create_integration_intervals_logarithmic

    SUBROUTINE create_integration_intervals_logarithmic_positiv(interval_end, &
                                                                n_partitions, &
                                                                int_intervals)
        REAL*8, INTENT(IN) :: interval_end
        INTEGER, INTENT(IN) :: n_partitions
        type(integration_intervals), INTENT(OUT) :: int_intervals

        REAL*8, DIMENSION(n_partitions) :: points_start, &
                                           points_end
        REAL*8 h, t_0
        INTEGER :: i

          t_0 = 0.15d0
          h = 1./real(n_partitions)*log ((interval_end)/t_0)

          points_start(1) = 0.d0

          do i=1, n_partitions - 1
            points_end(i) = t_0 * (exp(i*h)- 1.)
            points_start(i+1) = points_end(i)
          enddo

          points_end(n_partitions)= interval_end

          CALL create_integration_intervals(0.d0, interval_end, &
                                            n_partitions, &
                                            points_start, points_end, &
                                            int_intervals)

    END SUBROUTINE create_integration_intervals_logarithmic_positiv

    SUBROUTINE merge_integration_intervals(intervals1,intervals2, intervals_merged)
        type(integration_intervals), INTENT(IN) :: intervals1, intervals2
        type(integration_intervals), INTENT(OUT) :: intervals_merged

        INTEGER :: n_partitions
        REAL*8 :: interval_start, interval_end
        REAL*8, ALLOCATABLE, DIMENSION(:) :: points_start, points_end

        n_partitions = intervals1%n_partitions + intervals2%n_partitions

        interval_start = intervals1%total_interval_start
        interval_end =  intervals2%total_interval_end

!TODO check alloc
        ALLOCATE(points_start(n_partitions))
        ALLOCATE(points_end(n_partitions))


        points_start(1:intervals1%n_partitions) = intervals1%interval(:)%start
        points_start(intervals1%n_partitions+1:n_partitions) = intervals2%interval(:)%start

        points_end(1:intervals1%n_partitions) = intervals1%interval(:)%end
        points_end(intervals1%n_partitions+1:n_partitions) = intervals2%interval(:)%end


        CALL create_integration_intervals(interval_start, interval_end, &
                                          n_partitions, &
                                          points_start, points_end, &
                                          intervals_merged)

        DEALLOCATE(points_start)
        DEALLOCATE(points_end)

    END SUBROUTINE merge_integration_intervals

    SUBROUTINE mirror_integration_intervals(int_intervals)
        type(integration_intervals), INTENT(INOUT) :: int_intervals

        INTEGER :: i
        REAL*8 :: buffer

        buffer = int_intervals%total_interval_start
        int_intervals%total_interval_start = -int_intervals%total_interval_end
        int_intervals%total_interval_end   = -buffer

        do i = 1,int_intervals%n_partitions
            CALL mirror_integration_interval(int_intervals%interval(i))
        enddo
    END SUBROUTINE mirror_integration_intervals

    SUBROUTINE create_integration_grid_gauss_legendre(int_intervals,order, &
                                                     int_grid)
        type(integration_intervals), INTENT(INOUT) :: int_intervals
        INTEGER, INTENT(IN) :: order
        type(integration_grid), INTENT(OUT) :: int_grid

        INTEGER :: i, n_points_per_partion
        REAL*8, ALLOCATABLE, DIMENSION(:) :: points_per_partition, &
                                             weights_per_partition
        type(integration_interval) :: unit_interval

        if(order > 4) then
            stop "only order up to 4 is supported for gauss hermite integration"
        endif

        CALL create_integration_interval_unit(unit_interval)

        n_points_per_partion = order
!TODO check alloc
        ALLOCATE(points_per_partition(n_points_per_partion))
        ALLOCATE(weights_per_partition(n_points_per_partion))

        points_per_partition = GAUSS_LEGENDRE_ABSCISSAS(order,1:n_points_per_partion)
        weights_per_partition = GAUSS_LEGENDRE_WEIGHTS(order,1:n_points_per_partion)

        CALL set_integration_interval_points(n_points_per_partion, &
                                             points_per_partition, &
                                             weights_per_partition, &
                                             unit_interval)


        do i=1,int_intervals%n_partitions
            CALL linear_transform_unit_interval_to_interval(unit_interval, &
                                                            int_intervals%interval(i))
        enddo

        DEALLOCATE(points_per_partition)
        DEALLOCATE(weights_per_partition)
        CALL free_integration_interval(unit_interval)

        CALL create_integration_grid_from_intervals(int_intervals,int_grid)

    END SUBROUTINE create_integration_grid_gauss_legendre

    SUBROUTINE create_integration_grid_trapezodial(int_intervals, &
                                                     int_grid)
        type(integration_intervals), INTENT(INOUT) :: int_intervals
        type(integration_grid), INTENT(OUT) :: int_grid

        REAL*8, DIMENSION(int_intervals%n_partitions * 2) :: points, &
                                                             weights
        INTEGER :: i_point, n_points
        REAL*8 :: start, stepwidth

        n_points = int_intervals%n_partitions * 2
        start = int_intervals%total_interval_start
        stepwidth = (int_intervals%total_interval_end - &
                     int_intervals%total_interval_start) &
                    / &
                    dble(n_points)

        do i_point = 1,n_points
            points(i_point) = start + (i_point-1) * stepwidth
            if(i_point == 1 .OR. i_point == n_points) then
                weights(i_point) = 0.5d0
            else
                weights(i_point) = 1.d0
            endif
        enddo

        weights(:) = weights(:) * stepwidth

        CALL create_integration_grid(n_points, points, weights, &
                                     int_grid)

    END SUBROUTINE create_integration_grid_trapezodial

    SUBROUTINE create_integration_grid_from_intervals(int_intervals, int_grid)
        type(integration_intervals), INTENT(IN) :: int_intervals
        type(integration_grid), INTENT(OUT) :: int_grid

        INTEGER :: n_points, i,offset, n_points_per_partition
        REAL*8, ALLOCATABLE, DIMENSION(:) :: points, weights

        n_points = sum(int_intervals%interval(:)%n_points)

!TODO check alloc
        ALLOCATE(points(n_points))
        ALLOCATE(weights(n_points))

        offset = 1
        do i=1, int_intervals%n_partitions
            n_points_per_partition = int_intervals%interval(i)%n_points
            points(offset:offset+n_points_per_partition-1) = int_intervals%interval(i)%points(:)
            weights(offset:offset+n_points_per_partition-1) = int_intervals%interval(i)%weights(:)

            offset = offset + n_points_per_partition
        enddo

        CALL create_integration_grid(n_points, points, weights, int_grid)

        DEALLOCATE(points)
        DEALLOCATE(weights)
    END SUBROUTINE create_integration_grid_from_intervals

    SUBROUTINE create_integration_grid(n_points, points, weights, int_grid)
        INTEGER, INTENT(IN) :: n_points
        REAL*8, DIMENSION(n_points), INTENT(IN) :: points, weights
        type(integration_grid), INTENT(OUT) :: int_grid

        int_grid%n_points = n_points
!TODO check alloc
        ALLOCATE(int_grid%points(n_points))
        ALLOCATE(int_grid%weights(n_points))
        ALLOCATE(int_grid%weights_complex(n_points))

        int_grid%points(:) = points(:)
        int_grid%weights(:) = weights(:)
        int_grid%weights_complex(:) = DCMPLX(weights(:),0.d0)

    END SUBROUTINE create_integration_grid

    SUBROUTINE create_integration_grid_stepfunction(int_start, int_end, &
                                                    n_partitions, &
                                                    int_grid)
       REAL*8, INTENT(IN) :: int_start, int_end
       INTEGER, INTENT(IN) :: n_partitions
       type(integration_grid), INTENT(OUT) :: int_grid

       REAL*8, DIMENSION(n_partitions+1) :: points, weights
       INTEGER :: i

       do i = 0, n_partitions
            points(i+1) = int_start + (real(i) * (int_end-int_start)/n_partitions)
       enddo
       weights(:) = 1.d0
       weights(n_partitions+1) = 0.d0

       CALL create_integration_grid(n_partitions+1, points, weights, &
                                    int_grid)

    END SUBROUTINE create_integration_grid_stepfunction

    SUBROUTINE create_integration_grid_from_adaptive_grid(grid_adaptive, grid_int)
        type(adaptive_grid), INTENT(IN) :: grid_adaptive
        type(integration_grid), INTENT(OUT) :: grid_int

        REAL*8, DIMENSION(grid_adaptive%n_points) :: weights

        weights(:) = 1.d0
        weights(grid_adaptive%n_points) = 0.d0

        CALL create_integration_grid(grid_adaptive%n_points, &
                                     grid_adaptive%points(:), &
                                     weights, &
                                     grid_int)

    END SUBROUTINE create_integration_grid_from_adaptive_grid

    SUBROUTINE merge_integration_grids(grid1,grid2, grid_merged)
        type(integration_grid), INTENT(IN) :: grid1, grid2
        type(integration_grid), INTENT(OUT) :: grid_merged

        REAL*8, DIMENSION(grid1%n_points+grid2%n_points) :: merged_points, &
                                                            merged_weights
        INTEGER :: n_points

        n_points = grid1%n_points+grid2%n_points

        merged_points(1:grid1%n_points) = grid1%points(:)
        merged_points(grid1%n_points+1:n_points) = grid2%points(:)

        merged_weights(1:grid1%n_points) = grid1%weights(:)
        merged_weights(grid1%n_points+1:n_points) = grid2%weights(:)

        CALL create_integration_grid(n_points,merged_points,merged_weights, &
                                     grid_merged)
    END SUBROUTINE merge_integration_grids

    SUBROUTINE print_integration_grid(int_grid)
        type(integration_grid),INTENT(IN) :: int_grid

        INTEGER :: i

        write(use_unit,*) "Integration grid -  points/weights"

        do i=1, int_grid%n_points
            write(use_unit,*) int_grid%points(i), int_grid%weights(i)
        enddo
    END SUBROUTINE print_integration_grid

    SUBROUTINE free_integration_grid(int_grid)
        type(integration_grid), INTENT(INOUT) :: int_grid

        DEALLOCATE(int_grid%points)
        DEALLOCATE(int_grid%weights)
        DEALLOCATE(int_grid%weights_complex)

    END SUBROUTINE free_integration_grid

    SUBROUTINE create_integration_intervals(interval_start, &
                                            interval_end, &
                                            n_partitions, &
                                            points_start, &
                                            points_end, &
                                            int_intervals)
        REAL*8, INTENT(IN) :: interval_start, &
                              interval_end
        INTEGER, INTENT(IN) :: n_partitions
        REAL*8, DIMENSION(n_partitions),INTENT(IN) :: points_start, &
                                                      points_end
        type(integration_intervals), INTENT(OUT) :: int_intervals

        INTEGER :: i
        int_intervals%n_partitions = n_partitions
        int_intervals%total_interval_start = interval_start
        int_intervals%total_interval_end = interval_end

!TODO check allocation
        ALLOCATE(int_intervals%interval(n_partitions))

        do i=1, n_partitions
            CALL create_integration_interval(points_start(i), points_end(i), &
                                             int_intervals%interval(i))
        enddo

    END SUBROUTINE create_integration_intervals

    SUBROUTINE print_integration_intervals(int_intervals)
        type(integration_intervals), INTENT(IN) :: int_intervals
        INTEGER :: i

        write(use_unit,*) "Total interval start/stop", int_intervals%total_interval_start, &
                   int_intervals%total_interval_end
        do i=1, int_intervals%n_partitions
            write(use_unit,*) "Interval ",i
            CALL print_integration_interval(int_intervals%interval(i))
        enddo

    END SUBROUTINE print_integration_intervals

    SUBROUTINE free_integration_intervals(int_intervals)
        type(integration_intervals), INTENT(INOUT) :: int_intervals

        INTEGER :: i

        do i=1, int_intervals%n_partitions
            CALL free_integration_interval(int_intervals%interval(i))
        enddo

    END SUBROUTINE free_integration_intervals

    SUBROUTINE create_integration_interval_unit(an_unit_interval)
        type(integration_interval), INTENT(OUT) :: an_unit_interval

        CALL create_integration_interval(-1.d0,1.d0,an_unit_interval)
    END SUBROUTINE create_integration_interval_unit

    SUBROUTINE create_integration_interval(interval_start, interval_end, &
                                           an_interval)
        REAL*8,INTENT(IN) :: interval_start,interval_end
        type(integration_interval), INTENT(OUT) :: an_interval

        an_interval%start = interval_start
        an_interval%end = interval_end
        an_interval%n_points = 0

    END SUBROUTINE create_integration_interval

    SUBROUTINE mirror_integration_interval(int_interval)
        type(integration_interval), INTENT(INOUT) :: int_interval

        INTEGER :: i
        REAL *8 :: buffer

        buffer = int_interval%start
        int_interval%start = -int_interval%end
        int_interval%end = -buffer

    END SUBROUTINE mirror_integration_interval


    SUBROUTINE set_integration_interval_points(n_points,points, weights, &
                                               an_interval)
        INTEGER, INTENT(IN) :: n_points
        REAL*8, DIMENSION(n_points), INTENT(IN) :: points, weights
        type(integration_interval), INTENT(INOUT) :: an_interval

        if(an_interval%n_points/= n_points) then
            if(ALLOCATED(an_interval%points)) then
                DEALLOCATE(an_interval%points)
            endif

            if(ALLOCATED(an_interval%weights)) then
                DEALLOCATE(an_interval%weights)
            endif


!TODO check alloc
            ALLOCATE(an_interval%points(n_points))
            ALLOCATE(an_interval%weights(n_points))

            an_interval%n_points = n_points
        endif

        an_interval%points(:) = points(:)
        an_interval%weights(:) = weights(:)

    END SUBROUTINE set_integration_interval_points

    SUBROUTINE print_integration_interval(an_interval)
        type(integration_interval) :: an_interval

        write(use_unit,*) "Start ",an_interval%start," end ",an_interval%end, &
                   "points ", an_interval%n_points, an_interval%points
    END SUBROUTINE print_integration_interval

    SUBROUTINE free_integration_interval(an_interval)
        type(integration_interval), INTENT(INOUT) :: an_interval

        if(ALLOCATED(an_interval%points)) then
            DEALLOCATE(an_interval%points)
        endif

        if(ALLOCATED(an_interval%weights)) then
            DEALLOCATE(an_interval%weights)
        endif

    END SUBROUTINE free_integration_interval


    SUBROUTINE linear_transform_unit_interval_to_interval(unit_interval, &
                                                          transformed_interval)
        type(integration_interval), INTENT(IN) :: unit_interval
        type(integration_interval), INTENT(INOUT) :: transformed_interval

        INTEGER :: i
        REAL*8, DIMENSION(unit_interval%n_points) :: transformed_points, &
                                                     transformed_weights

        if(unit_interval%start /= -1 &
            .OR. &
           unit_interval%end /= 1) then
            stop "your input is not a unit interval!"
        endif

        transformed_points(:) = ( (transformed_interval%end - &
                                   transformed_interval%start) / 2.d0 ) &
                                * unit_interval%points(:)

        do i=1, unit_interval%n_points
            transformed_points(i) = transformed_points(i) + &
                                    ((transformed_interval%start + &
                                     transformed_interval%end) / 2.d0 )
        enddo


        transformed_weights(:) = ( (transformed_interval%end - &
                                   transformed_interval%start) / 2.d0 ) &
                                  * unit_interval%weights(:)


        CALL set_integration_interval_points(unit_interval%n_points, &
                                             transformed_points, &
                                             transformed_weights, &
                                             transformed_interval)
    END SUBROUTINE linear_transform_unit_interval_to_interval


    SUBROUTINE create_grid_old(x1,x2,n_points,int_grid)
        REAL*8, intent(in)  :: x1,x2       ! ranges
        INTEGER,          intent(in)  :: n_points    ! number of points
        type(integration_grid), INTENT(OUT) :: int_grid

        REAL*8 :: x(2*n_points+1), &
                  weight(2*n_points+1) ! abcsissa and weights
        ! internal vars
        real*8 h, a, t_0
        !real*8 t_0
        integer :: i,j,m

          t_0 = 0.00005
!          h = 1./real(n_points)*log ((x2-x1)/t_0)
          h = 1./real(n_points)*log ((x2)/t_0)

          do i=1, n_points, 1
            j = i+n_points + 1
            x(j) = t_0 * (exp(i*h)- 1.)
            x(n_points-i+1) = -x(j)
            weight(j) = h * t_0 * exp(i*h)
            weight(n_points-i+1) = weight(j)
          enddo

          x(n_points+1)=0
          weight(n_points+1)=t_0 * h

          h = 1./real(n_points)*log ((abs(x1))/t_0)

          do i=1, n_points, 1
            x(n_points-i+1) = -t_0 * (exp(i*h)- 1.)
            weight(n_points-i+1) = weight(j)
          enddo


          CALL create_integration_grid(2*n_points+1,x,weight,int_grid)
    END SUBROUTINE create_grid_old


    SUBROUTINE spline_extrapolate(tau_grid, &
                                  n_elements, f_data, &
                                  max_dist_between_points, &
                                  grid_splined, &
                                  splined_data)

    type(integration_grid), INTENT(IN) :: tau_grid
    INTEGER, INTENT(IN) :: n_elements
    REAL*8, DIMENSION(n_elements,tau_grid%n_points), INTENT(IN) :: f_data
    REAL*8, INTENT(IN) :: max_dist_between_points

    type(integration_grid), INTENT(OUT) :: grid_splined
    REAL*8, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: splined_data


    type(integration_intervals) :: int_intervals
    type(integration_grid) :: an_extra_grid
    real*8 spl_param (4, tau_grid%n_points)
    REAL*8, ALLOCATABLE,DIMENSION(:,:,:) :: spline_coeffs
    REAL*8 :: cur_dist_between_two_points, t_rel
    INTEGER :: i,j, i_tau, &
               n_order, n_partitions,  &
               i_splined_point, i_spline_pos, &
               n_splined_data
    INTEGER, DIMENSION(tau_grid%n_points-1) :: n_extra_partitions

    INTEGER, ALLOCATABLE, DIMENSION(:) :: i_spline_grid_to_tau
    REAL*8, ALLOCATABLE, DIMENSION(:) :: spline_points, spline_weights

    REAL*8, DIMENSION(n_elements) :: f_data_signs

!TODO check allocation
    n_order = 3

    n_splined_data = 0
    do i_tau = 1, tau_grid%n_points-1
        cur_dist_between_two_points = abs(tau_grid%points(i_tau) - &
                                                  tau_grid%points(i_tau+1))

        n_extra_partitions(i_tau) = int(ceiling((cur_dist_between_two_points / max_dist_between_points) &
                             / n_order))

        n_extra_partitions(i_tau) = max(1, n_extra_partitions(i_tau))

        if(n_extra_partitions(i_tau) > 0) then
            n_splined_data = n_splined_data + n_extra_partitions(i_tau)*n_order
        else
            n_splined_data = n_splined_data + 1
        endif
    enddo


!    CALL write_stdout("Points initially, splined:"//num2str(tau_grid%n_points)&
!                                                  //num2str(n_splined_data))

    ALLOCATE(i_spline_grid_to_tau(n_splined_data))
    ALLOCATE(spline_points(n_splined_data))
    ALLOCATE(spline_weights(n_splined_data))
    ALLOCATE(splined_data(n_elements,n_splined_data))
    ALLOCATE(spline_coeffs(4,n_elements,tau_grid%n_points))

    i_spline_pos = 1
    do i_tau = 1, tau_grid%n_points - 1

        if(n_extra_partitions(i_tau) > 0) then
            CALL create_integration_intervals_uniform(tau_grid%points(i_tau), &
                                                      tau_grid%points(i_tau+1), &
                                                      n_extra_partitions(i_tau), &
                                                      int_intervals)

            CALL create_integration_grid_gauss_legendre(int_intervals, n_order, &
                                                        an_extra_grid)

            spline_points(i_spline_pos:i_spline_pos + an_extra_grid%n_points-1) = &
                                                                    an_extra_grid%points(:)

            spline_weights(i_spline_pos:i_spline_pos + an_extra_grid%n_points-1) = &
                                                                    an_extra_grid%weights(:)

            i_spline_grid_to_tau(i_spline_pos:i_spline_pos + an_extra_grid%n_points-1) = &
                                                                                i_tau

            i_spline_pos = i_spline_pos + an_extra_grid%n_points


            CALL free_integration_intervals(int_intervals)
            CALL free_integration_grid(an_extra_grid)
        else
            spline_points(i_spline_pos) = tau_grid%points(i_tau)
            spline_weights(i_spline_pos) = abs(tau_grid%points(i_tau)-tau_grid%points(i_tau+1))
            i_spline_grid_to_tau(i_spline_pos) = i_tau

            i_spline_pos = i_spline_pos + 1
        endif
    enddo

    CALL create_integration_grid(n_splined_data, &
                                 spline_points, &
                                 spline_weights, &
                                 grid_splined)

    DEALLOCATE(spline_points)
    DEALLOCATE(spline_weights)

    do i=1, n_elements
       f_data_signs(i) = sign(1.d0, f_data(i,1))
       CALL cubic_spline(log(abs(f_data(i,:))), tau_grid%n_points, spline_coeffs(:,i,:) )

       do j = 1, tau_grid%n_points-1
!          if(abs(f_data(i,j)) <1.d-14) then
!            spline_coeffs(:,i,j) = 0.d0
!          endif

!          if( abs(f_data(i,j)-f_data(i,j+1)) <1.d-20) then
!write(use_unit,*) tau_grid%points(j),tau_grid%points(j+1),abs(f_data(i,j)-f_data(i,j+1))

!            spline_coeffs(:,i,j) = 0.d0
!          endif


!          if(tau_grid%points(j)*tau_grid%points(j+1) <0.d0 &
!             .OR. &
!             abs(tau_grid%points(j)*tau_grid%points(j+1)) < 1.d-30) then
!write(use_unit,*) tau_grid%points(j),tau_grid%points(j+1)
!            spline_coeffs(:,i,j) = 0.d0
!          endif

       enddo

    enddo

    do i_splined_point = 1, grid_splined%n_points
        i_tau = i_spline_grid_to_tau(i_splined_point)

        if(n_extra_partitions(i_tau) > 0) then
            t_rel = get_relative_point(i_splined_point, i_tau)
            CALL spline_value(t_rel, i_tau,splined_data(:,i_splined_point))
        else
            splined_data(:,i_splined_point) = f_data(:,i_tau)
        endif
    enddo

    do i=1, n_elements
        splined_data(i,:) = exp(splined_data(i,:)) * f_data_signs(i)
    enddo

    DEALLOCATE(i_spline_grid_to_tau)
    DEALLOCATE(spline_coeffs)

    CONTAINS
        REAL*8 FUNCTION get_relative_point(i_splined_point,i_tau)
            INTEGER, INTENT(IN) :: i_splined_point, i_tau
            get_relative_point = &
                    (grid_splined%points(i_splined_point) - tau_grid%points(i_tau)) &
                        / &
                    (tau_grid%points(i_tau+1) - tau_grid%points(i_tau))

        END FUNCTION get_relative_point

        SUBROUTINE spline_value(t,i_x, value)
            REAL*8, INTENT(IN) :: t
            INTEGER, INTENT(IN) :: i_x
            REAL*8,DIMENSION(n_elements), INTENT(OUT) :: value
            REAL*8, DIMENSION(4) :: t_vec

            t_vec(:) = (/1.d0, t, t**2, t**3 /)

            CALL dgemv("T", &
                       4, n_elements,&
                       1.d0, &
                       spline_coeffs(:,:,i_x), 4, &
                       t_vec, 1, &
                       0.d0, &
                       value, 1)


!            value(:) = spline_coeffs(1,i_x,:) +&
!                       t*spline_coeffs(2,i_x,:) +&
!                      (t**2)*spline_coeffs(3,i_x,:) +&
!                      (t**3)*spline_coeffs(4,i_x,:)
        END SUBROUTINE spline_value

    END SUBROUTINE spline_extrapolate

    SUBROUTINE linear_extrapolate_boundaries(max_interval_start, &
                                             min_interval_end, &
                                             max_dist_between_points, &
                                             n_elements, &
                                             grid_int, &
                                             data, &
                                             grid_extrapolated, &
                                             data_extrapolated)
        REAL*8,INTENT(IN) :: max_interval_start, &
                             min_interval_end, &
                             max_dist_between_points

        INTEGER, INTENT(IN) :: n_elements
        type(integration_grid), INTENT(IN) :: grid_int
        REAL*8, DIMENSION(n_elements, grid_int%n_points), INTENT(IN) :: data

        type(integration_grid), INTENT(OUT) :: grid_extrapolated
        REAL*8, ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: data_extrapolated

        REAL*8, ALLOCATABLE, DIMENSION(:) :: extrapolated_points, &
                                             extrapolated_weights
        REAL*8 :: cur_dist_between_two_points
        INTEGER :: i_size, n_order, i_elements, &
                   i_points_extra_lower, i_points_extra_upper, &
                   i_partitions_extra_lower, i_partitions_extra_upper
        LOGICAL :: extrapolate_bound_lower, &
                   extrapolate_bound_upper

        extrapolate_bound_lower = .FALSE.
        extrapolate_bound_upper = .FALSE.

        n_order = 3
        i_partitions_extra_lower = 0
        i_partitions_extra_upper = 0

        if(minval(grid_int%points(:))> max_interval_start) then

        stop
            extrapolate_bound_lower = .TRUE.
            cur_dist_between_two_points = minval(grid_int%points(:)) - max_interval_start
            CALL get_amount_partitions_to_add(cur_dist_between_two_points, &
                                              max_dist_between_points, &
                                              n_order, &
                                              i_partitions_extra_lower)
        endif

        if(maxval(grid_int%points(:)) < min_interval_end) then
            extrapolate_bound_upper = .TRUE.
            cur_dist_between_two_points = min_interval_end - maxval(grid_int%points(:))
            CALL get_amount_partitions_to_add(cur_dist_between_two_points, &
                                              max_dist_between_points, &
                                              n_order, &
                                              i_partitions_extra_upper)

        endif

        i_points_extra_lower = i_partitions_extra_lower * n_order
        i_points_extra_upper = i_partitions_extra_upper * n_order

        i_size = grid_int%n_points + &
                 i_points_extra_lower + &
                 i_points_extra_upper

!TODO check alloc
        ALLOCATE(data_extrapolated(n_elements,i_size))
        ALLOCATE(extrapolated_points(i_size))
        ALLOCATE(extrapolated_weights(i_size))

        data_extrapolated(:, i_points_extra_lower + 1 : &
                            i_size - i_points_extra_upper) = data(:,:)

        extrapolated_points(i_points_extra_lower + 1 : &
                            i_size - i_points_extra_upper) = grid_int%points(:)

        extrapolated_weights(i_points_extra_lower + 1 : &
                            i_size - i_points_extra_upper) = grid_int%weights(:)

        if(extrapolate_bound_lower) then
            CALL extend_bound(max_interval_start,minval(grid_int%points(:)), &
                              n_order, i_partitions_extra_lower, &
                              1)
            do i_elements = 1, n_elements
                data_extrapolated(i_elements,1:i_points_extra_lower) = data(i_elements,1)
            enddo
        endif

        if(extrapolate_bound_upper) then
            CALL extend_bound(maxval(grid_int%points(:)), min_interval_end, &
                              n_order, i_partitions_extra_upper, &
                              i_size - i_points_extra_upper +1)
            do i_elements = 1, n_elements
                data_extrapolated(i_elements,grid_int%n_points+i_partitions_extra_lower: &
                                  i_size) = data(i_elements,grid_int%n_points)
            enddo
        endif

        CALL create_integration_grid(i_size, &
                                     extrapolated_points, &
                                     extrapolated_weights, &
                                     grid_extrapolated)


        DEALLOCATE(extrapolated_points)
        DEALLOCATE(extrapolated_weights)
    CONTAINS
        SUBROUTINE extend_bound(interval_start, interval_end, n_order, &
                                n_extra_partitions, i_index_start)
            REAL*8, INTENT(IN) :: interval_start, interval_end
            INTEGER, INTENT(IN) :: n_order, n_extra_partitions, i_index_start

            type(integration_intervals) :: int_intervals
            type(integration_grid) :: an_extra_grid
            INTEGER :: i_index_end

            CALL create_integration_intervals_uniform(interval_start, &
                                                      interval_end, &
                                                      n_extra_partitions, &
                                                      int_intervals)

            CALL create_integration_grid_gauss_legendre(int_intervals, n_order, &
                                                        an_extra_grid)

            i_index_end = i_index_start + an_extra_grid%n_points - 1

            extrapolated_points(i_index_start : &
                                i_index_end) = an_extra_grid%points(:)

            extrapolated_weights(i_index_start : &
                                 i_index_end) = an_extra_grid%weights(:)


            CALL free_integration_intervals(int_intervals)
            CALL free_integration_grid(an_extra_grid)

        END SUBROUTINE

    END SUBROUTINE linear_extrapolate_boundaries


    SUBROUTINE spline_cos_integration(tau_grid, &
                                      n_elements, f_data, &
                                      min_interval_end, &
                                      omega, &
                                      cos_transform_data)

    type(integration_grid), INTENT(IN) :: tau_grid
    INTEGER, INTENT(IN) :: n_elements
    REAL*8, DIMENSION(n_elements,tau_grid%n_points), INTENT(IN) :: f_data
    REAL*8, INTENT(IN) :: min_interval_end,omega
    REAL*8, DIMENSION(n_elements), INTENT(OUT) :: cos_transform_data


    real*8 spl_param (4, tau_grid%n_points)
    REAL*8, ALLOCATABLE,DIMENSION(:,:) :: spline_coeffs
    REAL*8 :: point_last, int_lin,a,b
    INTEGER :: i,j, i_tau, &
               n_order, n_partitions,  &
               i_splined_point, i_spline_pos, &
               n_splined_data
    INTEGER, DIMENSION(tau_grid%n_points-1) :: n_extra_partitions

    ALLOCATE(spline_coeffs(4,tau_grid%n_points))


    point_last = maxval(tau_grid%points(:))

    do i=1, n_elements
       CALL cubic_spline(f_data(i,:), tau_grid%n_points, spline_coeffs(:,:) )

!BE CAREFuL!
!       do j = 1, tau_grid%n_points-1
!          if( abs(f_data(i,j)-f_data(i,j+1)) <1.d-20) then
!            spline_coeffs(:,j) = 0.d0
!          endif
!       enddo

       CALL integrate_spline_cos(tau_grid%n_points,spline_coeffs(:,:), &
                                 tau_grid%points, &
                                 omega, &
                                 cos_transform_data(i))

       if(min_interval_end>point_last) then
          a = point_last
          b = min_interval_end

          int_lin = f_data(i,tau_grid%n_points)/(a-b)* &
                    ( &
                    (-1.d0*a*omega*sin(a*omega) &
                     -1.d0*cos(a*omega) &
                     +b*omega*sin(b*omega) &
                     +cos(b*omega)) / omega**2 &
                    )&
                    - ((f_data(i,tau_grid%n_points)*b)/(a-b)) *&
                       (sin(b * omega) - sin(a * omega)) / omega


!          cos_transform_data(i)=cos_transform_data(i) + int_lin
       endif
    enddo

    DEALLOCATE(spline_coeffs)

    CONTAINS
        SUBROUTINE integrate_spline_cos(n_points, spline_coeffs, &
                                        points, &
                                        omega, &
                                        integral)
            INTEGER,INTENT(IN) :: n_points
            REAL*8, DIMENSION(4,n_points), INTENT(IN) :: spline_coeffs
            REAL*8, DIMENSION(n_points), INTENT(IN) :: points
            REAL*8, INTENT(IN) :: omega
            REAL*8, INTENT(OUT) :: integral

            REAL*8, DIMENSION(n_points-1) :: integral_parts
            INTEGER :: i
            REAL*8 :: a_0,b_0,c_0,d_0, &
                      a,b,c,d,k,w, s,t,&
                      i_a,i_b,i_c,i_d

            integral = 0.d0

            do i=1,n_points-1

                s = points(i)
                t = points(i+1)

                k = omega * (t-s)
                w = omega * s

                a_0 = spline_coeffs(1,i)
                b_0 = spline_coeffs(2,i)
                c_0 = spline_coeffs(3,i)
                d_0 = spline_coeffs(4,i)

                a = a_0 + b_0 * s + c_0 * s**2 + d_0 * s**3
                b = (b_0 + 2.d0 * c_0 * s + 3.d0 * d_0 * s**2) * (t - s)
                c = (c_0 + 3.d0 * d_0 * s) * (t - s)**2
                d = d_0 * (t - s)**3

                i_a = a*(sin(k+w)-sin(w)) / k
                i_b = b*(k*sin(k+w)+cos(k+w)-cos(w)) / k**2
                i_c = c*((k**2.d0-2.d0)*sin(k+w)+2.d0*k*cos(k+w)+2.d0*sin(w)) / k**3
                i_d = d*(k*(k**2-6.d0)*sin(k+w)+3.d0*(k**2-2.d0)*cos(k+w)+6.d0*cos(w)) / k**4

!if(i_a>0.d0) write(use_unit,*) i_a,i_b,i_c,i_d

                integral_parts(i) = (i_a + i_b + i_c + i_d)*(t-s)
!                integral = integral + &
!                           (i_a + i_b + i_c + i_d)*(t-s)
             enddo

            integral = sort_and_sum_vector(n_points-1, integral_parts)

        END SUBROUTINE integrate_spline_cos

    END SUBROUTINE spline_cos_integration


    SUBROUTINE get_amount_partitions_to_add(cur_dist_between_two_points, &
                                        max_dist_between_points, &
                                        n_order, &
                                        n_points_to_add)

        REAL*8, INTENT(IN) :: cur_dist_between_two_points, &
                              max_dist_between_points
        INTEGER, INTENT(IN) :: n_order
        INTEGER, INTENT(OUT) :: n_points_to_add

        n_points_to_add = int(ceiling((cur_dist_between_two_points / max_dist_between_points) &
                              / n_order))
    END SUBROUTINE get_amount_partitions_to_add

    REAL*8 FUNCTION sort_and_sum_vector(n_vec,vec)
        INTEGER, INTENT(IN) :: n_vec
        REAL*8, DIMENSION(n_vec), INTENT(INOUT) :: vec

        REAL*8 :: tmp
        INTEGER :: n,i
!bubblesort
        do n=n_vec, 1,-1
            do i=1,n-1
                if (abs(vec(i)) > abs(vec(i+1))) then
                   tmp = vec(i+1)
                   vec(i+1) = vec(i)
                   vec(i) = tmp
                endif
            enddo
        enddo

!        do n=1,n_vec
!            write(use_unit,*) n, vec(n)
!        enddo

!       sort_and_sum_vector = SUM(vec(:))
!       write(use_unit,*) "---",sort_and_sum_vector
       sort_and_sum_vector = sum_stable(n_vec, vec(:))
!       write(use_unit,*) "----",sort_and_sum_vector
    END FUNCTION sort_and_sum_vector

    REAL*8 FUNCTION sum_stable(n_vec, vec)
        INTEGER, INTENT(IN) :: n_vec
        REAL*8,DIMENSION(n_vec),INTENT(IN) :: vec

        REAL*8 :: c,t,y
        INTEGER :: i

        sum_stable = 0.d0
        c=0.d0

        do i = 1, n_vec
            y = vec(i) - c
            t = sum_stable + y
            c = (t - sum_stable) - y
            sum_stable = t
        enddo
    END FUNCTION sum_stable

    REAL*8 FUNCTION scalprod_stable(n_vec,vec_1,vec_2)
        INTEGER, INTENT(IN) :: n_vec
        REAL*8,DIMENSION(n_vec),INTENT(IN) :: vec_1, vec_2

        REAL*8 :: c,t,y
        INTEGER :: i

        scalprod_stable = 0.d0
        c=0.d0

        do i = 1, n_vec
            y = vec_1(i)*vec_2(i) - c
            t = scalprod_stable + y
            c = (t - scalprod_stable) - y
            scalprod_stable = t
        enddo
    END FUNCTION


    COMPLEX*16 FUNCTION scalprod_stable_comp(n_vec,vec_1,vec_2)
        INTEGER, INTENT(IN) :: n_vec
        COMPLEX*16,DIMENSION(n_vec),INTENT(IN) :: vec_1, vec_2

        REAL*8 :: part_real,part_img

        part_real = scalprod_stable(n_vec, dble(vec_1),dble(vec_2)) &
                    + &
                    scalprod_stable(n_vec, dimag(vec_1),dimag(vec_2))

        part_img =  scalprod_stable(n_vec, dimag(vec_1),dble(vec_2)) &
                    - &
                    scalprod_stable(n_vec, dble(vec_1),dimag(vec_2))

        scalprod_stable_comp = DCMPLX(part_real, part_img)
    END FUNCTION
END MODULE cRPA_calculation_integration
