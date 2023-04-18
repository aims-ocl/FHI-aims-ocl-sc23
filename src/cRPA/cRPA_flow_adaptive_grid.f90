!implements an adaptive grid
MODULE cRPA_flow_adaptive_grid
    USE cRPA_view
    implicit none

    ABSTRACT INTERFACE
        SUBROUTINE get_signal_template(point,signal)
            REAL*8, INTENT(IN) :: point
            REAL*8, INTENT(OUT) :: signal
        END SUBROUTINE get_signal_template
    END INTERFACE

    type adaptive_grid
        REAL*8, ALLOCATABLE, DIMENSION(:) :: points ,&
                                             values

        INTEGER :: n_points
    end type adaptive_grid

    INTEGER :: ADAPTIVE_GRID_POINT_NOT_FOUND
    PARAMETER(ADAPTIVE_GRID_POINT_NOT_FOUND = -1)

CONTAINS
    SUBROUTINE create_adaptive_grid_adaptivly(derivative_tolerance, &
                                              difference_tolerance, &
                                              get_signal, &
                                              a_grid)
        REAL*8, INTENT(IN) :: derivative_tolerance, &
                              difference_tolerance
        PROCEDURE(get_signal_template),POINTER :: get_signal
        type(adaptive_grid), INTENT(OUT) :: a_grid

        REAL*8 :: point_min, point_max
        INTEGER :: i_pos_added

        CALL create_adaptive_grid_adaptivly_rough(derivative_tolerance, &
                                                  difference_tolerance, &
                                                  get_signal, &
                                                  a_grid)


return

        CALL refine_adaptive_grid_adaptivly(derivative_tolerance, &
                                            difference_tolerance, &
                                            get_signal, &
                                            a_grid)


    END SUBROUTINE create_adaptive_grid_adaptivly

    SUBROUTINE create_adaptive_grid_adaptivly_rough(derivative_tolerance, &
                                                    difference_tolerance, &
                                                    get_signal, &
                                                    a_grid)
        REAL*8, INTENT(IN) :: derivative_tolerance, &
                              difference_tolerance

        PROCEDURE(get_signal_template),POINTER :: get_signal
        type(adaptive_grid), INTENT(OUT) :: a_grid

        real*8 h, t_0,s_0, point_next, signal, a_value
        INTEGER :: i_points, i_pos_added, i_slow
        type(adaptive_grid) :: tmp_grid


        !h = 0.1
        h = 0.3d0
        t_0 = 1.d-6!0.00001
        s_0 = 0.1

        CALL get_signal(1.d-20,a_value)
        CALL create_adaptive_grid(1,(/1.d-20/),(/a_value/), &
                                  a_grid)

        do i_points=1,2
            point_next = t_0 * dble(i_points)**5.d0!
            point_next = t_0 * (exp(i_points*h) -1.d0)
!            point_next = s_0 * dble(i_points)
            CALL add_point_next()
        enddo

        i_slow = 1

        do while(.NOT. is_last_derivation_smaller(derivative_tolerance,a_grid) &
                 .OR. &
                 .NOT. is_last_difference_smaller(difference_tolerance,a_grid) &
                 .OR. &
                 i_points<10 .OR. point_next<15 )

          i_points = i_points + 1

          point_next = t_0 * dble(i_points)**5.d0! (exp(i_points*h) -1.d0)
          point_next = t_0 * (exp(i_points*h) -1.d0)
 !         point_next = s_0 * dble(i_points)


          if(abs(point_next) > 5.d0) then
!             return
!             point_next=5.d0 + s_0 * dble(i_slow)**2.d0
             i_slow = i_slow + 1
          endif
!
!
          if(point_next > 15.d0) exit

          CALL add_point_next()
        enddo

    CONTAINS
        SUBROUTINE add_point_next()
          CALL get_signal(point_next,signal)
          CALL add_point_to_adaptive_grid(point_next ,signal, &
                                          a_grid,i_pos_added)
        END SUBROUTINE add_point_next

    END SUBROUTINE create_adaptive_grid_adaptivly_rough

    SUBROUTINE refine_adaptive_grid_adaptivly(derivative_tolerance, &
                                              difference_tolerance, &
                                              get_signal, &
                                              a_grid)

        REAL*8, INTENT(IN) :: derivative_tolerance, difference_tolerance
        PROCEDURE(get_signal_template),POINTER :: get_signal
        type(adaptive_grid), INTENT(INOUT) :: a_grid

        type(adaptive_grid) :: grid_rough
        INTEGER :: i_interval, i_recursion_max_deep

        CALL clone_adaptive_grid(a_grid, grid_rough)

        i_recursion_max_deep = 1

        do i_interval = 1, grid_rough%n_points-1
!CALL write_debug("Refining"//num2str(grid_rough%values(i_interval)) &
!                           //num2str(grid_rough%values(i_interval+1)))
            CALL refine_adaptive_grid_adaptivly_recursive(grid_rough%points(i_interval), &
                                                          grid_rough%values(i_interval), &
                                                          grid_rough%points(i_interval+1), &
                                                          grid_rough%values(i_interval+1), &
                                                          derivative_tolerance, &
                                                          difference_tolerance, &
                                                          get_signal, &
                                                          i_recursion_max_deep, &
                                                          a_grid)
        enddo
    END SUBROUTINE refine_adaptive_grid_adaptivly

    RECURSIVE SUBROUTINE refine_adaptive_grid_adaptivly_recursive(interval_start, &
                                                                  value_start, &
                                                                  interval_end, &
                                                                  value_end, &
                                                                  derivative_tolerance, &
                                                                  difference_tolerance, &
                                                                  get_signal, &
                                                                  i_recursions_left, &
                                                                  a_grid)

        REAL*8, INTENT(IN) :: interval_start, &
                              value_start, &
                              interval_end, &
                              value_end, &
                              derivative_tolerance, &
                              difference_tolerance
        PROCEDURE(get_signal_template),POINTER :: get_signal
        INTEGER, INTENT(IN) :: i_recursions_left
        type(adaptive_grid), INTENT(INOUT) :: a_grid

        INTEGER :: i_pos_added, i
        REAL*8 :: point_new, value_new, &
                  points(3), values(3)

        if(i_recursions_left == 0) return

        if(abs(value_start-value_end) <= difference_tolerance &
           .AND. &
           abs((value_start-value_end)/(interval_start-interval_end)) <= derivative_tolerance) then

            return
        endif

        point_new =  interval_start + &
                    (interval_end - interval_start) / 2.d0
!CALL write_debug("New point"//num2str(interval_start)//num2str(interval_end)//num2str(point_new))
        CALL get_signal(point_new, value_new)
        CALL add_point_to_adaptive_grid(point_new, value_new, &
                                        a_grid,i_pos_added)

        points(1) = interval_start!a_grid%points(i_pos_added-1)
        values(1) = value_start!a_grid%values(i_pos_added-1)

        points(2) = point_new
        values(2) = value_new

        points(3) = interval_end!a_grid%points(i_pos_added+1)
        values(3) = value_end!a_grid%values(i_pos_added+1)

        do i = 1,2
            if(abs(values(i)-values(i+1)) > difference_tolerance &
               .OR. &
               abs((values(i)-values(i+1))/(points(i)-points(i+1))) > derivative_tolerance) then

                CALL refine_adaptive_grid_adaptivly_recursive(points(i), &
                                                              values(i), &
                                                              points(i+1), &
                                                              values(i+1), &
                                                              derivative_tolerance, &
                                                              difference_tolerance, &
                                                              get_signal, &
                                                              i_recursions_left-1, &
                                                              a_grid)
            endif
        enddo


    END SUBROUTINE refine_adaptive_grid_adaptivly_recursive


    SUBROUTINE create_adaptive_grid(n_points,points,values, &
                                    a_grid)
        INTEGER, INTENT(IN) :: n_points
        REAL*8,DIMENSION(n_points), INTENT(IN) :: points, values
        type(adaptive_grid), INTENT(OUT) :: a_grid

!TODO check alloc
        ALLOCATE(a_grid%points(n_points))
        ALLOCATE(a_grid%values(n_points))
        a_grid%n_points = n_points

        a_grid%points(:) = points(:)
        a_grid%values(:) = values(:)

    END SUBROUTINE create_adaptive_grid

    SUBROUTINE clone_adaptive_grid(source,dest)
        type(adaptive_grid), INTENT(IN) :: source
        type(adaptive_grid), INTENT(OUT) :: dest

        CALL create_adaptive_grid(source%n_points, &
                                  source%points, &
                                  source%values, &
                                  dest)
    END SUBROUTINE clone_adaptive_grid

    SUBROUTINE add_point_to_adaptive_grid(point, val, a_grid, i_pos_added)
        REAL*8, INTENT(IN) :: point,val
        type(adaptive_grid), INTENT(INOUT) :: a_grid
        INTEGER, INTENT(OUT) :: i_pos_added

        REAL*8, DIMENSION(a_grid%n_points) :: points_old ,&
                                              values_old
        REAL*8, DIMENSION(a_grid%n_points+1) :: points_new ,&
                                                values_new

        INTEGER :: n_points_old, n_points_new, i_pos_first_greater

!CALL write_debug("Adding "//num2str(point)//num2str(val))

        points_old(:) = a_grid%points(:)
        values_old(:) = a_grid%values(:)
        n_points_old = a_grid%n_points

        CALL free_adaptive_grid(a_grid)

        do i_pos_first_greater = 1, a_grid%n_points
            if(points_old(i_pos_first_greater) > point) then
                exit
            endif
        enddo

        n_points_new = n_points_old + 1

        points_new(1:i_pos_first_greater-1) = points_old(1:i_pos_first_greater-1)
        values_new(1:i_pos_first_greater-1) = values_old(1:i_pos_first_greater-1)

        i_pos_added = i_pos_first_greater
        points_new(i_pos_added) = point
        values_new(i_pos_added) = val

        points_new(i_pos_added+1:n_points_new) = &
                                        points_old(i_pos_first_greater:n_points_old)
        values_new(i_pos_added+1:n_points_new) = &
                                        values_old(i_pos_first_greater:n_points_old)

        CALL create_adaptive_grid(n_points_new, points_new, values_new, &
                                  a_grid)
    END SUBROUTINE add_point_to_adaptive_grid

    SUBROUTINE find_point(a_grid, point, i_pos)
        type(adaptive_grid), INTENT(IN) :: a_grid
        REAL*8, INTENT(IN) :: point
        INTEGER, INTENT(OUT) :: i_pos

        INTEGER :: i_point

        do i_point = 1,a_grid%n_points
            if(abs(a_grid%points(i_point) - point) < 1.d-30) then
                i_pos = i_point
                return
            endif
        enddo

        i_pos = ADAPTIVE_GRID_POINT_NOT_FOUND
    END SUBROUTINE find_point

    LOGICAL FUNCTION has_point(a_grid, point)
        type(adaptive_grid), INTENT(IN) :: a_grid
        REAL*8, INTENT(IN) :: point

        INTEGER :: i_pos

        CALL find_point(a_grid, point, i_pos)

        has_point = ( i_pos/= ADAPTIVE_GRID_POINT_NOT_FOUND)

    END FUNCTION has_point

    SUBROUTINE get_discrete_derivation_of_adaptive_grid(a_grid, derivation)
        type(adaptive_grid), INTENT(IN) :: a_grid
        REAL*8, DIMENSION(a_grid%n_points-1), INTENT(OUT) :: derivation

        INTEGER :: i

        do i = 1, a_grid%n_points -1
            derivation(i) = discrete_derivation_of_interval(i,a_grid)
        enddo


    END SUBROUTINE get_discrete_derivation_of_adaptive_grid

    REAL*8 FUNCTION discrete_derivation_of_interval(i,a_grid)
        INTEGER, INTENT(IN) :: i
        type(adaptive_grid), INTENT(IN) :: a_grid

        discrete_derivation_of_interval = (a_grid%values(i)-a_grid%values(i+1)) / &
                                          (a_grid%points(i)-a_grid%points(i+1))
    END FUNCTION discrete_derivation_of_interval

    REAL*8 FUNCTION difference_of_interval(i,a_grid)
        INTEGER, INTENT(IN) :: i
        type(adaptive_grid), INTENT(IN) :: a_grid

        difference_of_interval = (a_grid%values(i)-a_grid%values(i+1))
    END FUNCTION difference_of_interval


    LOGICAL FUNCTION is_derivation_smaller(derivative_tolerance, i_interval, a_grid)
        REAL*8, INTENT(IN) :: derivative_tolerance
        INTEGER, INTENT(IN) :: i_interval
        type(adaptive_grid), INTENT(IN) :: a_grid

!CALL write_debug("values"//num2str(a_grid%values(i_interval))//&
!                 num2str(a_grid%values(i_interval+1)))

CALL write_stdout("derivation"//num2str(i_interval)//&
                 num2str(discrete_derivation_of_interval(i_interval, a_grid)))

        if(abs(discrete_derivation_of_interval(i_interval, a_grid)) < derivative_tolerance) then
            is_derivation_smaller = .TRUE.
            return
        endif

        is_derivation_smaller = .FALSE.

    END FUNCTION is_derivation_smaller

    LOGICAL FUNCTION is_last_derivation_smaller(derivative_tolerance, a_grid)
        REAL*8, INTENT(IN) :: derivative_tolerance
        type(adaptive_grid), INTENT(IN) :: a_grid

!TODO refactor
        if(a_grid%n_points < 2) then
            is_last_derivation_smaller = .FALSE.
            return
        endif

        is_last_derivation_smaller = is_derivation_smaller(derivative_tolerance, &
                                                           a_grid%n_points-1, &
                                                           a_grid)

    END FUNCTION is_last_derivation_smaller

    LOGICAL FUNCTION is_first_derivation_smaller(derivative_tolerance, a_grid)
        REAL*8, INTENT(IN) :: derivative_tolerance
        type(adaptive_grid), INTENT(IN) :: a_grid

        is_first_derivation_smaller = is_derivation_smaller(derivative_tolerance, &
                                                            1, &
                                                            a_grid)

    END FUNCTION is_first_derivation_smaller


    LOGICAL FUNCTION is_difference_smaller(difference_tolerance, i_interval, a_grid)
        REAL*8, INTENT(IN) :: difference_tolerance
        INTEGER, INTENT(IN) :: i_interval
        type(adaptive_grid), INTENT(IN) :: a_grid

!CALL write_debug("difference"//num2str(i_interval)//&
!                 num2str(difference_of_interval(i_interval, a_grid)))

        if(abs(difference_of_interval(i_interval, a_grid)) < difference_tolerance) then
            is_difference_smaller = .TRUE.
            return
        endif

        is_difference_smaller = .FALSE.

    END FUNCTION is_difference_smaller

    LOGICAL FUNCTION is_last_difference_smaller(difference_tolerance, a_grid)
        REAL*8, INTENT(IN) :: difference_tolerance
        type(adaptive_grid), INTENT(IN) :: a_grid

!TODO refactor
        if(a_grid%n_points < 2) then
            is_last_difference_smaller = .FALSE.
            return
        endif

        is_last_difference_smaller = is_difference_smaller(difference_tolerance, &
                                                           a_grid%n_points-1, &
                                                           a_grid)

    END FUNCTION is_last_difference_smaller

    LOGICAL FUNCTION is_first_difference_smaller(difference_tolerance, a_grid)
        REAL*8, INTENT(IN) :: difference_tolerance
        type(adaptive_grid), INTENT(IN) :: a_grid

        is_first_difference_smaller = is_difference_smaller(difference_tolerance, &
                                                            1, &
                                                            a_grid)

    END FUNCTION is_first_difference_smaller

    SUBROUTINE print_adaptive_grid(a_grid)
        type(adaptive_grid), INTENT(IN) :: a_grid

        INTEGER :: i

        CALL write_stdout("points values")

        do i = 1, a_grid%n_points
            CALL write_stdout(num2str(a_grid%points(i))// &
                              num2str(a_grid%values(i)))
        enddo

    END SUBROUTINE print_adaptive_grid

    SUBROUTINE free_adaptive_grid(a_grid)
        type(adaptive_grid), INTENT(INOUT) :: a_grid

        DEALLOCATE(a_grid%points)
        DEALLOCATE(a_grid%values)
    END SUBROUTINE free_adaptive_grid
END MODULE cRPA_flow_adaptive_grid
