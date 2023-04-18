MODULE cRPA_flow_adaptive_grid_test
    USE cRPA_flow_adaptive_grid
    implicit none

    REAL*8 :: decay_constant_positiv, &
              decay_constant_negativ

CONTAINS

    SUBROUTINE func_exp_decay(point,val)
        REAL*8, INTENT(IN) :: point
        REAL*8, INTENT(OUT) :: val

        REAL*8 :: decay_constant

        if(point >= 0.d0) then
            decay_constant = decay_constant_positiv
        else
            decay_constant = decay_constant_negativ
        endif

        val = exp(-decay_constant*abs(point)) + 1
    END SUBROUTINE func_exp_decay

    SUBROUTINE adaptive_grid_test()
        PROCEDURE(get_signal_template),POINTER :: a_signal
        type(adaptive_grid)  :: a_grid
        REAL*8 :: derivation_tolerance, &
                  difference_tolerance

        decay_constant_positiv = 1.d0
        decay_constant_negativ = 1.d1

        a_signal=>func_exp_decay
        derivation_tolerance = 1.d-3
        difference_tolerance = 1.d-1

        CALL create_adaptive_grid_adaptivly(derivation_tolerance, &
                                            difference_tolerance, &
                                            a_signal, &
                                            a_grid)
        CALL print_adaptive_grid(a_grid)
        CALL free_adaptive_grid(a_grid)
    END SUBROUTINE adaptive_grid_test
END MODULE cRPA_flow_adaptive_grid_test
