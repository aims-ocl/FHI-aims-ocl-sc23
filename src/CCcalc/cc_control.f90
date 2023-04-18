  Subroutine CC_main()

  Use runtime_choices
  Use dimensions

  if (n_periodic.eq.0) then
    ! Cluster calculation
    if ((n_spin.ne.1).or.(flag_cc_general)) then
      Call CC_Calculation()
    else
      !if (CC_solver.eq.'cluster') then
        Call CC_Calc_cluster()
      !else if (CC_solver.eq.'hybrid') then
      !  Call CC_Calc_hybrid()
      !else if (CC_solver.eq.'collective') then
      !  Call CC_Calc_cs()
      !else if (CC_solver.eq.'independent') then
      !  Call CC_Calc_is()
      !end if
    end if
  else
    if (CC_gamma_flag) then
      Call CC_Calc_cluster()
    !else if (CC_check_flag) then
    !  Call CC_Calc_periodic()
    else
      Call CC_calc_3d()
    end if
  end if

  End Subroutine CC_main
