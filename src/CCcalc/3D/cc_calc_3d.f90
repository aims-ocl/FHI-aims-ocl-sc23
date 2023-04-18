  Subroutine CC_Calc_3d()

! This is main subroutine of calculating CC correlation energy CC_E_corr

  Use timing
  Use CC_3d
  use localorb_io, only: use_unit

  Implicit None
  Integer :: i_tmp,j_tmp,k_tmp,rrr,sss,ttt,uuu,iii,jjj,aaa,bbb
  Integer :: s_tmp,i_dom,n_grp,i_grp,d_start,g_start
  Logical :: CC_flag_converge

  Double precision :: CC_time_stamp, CC_clock_stamp,CC_time_pre,CC_clock_pre
  Double precision :: CC_time_start, CC_clock_start
  Double precision :: CC_E_c_pre, CC_dE, CC_E_tot,int_test,ct1

! Get time stamp and clock time stamp
  Call get_timestamps(CC_time_stamp, CC_clock_stamp)

! Initialization
  Call CC_3d_initial()

! Get time stamp and clock time stamp
  CC_time_pre = CC_time_stamp
  CC_clock_pre = CC_clock_stamp
  Call get_timestamps(CC_time_stamp, CC_clock_stamp)

  if (myid.eq.0) then
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") 'Time of initialization', &
                     CC_time_stamp - CC_time_pre, 's (cpu)', &
                     CC_clock_stamp - CC_clock_pre, 's (wall clock)'
  end if

! Initial guess of wave function (MP2 calculation)
  if (CC_restart_point.eq.0) then
    Call CC_3d_ini_guess()
    CC_i_scf = 0
  end if

! Map t vector
  Call CC_3d_w2t()

! Calculate MP2 energy
  Call CC_3d_E_corr()

  CC_DIIS_ndc = 1
  CC_DIIS_m_bgn = 1
  Call CC_3d_DIIS_save('r',1)

  ! Calculate total energy
  CC_E_tot = E_HF + CC_E_corr

  if (CC_restart_point.eq.0) then
    if (myid.eq.0) then
      write(use_unit,"(2x,A)") 'MP2 calculation (initial guess of CC wave function)'
      write(use_unit,"(2x,A,F14.8)") 'Norm = ', CC_Norm
      write(use_unit,"(2x,A,F20.10,A,F16.10,A)") 'E(MP2) = ',CC_E_tot, &
                        ' Ha, Ec(MP2) = ', CC_E_corr, ' Ha'
    end if
  end if

! Get time stamp and clock time stamp
  CC_time_pre = CC_time_stamp
  CC_clock_pre = CC_clock_stamp
  CC_time_start = CC_time_stamp
  CC_clock_start = CC_clock_stamp
  Call get_timestamps(CC_time_stamp, CC_clock_stamp)

  if (myid.eq.0) then
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") 'Time of MP2 calculation', &
                     CC_time_stamp - CC_time_pre, 's (cpu)', &
                     CC_clock_stamp - CC_clock_pre, 's (wall clock)'
  end if

  CC_flag_converge = .false.

! CC SCF loops
  do while ((.not.(CC_flag_converge)).and.(CC_i_scf.lt.CC_max_cycle))
    CC_i_scf = CC_i_scf + 1

    !write(myid+80,*) 'CC_i_scf',CC_i_scf
    !write(myid+90,*) 'CC_i_scf',CC_i_scf

    CC_E_c_pre = CC_E_corr

    ! Evaluate w vector
    Call CC_3d_calc_w()

!    print*,'calc_w'

    ! Evaluate the best coefficient vector by DIIS algorithm
    if (CC_acc_method.eq.1) then
      Call CC_3d_Jacob()
    else
      Call CC_3d_DIIS()
    end if

    ! Update correlation energy
    Call CC_3d_w2t()
    Call CC_3d_E_corr()

    CC_dE = CC_E_corr - CC_E_c_pre
    if (abs(CC_dE).lt.CC_converg_thresh) then
      CC_flag_converge = .true.
    end if

    CC_E_tot = E_HF + CC_E_corr

    if (myid.eq.0) then
      write(use_unit,"(2x,A)") '----------------------'
      write(use_unit,"(2x,A,I3,A)") 'For ',CC_i_scf,'-th loop'
      write(use_unit,"(2x,A,F14.8)") 'Norm = ',CC_Norm
      write(use_unit,"(2x,A,F20.10,A,F16.10,A,F14.10,A)") 'E(CC) = ',CC_E_tot, &
                       ' Ha, E(corr) = ',CC_E_corr,' Ha, dE(CC) = ',CC_dE,' Ha'

! Get time stamp and clock time stamp
      CC_time_pre = CC_time_stamp
      CC_clock_pre = CC_clock_stamp

      Call get_timestamps(CC_time_stamp, CC_clock_stamp)
            
      if (myid.eq.0) then
        write(use_unit,"(2x,A,I3,A,F14.3,A,F14.3,A)") 'Time of ',CC_i_scf, &
                       'th iteration', CC_time_stamp - CC_time_pre, 's (cpu)', &
                       CC_clock_stamp - CC_clock_pre, 's (wall clock)'
      end if

    end if

  end do

! End of CC calculation
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") '------------------------'
    if (CC_flag_converge) then
      write(use_unit,"(2x,A,I4,A)") 'Convergence achieved in ',CC_i_scf, &
                                    ' iterations'
      write(use_unit,"(2x,A,F20.10,A)") 'The total energy is ',CC_E_tot,' Ha'
    else
      write(use_unit,"(2x,A,I4,A)") 'Warning: CC scf procedure failed in ', &
                                 CC_max_cycle,' iterations'
    end if
  end if

! Get time stamp and clock time stamp
  CC_time_pre = CC_time_stamp
  CC_clock_pre = CC_clock_stamp
  Call get_timestamps(CC_time_stamp, CC_clock_stamp)

  if (myid.eq.0) then
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") 'Time of CC calculation', &
                     CC_time_stamp - CC_time_start, 's (cpu)',           &
                     CC_clock_stamp - CC_clock_start, 's (wall clock)'
  end if

!  if (CC_use_disk) then
!    Call CC_3d_hdf5_cleanup()
!  end if

  End Subroutine CC_Calc_3d
