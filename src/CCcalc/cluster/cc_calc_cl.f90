  Subroutine CC_Calc_cluster()

! This is main subroutine of calculating CC correlation energy CC_E_corr

  Use timing
  Use CC_cl

  Implicit None
  Integer :: i_tmp,j_tmp,k_tmp,rrr,sss,ttt,uuu,iii,jjj,aaa,bbb
  Integer :: s_tmp,i_dom,n_grp,i_grp,d_start,g_start

  Double precision :: CC_time_stamp, CC_clock_stamp,CC_time_pre,CC_clock_pre
  Double precision :: CC_time_start, CC_clock_start
  Double precision :: CC_E_c_pre, CC_dE, CC_E_tot,int_test,ct1

! Get time stamp and clock time stamp
  Call get_timestamps(CC_time_stamp, CC_clock_stamp)

! Initialization
  Call CC_cl_initial()

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
  if (CC_res_pt.eq.0) then
    Call CC_cl_ini_guess()
  else
    Call CC_cl_read_t()
  end if

! Map t vector
  Call CC_cl_w2t()

! Calculate MP2 energy
  Call CC_cl_E_corr()

! save DIIS vectors
  CC_DIIS_ndc = 1
  CC_DIIS_m_bgn = 1
  Call CC_cl_DIIS_save('r',1)

  ! Calculate total energy
  CC_E_tot = E_HF + CC_E_corr

  if (CC_res_pt.eq.0) then
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
    write(use_unit,"(2x,A)") ' '
    write(use_unit,"(2x,A)") 'CCSD calculation starts..'
  end if

  if (CC_res_pt.ne.-1) then
    CC_flag_converge = .false.
  else
    CC_flag_converge = .true.
  end if

  CC_i_scf = 0

! CC SCF loops
  do while ((.not.(CC_flag_converge)).and.(CC_i_scf.lt.CC_max_cycle))
    CC_i_scf = CC_i_scf + 1

    CC_E_c_pre = CC_E_corr

!    write(70+myid,*) 'i_scf',CC_i_scf
!    write(80+myid,*) 'i_scf',CC_i_scf

    ! Evaluate w vector
    Call CC_cl_calc_w()

    ! Evaluate the best coefficient vector by RLE(or DIIS) algorithm
    Call CC_cl_DIIS()

    !Call CC_cl_Jacob()

    ! Update correlation energy
    Call CC_cl_w2t()
    Call CC_cl_E_corr()

!    print*,'E_corr'

    ! Calculate normaliztion coefficient

    CC_dE = CC_E_corr - CC_E_c_pre
    if (abs(CC_dE).lt.CC_converg_thresh) then
      CC_flag_converge = .true.
    end if

    CC_E_tot = E_HF + CC_E_corr

    if (myid.eq.0) then
      write(use_unit,"(2x,A)") '----------------------'
      write(use_unit,"(2x,A,I3,A)") 'For ',CC_i_scf+CC_res_pt,'-th loop'
      write(use_unit,"(2x,A,F14.8)") 'Norm = ',CC_Norm
      write(use_unit,"(2x,A,F20.10,A,F16.10,A,F14.10,A)") 'E(CC) = ',CC_E_tot, &
                       ' Ha, E(corr) = ',CC_E_corr,' Ha, dE(CC) = ',CC_dE,' Ha'
      if (CC_RLE_flag_singular) then
        write(use_unit,"(2X,A)") 'Warnning: R matrix in RLE procedure is & 
                                  &singular, Jacob method is used'
      end if

      if (CC_RLE_flag_sv_overflow) then
        write(use_unit,"(2X,A)") 'Warnning: The RLE saving vector is &
                                  &overflowed, Jacob method is used'
      end if

! Get time stamp and clock time stamp
      CC_time_pre = CC_time_stamp
      CC_clock_pre = CC_clock_stamp

      Call get_timestamps(CC_time_stamp, CC_clock_stamp)
            
      if (myid.eq.0) then
        write(use_unit,"(2x,A,I3,A,F14.3,A,F14.3,A)") 'Time of ',CC_i_scf+CC_res_pt, &
                       'th iteration', CC_time_stamp - CC_time_pre, 's (cpu)', &
                       CC_clock_stamp - CC_clock_pre, 's (wall clock)'
        write(use_unit,"(2x,A)") '----------------------'
      end if
 
    end if

    if ((CC_restart_flag).and.(mod(CC_i_scf,CC_restart_step).eq.0)) then
      Call CC_cl_write_t(CC_i_scf+CC_res_pt)
    end if

  end do

! End of CC calculation
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") '------------------------'
    if (CC_flag_converge) then
      write(use_unit,"(2x,A,I4,A)") 'Convergence achieved in ',CC_i_scf+CC_res_pt, &
                                    ' iterations'
      write(use_unit,"(2x,A,F20.10,A)") 'The CCSD energy is ',CC_E_tot,' Ha'
    else
      write(use_unit,"(2x,A,I4,A)") 'Warning: CCSD scf procedure failed in ', &
                                 CC_max_cycle,' iterations'
    end if
  end if

! Get time stamp and clock time stamp
  CC_time_pre = CC_time_stamp
  CC_clock_pre = CC_clock_stamp
  Call get_timestamps(CC_time_stamp, CC_clock_stamp)

  if (myid.eq.0) then
    write(use_unit,"(2x,A,F14.3,A,F14.3,A)") 'Time of CCSD calculation', &
                     CC_time_stamp - CC_time_start, 's (cpu)',           &
                     CC_clock_stamp - CC_clock_start, 's (wall clock)'
  end if

  if ((CC_restart_flag).and.(CC_res_pt.ne.-1)) then
    if (CC_flag_converge) then
      Call CC_cl_write_t(-1)
    else
      Call CC_cl_write_t(CC_i_scf+CC_res_pt)
    end if
  end if

  if ((CC_flag_converge).and.(CC_calc_method.eq.2)) then

    if (myid.eq.0) then
      write(use_unit,"(2x,A)") 'Perturbative triples calculation starts...'
    end if

    if (CC_check_flag) then
!      Call CC_cl_PT_check()
    else
      Call CC_cl_PT()
    end if

    CC_time_pre = CC_time_stamp
    CC_clock_pre = CC_clock_stamp
    Call get_timestamps(CC_time_stamp, CC_clock_stamp)
  
    if (myid.eq.0) then
      write(use_unit,"(2x,A)") 'Perturbative triples calculation completed.'
      write(use_unit,"(2x,A,F20.10,A)") 'Perturbative triples correction: ', &
                                        CC_E_PT,' Ha.'

      write(use_unit,"(2x,A,F20.10,A)") 'CCSD(T) correlation energy = ', &
                                        CC_E_corr + CC_E_PT,' Ha.'

      write(use_unit,"(2x,A,F20.10,A)") 'CCSD(T) Total energy = ', &
                                        E_HF + CC_E_corr + CC_E_PT,' Ha.'

      write(use_unit,"(2x,A,F14.3,A,F14.3,A)") 'Time of perturbative triples calculation', &
                       CC_time_stamp - CC_time_pre, 's (cpu)',           &
                       CC_clock_stamp - CC_clock_pre, 's (wall clock)'
    end if

  end if

  End Subroutine CC_Calc_cluster
