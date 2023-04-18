  Subroutine CC_Calculation()

! This is main subroutine of calculating CC correlation energy CC_E_corr

  use dimensions
  use runtime_choices
  use physics
  use synchronize_mpi
  use basis
  use prodbas
  use hartree_fock
  use mpi_tasks
  use timing
  Use CC_corr

  Implicit None
  Integer :: i_tmp,j_tmp,k_tmp,CC_i_scf,rrr,sss,ttt,uuu
  Logical :: CC_flag_converge

  Double precision :: CC_time_stamp, CC_clock_stamp,CC_time_pre,CC_clock_pre
  Double precision :: CC_E_c_pre, CC_Norm, CC_dE, CC_E_tot,int_test,ct1,rlt,ddot

! Get time stamp and clock time stamp
  Call get_timestamps(CC_time_stamp, CC_clock_stamp)

! Initialization
  Call cc_initialization()

  CC_i_scf = 0
  CC_flag_converge = .false.

! Get time stamp and clock time stamp
  CC_time_pre = CC_time_stamp
  CC_clock_pre = CC_clock_stamp
  Call get_timestamps(CC_time_stamp, CC_clock_stamp)

  if (myid.eq.0) then
    write(use_unit,"(2x,A,F10.3,A,F10.3,A)") 'Time of initialization', &
                     CC_time_stamp - CC_time_pre, 's (cpu)', &
                     CC_clock_stamp - CC_clock_pre, 's (wall clock)'
  end if

! Initial guess of wave function (MP2 calculation)
  Call CC_ini_guess()

!  Call CC_jacob()
!  Call CC_Calc_E_corr()
  CC_E_tot = E_HF + CC_E_corr 

!  if (myid.eq.0) then
!    print*,'CC_t_1a'
!    do rrr = 1,CC_n_elec(1)
!      do sss = 1,CC_n_vir(1)
!        print*,CC_t_1a(rrr,sss)
!      end do
!    end do
!    print*,'CC_t_ab'
!    i_tmp = 0
!    do rrr = 1,CC_n_elec(1)
!      do sss = 1,CC_n_elec(2)
!        do ttt = CC_n_elec(1) + 1, n_states
!          do uuu = CC_n_elec(2) + 1, n_states
!            i_tmp = i_tmp + 1
!            print*,rrr,sss,ttt,uuu,CC_t_ab(i_tmp)
!          end do
!        end do
!      end do
!    end do

!    print*,'CC_t_2a'
!    do i_tmp = 1, CC_n_2a
!      Call CC_decode_index_xx(i_tmp,1,rrr,sss,ttt,uuu)
!      print*,rrr,sss,ttt,uuu,CC_t_2a(i_tmp)
!    end do

!  end if

    write(42,*) n_basbas
    do rrr = 1,n_states
      do sss = 1,n_states
        do ttt = 1, n_states
          do uuu = 1, n_states
            rlt = ddot(n_basbas,ovlp_3ks(:,rrr,sss,1),1, &
                                ovlp_3ks(:,ttt,uuu,1),1)
            write(42,"(4I8,F16.10)") rrr,sss,ttt,uuu,rlt
          end do
        end do
      end do
    end do


!  Call CC_check()

  Call CC_Calc_Norm(CC_Norm)
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'MP2 calculation (initial guess of CC wave function)'
    write(use_unit,"(2x,A,F14.8)") 'Norm = ', CC_Norm
    write(use_unit,"(2x,A,F20.10,A,F16.10,A)") 'E(MP2) = ',CC_E_tot, &
                      ' Ha, Ec(MP2) = ', CC_E_corr, ' Ha'
  end if

! Get time stamp and clock time stamp
  CC_time_pre = CC_time_stamp
  CC_clock_pre = CC_clock_stamp
  Call get_timestamps(CC_time_stamp, CC_clock_stamp)

  if (myid.eq.0) then
    write(use_unit,"(2x,A,F10.3,A,F10.3,A)") 'Time of MP2 calculation', &
                     CC_time_stamp - CC_time_pre, 's (cpu)', &
                     CC_clock_stamp - CC_clock_pre, 's (wall clock)'
  end if


! CC SCF loops
  do while ((.not.(CC_flag_converge)).and.(CC_i_scf.le.CC_max_cycle))
    CC_i_scf = CC_i_scf + 1

    CC_E_c_pre = CC_E_corr

    ! Calculate auxiliary vectors
    Call CC_Calc_aux_vect()

    ! Evaluate w vector
    Call CC_Calc_w()

    ! Evaluate the best coefficient vector by RLE algorithm
!    if (CC_i_scf.eq.1) then
      Call CC_Jacob()
!    else
!      Call CC_RLE()
!    end if

    ! Update correlation energy
    Call CC_Calc_E_corr()

    ! Calculate normalization factor
    Call CC_Calc_Norm(CC_Norm)

!    if (myid.eq.0) then
!    print*,'CC_t_1a'
!    do rrr = 1,CC_n_elec(1)
!      do sss = 1,CC_n_vir(1)
!        print*,CC_t_1a(rrr,sss)
!      end do
!    end do
!    print*,'CC_t_ab'
!    i_tmp = 0
!    do rrr = 1,CC_n_elec(1)
!      do sss = 1,CC_n_elec(2)
!        do ttt = CC_n_elec(1) + 1, n_states
!          do uuu = CC_n_elec(2) + 1, n_states
!            i_tmp = i_tmp + 1
!            print*,rrr,sss,ttt,uuu,CC_t_ab(i_tmp)
!          end do
!        end do
!      end do
!    end do
!    end if
!    print*,'CC_t_2a'
!    do i_tmp = 1, CC_n_2a
!      Call CC_decode_index_xx(i_tmp,1,rrr,sss,ttt,uuu)
!      print*,rrr,sss,ttt,uuu,CC_t_2a(i_tmp)
!    end do

    ! Convergence test
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
      if (CC_RLE_flag_singular) then
        write(use_unit,"(2X,A)") 'Warnning: R matrix in RLE procedure is & 
                                  &singular, Jacob method is used'
      end if

      if (CC_RLE_flag_sv_overflow) then
        write(use_unit,"(2X,A)") 'Warnning: The RLE saving vector is &
                                  &overflowed, Jacob method is used'
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
    write(use_unit,"(2x,A,F10.3,A,F10.3,A)") 'Time of CC calculation', &
                     CC_time_stamp - CC_time_pre, 's (cpu)',           &
                     CC_clock_stamp - CC_clock_pre, 's (wall clock)'
  end if

  End Subroutine CC_Calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_check()

  Use CC_corr
  USe dimensions

  Implicit None
  Integer :: i,j,a,b,k,l,c,d,i_run,i_state,i_spin

  print*,'CC_check'
  do i_spin = 1, n_spin
    print*,'i_spin',i_spin
    do i = CC_valence, CC_n_elec(i_spin)
      do a = CC_n_elec(i_spin) + 1, n_states
        Call CC_code_index_s(i_state,i_spin,i,a)
        print*,'code_s---i,a,i_state',i,a,i_state
        Call CC_decode_index_s(i_state,i_spin,k,l)
        print*,'decode_s---i,a,i_state',k,l,i_state
      end do
    end do

    do i_state = CC_index_h_ij(myid+1,i_spin), CC_index_h_ij(myid+1,i_spin) + CC_mem_h_ij(myid+1,i_spin) - 1
      Call CC_decode_ij4h(i_state,i,j,i_spin)
      print*,'decode_ij4h---i,j,i_state',i,j,i_state
    end do

    do i_state = CC_index_h_ba(myid+1,i_spin), CC_index_h_ba(myid+1,i_spin) + CC_mem_h_ba(myid+1,i_spin) - 1
      Call CC_decode_ab4h(i_state,a,b,i_spin)
      print*,'decode_ab4h---a,b,i_state',a,b,i_state
    end do

    do i_state = CC_index_h_bj(myid+1,i_spin), CC_index_h_bj(myid+1,i_spin) + CC_mem_h_bj(myid+1,i_spin) - 1
      Call CC_decode_bj4h(i_state,b,j,i_spin)
      print*,'decode_bj4h---b,j,i_state',b,j,i_state
    end do


    do i = CC_valence, CC_n_elec(i_spin) - 1
      do j = i+1, CC_n_elec(i_spin)
        do a = CC_n_elec(i_spin) + 1, n_states - 1
          do b = a+1, n_states
            Call CC_code_index_xx(i_state,i_spin,i,j,a,b)
            print*,'code_index_xx,i_state,i,j,a,b',i_state,i,j,a,b
            Call CC_decode_index_xx(i_state,i_spin,k,l,c,d)
            print*,'decode_index_xx,i_state,i,j,a,b',i_state,k,l,c,d
          end do
        end do
      end do
    end do
  end do

  End Subroutine CC_check
