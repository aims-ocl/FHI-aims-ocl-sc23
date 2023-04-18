  Subroutine CC_Calc_E_corr()

! This subroutine calculates the correlation energy CC_E_corr

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_run
  Integer :: i_spin,iii,jjj,aaa,bbb
  Double precision :: int_rlt,tau,ct1,ct2,ct3,ct4,ct5

  CC_E_corr = 0.0D0

!  write(myid+40,*) 'E_corr',CC_E_corr
  do i_spin = 1, n_spin
!    write(myid+40,*) 'i_spin',i_spin
    i_state = CC_index_w(myid+1,3+i_spin) - 1
!    write(myid+40,*) 'i_state',i_state

    do i_run = 1, CC_mem_w(myid+1,3+i_spin)

      i_state = i_state + 1
      Call CC_decode_index_xx(i_state,i_spin,iii,jjj,aaa,bbb)

!    write(myid+40,*) 'i_state',i_state
!    write(myid+40,*) 'iii,jjj,aaa,bbb',iii,jjj,aaa,bbb


      Call CC_Calc_integral(int_rlt,iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin)

      Call CC_get_coeff_single(i_spin,iii,aaa,ct1)
      Call CC_get_coeff_single(i_spin,jjj,bbb,ct2)
      Call CC_get_coeff_single(i_spin,jjj,aaa,ct3)
      Call CC_get_coeff_single(i_spin,iii,bbb,ct4)
      Call CC_get_coeff_double(iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin,ct5)

      tau = ct5 + ct1 * ct2 - ct3 * ct4

      CC_E_corr = CC_E_corr + int_rlt * tau

    end do
  end do

  if (n_spin.eq.1) then
    CC_E_corr = CC_E_corr * 2.0D0
  end if

  i_state = CC_index_w(myid+1,3) - 1

  Call CC_decode_index_ab(i_state,iii,jjj,aaa,bbb)

  do i_run = 1, CC_mem_w(myid+1,3)
  
    Call CC_inc_index_ab(iii,jjj,aaa,bbb)

    Call CC_Calc_integral(int_rlt,iii,jjj,aaa,bbb,1,2,1,2)
    
    Call CC_get_coeff_single(1,iii,aaa,ct1)
    Call CC_get_coeff_single(2,jjj,bbb,ct2)
    Call CC_get_coeff_double(iii,jjj,aaa,bbb,1,2,1,2,ct5)

    tau = ct5 + ct1 * ct2 
    
    CC_E_corr = CC_E_corr + int_rlt * tau

  end do

  Call sync_real_number(CC_E_corr)

!  write(40+myid,*) 'end of E_corr'
  End Subroutine CC_Calc_E_corr
