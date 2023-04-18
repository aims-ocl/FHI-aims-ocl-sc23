  Subroutine CC_calc_w()
! This subroutine calculates the working vector w

  Use dimensions
  Use CC_corr

  Implicit None

  Call CC_calc_w_1a(1)
  Call CC_calc_w_2a(1)
  Call CC_calc_w_ab()

!print*, 'CC_w_1a',CC_w_1a

  if (n_spin.ne.1) then
    Call CC_calc_w_1a(2)
    Call CC_calc_w_2a(2)
  end if

  End Subroutine cc_calc_w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_calc_w_1a(i_spin)

! This subroutine calculates w_1a(i_spin = 1) and w_1b(i_spin = 2)
! (iii,ainx)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None

  Integer :: i_spin,o_spin,iii,aaa,jjj,bbb,kkk,ccc,ainx,binx,cinx,iinx,jinx
  Integer (kind = CC_ip) :: i_state,i_run
  Double precision :: rlt,int_rlt,tau,ct1,ct2,ct3,ct4,ct5,rlt2

  o_spin = 3 - i_spin

  i_state = CC_index_w(myid+1,i_spin) - 1
  Call CC_decode_index_s(i_state,i_spin,iii,aaa)

!  print*,'w_1a'
  do i_run = 1, CC_mem_w(myid+1,i_spin)
    Call CC_inc_index_s(i_spin,iii,aaa)
    ainx = aaa - CC_n_elec(i_spin)
    iinx = iii - CC_valence + 1
!    print*,'iii,aaa',iii,aaa

    rlt = 0.0D0

    ! First term
    do binx = 1, CC_n_vir(i_spin)
      if (i_spin.eq.1) then
        rlt = rlt + CC_h_ba_a(binx,ainx) * CC_t_1a(iinx,binx)
      else
        rlt = rlt + CC_h_ba_b(binx,ainx) * CC_t_1b(iinx,binx)
      end if
    end do

!    print*,'sum hca',rlt
!    rlt2 = rlt

    ! Second term
    do jinx = 1, CC_n_occ(i_spin)
      if (i_spin.eq.1) then
        rlt = rlt - CC_h_ij_a(iinx,jinx) * CC_t_1a(jinx,ainx)
      else
        rlt = rlt - CC_h_ij_b(iinx,jinx) * CC_t_1b(jinx,ainx)
      end if
    end do

!    print*,'sum hik', rlt2 - rlt
!    rlt2 = rlt

    ! Third term
    do jjj = CC_valence, CC_n_elec(i_spin)
      do binx = 1, CC_n_vir(i_spin)
        bbb = binx + CC_n_elec(i_spin)
        jinx = jjj - CC_valence + 1
        Call CC_get_coeff_single(i_spin,iii,bbb,ct1)
        Call CC_get_coeff_single(i_spin,jjj,aaa,ct2)
        if ((iii.eq.jjj).or.(ainx.eq.binx)) then
        ! epv term
          tau = ct1 * ct2
        else
        ! non-epv term
          Call CC_get_coeff_double(iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin,ct3)
          tau = ct3 + ct1 * ct2
        end if

        if (i_spin.eq.1) then
          rlt = rlt + CC_h_bj_a(binx,jinx) * tau
        else
          rlt = rlt + CC_h_bj_b(binx,jinx) * tau
        end if

      end do
    end do

    ! Forth term
    do jjj = CC_valence, CC_n_elec(o_spin)
      do binx = 1, CC_n_vir(o_spin)
        bbb = binx + CC_n_elec(o_spin)
        jinx = jjj - CC_valence + 1
        Call CC_get_coeff_double(iii,jjj,aaa,bbb,i_spin,o_spin,i_spin,o_spin,ct3)
        if ((o_spin.eq.1).or.(n_spin.eq.1)) then
          rlt = rlt + CC_h_bj_a(binx,jinx) * ct3
        else
          rlt = rlt + CC_h_bj_b(binx,jinx) * ct3
        end if
      end do
    end do

!    print*,'sum hck',rlt - rlt2
!    rlt2 = rlt

    ! 5-th term
    do jjj = CC_valence, CC_n_elec(i_spin)
      do bbb = CC_n_elec(i_spin) + 1, n_states
        Call CC_Calc_integral(int_rlt,iii,bbb,aaa,jjj,i_spin,i_spin,i_spin,i_spin)
        Call CC_get_coeff_single(i_spin,jjj,bbb,ct1)
        rlt = rlt + int_rlt * ct1
      end do
    end do

    ! 6-th term
    do jjj = CC_valence, CC_n_elec(o_spin)
      do bbb = CC_n_elec(o_spin) + 1, n_states
        Call CC_Calc_integral(int_rlt,iii,bbb,aaa,jjj,i_spin,o_spin,i_spin,o_spin)
        Call CC_get_coeff_single(o_spin,jjj,bbb,ct1)
        rlt = rlt + int_rlt * ct1
      end do
    end do

!    print*,'sum icak',rlt - rlt2
!    rlt2 = rlt

    ! 7-th term
    do jjj = CC_valence, CC_n_elec(i_spin)
      do bbb = CC_n_elec(i_spin) + 1, n_states - 1
        do ccc = bbb + 1, n_states
          Call CC_Calc_integral(int_rlt,jjj,aaa,bbb,ccc,i_spin,i_spin,i_spin,i_spin)
          Call CC_get_coeff_single(i_spin,iii,bbb,ct1)
          Call CC_get_coeff_single(i_spin,jjj,ccc,ct2)
          Call CC_get_coeff_single(i_spin,iii,ccc,ct3)
          Call CC_get_coeff_single(i_spin,jjj,bbb,ct4)

          if (jjj.eq.iii) then
          ! epv term
            tau = ct1 * ct2 - ct3 * ct4
          else
          ! non-epv term
            Call CC_get_coeff_double(iii,jjj,bbb,ccc,i_spin,i_spin,i_spin,i_spin,ct5)
            tau = ct5 + ct1 * ct2 - ct3 * ct4
          end if

          rlt = rlt - int_rlt * tau

        end do
      end do
    end do

    ! 8-th term
    do jjj = CC_valence, CC_n_elec(o_spin)
      do bbb = CC_n_elec(i_spin) + 1, n_states
        do ccc = CC_n_elec(o_spin) + 1, n_states

          Call CC_Calc_integral(int_rlt,jjj,aaa,bbb,ccc,o_spin,i_spin,i_spin,o_spin)
          Call CC_get_coeff_single(i_spin,iii,bbb,ct1)
          Call CC_get_coeff_single(o_spin,jjj,ccc,ct2)
          Call CC_get_coeff_double(iii,jjj,bbb,ccc,i_spin,o_spin,i_spin,o_spin,ct5)

          tau = ct5 + ct1 * ct2
          rlt = rlt - int_rlt * tau
        end do
      end do
    end do

!    print*,'sum kacd',rlt - rlt2
!    rlt2 = rlt

    ! 9-th term
    do jjj = CC_valence, CC_n_elec(i_spin) - 1
      do kkk = jjj + 1, CC_n_elec(i_spin)
        do bbb = CC_n_elec(i_spin) + 1, n_states

          Call CC_Calc_integral(int_rlt,jjj,kkk,iii,bbb,i_spin,i_spin,i_spin,i_spin)
          Call CC_get_coeff_single(i_spin,jjj,aaa,ct1)
          Call CC_get_coeff_single(i_spin,kkk,bbb,ct2)
          Call CC_get_coeff_single(i_spin,kkk,aaa,ct3)
          Call CC_get_coeff_single(i_spin,jjj,bbb,ct4)

          if (bbb.eq.aaa) then
          ! epv term
            tau = ct1 * ct2 - ct3 * ct4
          else
          ! non-epv term
            Call CC_get_coeff_double(jjj,kkk,aaa,bbb,i_spin,i_spin,i_spin,i_spin,ct5)
            tau = ct5 + ct1 * ct2 - ct3 * ct4
          end if

          rlt = rlt - int_rlt * tau
        end do
      end do
    end do

    ! 10-th term
    do jjj = CC_valence, CC_n_elec(i_spin)
      do kkk = CC_valence, CC_n_elec(o_spin)
        do bbb = CC_n_elec(o_spin) + 1, n_states

          Call CC_Calc_integral(int_rlt,jjj,kkk,iii,bbb,i_spin,o_spin,i_spin,o_spin)
          Call CC_get_coeff_single(i_spin,jjj,aaa,ct1)
          Call CC_get_coeff_single(o_spin,kkk,bbb,ct2)
          Call CC_get_coeff_double(jjj,kkk,aaa,bbb,i_spin,o_spin,i_spin,o_spin,ct5)
          tau = ct5 + ct1 * ct2 

          rlt = rlt - int_rlt * tau
        end do
      end do
    end do

!    print*,'sum klci', rlt2 - rlt
!    rlt2 = rlt


    if (i_spin.eq.1) then
      CC_w_1a(i_run) = rlt
    else
      CC_w_1b(i_run) = rlt
    end if

  end do


  End Subroutine CC_calc_w_1a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_calc_w_2a(i_spin)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_run,task_state
  Integer :: i_spin,o_spin,iii,jjj,aaa,bbb
  Integer :: i_task,errnum
  Double precision :: rlt,h_rlt,a_rlt,b_rlt
  Double precision , dimension(:), allocatable :: CC_w_local

  o_spin = 3 - i_spin

! First step, calculate the terms including auxiliary vectors which are saved separately
  do i_task = 1, n_tasks
    task_state = CC_mem_w(i_task,i_spin+3)

    Allocate(CC_w_local(task_state),stat=errnum)
    Call check_allocation(errnum, 'CC_w_local')

    i_state = CC_index_w(i_task,i_spin+3) - 1

!   Start calculation
    do i_run = 1, task_state
      i_state = i_state + 1

      Call CC_decode_index_xx(i_state,i_spin,iii,jjj,aaa,bbb)

      rlt = 0.0D0

!     Calculate first term (a_ijkl term)
      Call CC_sum_aux_a_xx(a_rlt,i_spin,i_task,iii,jjj,aaa,bbb)
      rlt = rlt + a_rlt
  
!     Calculate second term (b_cdab term)
      Call CC_sum_aux_b_xx(b_rlt,i_spin,i_task,iii,jjj,aaa,bbb)
      rlt = rlt + b_rlt

!     Calculate third term (h_icak_aaaa for i_spin = 1, and h_icak_bbbb for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,iii,jjj,aaa,bbb, &
                        i_spin,i_spin,i_spin,i_spin,i_spin,i_spin)
      rlt = rlt + h_rlt

!     Calculate 4-th term (h_icak_abab for i_spin = 1, and h_icak_baba for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,iii,jjj,aaa,bbb, &
                        i_spin,i_spin,i_spin,i_spin,o_spin,o_spin)
      rlt = rlt + h_rlt

!     Calculate 5-th term (h_icak_aaaa for i_spin = 1, and h_icak_bbbb for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,iii,jjj,bbb,aaa, &
                        i_spin,i_spin,i_spin,i_spin,i_spin,i_spin)
      rlt = rlt - h_rlt

!     Calculate 6-th term (h_icak_abab for i_spin = 1, and h_icak_baba for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,iii,jjj,bbb,aaa, &
                        i_spin,i_spin,i_spin,i_spin,o_spin,o_spin)
      rlt = rlt - h_rlt

!     Calculate 7-th term (h_icak_aaaa for i_spin = 1, and h_icak_bbbb for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,jjj,iii,bbb,aaa, &
                        i_spin,i_spin,i_spin,i_spin,i_spin,i_spin)
      rlt = rlt + h_rlt

!     Calculate 8-th term (h_icak_abab for i_spin = 1, and h_icak_baba for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,jjj,iii,bbb,aaa, &
                        i_spin,i_spin,i_spin,i_spin,o_spin,o_spin)
      rlt = rlt + h_rlt

!     Calculate 9-th term (h_icak_aaaa for i_spin = 1, and h_icak_bbbb for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,jjj,iii,aaa,bbb, &
                        i_spin,i_spin,i_spin,i_spin,i_spin,i_spin)
      rlt = rlt - h_rlt

!     Calculate 10-th term (h_icak_abab for i_spin = 1, and h_icak_baba for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,jjj,iii,aaa,bbb, &
                        i_spin,i_spin,i_spin,i_spin,o_spin,o_spin)
      rlt = rlt - h_rlt

      CC_w_local(i_run) = rlt

    end do

    Call sync_vector(CC_w_local,task_state)

    if (i_task.eq.myid+1) then
      if (i_spin.eq.1) then
        CC_w_2a = CC_w_local
      else
        CC_w_2b = CC_w_local
      end if
    end if

    Deallocate(CC_w_local)
  end do

! Second step, calculate the terms that saved in all cores

  task_state = CC_mem_w(myid+1,i_spin+3)

  i_state = CC_index_w(myid+1,i_spin+3) - 1

! Start calculation
  do i_run = 1, task_state
    i_state = i_state + 1

    Call CC_decode_index_xx(i_state,i_spin,iii,jjj,aaa,bbb)

    rlt = 0.0D0

    ! First term, sum(c) [g_ca * t(i,j,c,b)]
    Call CC_sum_aux_g_ca(h_rlt,iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin)
    rlt = rlt + h_rlt

    ! Second term, sum(c) [g_cb * t(i,j,c,a)]
    Call CC_sum_aux_g_ca(h_rlt,iii,jjj,bbb,aaa,i_spin,i_spin,i_spin,i_spin)
    rlt = rlt - h_rlt

    ! Third term, sum(k) [<ka||ij> * t(k,b)]
    Call CC_sum_intl_kaij(h_rlt,iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin)
    rlt = rlt + h_rlt

    ! 4-th term, sum(k) [<kb||ij> * t(k,a)]
    Call CC_sum_intl_kaij(h_rlt,iii,jjj,bbb,aaa,i_spin,i_spin,i_spin,i_spin)
    rlt = rlt - h_rlt

    ! 5-th term, sum(k) [g_ik * t(k,j,a,b)]
    Call CC_sum_aux_g_ik(h_rlt,iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin)
    rlt = rlt - h_rlt

    ! 6-th term, sum(k) [g_jk * t(k,i,a,b)]
    Call CC_sum_aux_g_ik(h_rlt,jjj,iii,aaa,bbb,i_spin,i_spin,i_spin,i_spin)
    rlt = rlt + h_rlt

    ! 7-th term, sum(c) [<ab||cj> * t(i,c)]
    Call CC_sum_intl_abcj(h_rlt,iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin)
    rlt = rlt + h_rlt

    ! 8-th term, sum(c) [<ab||ci> * t(j,c)]
    Call CC_sum_intl_abcj(h_rlt,jjj,iii,aaa,bbb,i_spin,i_spin,i_spin,i_spin)
    rlt = rlt - h_rlt

    ! 9-th term, sum(k,c) [<ic||ak> * t(j,c) * t(k,b)]
    Call CC_sum_intl_icak(h_rlt,iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin)
    rlt = rlt - h_rlt

    ! 10-th term, sum(k,c) [<ic||bk> * t(j,c) * t(k,a)]
    Call CC_sum_intl_icak(h_rlt,iii,jjj,bbb,aaa,i_spin,i_spin,i_spin,i_spin)
    rlt = rlt + h_rlt

    ! 11-th term, sum(k,c) [<jc||ak> * t(i,c) * t(k,b)]
    Call CC_sum_intl_icak(h_rlt,jjj,iii,aaa,bbb,i_spin,i_spin,i_spin,i_spin)
    rlt = rlt + h_rlt

    ! 12-th term, sum(k,c) [<jc||bk> * t(i,c) * t(k,a)]
    Call CC_sum_intl_icak(h_rlt,jjj,iii,bbb,aaa,i_spin,i_spin,i_spin,i_spin)
    rlt = rlt - h_rlt

    if (i_spin.eq.1) then
      CC_w_2a(i_run) = CC_w_2a(i_run) + rlt
    else
      CC_w_2b(i_run) = CC_w_2b(i_run) + rlt
    end if

  end do


  End Subroutine CC_calc_w_2a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_calc_w_ab()

! This Subroutine Calculates w vector in alpha-beta space

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None

  Integer (kind = CC_ip) :: i_state,task_state,i_run
  Integer :: i_spin,o_spin,iii,jjj,aaa,bbb,i_task,errnum
  Double precision :: rlt,h_rlt,a_rlt,b_rlt
  Double precision , dimension(:), allocatable :: CC_w_local

  i_spin = 1  
  o_spin = 2

!  print*,'first part in w_ab'
! First step, calculate the terms including auxiliary vectors which are saved separately
  do i_task = 1, n_tasks
    task_state = CC_mem_w(i_task,3)

    Allocate(CC_w_local(task_state),stat=errnum)
    Call check_allocation(errnum, 'CC_w_local')

    CC_w_local = 0.0D0
    i_state = CC_index_w(i_task,3) - 1

!   Start calculation
    do i_run = 1, task_state
      i_state = i_state + 1

      Call CC_decode_index_ab(i_state,iii,jjj,aaa,bbb)

!      print*,'iii,jjj,aaa,bbb',iii,jjj,aaa,bbb
      rlt = 0.0D0

!     Calculate first term (a_ijkl term)
      Call CC_sum_aux_a_ab(a_rlt,i_task,iii,jjj,aaa,bbb)
      rlt = rlt + a_rlt

!      print*,'sum_a',a_rlt  
!     Calculate second term (b_cdab term)
      Call CC_sum_aux_b_ab(b_rlt,i_task,iii,jjj,aaa,bbb)
      rlt = rlt + b_rlt

!      print*,'sum_b',b_rlt
!     Calculate third term (h_icak_aaaa for i_spin = 1, and h_icak_bbbb for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,iii,jjj,aaa,bbb, &
                        i_spin,o_spin,i_spin,o_spin,i_spin,i_spin)
      rlt = rlt + h_rlt 

!     Calculate 4-th term (h_icak_abab for i_spin = 1, and h_icak_baba for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,iii,jjj,aaa,bbb, &
                        i_spin,o_spin,i_spin,o_spin,o_spin,o_spin)
      rlt = rlt + h_rlt

!     Calculate 5-th term (h_icak_abba for i_spin = 1, and h_icak_baab for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,iii,jjj,bbb,aaa, &
                        i_spin,o_spin,o_spin,i_spin,i_spin,o_spin)
      rlt = rlt - h_rlt

!     Calculate 6-th term (h_icak_abab for i_spin = 1, and h_icak_baba for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,jjj,iii,bbb,aaa, &
                        o_spin,i_spin,o_spin,i_spin,i_spin,i_spin)
      rlt = rlt + h_rlt

!     Calculate 7-th term (h_icak_aaaa for i_spin = 1, and h_icak_bbbb for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,jjj,iii,bbb,aaa, &
                        o_spin,i_spin,o_spin,i_spin,o_spin,o_spin)
      rlt = rlt + h_rlt

!     Calculate 8-th term (h_icak_abab for i_spin = 1, and h_icak_baba for i_spin = 2)
      Call CC_sum_aux_h(h_rlt,i_task,jjj,iii,aaa,bbb, &
                        o_spin,i_spin,i_spin,o_spin,o_spin,i_spin)
      rlt = rlt - h_rlt

      CC_w_local(i_run) = rlt

!      print*,'i_run',i_run
!      print*,'rlt',rlt
    end do

    Call sync_vector(CC_w_local,task_state)

    if (i_task.eq.myid+1) then
      CC_w_ab = CC_w_local
    end if

    Deallocate(CC_w_local)
  end do

! Second step, calculate the terms that saved in all cores

  task_state = CC_mem_w(myid+1,3)

  i_state = CC_index_w(myid+1,3) - 1

!  print*,'second part in w_ab'
! Start calculation
  do i_run = 1, task_state
    i_state = i_state + 1

    Call CC_decode_index_ab(i_state,iii,jjj,aaa,bbb)

!    print*,'iii,jjj,aaa,bbb',iii,jjj,aaa,bbb
    rlt = 0.0D0

    ! First term, sum(c) [g_ca * t(i,j,c,b)]
    Call CC_sum_aux_g_ca(h_rlt,iii,jjj,aaa,bbb,i_spin,o_spin,i_spin,o_spin)
    rlt = rlt + h_rlt

!    print*,'sum gca',h_rlt
    ! Second term, sum(c) [g_cb * t(i,j,c,a)]
    Call CC_sum_aux_g_ca(h_rlt,iii,jjj,bbb,aaa,i_spin,o_spin,o_spin,i_spin)
    rlt = rlt - h_rlt
!    print*,'sum gca',h_rlt

    ! Third term, sum(k) [<ka||ij> * t(k,b)]
    Call CC_sum_intl_kaij(h_rlt,iii,jjj,aaa,bbb,i_spin,o_spin,i_spin,o_spin)
    rlt = rlt + h_rlt

    ! 4-th term, sum(k) [<kb||ij> * t(k,a)]
    Call CC_sum_intl_kaij(h_rlt,iii,jjj,bbb,aaa,i_spin,o_spin,o_spin,i_spin)
    rlt = rlt - h_rlt

    ! 5-th term, sum(k) [g_ik * t(k,j,a,b)]
    Call CC_sum_aux_g_ik(h_rlt,iii,jjj,aaa,bbb,i_spin,o_spin,i_spin,o_spin)
    rlt = rlt - h_rlt

!    print*,'sum_gik',h_rlt
    ! 6-th term, sum(k) [g_jk * t(k,i,a,b)]
    Call CC_sum_aux_g_ik(h_rlt,jjj,iii,aaa,bbb,o_spin,i_spin,i_spin,o_spin)
    rlt = rlt + h_rlt
!    print*,'sum_gik',h_rlt

    ! 7-th term, sum(c) [<ab||cj> * t(i,c)]
    Call CC_sum_intl_abcj(h_rlt,iii,jjj,aaa,bbb,i_spin,o_spin,i_spin,o_spin)
    rlt = rlt + h_rlt

    ! 8-th term, sum(c) [<ab||ci> * t(j,c)]
    Call CC_sum_intl_abcj(h_rlt,jjj,iii,aaa,bbb,o_spin,i_spin,i_spin,o_spin)
    rlt = rlt - h_rlt

    ! 9-th term, sum(k,c) [<ic||ak> * t(j,c) * t(k,b)]
    Call CC_sum_intl_icak(h_rlt,iii,jjj,aaa,bbb,i_spin,o_spin,i_spin,o_spin)
    rlt = rlt - h_rlt

    ! 10-th term, sum(k,c) [<ic||bk> * t(j,c) * t(k,a)]
    Call CC_sum_intl_icak(h_rlt,iii,jjj,bbb,aaa,i_spin,o_spin,o_spin,i_spin)
    rlt = rlt + h_rlt

    ! 11-th term, sum(k,c) [<jc||ak> * t(i,c) * t(k,b)]
    Call CC_sum_intl_icak(h_rlt,jjj,iii,aaa,bbb,o_spin,i_spin,i_spin,o_spin)
    rlt = rlt + h_rlt

    ! 12-th term, sum(k,c) [<jc||bk> * t(i,c) * t(k,a)]
    Call CC_sum_intl_icak(h_rlt,jjj,iii,bbb,aaa,o_spin,i_spin,o_spin,i_spin)
    rlt = rlt - h_rlt

    CC_w_ab(i_run) = CC_w_ab(i_run) + rlt
!    print*,'i_run',i_run
!    print*,'rlt',rlt
!    print*,'CC_w_ab(i_run)'
!    Call CC_Calc_integral(rlt,iii,jjj,aaa,bbb,1,2,1,2)
!    print*,iii,jjj,aaa,bbb,CC_w_ab(i_run)+rlt,rlt
  end do



  End Subroutine CC_calc_w_ab


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_sum_aux_a_xx(rlt,i_spin,i_task,iii,jjj,aaa,bbb)

! This subroutine calculates 

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None

  Integer (kind = CC_ip) :: i_step,task_aux,i_aux,i_local,code_kl,code_tp,code_ij
  Integer (kind = CC_ip) :: i_ijkl
  Integer :: i_spin,i_task,iii,jjj,kkk,lll,aaa,bbb,i_tmp,j_tmp
  Double precision :: rlt,ct1,ct2,ct3,ct4,ct5,tau
  Logical :: flag_circle

  rlt = 0.0D0
  i_step = (CC_n_occ(i_spin) - 1) * CC_n_occ(i_spin) / 2
  task_aux = CC_mem_a(myid+1,i_spin+1)
  i_aux = CC_index_a(myid+1,i_spin+1)

  Call CC_aux_decode_a_xx(i_aux,i_spin,kkk,lll,i_tmp,j_tmp)

  ! map kkk and lll to code_kl
!  code_kl = (lll - 2) * (lll - 1) / 2 + kkk
  Call CC_code_ij(code_kl,kkk,lll,i_spin,i_spin)

  ! map iii and jjj to code_ij
!  code_ij = (jjj - 2) * (jjj - 1) / 2 + iii
  Call CC_code_ij(code_ij,iii,jjj,i_spin,i_spin)

  ! map i_tmp and j_tmp to code_tp
!  code_tp = (j_tmp - 2) * (j_tmp - 1) / 2 + i_tmp
  Call CC_code_ij(code_tp,i_tmp,j_tmp,i_spin,i_spin)

  if (code_tp.gt.code_ij) then
    code_kl = code_kl + 1
    kkk = kkk + 1
    if (kkk.ge.lll) then
      lll = lll + 1
      kkk = CC_valence
    end if
  end if

  i_ijkl = (code_kl - 1) * i_step + code_ij
  i_local = i_ijkl + 1 - i_aux

  if (i_ijkl.le.(i_aux + task_aux - 1)) then
    flag_circle = .true.
  else
    flag_circle = .false.
  end if

!  print*,'sum a'
!  print*,'i_aux',i_aux
!  print*,'task_aux',task_aux
!  print*,'i_step',i_step

  do while (flag_circle)

!    print*,'i_ijkl',i_ijkl
!    print*,'iii,jjj,kkk,lll',iii,jjj,kkk,lll
    Call CC_get_coeff_single(i_spin,kkk,aaa,ct1)
    Call CC_get_coeff_single(i_spin,lll,bbb,ct2)
    Call CC_get_coeff_single(i_spin,lll,aaa,ct3)
    Call CC_get_coeff_single(i_spin,kkk,bbb,ct4)
    Call CC_get_coeff_double(kkk,lll,aaa,bbb,i_spin,i_spin,i_spin,i_spin,ct5)
    tau = ct5 + ct1 * ct2 - ct3 * ct4

    if (i_spin.eq.1) then
      rlt = rlt + CC_aux_a_aa(i_local) * tau
    else
      rlt = rlt + CC_aux_a_bb(i_local) * tau
    end if

    i_local = i_local + i_step
    kkk = kkk + 1
    if (kkk.ge.lll) then
      lll = lll + 1
      kkk = CC_valence
    end if

    if (i_local.gt.task_aux) then
      flag_circle = .false.
    end if
  end do

  End Subroutine CC_sum_aux_a_xx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_sum_aux_a_ab(rlt,i_task,iii,jjj,aaa,bbb)

! This subroutine calculates 

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None

  Integer (kind = CC_ip) :: i_step,task_aux,i_aux,i_local,code_kl,code_tp,code_ij
  INteger (kind = CC_ip) :: i_ijkl
  Integer :: i_spin,o_spin,i_task,aaa,bbb,iii,jjj,kkk,lll,i_tmp,j_tmp
  Double precision :: rlt,ct1,ct2,ct3,ct4,ct5,tau
  Logical :: flag_circle

  i_spin = 1
  o_spin = 2
  i_step = CC_n_occ(i_spin) * CC_n_occ(o_spin)
  task_aux = CC_mem_a(myid+1,1)
  i_aux = CC_index_a(myid+1,1)

  Call CC_aux_decode_a_ab(i_aux,kkk,lll,i_tmp,j_tmp)

  ! map kkk and lll to code_kl
!  code_kl = (kkk - 1) * CC_n_elec(o_spin) + lll
  Call CC_code_ij(code_kl,kkk,lll,i_spin,o_spin)

  ! map iii and jjj to code_ij
!  code_ij = (iii - 1) * CC_n_elec(o_spin) + jjj
  Call CC_code_ij(code_ij,iii,jjj,i_spin,o_spin)

  ! map i_tmp and j_tmp to code_tp
!  code_tp = (i_tmp - 1) * CC_n_elec(o_spin) + j_tmp
  Call CC_code_ij(code_tp,i_tmp,j_tmp,i_spin,o_spin)

  if (code_tp.gt.code_ij) then
    code_kl = code_kl + 1
    lll = lll + 1
    if (lll.gt.CC_n_elec(o_spin)) then
      lll = CC_valence
      kkk = kkk + 1
    end if
  end if

  i_ijkl = (code_kl - 1) * i_step + code_ij
  i_local = i_ijkl + 1 - i_aux

  if (i_ijkl.le.(i_aux + task_aux - 1)) then
    flag_circle = .true.
  else
    flag_circle = .false.
  end if

!  write(40+myid,*) 'sum_a_ab'
!  write(40+myid,*) 'i_task',myid+1
!  write(40+myid,*) 'task_aux',task_aux
!  write(40+myid,*) 'i_aux',i_aux
!  write(40+myid,*) 'code_ij',code_ij
!  write(40+myid,*) 'code_kl',code_kl
!  write(40+myid,*) 'code_tp',code_tp
!  print*,'sum_a_ab'
!  print*,'i_ijkl',i_ijkl
!  print*,'iii,jjj,aaa,bbb',iii,jjj,aaa,bbb
!  print*,'code_ij',code_ij
!  print*,'code_kl',code_kl
!  print*,'code_tp',code_tp

  rlt = 0.0D0
  do while (flag_circle)

!    print*,'kkk,lll,i_tmp,j_tmp',kkk,lll,i_tmp,j_tmp
!    print*,'i_ijkl',i_ijkl

    Call CC_get_coeff_single(i_spin,kkk,aaa,ct1)
    Call CC_get_coeff_single(o_spin,lll,bbb,ct2)
    Call CC_get_coeff_double(kkk,lll,aaa,bbb,i_spin,o_spin,i_spin,o_spin,ct5)
    tau = ct5 + ct1 * ct2

    rlt = rlt + CC_aux_a_ab(i_local) * tau

    i_local = i_local + i_step
    lll = lll + 1
    if (lll.gt.CC_n_elec(o_spin)) then
      lll = CC_valence
      kkk = kkk + 1
    end if

    if (i_local.gt.task_aux) then
      flag_circle = .false.
    end if
  end do

  End Subroutine CC_sum_aux_a_ab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_sum_aux_b_xx(rlt,i_spin,i_task,iii,jjj,aaa,bbb)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_step,task_aux,i_aux,i_local
  Integer (kind = CC_ip) :: code_tp,code_cd,code_ab,i_abcd
  Integer :: i_spin,i_task
  Integer :: iii,jjj,aaa,bbb,ccc,ddd,ainx,binx,cinx,dinx,a_tmp,b_tmp 
  Double precision :: rlt,ct1,ct2,ct3,ct4,ct5,tau
  Logical :: flag_circle

  rlt = 0.0D0
!  ainx = aaa - CC_n_elec(i_spin)
!  binx = bbb - CC_n_elec(i_spin)

  i_step = (CC_n_vir(i_spin) - 1) * CC_n_vir(i_spin) / 2
  task_aux = CC_mem_b(myid+1,i_spin+1)
  i_aux = CC_index_b(myid+1,i_spin+1)

  Call CC_aux_decode_b_xx(i_aux,i_spin,ccc,ddd,a_tmp,b_tmp)
!  cinx = ccc - CC_n_elec(i_spin)
!  dinx = ddd - CC_n_elec(i_spin)
!  a_tmp = a_tmp - CC_n_elec(i_spin)
!  b_tmp = b_tmp - CC_n_elec(i_spin)

  ! map ccc and ddd to code_cd
!  code_cd = (dinx - 2) * (dinx - 1) / 2 + cinx
  Call CC_code_ab(code_cd,ccc,ddd,i_spin,i_spin)

  ! map aaa and bbb to code_ab
!  code_ab = (binx - 2) * (binx - 1) / 2 + ainx
  Call CC_code_ab(code_ab,aaa,bbb,i_spin,i_spin)

  ! map a_tmp and b_tmp to code_tp
!  code_tp = (b_tmp - 2) * (b_tmp - 1) / 2 + a_tmp
  Call CC_code_ab(code_tp,a_tmp,b_tmp,i_spin,i_spin)

  if (code_tp.gt.code_ab) then
    code_cd = code_cd + 1
    ccc = ccc + 1
    if (ccc.ge.ddd) then
      ddd = ddd + 1
      ccc = CC_n_elec(i_spin) + 1
    end if
  end if

  i_abcd = (code_cd - 1) * i_step + code_ab
  i_local = i_abcd + 1 - i_aux

  if (i_abcd.le.(i_aux + task_aux - 1)) then
    flag_circle = .true.
  else
    flag_circle = .false.
  end if

!  ccc = cinx + CC_n_elec(i_spin)
!  ddd = dinx + CC_n_elec(i_spin)

  do while (flag_circle)

    Call CC_get_coeff_single(i_spin,iii,ccc,ct1)
    Call CC_get_coeff_single(i_spin,jjj,ddd,ct2)
    Call CC_get_coeff_single(i_spin,iii,ddd,ct3)
    Call CC_get_coeff_single(i_spin,jjj,ccc,ct4)
    Call CC_get_coeff_double(iii,jjj,ccc,ddd,i_spin,i_spin,i_spin,i_spin,ct5)
    tau = ct5 + ct1 * ct2 - ct3 * ct4

    if (i_spin.eq.1) then
      rlt = rlt + CC_aux_b_aa(i_local) * tau
    else
      rlt = rlt + CC_aux_b_bb(i_local) * tau
    end if

    i_local = i_local + i_step
    ccc = ccc + 1
    if (ccc.ge.ddd) then
      ddd = ddd + 1
      ccc = CC_n_elec(i_spin) + 1
    end if

    if (i_local.gt.task_aux) then
      flag_circle = .false.
    end if
  end do
 
  End Subroutine CC_sum_aux_b_xx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_sum_aux_b_ab(rlt,i_task,iii,jjj,aaa,bbb)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer :: i_spin,o_spin,i_task
  Integer (kind = CC_ip) :: i_step,task_aux,i_aux,i_local
  Integer :: iii,jjj,aaa,bbb,ccc,ddd,ainx,binx,cinx,dinx,a_tmp,b_tmp 
  Integer (kind = CC_ip) :: code_tp,code_cd,code_ab,i_abcd
  Double precision :: rlt,ct1,ct2,ct3,ct4,ct5,tau
  Logical :: flag_circle

  i_spin = 1
  o_spin = 2

!  ainx = aaa - CC_n_elec(i_spin)
!  binx = bbb - CC_n_elec(o_spin)

  i_step = CC_n_vir(i_spin) * CC_n_vir(o_spin)
  task_aux = CC_mem_b(myid+1,1)
  i_aux = CC_index_b(myid+1,1)

  Call CC_aux_decode_b_ab(i_aux,ccc,ddd,a_tmp,b_tmp)

  ! map ccc and ddd to code_cd
!  code_cd = (cinx - 1) * CC_n_vir(o_spin) + dinx
  Call CC_code_ab(code_cd,ccc,ddd,i_spin,o_spin)

  ! map aaa and bbb to code_ab
!  code_ab = (ainx - 1) * CC_n_vir(o_spin) + binx
  Call CC_code_ab(code_ab,aaa,bbb,i_spin,o_spin)

  ! map a_tmp and b_tmp to code_tp
!  code_tp = (a_tmp - 1) * CC_n_vir(o_spin) + b_tmp
  Call CC_code_ab(code_tp,a_tmp,b_tmp,i_spin,o_spin)

!   write(40+myid,*) 'sum_b_ab'
!   write(40+myid,*) 'i_task',myid+1
!   write(40+myid,*) 'task_aux',task_aux
!   write(40+myid,*) 'i_aux',i_aux
!   write(40+myid,*) 'code_ab',code_ab
!   write(40+myid,*) 'code_cd',code_cd
!   write(40+myid,*) 'code_tp',code_tp
!   write(40+myid,*) 'i_step',i_step
!  print*, 'sum_b'
!   write(40+myid,*) 'iii,jjj,aaa,bbb',iii,jjj,aaa,bbb
!  print*, 'i_step',i_step
!  print*, 'code_ab',code_ab
!  print*, 'code_cd',code_cd
!  print*, 'code_tp',code_tp

  if (code_tp.gt.code_ab) then
    code_cd = code_cd + 1
    ddd = ddd + 1
    if (ddd.gt.n_states) then
      ddd = CC_n_elec(o_spin) + 1
      ccc = ccc + 1
    end if
  end if

  i_abcd = (code_cd - 1) * i_step + code_ab
  i_local = i_abcd + 1 - i_aux
 
  if (i_abcd.le.(i_aux + task_aux - 1)) then
    flag_circle = .true.
  else
    flag_circle = .false.
  end if

!  ccc = cinx + CC_n_elec(i_spin)
!  ddd = dinx + CC_n_elec(o_spin)

  rlt = 0.0D0
  do while (flag_circle)

!    write(40+myid,*) 'ccc,ddd,aaa,bbb',ccc,ddd,aaa,bbb
!    write(40+myid,*) 'i_abcd',i_local

    Call CC_get_coeff_single(i_spin,iii,ccc,ct1)
    Call CC_get_coeff_single(o_spin,jjj,ddd,ct2)
    Call CC_get_coeff_double(iii,jjj,ccc,ddd,i_spin,o_spin,i_spin,o_spin,ct5)
    tau = ct5 + ct1 * ct2

    rlt = rlt + CC_aux_b_ab(i_local) * tau

    i_local = i_local + i_step
    ddd = ddd + 1
    if (ddd.gt.n_states) then
      ddd = CC_n_elec(o_spin) + 1
      ccc = ccc + 1
    end if

    if (i_local.gt.task_aux) then
      flag_circle = .false.
    end if
  end do
 
  End Subroutine CC_sum_aux_b_ab



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_sum_aux_h(rlt,i_task,iii,jjj,aaa,bbb, &
                          i_spin,j_spin,a_spin,b_spin,k_spin,c_spin)

! This subroutine calculates the summation of h_icak term
! rlt = sum(k,c) [h(i,c,a,k) * t(j,k,b,c)]

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None

  Integer (kind = CC_ip) :: i_step,task_aux,i_aux,i_icak,i_local
  Integer :: i_spin,j_spin,a_spin,b_spin,k_spin,c_spin,h_type
  Integer :: i_task,iii,jjj,kkk,aaa,bbb,ccc
  Double precision :: rlt,ct5
  Logical :: flag_circle

  i_step = CC_n_occ(i_spin) * CC_n_vir(a_spin)

  ! Find out the h vector involved in summation according to spin
  if (i_spin.eq.c_spin) then
    ! i_spin = c_spin = a_spin = k_spin, for h_aaaa or h_bbbb vector
    if ((i_spin.eq.1).or.(n_spin.eq.1)) then
      ! i_spin = c_spin = a_spin = k_spin = 1, or n_spin = 1, h_aaaa
      h_type = 1
    else
      ! i_spin = c_spin = a_spin = k_spin = 2, and n_spin = 2, h_bbbb
      h_type = 6
    end if
  else if (i_spin.eq.a_spin) then
    ! i_spin = a_spin /= c_spin = k_spin, for h_abab or h_baba vector
    if ((i_spin.eq.1).or.(n_spin.eq.1)) then
      ! i_spin = a_spin = 1, c_spin = k_spin = 2, or n_spin = 1, h_abab
      h_type = 2
    else
      ! i_spin = a_spin = 2, c_spin = k_spin = 1, and n_spin = 2, h_baba
      h_type = 5
    end if
  else
    ! i_spin = k_spin /= c_spin = a_spin, for h_abba or h_baab vector
    if ((i_spin.eq.1).or.(n_spin.eq.1)) then
      ! i_spin = k_spin = 1, c_spin = a_spin = 2, or n_spin = 1, h_abba
      h_type = 3
    else
      ! i_spin = k_spin = 2, c_spin = a_spin = 1, and n_spin = 2, h_baab
      h_type = 4
    end if
  end if

  task_aux = CC_mem_h(myid+1,h_type)
  i_aux = CC_index_h(myid+1,h_type)

  Call CC_aux_h_pointer(i_aux,iii,ccc,aaa,kkk,i_spin,c_spin,a_spin,k_spin,i_icak)
  i_local = i_icak + 1 - i_aux

  if (i_icak.le.(i_aux + task_aux - 1)) then
    flag_circle = .true.
  else
    flag_circle = .false.
  end if

!  print*, 'i_aux',i_aux
!  print*, 'task_aux',task_aux
!  write(32,*) 'sum_h'
!  write(32,*) 'icak_spin',i_spin,c_spin,a_spin,k_spin
!  write(32,*) 'h_type',h_type
!  print*, 'i_step',i_step
!  print*, 'iii,jjj,aaa,bbb',iii,jjj,aaa,bbb
!  print*, 'ijab_spin',i_spin,j_spin,a_spin,b_spin
!  print*, 'jkbc_spin',j_spin,k_spin,b_spin,c_spin
!  print*, 'i_icak',i_icak

  rlt = 0.0D0
  do while (flag_circle)
    if ((ccc.ne.bbb).or.(c_spin.ne.b_spin)) then
      if ((kkk.ne.jjj).or.(k_spin.ne.j_spin)) then
!        write(32,*) 'i_icak',i_icak
!        write(32,*) 'iii,ccc,aaa,kkk',iii,ccc,aaa,kkk
!        write(32,*) 'jjj,kkk,bbb,ccc',jjj,kkk,bbb,ccc

        Call CC_get_coeff_double(jjj,kkk,bbb,ccc,j_spin,k_spin,b_spin,c_spin,ct5)

        if (h_type.eq.1) then
          rlt = rlt + CC_aux_h_aaaa(i_local) * ct5
        else if (h_type.eq.2) then
          rlt = rlt + CC_aux_h_abab(i_local) * ct5
        else if (h_type.eq.3) then
          rlt = rlt + CC_aux_h_abba(i_local) * ct5
        else if (h_type.eq.4) then
          rlt = rlt + CC_aux_h_baab(i_local) * ct5
        else if (h_type.eq.5) then
          rlt = rlt + CC_aux_h_baba(i_local) * ct5
        else if (h_type.eq.6) then
          rlt = rlt + CC_aux_h_bbbb(i_local) * ct5
        end if
      end if
    end if

    i_local = i_local + i_step
    ccc = ccc + 1
    if (ccc.gt.n_states) then
      ccc = CC_n_elec(c_spin) + 1
      kkk = kkk + 1
    end if

    if (i_local.gt.task_aux) then
      flag_circle = .false.
    end if
  end do

  End Subroutine CC_sum_aux_h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_sum_aux_g_ca(rlt,iii,jjj,aaa,bbb,i_spin,j_spin,a_spin,b_spin)

! This subroutine calculates g_ca term
! rlt = sum(c) [g(c,a) * t(i,j,c,b)]

  Use dimensions
  Use CC_corr

  Implicit None

  Integer :: iii,jjj,aaa,bbb,ccc,ainx,binx,cinx,i_spin,j_spin,a_spin,b_spin
  Integer (kind = CC_ip) :: i_run
  Double precision :: rlt,ct1

  rlt = 0.0D0
  ainx = aaa - CC_n_elec(a_spin)

  do cinx = 1, CC_n_vir(a_spin)
    ccc = cinx + CC_n_elec(a_spin)
    if ((ccc.ne.bbb).or.(a_spin.ne.b_spin)) then
      Call CC_get_coeff_double(iii,jjj,ccc,bbb,i_spin,j_spin,a_spin,b_spin,ct1)
      if ((a_spin.eq.1).or.(n_spin.eq.1)) then
        rlt = rlt + CC_g_ca_a(cinx,ainx) * ct1
      else
        rlt = rlt + CC_g_ca_b(cinx,ainx) * ct1
      end if
    end if
  end do

  End Subroutine CC_sum_aux_g_ca

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_sum_intl_kaij(rlt,iii,jjj,aaa,bbb,i_spin,j_spin,a_spin,b_spin)
! This subroutine calculates the summation of <ka||ij> * t(k,b)
! rlt = sum(k) [<ka||ij> * t(k,b)]

  Use dimensions
  Use CC_corr

  Implicit None

  Integer :: iii,jjj,aaa,bbb,kkk,i_spin,j_spin,a_spin,b_spin
  Double precision :: rlt,int_rlt,ct1

  rlt = 0.0D0
  do kkk = CC_valence, CC_n_elec(b_spin)

    Call CC_Calc_integral(int_rlt,kkk,aaa,iii,jjj,b_spin,a_spin,i_spin,j_spin)
    Call CC_get_coeff_single(b_spin,kkk,bbb,ct1)

    rlt = rlt + int_rlt * ct1

  end do

  End Subroutine CC_sum_intl_kaij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_sum_aux_g_ik(rlt,iii,jjj,aaa,bbb,i_spin,j_spin,a_spin,b_spin)
! This subroutine calculates the summation of g_ik term
! rlt = sum(k) [g(i,k) * t(j,k,a,b)]

  Use dimensions
  Use CC_corr

  Implicit None

  Integer :: iii,jjj,kkk,aaa,bbb,i_spin,j_spin,a_spin,b_spin,iinx,kinx
  Double precision :: rlt,ct1

  rlt = 0.0D0
  iinx = iii - CC_valence + 1
  do kkk = CC_valence, CC_n_elec(i_spin)
    if ((kkk.ne.jjj).or.(i_spin.ne.j_spin)) then
      kinx = kkk - CC_valence + 1
      Call CC_get_coeff_double(kkk,jjj,aaa,bbb,i_spin,j_spin,a_spin,b_spin,ct1)
      if ((i_spin.eq.1).or.(n_spin.eq.1)) then
        rlt = rlt + CC_g_ik_a(iinx,kinx) * ct1
      else
        rlt = rlt + CC_g_ik_b(iinx,kinx) * ct1
      end if
    end if
  end do

  End Subroutine CC_sum_aux_g_ik

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_sum_intl_abcj(rlt,iii,jjj,aaa,bbb,i_spin,j_spin,a_spin,b_spin)

! This subroutine calculates the term <ab||cj> * t(i,c)
! rlt = sum(c) [<ab||cj> * t(i,c)]

  Use dimensions
  Use CC_corr

  Implicit None

  Integer :: iii,jjj,aaa,bbb,ccc,cinx,i_spin,j_spin,a_spin,b_spin
  Double precision :: rlt,ct1,int_rlt

  rlt = 0.0D0
  do ccc = CC_n_elec(i_spin) + 1, n_states

    Call CC_Calc_integral(int_rlt,aaa,bbb,ccc,jjj,a_spin,b_spin,i_spin,j_spin)
    Call CC_get_coeff_single(i_spin,iii,ccc,ct1)

    rlt = rlt + int_rlt * ct1

  end do

  End Subroutine CC_sum_intl_abcj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_sum_intl_icak(rlt,iii,jjj,aaa,bbb,i_spin,j_spin,a_spin,b_spin)

! This subroutine calculates the sum of <ic||ak> * t(j,c) * t(k,b)
! rlt = sum(k,c) [<ic||ak> * t(j,c) * t(k,b)]

  Use dimensions
  Use CC_corr

  Implicit None

  Integer :: iii,jjj,kkk,aaa,bbb,ccc,i_spin,j_spin,a_spin,b_spin
  Double precision :: rlt,ct1,ct2,int_rlt

  rlt = 0.0D0
  do kkk = CC_valence, CC_n_elec(b_spin)
    do ccc = CC_n_elec(j_spin) + 1, n_states

      Call CC_Calc_integral(int_rlt,iii,ccc,aaa,kkk,i_spin,j_spin,a_spin,b_spin)
      Call CC_get_coeff_single(j_spin,jjj,ccc,ct1)
      Call CC_get_coeff_single(b_spin,kkk,bbb,ct2)

      rlt = rlt + int_rlt * ct1 * ct2

    end do
  end do

  End Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

