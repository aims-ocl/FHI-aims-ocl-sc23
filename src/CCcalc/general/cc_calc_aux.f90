!  PURPOSE
!
!  This file contains the subroutines for calculating auxiliary vetors needed 
!  in the Coupled-Cluster (CC) FCI calculations
!
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!    This file was written by Tonghao Shen
!  SOURCE

  Subroutine CC_Calc_aux_vect()

! This subroutine calculates the auxiliary vectors

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_tmp

  ! Clear auxiliary vectors
  CC_h_ij_a = 0.0D0
  CC_h_ba_a = 0.0D0
  CC_h_bj_a = 0.0D0
  CC_g_ca_a = 0.0D0
  CC_g_ik_a = 0.0D0
  CC_aux_a_aa = 0.0D0
  CC_aux_a_ab = 0.0D0
  CC_aux_b_aa = 0.0D0
  CC_aux_b_ab = 0.0D0
  CC_aux_h_aaaa = 0.0D0
  CC_aux_h_abab = 0.0D0
  CC_aux_h_abba = 0.0D0

  if (n_spin.ne.1) then
    CC_h_ij_b = 0.0D0
    CC_h_ba_b = 0.0D0
    CC_h_bj_b = 0.0D0
    CC_g_ca_b = 0.0D0
    CC_g_ik_b = 0.0D0
    CC_aux_a_bb = 0.0D0
    CC_aux_b_bb = 0.0D0
    CC_aux_h_baab = 0.0D0
    CC_aux_h_baba = 0.0D0
    CC_aux_h_bbbb = 0.0D0
  end if

! Calculate auxiliary vectors
  Call CC_Calc_aux_h_ij(1)
  Call CC_Calc_aux_h_ba(1)
  Call CC_Calc_aux_h_bj(1)

  i_tmp = CC_n_occ(1) * CC_n_occ(1)
  Call sync_vector(CC_h_ij_a,i_tmp)

  i_tmp = CC_n_vir(1) * CC_n_vir(1)
  Call sync_vector(CC_h_ba_a,i_tmp)

  i_tmp = CC_n_vir(1) * CC_n_occ(1)
  Call sync_vector(CC_h_bj_a,i_tmp)

  Call CC_Calc_aux_g_ca(1)
  Call CC_Calc_aux_g_ik(1)

  Call CC_Calc_aux_a_xx(1)
  Call CC_Calc_aux_a_ab()
  Call CC_Calc_aux_b_xx(1)
  Call CC_Calc_aux_b_ab()

  Call CC_Calc_aux_h_icak_aaaa(1)
  Call CC_Calc_aux_h_icak_abab(1)
  Call CC_Calc_aux_h_icak_abba(1)

! Synchronize vectors
  i_tmp = CC_n_vir(1) * CC_n_vir(1)
  Call sync_vector(CC_g_ca_a,i_tmp)

  i_tmp = CC_n_occ(1) * CC_n_occ(1)
  Call sync_vector(CC_g_ik_a,i_tmp)

  if (n_spin.ne.1) then
! Calculate beta part for open shell molecules
    Call CC_Calc_aux_h_ij(2)
    Call CC_Calc_aux_h_ba(2)
    Call CC_Calc_aux_h_bj(2)

    i_tmp = CC_n_occ(2) * CC_n_occ(2)
    Call sync_vector(CC_h_ij_b,i_tmp)

    i_tmp = CC_n_vir(2) * CC_n_vir(2)
    Call sync_vector(CC_h_ba_b,i_tmp)

    i_tmp = CC_n_vir(2) * CC_n_occ(2)
    Call sync_vector(CC_h_bj_b,i_tmp)

    Call CC_Calc_aux_g_ca(2)
    Call CC_Calc_aux_g_ik(2)

    Call CC_Calc_aux_a_xx(2)
    Call CC_Calc_aux_b_xx(2)

    Call CC_Calc_aux_h_icak_aaaa(2)
    Call CC_Calc_aux_h_icak_abab(2)
    Call CC_Calc_aux_h_icak_abba(2)

! Synchronize vectors
    i_tmp = CC_n_vir(2) * CC_n_vir(2)
    Call sync_vector(CC_g_ca_b,i_tmp)

    i_tmp = CC_n_occ(2) * CC_n_occ(2)
    Call sync_vector(CC_g_ik_b,i_tmp)

  end if


  End Subroutine CC_Calc_aux_vect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_aux_h_ij(i_spin)

! This subroutine calculates h_ij vector
! (iii,jjj)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer :: i_spin,o_spin
  Integer (kind = CC_ip) :: i_state,i_tmp
  Integer :: iii,jjj,kkk,ccc,ddd,lll,iinx,jinx
  Double precision :: int_rlt,rlt,tau,ct1,ct2,ct3,ct4,ct5

  o_spin = 3 - i_spin
  i_state = CC_index_h_ij(myid+1,i_spin)

!  jjj = Mod(i_state,CC_n_occ(i_spin))
!  if ((jjj.eq.0).and.(i_state.ne.0)) then
!    jjj = CC_n_occ(i_spin)
!  end if
!  iii = Int((i_state - jjj) / CC_n_occ(i_spin)) + CC_valence
!  jjj = jjj + CC_valence - 1
   Call CC_decode_ij4h(i_state,iii,jjj,i_spin)
   jjj = jjj - 1

!  print*,'h_ij'
  do i_tmp = 1, CC_mem_h_ij(myid+1,i_spin)
!    jjj = jjj + 1
!    if (jjj.gt.CC_n_elec(i_spin)) then
!      iii = iii + 1
!      jjj = 1
!    end if

    Call CC_inc_ij4h(iii,jjj,i_spin)

    rlt = 0.0D0
    ! Calculate h_ij(iii,jjj)
    ! first term, (aaaa) for i_spin = 1, (bbbb) for i_spin = 2
    do lll = CC_valence, CC_n_elec(i_spin)
      if (jjj.ne.lll) then

        do ccc = CC_n_elec(i_spin) + 1, n_states - 1
          do ddd = ccc + 1, n_states
            Call CC_Calc_integral(int_rlt,jjj,lll,ccc,ddd,i_spin,i_spin,i_spin,i_spin)

            Call CC_get_coeff_single(i_spin,iii,ccc,ct1)
            Call CC_get_coeff_single(i_spin,lll,ddd,ct2)
            Call CC_get_coeff_single(i_spin,iii,ddd,ct3)
            Call CC_get_coeff_single(i_spin,lll,ccc,ct4)

            if (iii.eq.lll) then
              ! iii = lll, epv term
              tau = ct1 * ct2 - ct3 * ct4
            else
              ! iii /= lll, non-epv term
              Call CC_get_coeff_double(iii,lll,ccc,ddd,i_spin,i_spin,i_spin,i_spin,ct5)
              tau = ct5 + ct1 * ct2 - ct3 * ct4
            end if

            rlt = rlt + int_rlt * tau
          end do
        end do

      end if
    end do

    ! second term (abab)
    do lll = CC_valence, CC_n_elec(o_spin)
      do ccc = CC_n_elec(i_spin) + 1, n_states
        do ddd = CC_n_elec(o_spin) + 1, n_states

          Call CC_Calc_integral(int_rlt,jjj,lll,ccc,ddd,i_spin,o_spin,i_spin,o_spin)

          Call CC_get_coeff_single(i_spin,iii,ccc,ct1)
          Call CC_get_coeff_single(o_spin,lll,ddd,ct2)
          Call CC_get_coeff_double(iii,lll,ccc,ddd,i_spin,o_spin,i_spin,o_spin,ct3)

          tau = ct3 + ct1 * ct2

          rlt = rlt + int_rlt * tau
        end do
      end do
    end do

    iinx = iii - CC_valence + 1
    jinx = jjj - CC_valence + 1

!    print*,iii,jjj,rlt
    if (i_spin.eq.1) then
      CC_h_ij_a(iinx,jinx) = rlt
    else
      CC_h_ij_b(iinx,jinx) = rlt
    end if

  end do

  End Subroutine CC_Calc_aux_h_ij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_aux_h_ba(i_spin)

! This subroutine calculates the auxiliary vector h_ba
! (bbb,aaa)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None

  Integer (kind = CC_ip) :: i_state,i_run
  Integer :: i_spin,o_spin,aaa,bbb,kkk,lll,ddd,ainx,binx
  Double precision :: rlt,int_rlt,tau,ct1,ct2,ct3,ct4,ct5

  o_spin = 3 - i_spin
  i_state = CC_index_h_ba(myid+1,i_spin)

!  ainx = Mod(i_state,CC_n_vir(i_spin))
!  if (ainx.eq.0) then
!    ainx = CC_n_vir(i_spin)
!  end if
!  binx = Int((i_state - ainx) / CC_n_vir(i_spin)) + 1
!  ainx = ainx - 1

  Call CC_decode_ab4h(i_state,aaa,bbb,i_spin)
  aaa = aaa - 1

!  print*,'h_ba'
  do i_run = 1, CC_mem_h_ba(myid+1,i_spin)
!    ainx = ainx + 1
!    if (ainx.gt.CC_n_vir(i_spin)) then
!      ainx = 1
!      binx = binx + 1
!    end if

    Call CC_inc_ab4h(aaa,bbb,i_spin)

    ainx = aaa - CC_n_elec(i_spin)
    binx = bbb - CC_n_elec(i_spin)
!    aaa = ainx + CC_n_elec(i_spin)
!    bbb = binx + CC_n_elec(i_spin)

!    print*,'binx,ainx',binx,ainx

    rlt = 0.0D0

    ! First term, (aaaa) for i_spin = 1, (bbbb) for i_spin = 2
    do ddd = CC_n_elec(i_spin) + 1, n_states
      if (ddd.ne.bbb) then

        do kkk = CC_valence, CC_n_elec(i_spin) - 1
          do lll = kkk + 1, CC_n_elec(i_spin)
            Call CC_Calc_integral(int_rlt,kkk,lll,bbb,ddd,i_spin,i_spin,i_spin,i_spin)
            Call CC_get_coeff_single(i_spin,kkk,aaa,ct1)
            Call CC_get_coeff_single(i_spin,lll,ddd,ct2)
            Call CC_get_coeff_single(i_spin,lll,aaa,ct3)
            Call CC_get_coeff_single(i_spin,kkk,ddd,ct4)

            if (ddd.eq.aaa) then
              ! ddd = aaa, epv term
              tau = ct1 * ct2 - ct3 * ct4
            else 
              ! non-epv term
              Call CC_get_coeff_double(kkk,lll,aaa,ddd,i_spin,i_spin,i_spin,i_spin,ct5)
              tau = ct5 + ct1 * ct2 - ct3 * ct4
            end if

            rlt = rlt - int_rlt * tau
          end do
        end do
      end if
    end do

    ! Second term, (abab)
    do ddd = CC_n_elec(o_spin) + 1, n_states
      do kkk = CC_valence, CC_n_elec(i_spin)
        do lll = CC_valence, CC_n_elec(o_spin)
          Call CC_Calc_integral(int_rlt,kkk,lll,bbb,ddd,i_spin,o_spin,i_spin,o_spin)

          Call CC_get_coeff_single(i_spin,kkk,aaa,ct1)
          Call CC_get_coeff_single(o_spin,lll,ddd,ct2)
          Call CC_get_coeff_double(kkk,lll,aaa,ddd,i_spin,o_spin,i_spin,o_spin,ct5)
          tau = ct5 + ct1 * ct2

          rlt = rlt - int_rlt * tau

        end do
      end do
    end do

!    print*,bbb,aaa,rlt
    if (i_spin.eq.1) then
      CC_h_ba_a(binx,ainx) = rlt
    else
      CC_h_ba_b(binx,ainx) = rlt
    end if

  end do

  End Subroutine CC_Calc_aux_h_ba

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_aux_h_bj(i_spin)
! This subrtoutine calculates auxiliary vector h_bj
! (bbb,jjj)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp
  Integer :: i_spin,o_spin,jjj,kkk,bbb,ccc,binx,jinx
  Double precision :: rlt,int_rlt,ct1,ct2,ct3,ct4,ct5

  o_spin = 3 - i_spin
  i_state = CC_index_h_bj(myid+1,i_spin)

  Call CC_decode_bj4h(i_state,bbb,jjj,i_spin)

!  jjj = Mod(i_state,CC_n_elec(i_spin))
!  if (jjj.eq.0) then
!    jjj = CC_n_elec(i_spin)
!  end if
!  binx = Int((i_state - jjj) / CC_n_elec(i_spin)) + 1
  jjj = jjj - 1

!  print*, 'h_bj'
  do i_tmp = 1, CC_mem_h_bj(myid+1,i_spin)

!    jjj = jjj + 1
!    if (jjj.gt.CC_n_elec(i_spin)) then
!      jjj = 1
!      binx = binx + 1
!    end if

!    bbb = binx + CC_n_elec(i_spin)

    Call CC_inc_bj4h(bbb,jjj,i_spin)

    binx = bbb - CC_n_elec(i_spin)
    jinx = jjj - CC_valence + 1

!    print*,'binx,jjj',binx,jjj
    rlt = 0.0D0
    ! First term (same spin part)
    do kkk = CC_valence, CC_n_elec(i_spin)
      if (kkk.ne.jjj) then
        do ccc = CC_n_elec(i_spin) + 1, n_states
          if (ccc.ne.bbb) then
            Call CC_Calc_integral(int_rlt,jjj,kkk,bbb,ccc,i_spin,i_spin,i_spin,i_spin)
            Call CC_get_coeff_single(i_spin,kkk,ccc,ct1)
            rlt = rlt + int_rlt * ct1

          end if
        end do
      end if
    end do

    ! Second term (opposite spin part)
    do kkk = CC_valence, CC_n_elec(o_spin)
      do ccc = CC_n_elec(o_spin) + 1, n_states
        Call CC_Calc_integral(int_rlt,jjj,kkk,bbb,ccc,i_spin,o_spin,i_spin,o_spin)
        Call CC_get_coeff_single(o_spin,kkk,ccc,ct1)
        rlt = rlt + int_rlt * ct1
      end do
    end do

    if (i_spin.eq.1) then
      CC_h_bj_a(binx,jinx) = rlt
    else
      CC_h_bj_b(binx,jinx) = rlt
    end if
!    print*,bbb,jjj,rlt
  end do

  End Subroutine CC_Calc_aux_h_bj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_aux_g_ca(i_spin)
! This subroutine calculates auxiliary vector g_ca
! (ccc,aaa)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp
  Integer :: i_spin,o_spin,aaa,kkk,ccc,ddd,ainx,cinx
  Double precision :: rlt,int_rlt,ct1,ct2,ct3,ct4,ct5

  o_spin = 3 - i_spin
  i_state = CC_index_g_ca(myid+1,i_spin)

  Call CC_decode_ab4h(i_state,aaa,ccc,i_spin)
  aaa = aaa - 1

!  ainx = Mod(i_state,CC_n_vir(i_spin))
!  if (ainx.eq.0) then
!    ainx = CC_n_vir(i_spin)
!  end if

!  cinx = Int((i_state - ainx) / CC_n_vir(i_spin)) + 1
!  ainx = ainx - 1

!  print*,'g_ca'
  do i_tmp = 1, CC_mem_g_ca(myid+1,i_spin)

!    ainx = ainx + 1
!    if (ainx.gt.CC_n_vir(i_spin)) then
!      ainx = 1
!      cinx = cinx + 1
!    end if

    Call CC_inc_ab4h(aaa,ccc,i_spin)

    cinx = ccc - CC_n_elec(i_spin)
    ainx = aaa - CC_n_elec(i_spin) 

!    aaa = ainx + CC_n_elec(i_spin)
!    ccc = cinx + CC_n_elec(i_spin)
!    print*,'cinx,ainx',cinx,ainx

    rlt = 0.0D0
    ! First term (same spin part)
    do kkk = CC_valence, CC_n_elec(i_spin)
      do ddd = CC_n_elec(i_spin) + 1, n_states
        if (ddd.ne.ccc) then
          Call CC_Calc_integral(int_rlt,aaa,kkk,ccc,ddd,i_spin,i_spin,i_spin,i_spin)
          Call CC_get_coeff_single(i_spin,kkk,ddd,ct1)

          rlt = rlt + int_rlt * ct1

        end if
      end do
    end do

    ! Second term (opposite spin part)
    do kkk = CC_valence, CC_n_elec(o_spin)
      do ddd = CC_n_elec(o_spin) + 1, n_states
        Call CC_Calc_integral(int_rlt,aaa,kkk,ccc,ddd,i_spin,o_spin,i_spin,o_spin)
        Call CC_get_coeff_single(o_spin,kkk,ddd,ct1)

        rlt = rlt + int_rlt * ct1
      end do
    end do

!    print*,ccc,aaa,CC_h_ba_a(cinx,ainx) + rlt
!    print*,'h_ba(cinx,ainx)',CC_h_ba_a(cinx,ainx)
    if (i_spin.eq.1) then
      CC_g_ca_a(cinx,ainx) = CC_h_ba_a(cinx,ainx) + rlt
    else
      CC_g_ca_b(cinx,ainx) = CC_h_ba_b(cinx,ainx) + rlt
    end if

  end do

  End Subroutine CC_Calc_aux_g_ca

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_aux_g_ik(i_spin)
! This subroutine calculates auxiliary vector g_ik
! (iii,kkk)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp
  Integer :: i_spin,o_spin,aaa,iii,kkk,ccc,ddd,lll,iinx,kinx
  Double precision :: rlt,int_rlt,ct1,ct2,ct3,ct4,ct5

  o_spin = 3 - i_spin
  i_state = CC_index_g_ik(myid+1,i_spin)

!  kkk = Mod(i_state,CC_n_elec(i_spin))
!  if (kkk.eq.0) then
!    kkk = CC_n_elec(i_spin)
!  end if
!  iii = Int((i_state - kkk) / CC_n_elec(i_spin)) + 1
!  kkk = kkk - 1

  Call CC_decode_ij4h(i_state,iii,kkk,i_spin)
  kkk = kkk - 1

!  print*,'g_ik'
  do i_tmp = 1, CC_mem_g_ik(myid+1,i_spin)
!    kkk = kkk + 1
!    if (kkk.gt.CC_n_elec(i_spin)) then
!      kkk = 1
!      iii = iii + 1
!    end if

    Call CC_inc_ij4h(iii,kkk,i_spin)

    iinx = iii - CC_valence + 1
    kinx = kkk - CC_valence + 1


!    print*,'iii,kkk',iii,kkk
    rlt = 0.0D0

    ! First term (same spin part)
    do lll = CC_valence, CC_n_elec(i_spin)
      if (lll.ne.kkk) then
        do ccc = CC_n_elec(i_spin) + 1, n_states
          Call CC_Calc_integral(int_rlt,kkk,lll,iii,ccc,i_spin,i_spin,i_spin,i_spin)
          Call CC_get_coeff_single(i_spin,lll,ccc,ct1)
          rlt = rlt + int_rlt * ct1
        end do
      end if
    end do

    ! Second term (opposite spin part)
    do lll = CC_valence, CC_n_elec(o_spin)
      do ccc = CC_n_elec(o_spin) + 1, n_states
        Call CC_Calc_integral(int_rlt,kkk,lll,iii,ccc,i_spin,o_spin,i_spin,o_spin)
        Call CC_get_coeff_single(o_spin,lll,ccc,ct1)
        rlt = rlt + int_rlt * ct1
      end do
    end do

!    print*,iii,kkk,CC_h_ij_a(iii,kkk) + rlt
!    print*,'h_ij(iii,kkk)',CC_h_ij_a(iii,kkk)
    if (i_spin.eq.1) then
      CC_g_ik_a(iinx,kinx) = CC_h_ij_a(iinx,kinx) + rlt
    else
      CC_g_ik_b(iinx,kinx) = CC_h_ij_b(iinx,kinx) + rlt
    end if

  end do

  End Subroutine CC_Calc_aux_g_ik

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_aux_a_xx(i_spin)
! This subroutine calculates the aa(bb) part of auxiliary vector a_ijkl
! (kkk,lll,iii,jjj)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp
  Integer :: i_spin,o_spin,iii,jjj,kkk,lll,ccc,ddd,cinx,dinx
  Double precision :: rlt,int_rlt,tau,ct1,ct2,ct3,ct4,ct5

  o_spin = 3 - i_spin
  ! for i_spin = 1, calculations start at CC_index_a(myid+1,2)
  ! for i_spin = 2, calculations start at CC_index_a(myid+1,3)
  ! also see subroutine CC_initialization() for definitions of these indices
  i_state = CC_index_a(myid+1,i_spin+1) - 1
 
  do i_tmp = 1, CC_mem_a(myid+1,i_spin+1) 

    i_state = i_state + 1
    Call CC_aux_decode_a_xx(i_state,i_spin,kkk,lll,iii,jjj)
    rlt = 0.0D0

    Call CC_Calc_integral(rlt,iii,jjj,kkk,lll,i_spin,i_spin,i_spin,i_spin)

    ! Firt term
    do ccc = CC_n_elec(i_spin) + 1, n_states
      Call CC_Calc_integral(int_rlt,iii,ccc,kkk,lll,i_spin,i_spin,i_spin,i_spin)
      Call CC_get_coeff_single(i_spin,jjj,ccc,ct1)
      rlt = rlt + int_rlt * ct1
    end do

    ! Second term
    do ccc = CC_n_elec(i_spin) + 1, n_states
      Call CC_Calc_integral(int_rlt,jjj,ccc,kkk,lll,i_spin,i_spin,i_spin,i_spin)
      Call CC_get_coeff_single(i_spin,iii,ccc,ct1)
      rlt = rlt - int_rlt * ct1
    end do

    ! Third term
    do ccc = CC_n_elec(i_spin) + 1, n_states - 1
      do ddd = ccc + 1, n_states
        Call CC_Calc_integral(int_rlt,kkk,lll,ccc,ddd,i_spin,i_spin,i_spin,i_spin)
        Call CC_get_coeff_single(i_spin,iii,ccc,ct1)
        Call CC_get_coeff_single(i_spin,jjj,ddd,ct2)
        Call CC_get_coeff_single(i_spin,jjj,ccc,ct3)
        Call CC_get_coeff_single(i_spin,iii,ddd,ct4)
        Call CC_get_coeff_double(iii,jjj,ccc,ddd,i_spin,i_spin,i_spin,i_spin,ct5)
        tau = ct5 + ct1 * ct2 - ct3 * ct4
        rlt = rlt + int_rlt * tau
      end do
    end do

    if (i_spin.eq.1) then
      CC_aux_a_aa(i_tmp) = rlt
    else
      CC_aux_a_bb(i_tmp) = rlt
    end if

  end do

  End Subroutine CC_Calc_aux_a_xx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_aux_a_ab()
! This subroutine calculates auxiliary vector a in alpha-beta space
! (kkk,lll,iii,jjj)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp
  Integer :: iii,jjj,kkk,lll,ccc,ddd,cinx,dinx
  Double precision :: rlt,int_rlt,tau,ct1,ct2,ct3,ct4,ct5

  i_state = CC_index_a(myid+1,1) - 1
!  Call CC_aux_decode_a_ab(i_state,kkk,lll,iii,jjj)

!  print*,'a_ijkl'
  do i_tmp = 1, CC_mem_a(myid+1,1)
!    jjj = jjj + 1
!    if (jjj.gt.CC_n_elec(2)) then
!      jjj = 1
!      iii = iii + 1
!    end if

!    if (iii.gt.CC_n_elec(1)) then
!      iii = 1
!      lll = lll + 1
!    end if

!    if (lll.gt.CC_n_elec(2)) then
!      lll = 1
!      kkk = kkk + 1
!    end if
    i_state = i_state + 1
    Call CC_aux_decode_a_ab(i_state,kkk,lll,iii,jjj)

!    print*,'iii,jjj,kkk,lll',iii,jjj,kkk,lll

    rlt = 0.0D0

    Call CC_Calc_integral(rlt,iii,jjj,kkk,lll,1,2,1,2)

    ! First term
    do ccc = CC_n_elec(2) + 1, n_states
      Call CC_Calc_integral(int_rlt,iii,ccc,kkk,lll,1,2,1,2)
      Call CC_get_coeff_single(2,jjj,ccc,ct1)
      rlt = rlt + int_rlt * ct1
    end do

    ! Second term
    do ccc = CC_n_elec(1) + 1, n_states
      Call CC_Calc_integral(int_rlt,jjj,ccc,kkk,lll,2,1,1,2)
      Call CC_get_coeff_single(1,iii,ccc,ct1)
      rlt = rlt - int_rlt * ct1
    end do

    ! Third term
    do ccc = CC_n_elec(1) + 1, n_states
      do ddd = CC_n_elec(2) + 1, n_states
        Call CC_Calc_integral(int_rlt,kkk,lll,ccc,ddd,1,2,1,2)
        Call CC_get_coeff_single(1,iii,ccc,ct1)
        Call CC_get_coeff_single(2,jjj,ddd,ct2)
        Call CC_get_coeff_double(iii,jjj,ccc,ddd,1,2,1,2,ct5)
        tau = ct5 + ct1 * ct2
        rlt = rlt + int_rlt * tau
      end do
    end do
!    print*,iii,jjj,kkk,lll,i_tmp,rlt
    CC_aux_a_ab(i_tmp) = rlt
  end do


  End Subroutine CC_Calc_aux_a_ab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_aux_b_xx(i_spin)
! This subroutine calculates auxiliary vector b in alpha-alpha(beta-beta) space.
! (ccc,ddd,aaa,bbb)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp
  Integer :: i_spin,o_spin,aaa,bbb,ccc,ddd,kkk
  Double precision :: rlt,int_rlt,ct1,ct2,ct3,ct4,ct5

  o_spin = 3 - i_spin
  i_state = CC_index_b(myid+1,i_spin+1) - 1

  do i_tmp = 1, CC_mem_b(myid+1,i_spin+1)
    i_state = i_state + 1
    Call CC_aux_decode_b_xx(i_state,i_spin,ccc,ddd,aaa,bbb)

    rlt = 0.0D0
    Call CC_Calc_integral(rlt,aaa,bbb,ccc,ddd,i_spin,i_spin,i_spin,i_spin)

    ! First term
    do kkk = CC_valence, CC_n_elec(i_spin)
      Call CC_Calc_integral(int_rlt,aaa,kkk,ccc,ddd,i_spin,i_spin,i_spin,i_spin)
      Call CC_get_coeff_single(i_spin,kkk,bbb,ct1)
      rlt = rlt - int_rlt * ct1
    end do

    ! Second term
    do kkk = CC_valence, CC_n_elec(i_spin)
      Call CC_Calc_integral(int_rlt,bbb,kkk,ccc,ddd,i_spin,i_spin,i_spin,i_spin)
      Call CC_get_coeff_single(i_spin,kkk,aaa,ct1)
      rlt = rlt + int_rlt * ct1
    end do

    if (i_spin.eq.1) then
      CC_aux_b_aa(i_tmp) = rlt
    else
      CC_aux_b_bb(i_tmp) = rlt
    end if
  end do

  End Subroutine CC_Calc_aux_b_xx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_aux_b_ab()
! This subroutine calculates auxiliary vector b in alpha-beta space
! (ccc,ddd,aaa,bbb)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp
  Integer :: aaa,bbb,ccc,ddd,kkk,ainx,binx,cinx,dinx
  Double precision :: rlt,int_rlt,ct1

  i_state = CC_index_b(myid+1,1) - 1
!  Call CC_aux_decode_b_ab(i_state,ccc,ddd,aaa,bbb)

!  print*, 'b_cdab'
  do i_tmp = 1, CC_mem_b(myid+1,1)
!    bbb = bbb + 1
!    if (bbb.gt.n_states) then
!      bbb = CC_n_elec(2) + 1
!      aaa = aaaa + 1
!    end if

!    if (aaa.gt.n_states) then
!      aaa = CC_n_elec(1) + 1
!      ddd = ddd + 1
!    end if

!    if (ddd.gt.n_states) then
!      ddd = CC_n_elec(2) + 1
!      ccc = ccc + 1
!    end if
    i_state = i_state + 1
    Call CC_aux_decode_b_ab(i_state,ccc,ddd,aaa,bbb)

!    write(31,*) 'ccc,ddd,aaa,bbb',ccc,ddd,aaa,bbb
    rlt = 0.0D0
    Call CC_Calc_integral(rlt,aaa,bbb,ccc,ddd,1,2,1,2)

    ! First term
    do kkk = CC_valence, CC_n_elec(2)
      Call CC_Calc_integral(int_rlt,aaa,kkk,ccc,ddd,1,2,1,2)
      Call CC_get_coeff_single(2,kkk,bbb,ct1)
      rlt = rlt - int_rlt * ct1
    end do

    ! Second term
    do kkk = CC_valence, CC_n_elec(1)
      Call CC_Calc_integral(int_rlt,bbb,kkk,ccc,ddd,2,1,1,2)
      Call CC_get_coeff_single(1,kkk,aaa,ct1)
      rlt = rlt + int_rlt * ct1
    end do

!    print*, aaa,bbb,ccc,ddd,i_tmp,rlt
!    print*,'rlt',rlt
    CC_aux_b_ab(i_tmp) = rlt
  end do

  End Subroutine CC_Calc_aux_b_ab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_aux_h_icak_aaaa(i_spin)
! This subroutine calculates auxiliary vector h_icak in alpha-alpha-alpha-alpha(aaaa)
! or in beta-beta-beta-beta(bbbb) space
! (kkk,ccc,iii,aaa)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp
  Integer :: i_spin,o_spin,iii,kkk,lll,aaa,ccc,ddd,ainx,cinx,dinx,j_tmp
  Double precision :: rlt,int_rlt,tau,ct1,ct2,ct3,ct4,ct5

  o_spin = 3 - i_spin
  if (i_spin.eq.1) then
  ! in aaaa space
    j_tmp = 1
  else
  ! in bbbb space
    j_tmp = 6
  end if

  i_state = CC_index_h(myid+1,j_tmp)
  Call CC_aux_decode_h(i_state,iii,ccc,aaa,kkk, & 
                       i_spin,i_spin,i_spin,i_spin,CC_n_elec,CC_n_vir)

  aaa = aaa - 1
!  print*, 'h_aaaa'
  do i_tmp = 1, CC_mem_h(myid+1,j_tmp)
!    aaa = aaa + 1
!    if (aaa.gt.n_states) then
!      aaa = CC_n_elec(i_spin) + 1
!      iii = iii + 1
!    end if

!    if (iii.gt.CC_n_elec(i_spin)) then
!      iii = 1
!      ccc = ccc + 1
!    end if

!    if (ccc.gt.n_states) then
!      ccc = CC_n_elec(i_spin) + 1
!      kkk = kkk + 1
!    end if

    Call CC_inc_icak4h(iii,ccc,aaa,kkk,i_spin,i_spin,i_spin,i_spin)

!    ainx = aaa - CC_n_elec(i_spin)
!    cinx = ccc - CC_n_elec(i_spin)

!    write(31,*) 'iii,ccc,aaa,kkk',iii,ccc,aaa,kkk

! Calculation starts
    rlt = 0.0D0
    Call CC_Calc_integral(rlt,iii,ccc,aaa,kkk,i_spin,i_spin,i_spin,i_spin)

    ! First term
    do lll = CC_valence, CC_n_elec(i_spin)
      if (lll.ne.kkk) then
        Call CC_Calc_integral(int_rlt,iii,ccc,lll,kkk,i_spin,i_spin,i_spin,i_spin)
        Call CC_get_coeff_single(i_spin,lll,aaa,ct1)
        rlt = rlt - int_rlt * ct1
      end if
    end do

    ! Second term
    do ddd = CC_n_elec(i_spin) + 1, n_states
      if (ddd.ne.ccc) then
        Call CC_Calc_integral(int_rlt,ddd,ccc,aaa,kkk,i_spin,i_spin,i_spin,i_spin)
        Call CC_get_coeff_single(i_spin,iii,ddd,ct1)
        rlt = rlt + int_rlt * ct1
      end if
    end do

    ! Third term
    do lll = CC_valence, CC_n_elec(i_spin)
      if (lll.ne.kkk) then
        do ddd = CC_n_elec(i_spin) + 1, n_states
          if (ddd.ne.ccc) then
            Call CC_Calc_integral(int_rlt,kkk,lll,ccc,ddd,i_spin,i_spin,i_spin,i_spin)
            Call CC_get_coeff_single(i_spin,iii,ddd,ct1)
            Call CC_get_coeff_single(i_spin,lll,aaa,ct2)

            if ((iii.eq.lll).or.(ddd.eq.aaa)) then
              ! epv term
              tau = ct1 * ct2
            else
              ! non-epv term
              Call CC_get_coeff_double(iii,lll,ddd,aaa,i_spin,i_spin,i_spin,i_spin,ct5)
              tau = 0.5D0 * ct5 + ct1 * ct2
            end if

            rlt = rlt - int_rlt * tau
          end if
        end do
      end if
    end do

    ! Forth term
    do lll = CC_valence, CC_n_elec(o_spin)
      do ddd = CC_n_elec(o_spin)+1, n_states
        Call CC_Calc_integral(int_rlt,kkk,lll,ccc,ddd,i_spin,o_spin,i_spin,o_spin)
        Call CC_get_coeff_double(iii,lll,ddd,aaa,i_spin,o_spin,o_spin,i_spin,ct5)
        tau = 0.5D0 * ct5
        rlt = rlt - int_rlt * tau

      end do
    end do

!    write(31,*) 'i_tmp',i_tmp
!    print*,iii,ccc,aaa,kkk,i_tmp,rlt
    if (i_spin.eq.1) then
      CC_aux_h_aaaa(i_tmp) = rlt
    else
      CC_aux_h_bbbb(i_tmp) = rlt
    end if

  end do

  End Subroutine CC_Calc_aux_h_icak_aaaa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_aux_h_icak_abab(i_spin)
! This subroutine calculates alpha-beta-alpha-beta part (and baba part) of h_icak aux array
! (kkk,ccc,iii,aaa)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp
  Integer :: i_spin,o_spin,iii,kkk,lll,aaa,ccc,ddd,j_tmp
  Double precision :: rlt,int_rlt,tau,coef_sign,ct1,ct2,ct3,ct4,ct5

  o_spin = 3 - i_spin
  if (i_spin.eq.1) then
  ! in aaaa space
    j_tmp = 2
  else
  ! in bbbb space
    j_tmp = 5
  end if

!  print*, 'h_abab'

  i_state = CC_index_h(myid+1,j_tmp)
  Call CC_aux_decode_h(i_state,iii,ccc,aaa,kkk, & 
                       i_spin,o_spin,i_spin,o_spin,CC_n_elec,CC_n_vir)

  aaa = aaa - 1
  do i_tmp = 1, CC_mem_h(myid+1,j_tmp)
!    aaa = aaa + 1
!    if (aaa.gt.n_states) then
!      aaa = CC_n_elec(i_spin) + 1
!      iii = iii + 1
!    end if

!    if (iii.gt.CC_n_elec(i_spin)) then
!      iii = 1
!      ccc = ccc + 1
!    end if

!    if (ccc.gt.n_states) then
!      ccc = CC_n_elec(o_spin) + 1
!      kkk = kkk + 1
!    end if

!    write(31,*) 'iii,ccc,aaa,kkk',iii,ccc,aaa,kkk

    Call CC_inc_icak4h(iii,ccc,aaa,kkk,i_spin,o_spin,i_spin,o_spin)

! Calculation starts
    rlt = 0.0D0
    Call CC_Calc_integral(rlt,iii,ccc,aaa,kkk,i_spin,o_spin,i_spin,o_spin)

    ! First term
    do lll = CC_valence, CC_n_elec(i_spin)
      Call CC_Calc_integral(int_rlt,iii,ccc,lll,kkk,i_spin,o_spin,i_spin,o_spin)
      Call CC_get_coeff_single(i_spin,lll,aaa,ct1)
      rlt = rlt - int_rlt * ct1
    end do

    ! Second term
    do ddd = CC_n_elec(i_spin) + 1, n_states
      Call CC_Calc_integral(int_rlt,ddd,ccc,aaa,kkk,i_spin,o_spin,i_spin,o_spin)
      Call CC_get_coeff_single(i_spin,iii,ddd,ct1)
      rlt = rlt + int_rlt * ct1
    end do

    ! Third term
    do lll = CC_valence, CC_n_elec(i_spin)
      do ddd = CC_n_elec(i_spin) + 1, n_states
        Call CC_Calc_integral(int_rlt,kkk,lll,ccc,ddd,o_spin,i_spin,o_spin,i_spin)
        Call CC_get_coeff_single(i_spin,iii,ddd,ct1)
        Call CC_get_coeff_single(i_spin,lll,aaa,ct2)

        if ((iii.eq.lll).or.(ddd.eq.aaa)) then
        ! epv term
          tau = ct1 * ct2
        else
        ! non-epv term
          Call CC_get_coeff_double(iii,lll,ddd,aaa,i_spin,i_spin,i_spin,i_spin,ct5)
          tau = 0.5D0 * ct5 + ct1 * ct2
        end if

        rlt = rlt - int_rlt * tau
      end do
    end do

    ! Forth term
    do lll = CC_valence, CC_n_elec(o_spin)
      if (lll.ne.kkk) then
        do ddd = CC_n_elec(o_spin)+1, n_states
          if (ddd.ne.ccc) then
            Call CC_Calc_integral(int_rlt,kkk,lll,ccc,ddd,o_spin,o_spin,o_spin,o_spin)
            Call CC_get_coeff_double(iii,lll,ddd,aaa,i_spin,o_spin,o_spin,i_spin,ct5)

            tau = 0.5D0 * ct5
            rlt = rlt - int_rlt * tau

          end if
        end do
      end if
    end do

!    write(31,*) 'i_tmp',i_tmp
!    print*,iii,ccc,aaa,kkk,i_tmp,rlt

    if (i_spin.eq.1) then
      CC_aux_h_abab(i_tmp) = rlt
    else
      CC_aux_h_baba(i_tmp) = rlt
    end if

  end do

  End Subroutine CC_Calc_aux_h_icak_abab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine  CC_Calc_aux_h_icak_abba(i_spin)

! This subroutine calculates alpha-beta-beta-alpha part (and baab part) of h_icak aux array
! (kkk,ccc,iii,aaa)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp
  Integer :: i_spin,o_spin,iii,kkk,lll,aaa,ccc,ddd,j_tmp
  Double precision :: rlt,int_rlt,tau,ct1,ct2,ct3,ct4,ct5

  o_spin = 3 - i_spin
  if (i_spin.eq.1) then
  ! in aaaa space
    j_tmp = 3
  else
  ! in bbbb space
    j_tmp = 4
  end if

!  print*,'h_abba'
  i_state = CC_index_h(myid+1,j_tmp)
  Call CC_aux_decode_h(i_state,iii,ccc,aaa,kkk, & 
                       i_spin,o_spin,o_spin,i_spin,CC_n_elec,CC_n_vir)

  aaa = aaa - 1
  do i_tmp = 1, CC_mem_h(myid+1,j_tmp)
!    aaa = aaa + 1
!    if (aaa.gt.n_states) then
!      aaa = CC_n_elec(o_spin) + 1
!      iii = iii + 1
!    end if

!    if (iii.gt.CC_n_elec(i_spin)) then
!      iii = 1
!      ccc = ccc + 1
!    end if

!    if (ccc.gt.n_states) then
!      ccc = CC_n_elec(o_spin) + 1
!      kkk = kkk + 1
!    end if

    Call CC_inc_icak4h(iii,ccc,aaa,kkk,i_spin,o_spin,o_spin,i_spin)

! Calculation starts
    rlt = 0.0D0
    Call CC_Calc_integral(rlt,iii,ccc,aaa,kkk,i_spin,o_spin,o_spin,i_spin)

    ! First term
    do lll = CC_valence, CC_n_elec(o_spin)
      Call CC_Calc_integral(int_rlt,iii,ccc,lll,kkk,i_spin,o_spin,o_spin,i_spin)
      Call CC_get_coeff_single(o_spin,lll,aaa,ct1)
      rlt = rlt - int_rlt * ct1
    end do

    ! Second term
    do ddd = CC_n_elec(i_spin) + 1, n_states
      Call CC_Calc_integral(int_rlt,ddd,ccc,aaa,kkk,i_spin,o_spin,o_spin,i_spin)
      Call CC_get_coeff_single(i_spin,iii,ddd,ct1)
      rlt = rlt + int_rlt * ct1
    end do

    ! Third term
    do lll = CC_valence, CC_n_elec(o_spin)
      do ddd = CC_n_elec(i_spin) + 1, n_states
        Call CC_Calc_integral(int_rlt,kkk,lll,ccc,ddd,i_spin,o_spin,o_spin,i_spin)
        Call CC_get_coeff_single(i_spin,iii,ddd,ct1)
        Call CC_get_coeff_single(o_spin,lll,aaa,ct2)
        Call CC_get_coeff_double(iii,lll,ddd,aaa,i_spin,o_spin,i_spin,o_spin,ct5)

        tau = 0.5D0 * ct5 + ct1 * ct2
        rlt = rlt - int_rlt * tau
      end do
    end do

!    write(31,*) 'i_tmp',i_tmp
!    print*,iii,ccc,aaa,kkk,i_tmp,rlt

    if (i_spin.eq.1) then
      CC_aux_h_abba(i_tmp) = rlt
    else
      CC_aux_h_baab(i_tmp) = rlt
    end if

  end do


  End Subroutine CC_Calc_aux_h_icak_abba
