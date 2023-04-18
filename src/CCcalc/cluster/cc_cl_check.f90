  Subroutine CC_cl_PT_check()

  Use CC_cl

  Implicit None
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,i_len,j_len,errnum
  Integer :: i_a,i_b,i_c,i_d
  Double precision :: sum_rlt,E_pre,D_ijk,delta_ijk,rlt

  Call CC_cl_w2t(2)

  Allocate(CC_PT_Mat_trans(CC_n_vir,CC_n_vir), stat=errnum)
  Call check_allocation(errnum,'CC_PT_Mat_trans in CC')

  Allocate(CC_t_check(CC_n_vir,CC_n_vir,CC_n_occ,CC_n_occ), stat=errnum)
  Call check_allocation(errnum,'CC_t_tmp in CC')

  CC_E_PT = 0.0D0

  i_len = CC_n_vir * CC_n_vir
  j_len = CC_n_vir * CC_n_vir

  print*,'CC_PT_check'
  print*,'t'
  do iii = 1, CC_n_occ
    do jjj = 1, CC_n_occ
      CC_PT_Mat_trans = 0.0D0
      Call CC_Mat_regular(i_len,CC_PT_Mat_trans, &
                          j_len,CC_t_d_T(:,iii,jjj,1),1,i_len)
      CC_t_check(:,:,iii,jjj) = transpose(CC_PT_Mat_trans)
    end do
  end do

  do aaa = 1, CC_n_vir
    do bbb = 1, CC_n_vir
      do iii = 1, CC_n_occ
        do jjj = 1, CC_n_occ
          print*,iii,jjj,aaa,bbb,CC_t_check(aaa,bbb,iii,jjj)
        end do
      end do
    end do
  end do

  Allocate(CC_intl_bdai_check(CC_n_vir,CC_n_vir,CC_n_vir,CC_n_occ), stat=errnum)
  Call check_allocation(errnum,'CC_t_tmp in CC')

  print*,'bdai'
  do ddd = 1, CC_n_vir
    do iii = 1, CC_n_occ
      do bbb = 1, CC_n_vir
        do aaa = 1, CC_n_vir
          i_d = ddd + CC_n_occ
          i_b = bbb + CC_n_occ
          i_a = aaa + CC_n_occ
          CC_intl_bdai_check(bbb,ddd,aaa,iii) = dot_product(CC_RI(:,i_b,i_d,1), &
                                                            CC_RI(:,i_a,iii,1))
          print*,ddd,iii,bbb,aaa,CC_intl_bdai_check(bbb,ddd,aaa,iii)
        end do
      end do
    end do
  end do

  Allocate(CC_intl_ckjl_check(CC_n_vir,CC_n_occ,CC_n_occ,CC_n_occ), stat=errnum)
  Call check_allocation(errnum,'CC_t_tmp in CC')

  print*,'ckjl'
  do lll = 1, CC_n_occ
    do jjj = 1, CC_n_occ
      do kkk = 1, CC_n_occ
        do ccc = 1, CC_n_vir
          i_c = ccc + CC_n_occ
          CC_intl_ckjl_check(ccc,kkk,jjj,lll) = dot_product(CC_RI(:,i_c,kkk,1), &
                                                            CC_RI(:,lll,jjj,1))
          print*,lll,jjj,kkk,ccc,CC_intl_ckjl_check(ccc,kkk,jjj,lll)
        end do
      end do
    end do
  end do

  Allocate(CC_intl_bjck_check(CC_n_vir,CC_n_occ,CC_n_vir,CC_n_occ), stat=errnum)
  Call check_allocation(errnum,'CC_t_tmp in CC')

  print*,'bjck'
  do jjj = 1, CC_n_occ
    do kkk = 1, CC_n_occ
      do bbb = 1, CC_n_vir
        do ccc = 1, CC_n_vir
          i_b = bbb + CC_n_occ
          i_c = ccc + CC_n_occ
          CC_intl_bjck_check(bbb,jjj,ccc,kkk) = dot_product(CC_RI(:,i_b,jjj,1), &
                                                            CC_RI(:,i_c,kkk,1))
          print*,jjj,kkk,bbb,ccc,CC_intl_bjck_check(bbb,jjj,ccc,kkk)
        end do
      end do
    end do
  end do


!  do iii = 1, CC_n_occ
!    do jjj = 1, CC_n_occ
!      do kkk = 1, CC_n_occ
!        E_pre = CC_E_PT
!        Call CC_cl_D_ijk(iii,jjj,kkk,delta_ijk,D_ijk)
!        do aaa = 1, CC_n_vir
!          do bbb = 1, CC_n_vir
!            do ccc = 1, CC_n_vir
!              Call CC_cl_check_calc_W(iii,jjj,kkk,aaa,bbb,ccc,sum_rlt)
!              CC_E_PT = CC_E_PT + sum_rlt
!            end do
!          end do
!        end do
!        print*,iii,jjj,kkk,CC_E_PT - E_pre,D_ijk
!      end do
!    end do
!  end do
!
!  CC_E_PT = CC_E_PT / 3.0D0

!  do iii = 1, CC_n_occ
!    do jjj = 1, iii
!      do kkk = 1, jjj
!        E_pre = CC_E_PT
!        Call CC_cl_D_ijk(iii,jjj,kkk,delta_ijk,D_ijk)
!        do aaa = 1, CC_n_vir
!          do bbb = 1, aaa
!            do ccc = 1, bbb
!              Call CC_cl_check_calc_W_2(iii,jjj,kkk,aaa,bbb,ccc,rlt)
!              CC_E_PT = CC_E_PT + rlt
!            end do
!          end do
!        end do
!        
!        !print*,iii,jjj,kkk,CC_E_PT - E_pre,D_ijk
!      end do
!    end do
!  end do

  do aaa = 1, CC_n_vir
    do bbb = 1, aaa
      do ccc = 1, bbb
        do iii = 1, CC_n_occ
          do jjj = 1, iii
            do kkk = 1, jjj
              E_pre = CC_E_PT
              Call CC_cl_D_ijk(iii,jjj,kkk,delta_ijk,D_ijk)
              Call CC_cl_check_calc_W_2(iii,jjj,kkk,aaa,bbb,ccc,rlt)
              CC_E_PT = CC_E_PT + rlt
              print*,'E_PT',CC_E_PT
            end do
          end do
        end do
        
        !print*,iii,jjj,kkk,CC_E_PT - E_pre,D_ijk
      end do
    end do
  end do

  Deallocate(CC_PT_Mat_trans,CC_t_check,CC_intl_bdai_check,CC_intl_ckjl_check,CC_intl_bjck_check)

  End Subroutine CC_cl_PT_check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_check_calc_W(iii,jjj,kkk,aaa,bbb,ccc,sum_rlt)

  Use CC_cl

  Implicit None
  Integer , intent(in) :: iii,jjj,kkk,aaa,bbb,ccc
  Double precision , intent(out) :: sum_rlt
  Double precision :: W_abc, W_bca, W_cab, W_cba, V_abc, V_cba
  Double precision :: D_abc, delta_ijk, D_ijk

  Call CC_cl_D_ijk(iii,jjj,kkk,delta_ijk,D_ijk)
  Call CC_cl_D_abc(aaa,bbb,ccc,D_abc)

  Call CC_cl_check_W_abc(iii,jjj,kkk,aaa,bbb,ccc,W_abc)
  Call CC_cl_check_W_abc(iii,jjj,kkk,bbb,ccc,aaa,W_bca)
  Call CC_cl_check_W_abc(iii,jjj,kkk,ccc,aaa,bbb,W_cab)
  Call CC_cl_check_W_abc(iii,jjj,kkk,ccc,bbb,aaa,W_cba)

  Call CC_cl_check_V_abc(iii,jjj,kkk,aaa,bbb,ccc,V_abc)
  Call CC_cl_check_V_abc(iii,jjj,kkk,ccc,bbb,aaa,V_cba)

  V_abc = V_abc + W_abc
  V_cba = V_cba + W_cba

  sum_rlt = (4.0D0 * W_abc + W_bca + W_cab) * (V_abc - V_cba) / (D_ijk - D_abc)

!  print*,iii,jjj,kkk,aaa,bbb,ccc
!  print*,'W_abc,W_bca,W_cab',W_abc,W_bca,W_cab
!  print*,'V_abc,V_cba',V_abc,V_cba
!  print*,'D_ijkabc',D_ijk - D_abc

  !if ((abs(W_abc).gt.1.0E-8).and.(iii.ge.jjj).and.(jjj.ge.kkk)) then
  !  write(60,*) iii,jjj,kkk,aaa,bbb,ccc,W_abc
  !!  print*,'V_abc',V_abc
  !!  print*,'D_ijkabc',D_ijk - D_abc
  !end if

  End Subroutine CC_cl_check_calc_W

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_check_calc_W_2(iii,jjj,kkk,aaa,bbb,ccc,sum_rlt)

  Use CC_cl

  Implicit None
  Integer , intent(in) :: iii,jjj,kkk,aaa,bbb,ccc
  Double precision , intent(out) :: sum_rlt
  Double precision :: W_abc, W_bca, W_cab, W_acb, W_bac, W_cba
  Double precision :: V_abc, V_bca, V_cab, V_acb, V_bac, V_cba
  Double precision :: D_abc, delta_ijk, D_ijk, delta_abc
  Double precision :: X,Y,Z,W1,W2

  Call CC_cl_D_ijk(iii,jjj,kkk,delta_ijk,D_ijk)
  Call CC_cl_D_abc(aaa,bbb,ccc,D_abc)

  delta_abc = 1.0D0
  if (aaa.eq.bbb) then
    delta_abc = delta_abc + 1.0D0
  end if

  if (bbb.eq.ccc) then
    delta_abc = delta_abc + 1.0D0
  end if

  print*,'term',iii,jjj,kkk,aaa,bbb,ccc

  Call CC_cl_check_W_abc(iii,jjj,kkk,aaa,bbb,ccc,W_abc)
  Call CC_cl_check_W_abc(iii,jjj,kkk,bbb,ccc,aaa,W_bca)
  Call CC_cl_check_W_abc(iii,jjj,kkk,ccc,aaa,bbb,W_cab)
  Call CC_cl_check_W_abc(iii,jjj,kkk,aaa,ccc,bbb,W_acb)
  Call CC_cl_check_W_abc(iii,jjj,kkk,bbb,aaa,ccc,W_bac)
  Call CC_cl_check_W_abc(iii,jjj,kkk,ccc,bbb,aaa,W_cba)

  !print*,'W'
  !print*,iii,jjj,kkk,aaa,bbb,ccc,W_abc
  !print*,kkk,iii,jjj,aaa,bbb,ccc,W_bca
  !print*,jjj,kkk,iii,aaa,bbb,ccc,W_cab
  !print*,iii,kkk,jjj,aaa,bbb,ccc,W_acb
  !print*,jjj,iii,kkk,aaa,bbb,ccc,W_bac
  !print*,kkk,jjj,iii,aaa,bbb,ccc,W_cba

  Call CC_cl_check_V_abc(iii,jjj,kkk,aaa,bbb,ccc,V_abc)
  Call CC_cl_check_V_abc(iii,jjj,kkk,bbb,ccc,aaa,V_bca)
  Call CC_cl_check_V_abc(iii,jjj,kkk,ccc,aaa,bbb,V_cab)
  Call CC_cl_check_V_abc(iii,jjj,kkk,aaa,ccc,bbb,V_acb)
  Call CC_cl_check_V_abc(iii,jjj,kkk,bbb,aaa,ccc,V_bac)
  Call CC_cl_check_V_abc(iii,jjj,kkk,ccc,bbb,aaa,V_cba)

  V_abc = (V_abc + W_abc) / delta_abc 
  V_bca = (V_bca + W_bca) / delta_abc 
  V_cab = (V_cab + W_cab) / delta_abc 
  V_acb = (V_acb + W_acb) / delta_abc 
  V_bac = (V_bac + W_bac) / delta_abc 
  V_cba = (V_cba + W_cba) / delta_abc 

  !print*,'V'
  !print*,iii,jjj,kkk,aaa,bbb,ccc,V_abc
  !print*,kkk,iii,jjj,aaa,bbb,ccc,V_bca
  !print*,jjj,kkk,iii,aaa,bbb,ccc,V_cab
  !print*,iii,kkk,jjj,aaa,bbb,ccc,V_acb
  !print*,jjj,iii,kkk,aaa,bbb,ccc,V_bac
  !print*,kkk,jjj,iii,aaa,bbb,ccc,V_cba

  X = W_abc * V_abc + W_bca * V_bca + W_cab * V_cab &
    + W_acb * V_acb + W_bac * V_bac + W_cba * V_cba

  Y = V_abc + V_bca + V_cab
  Z = V_acb + V_bac + V_cba

  !print*,'X Y Z',X,Y,Z

  W1 = W_abc + W_bca + W_cab
  W2 = W_acb + W_bac + W_cba

  !print*,'W1',W1,'W2',W2
  print*,'delta_ijk',delta_ijk
  print*,'Dijkabc',D_ijk - D_abc

  sum_rlt = ((Y - 2.0D0 * Z) * W1 + (Z - 2.0D0 * Y) * W2 + 3.0D0 * X) &
            * delta_ijk / (D_ijk - D_abc)

  print*,'rlt',sum_rlt

  !write(61,*) iii,jjj,kkk,aaa,bbb,ccc,W1
  !write(62,*) iii,jjj,kkk,aaa,bbb,ccc,W2
  !write(63,*) iii,jjj,kkk,aaa,bbb,ccc,X
  !write(64,*) iii,jjj,kkk,aaa,bbb,ccc,Y
  !write(65,*) iii,jjj,kkk,aaa,bbb,ccc,Z
  !write(66,*) iii,jjj,kkk,aaa,bbb,ccc,V_abc,delta_abc

  End Subroutine CC_cl_check_calc_W_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_check_W_abc(iii,jjj,kkk,aaa,bbb,ccc,W_abc)

  Use CC_cl

  Implicit None
  Integer , intent(in) :: iii,jjj,kkk,aaa,bbb,ccc
  Double precision , intent(out) :: W_abc
  Double precision :: rlt
  Double precision , dimension(6) :: rlt_sv

  W_abc = 0.0D0

  Call CC_cl_check_W_PM(iii,jjj,kkk,aaa,bbb,ccc,rlt)
  W_abc = W_abc + rlt
  rlt_sv(1) = rlt

  Call CC_cl_check_W_PM(iii,kkk,jjj,aaa,ccc,bbb,rlt)
  W_abc = W_abc + rlt
  rlt_sv(2) = rlt

  Call CC_cl_check_W_PM(kkk,iii,jjj,ccc,aaa,bbb,rlt)
  W_abc = W_abc + rlt
  rlt_sv(3) = rlt

  Call CC_cl_check_W_PM(kkk,jjj,iii,ccc,bbb,aaa,rlt)
  W_abc = W_abc + rlt
  rlt_sv(4) = rlt

  Call CC_cl_check_W_PM(jjj,kkk,iii,bbb,ccc,aaa,rlt)
  W_abc = W_abc + rlt
  rlt_sv(5) = rlt

  Call CC_cl_check_W_PM(jjj,iii,kkk,bbb,aaa,ccc,rlt)
  W_abc = W_abc + rlt
  rlt_sv(6) = rlt

!  print*,'U'
!  print*,iii,jjj,kkk,aaa,bbb,ccc
!  print*,rlt_sv(1),rlt_sv(2),rlt_sv(3),rlt_sv(4),rlt_sv(5),rlt_sv(6)

  End Subroutine CC_cl_check_W_abc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_check_W_PM(iii,jjj,kkk,aaa,bbb,ccc,rlt)

  Use CC_cl

  Implicit None
  Integer , intent(in) :: iii,jjj,kkk,aaa,bbb,ccc
  Double precision , intent(out) :: rlt
  Double precision :: rlt1, rlt2

  rlt1 = dot_product(CC_intl_bdai_check(bbb,:,aaa,iii),CC_t_check(ccc,:,kkk,jjj))
  rlt2 = dot_product(CC_intl_ckjl_check(ccc,kkk,jjj,:),CC_t_check(aaa,bbb,iii,:))
  rlt = rlt1 - rlt2

  End Subroutine CC_cl_check_W_PM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_check_V_abc(iii,jjj,kkk,aaa,bbb,ccc,V_abc)

  Use CC_cl

  Implicit None
  Integer , intent(in) :: iii,jjj,kkk,aaa,bbb,ccc
  Double precision , intent(out) :: V_abc

  V_abc = CC_intl_bjck_check(bbb,jjj,ccc,kkk) * CC_t_s(iii,aaa) &
        + CC_intl_bjck_check(aaa,iii,ccc,kkk) * CC_t_s(jjj,bbb) &
        + CC_intl_bjck_check(aaa,iii,bbb,jjj) * CC_t_s(kkk,ccc)

  End Subroutine CC_cl_check_V_abc
