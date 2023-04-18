  Subroutine CC_mem_distribution(total_distr,n_tasks,CC_mem_distr,CC_mem_index)

! this subroutine determines the memory distribution for a vector/matrix on each core

  Implicit None

  Integer , Parameter :: ip = Selected_Int_kind(8)

  Integer (kind = ip) :: total_distr, i_tmp, j_tmp
  Integer :: n_tasks, i_task
  Integer (kind = ip) , dimension(n_tasks) :: CC_mem_distr,CC_mem_index

  i_tmp = Int(total_distr/n_tasks)
  j_tmp = Mod(total_distr,n_tasks)

  do i_task = 1, n_tasks
    if (i_task.le.j_tmp) then
      CC_mem_distr(i_task) = i_tmp + 1
    else
      CC_mem_distr(i_task) = i_tmp
    end if
  end do

  CC_mem_index(1) = 1

  do i_task = 2, n_tasks
    CC_mem_index(i_task) = CC_mem_index(i_task - 1) + CC_mem_distr(i_task - 1)
  end do

  End Subroutine CC_mem_distribution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_clear_w()

! this subroutine is used to clean up w vector

  Use dimensions
  Use CC_corr

  Implicit None

  CC_w_1a = 0.0D0
  CC_w_2a = 0.0D0
  CC_w_ab = 0.0D0

  if (n_spin.ne.1) then
    CC_w_1b = 0.0D0
    CC_w_2b = 0.0D0
  end if

  End Subroutine CC_clear_w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_clear_t()

! this subroutine is used to clean up w vector

  Use dimensions
  Use CC_corr

  Implicit None

  CC_t_1a = 0.0D0
  CC_t_2a = 0.0D0
  CC_t_ab = 0.0D0

  if (n_spin.ne.1) then
    CC_t_1b = 0.0D0
    CC_t_2b = 0.0D0
  end if

  End Subroutine CC_clear_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_code_index_s(i_state,i_spin,i,a)
! This subroutine calculates the global index of w_1a(w_1b) for a given excitation
! state (i,a)
! i_state = (i - valence) * nv(i_spin) + (a - ne(i_spin))

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state
  Integer :: i_spin,i,a

  i_state = (i - CC_valence) * CC_n_vir(i_spin) + (a - CC_n_elec(i_spin))

  End Subroutine CC_code_index_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_decode_index_s(i_state,i_spin,i,a)
! This subroutine calculates the single-exciatation state (i,a) of w_1a(w_1b) 
! for a given global index i_state
! i_state = (i - valence) * nv(i_spin) + (a - ne(i_spin))

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state
  Integer :: i_spin,i,a

  a = Mod(i_state,CC_n_vir(i_spin))

  if ((a.eq.0).and.(i_state.ne.0)) then
    a = CC_n_vir(i_spin)
  end if

  i = Int((i_state - a) / CC_n_vir(i_spin)) + CC_valence

  a = a + CC_n_elec(i_spin)

  End Subroutine CC_decode_index_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_inc_index_s(i_spin,iii,aaa)

  Use dimensions
  Use CC_corr

  Implicit None
  Integer :: i_spin,iii,aaa

  aaa = aaa + 1
  if (aaa.gt.n_states) then
    aaa = CC_n_elec(i_spin) + 1
    iii = iii + 1
  end if

  End Subroutine CC_inc_index_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_code_ij(code_ij,iii,jjj,i_spin,j_spin)

! This subroutine code index i and j into code_ij
! if i_spin.eq.j_spin, code_ij = (j - valence - 1) * (j - valcence) / 2 + i - valence + 1
! if i_spin.ne.j_spin, code_ij = (i - valence) * no(j_spin) + j - valence + 1

  Use CC_corr

  Implicit none
  Integer (kind = CC_ip) :: code_ij
  Integer :: iii,jjj,i_spin,j_spin,f,g

  if (i_spin.eq.j_spin) then
    f = jjj - CC_valence
    g = iii - CC_valence + 1
    code_ij = (f - 1) * f / 2 + g
  else
    f = iii - CC_valence
    g = jjj - CC_valence + 1
    code_ij = f * CC_n_occ(j_spin) + g
  end if

  End Subroutine CC_code_ij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_decode_ij(code_ij,iii,jjj,i_spin,j_spin)

! This subroutine decode index code_ij into i and j
! if i_spin.eq.j_spin, code_ij = (j - valence - 1) * (j - valence) / 2 + i - valence + 1
!    code_ij = (j - valence - 1) * (j - valence) / 2 + i - valence + 1
!    assume that f = j - valence, and g = i - valence + 1, we have
!    code_ij = (f - 1) * f / 2 + g
!    note that i < j, we have g <= f, therefore, for 
!    c = 2 * code_ij = f**2 - f + 2 * g
!    we have (f - 1)**2 < c < (f + 1)**2
!    and thus int(sqrt(c)) = f or int(sqrt(c)) = f - 1
!    if g > int(sqrt(c)), f = int(sqrt(c)) + 1
!    else f = int(sqrt(c)) 

! if i_spin.ne.j_spin, code_ij = (i - valence) * no(j_spin) + j - valence + 1
! j = Mod(code_ij,no(j_spin)) + valence - 1
! i = Int(code_ij/no(j_spin)) + valence

  Use CC_corr

  Implicit none
  Integer (kind = CC_ip) :: code_ij
  Integer :: iii,jjj,i_spin,j_spin,f,g

  if (i_spin.eq.j_spin) then
    f = Int(sqrt(2.0D0 * dble(code_ij)))
    g = code_ij - f * (f - 1) / 2
    if (g.gt.f) then
      f = f + 1
      g = code_ij - f * (f - 1) / 2
    else if (g.eq.0) then
      f = f - 1
      g = code_ij - f * (f - 1) / 2
    end if

    iii = g + CC_valence - 1
    jjj = f + CC_valence
  else
    g = Mod(code_ij,CC_n_occ(j_spin))

    if ((g.eq.0).and.(code_ij.ne.0)) then
      g = CC_n_occ(j_spin)
    end if
    f = Int((code_ij - g) / CC_n_occ(j_spin))

    iii = f + CC_valence
    jjj = g + CC_valence - 1
  end if

  End Subroutine CC_decode_ij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_decode_ij4h(i_state,iii,jjj,i_spin)

! (iii,jjj)

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state
  Integer :: iii,jjj,i_spin

  jjj = Mod(i_state,CC_n_occ(i_spin))
  if ((jjj.eq.0).and.(i_state.ne.0)) then
    jjj = CC_n_occ(i_spin)
  end if
  iii = Int((i_state - jjj) / CC_n_occ(i_spin)) + CC_valence
  jjj = jjj + CC_valence - 1

  End Subroutine CC_decode_ij4h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_code_ab(code_ab,aaa,bbb,a_spin,b_spin)

! This subroutine code index a and b into code_ab
! if a_spin.eq.b_spin, code_ab = (b - ne(b_spin) - 2) * (b - ne(b_spin) - 1) / 2 + a - ne(a_spin)
! if a_spin.ne.b_spin, code_ab = (a - ne(a_spin) - 1) * nv(b_spin) + b - ne(b_spin)

  Use CC_corr

  Implicit none
  Integer (kind = CC_ip) :: code_ab
  Integer :: aaa,bbb,a_spin,b_spin,f,g

  if (a_spin.eq.b_spin) then
    f = bbb - CC_n_elec(b_spin) - 1
    g = aaa - CC_n_elec(a_spin)
    code_ab = (f - 1) * f / 2 + g
  else
    f = aaa - CC_n_elec(a_spin) - 1
    g = bbb - CC_n_elec(b_spin)
    code_ab = f * CC_n_vir(b_spin) + g
  end if

  End Subroutine CC_code_ab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_decode_ab(code_ab,aaa,bbb,a_spin,b_spin)

! This subroutine decode index code_ab into a and b
! if a_spin.eq.b_spin, code_ab = (b - ne(b_spin) - 2) * (b - ne(b_spin) - 1) / 2 + a - ne(a_spin)
!    code_ab = (b - ne(b_spin) - 2) * (b - ne(b_spin) - 1) / 2 + a - ne(a_spin)
!    assume that f = b - ne(b_spin) - 1, and g = a - ne(a_spin), we have
!    code_ab = (f - 1) * f / 2 + g
!    note that a < b, we have g <= f, therefore, for 
!    c = 2 * code_ab = f**2 - f + 2 * g
!    we have (f - 1)**2 < c < (f + 1)**2
!    and thus int(sqrt(c)) = f or int(sqrt(c)) = f - 1
!    if g > int(sqrt(c)), f = int(sqrt(c)) + 1
!    else f = int(sqrt(c)) 

! if a_spin.ne.b_spin, code_ab = (a - ne(a_spin) - 1) * nv(b_spin) + b - ne(b_spin)
! b = Mod(code_ab,nv(b_spin)) + ne(b_spin)
! a = Int(code_ab/nv(b_spin)) + ne(a_spin)

  Use CC_corr

  Implicit none
  Integer (kind = CC_ip) :: code_ab
  Integer :: aaa,bbb,a_spin,b_spin,f,g

  if (a_spin.eq.b_spin) then
    f = Int(sqrt(2.0D0 * dble(code_ab)))
    g = code_ab - f * (f - 1) / 2
    if (g.gt.f) then
      f = f + 1
      g = code_ab - f * (f - 1) / 2
    else if (g.eq.0) then
      f = f - 1
      g = code_ab - f * (f - 1) / 2
    end if

    aaa = g + CC_n_elec(a_spin)
    bbb = f + 1 + CC_n_elec(b_spin)
  else
    g = Mod(code_ab,CC_n_vir(b_spin))

    if ((g.eq.0).and.(code_ab.ne.0)) then
      g = CC_n_vir(b_spin)
    end if
    f = Int((code_ab - g) / CC_n_vir(b_spin))

    aaa = f + CC_n_elec(a_spin) + 1
    bbb = g + CC_n_elec(b_spin)
  end if

  End Subroutine CC_decode_ab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_decode_ab4h(i_state,aaa,bbb,i_spin)

! (bbb,aaa)

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state
  Integer :: aaa,bbb,i_spin

  aaa = Mod(i_state,CC_n_vir(i_spin))
  if (aaa.eq.0) then
    aaa = CC_n_vir(i_spin)
  end if
  bbb = Int((i_state - aaa) / CC_n_vir(i_spin)) + CC_n_elec(i_spin) + 1
  aaa = aaa + CC_n_elec(i_spin)

  End Subroutine CC_decode_ab4h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_decode_bj4h(i_state,bbb,jjj,i_spin)

! (bbb,jjj)

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state
  Integer :: bbb,jjj,i_spin

  jjj = Mod(i_state,CC_n_occ(i_spin))
  if ((jjj.eq.0).and.(i_state.ne.0)) then
    jjj = CC_n_occ(i_spin)
  end if
  bbb = Int((i_state - jjj) / CC_n_occ(i_spin)) + 1 + CC_n_elec(i_spin)
  jjj = jjj + CC_valence - 1

  End Subroutine CC_decode_bj4h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_code_index_ab(i_state,i,j,a,b)

! This subroutine determines the global index of saving w_ab vetor for a given 
! excitaition state (i,j,a,b)
! i->a (alpha), j->b (beta)
! i_state = (((i-valence)*no(2)+(j-valence))*nv(1)+(a-ne(1)-1))*nv(2)+(b-ne(2))

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,code_ij,code_ab
  Integer :: i,j,a,b

  Call CC_code_ij(code_ij,i,j,1,2)
  Call CC_code_ab(code_ab,a,b,1,2)

  i_state = (code_ij - 1) * CC_n_vir(1) * CC_n_vir(2) + code_ab

!  i_state = (((i - CC_valence) * CC_n_occ(2) + (j - CC_valence)) * CC_n_vir(1) &
!            + (a - CC_n_elec(1) - 1)) * CC_n_vir(2) + b - CC_n_elec(2)

  End Subroutine CC_code_index_ab


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_decode_index_ab(i_state,i,j,a,b)

! this subroutine determines the 4 index (ijab) for a given excitation state (i_state)
! i->a (alpha), j->b (beta)
! i_state = (((i-valence)*no(2)+(j-valence))*nv(1)+(a-ne(1)-1))*nv(2)+(b-ne(2))

  Use CC_corr

  Implicit None

  Integer (kind = CC_ip) :: i_state,i_tmp,j_tmp,k_tmp,l_tmp,code_ij,code_ab
  Integer :: i,j,a,b

!  i_tmp = mod(i_state,CC_n_vir(2))

!  if ((i_tmp.eq.0).and.(i_state.ne.0)) then
    ! b should not be 0 if i_state is non-zero
!    i_tmp = CC_n_vir(2)
!  end if

!  j_tmp = int((i_state - i_tmp) / CC_n_vir(2))
!  b = i_tmp + CC_n_elec(2)

!  i_tmp = mod(j_tmp,CC_n_vir(1))
!  j_tmp = int((j_tmp - i_tmp) / CC_n_vir(1))
!  a = i_tmp + CC_n_elec(1) + 1

!  i_tmp = mod(j_tmp,CC_n_occ(2))
!  j_tmp = int((j_tmp - i_tmp) / CC_n_occ(2))
!  j = i_tmp + valence
!  i = j_tmp + valence

  i_tmp = CC_n_vir(1) * CC_n_vir(2)

  code_ab = Mod(i_state,i_tmp)
  if ((code_ab.eq.0).and.(i_state.ne.0)) then
    code_ab = i_tmp
  end if

  code_ij = Int((i_state - code_ab) / i_tmp) + 1

  Call CC_decode_ab(code_ab,a,b,1,2)
  Call CC_decode_ij(code_ij,i,j,1,2)

!  print*,'i_state',i_state
!  print*,'code_ij',code_ij
!  print*,'code_ab',code_ab
!  print*,'ijab',i,j,a,b

  End Subroutine CC_decode_index_ab


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_inc_index_ab(iii,jjj,aaa,bbb)

  Use dimensions
  Use CC_corr

  Implicit None
  Integer :: iii,jjj,aaa,bbb

  bbb = bbb + 1
  if (bbb.gt.n_states) then
    bbb = CC_n_elec(2) + 1
    aaa = aaa + 1
  end if

  if (aaa.gt.n_states) then
    aaa = CC_n_elec(1) + 1
    jjj = jjj + 1
  end if

  if (jjj.gt.CC_n_elec(2)) then
    jjj = CC_valence
    iii = iii + 1
  end if

  End Subroutine CC_inc_index_ab


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_inc_ij4h(iii,jjj,i_spin)

  Use dimensions
  Use CC_corr

  Implicit None
  Integer :: iii,jjj,i_spin

  jjj = jjj + 1

  if (jjj.gt.CC_n_elec(i_spin)) then
    jjj = CC_valence
    iii = iii + 1
  end if

  End Subroutine CC_inc_ij4h


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_inc_ab4h(aaa,bbb,i_spin)

! (bbb,aaa)

  Use dimensions
  Use CC_corr

  Implicit None
  Integer :: aaa,bbb,i_spin

  aaa = aaa + 1

  if (aaa.gt.n_states) then
    aaa = CC_n_elec(i_spin) + 1
    bbb = bbb + 1
  end if

  End Subroutine CC_inc_ab4h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_inc_bj4h(bbb,jjj,i_spin)

  Use dimensions
  Use CC_corr

  Implicit None
  Integer :: bbb,jjj,i_spin

  jjj = jjj + 1

  if (jjj.gt.CC_n_elec(i_spin)) then
    jjj = CC_valence
    bbb = bbb + 1
  end if

  End Subroutine CC_inc_bj4h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_inc_icak4h(iii,ccc,aaa,kkk,i_spin,c_spin,a_spin,k_spin)

! (kkk,ccc,iii,aaa)

  Use dimensions
  Use CC_corr

  Implicit None
  Integer :: iii,ccc,aaa,kkk,i_spin,c_spin,a_spin,k_spin

  aaa = aaa + 1
  if (aaa.gt.n_states) then
    aaa = CC_n_elec(a_spin) + 1
    iii = iii + 1
  end if

  if (iii.gt.CC_n_elec(i_spin)) then
    iii = CC_valence
    ccc = ccc + 1
  end if

  if (ccc.gt.n_states) then
    ccc = CC_n_elec(c_spin) + 1
    kkk = kkk + 1
  end if

  End Subroutine CC_inc_icak4h


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_code_index_xx(i_state,i_spin,i,j,a,b)

! This subroutine calculate the global index (i_state) of for variables i<j,a<b
! code_ij = (j - valence - 1) * (j - valence) / 2 + i - valence + 1
! code_ab = (b - 2 - ne) * (b - 1 - ne) / 2 + a - ne
! i_state = (code_ij - 1) * (nv - 1) * nv / 2 + code_ab

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,code_ij,code_ab
  Integer :: i_spin,i,j,a,b
  Integer :: f,g,p,q,r,s


! Calculate code_ij
  Call CC_code_ij(code_ij,i,j,i_spin,i_spin)

!  f = j - CC_valence + 1
!  g = i - CC_valence + 1
!  code_ij = (f - 2) * (f - 1) / 2 + g

! Calculate code_ab
!  f = b - CC_n_elec(i_spin)
!  g = a - CC_n_elec(i_spin)
!  code_ab = (f - 2) * (f - 1) / 2 + g
  Call CC_code_ab(code_ab,a,b,i_spin,i_spin)

! Calculate i_state
  i_state = (code_ij - 1) * (CC_n_vir(i_spin) - 1) * CC_n_vir(i_spin) / 2 + code_ab

!  print*,'CC_code_index_xx'
!  print*,'i,j,a,b',i,j,a,b
!  print*,i_state

  End Subroutine CC_code_index_xx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_decode_index_xx(i_state,i_spin,i,j,a,b)

! This subroutine calculates four-index (i<j,a<b) for a given global index (i_state)
! First: decode i_state into code_ij and code_ab
!        i_state = (code_ij - 1) * (nv - 1) * nv / 2 + code_ab
! Second: decode code_ij and code_ab
!        code_ij = (j - valence - 1) * (j - valence) / 2 + i - valence + 1
!        assume that f = j - valence, and g = i - valence + 1, we have
!        code_ij = (f - 1) * f / 2 + i
!        note that i < j, we have g <= f, therefore, for 
!        c = 2 * code_ij = f**2 - f + 2 * g
!        we have (f - 1)**2 < c < (f + 1)**2
!        and thus int(sqrt(c)) = f or int(sqrt(c)) = f - 1
!        if g > int(sqrt(c)), f = int(sqrt(c)) + 1
!        else f = int(sqrt(c)) 
! the decoding procedure for code_ab is much similar

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,code_ij,code_ab
  Integer :: i_spin,i,j,a,b
  Integer :: f,g,p,q,r,s

! Calculate code_ij and code_ab
  s = (CC_n_vir(i_spin) - 1) * CC_n_vir(i_spin) / 2

  code_ab = mod(i_state,s)

  if ((code_ab.eq.0).and.(i_state.ne.0)) then
    code_ab = s
  end if

  code_ij = Int((i_state - code_ab) / s) + 1

! Decode code_ij
!  f = Int(sqrt(2.0D0 * dble(code_ij)))
!  g = code_ij - f * (f - 1) / 2
!  if (g.gt.f) then
!    f = f + 1
!    g = code_ij - f * (f - 1) / 2
!  else if (g.eq.0) then
!    f = f - 1
!    g = code_ij - f * (f - 1) / 2
!  end if

!  i = g + CC_valence - 1
!  j = f + CC_valence

  Call CC_decode_ij(code_ij,i,j,i_spin,i_spin)

! Decode code_ab
!  f = Int(sqrt(2.0D0 * dble(code_ab)))
!  g = code_ab - f * (f - 1) / 2
!  if (g.gt.f) then
!    f = f + 1
!    g = code_ab - f * (f - 1) / 2
!  else if (g.eq.0) then
!    f = f - 1
!    g = code_ab - f * (f - 1) / 2
!  end if

!  a = g + CC_n_elec(i_spin)
!  b = f + 1 + CC_n_elec(i_spin)

  Call CC_decode_ab(code_ab,a,b,i_spin,i_spin)

  End Subroutine CC_decode_index_xx



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  Subroutine CC_aux_decode_a_xx(i_state,i_spin,k,l,i,j)

! This subroutine calculates four-index (k<l,i<j) for a given index of auxiliary vector a_xx
! First: decode i_state into code_ij and code_kl
!        i_state = (code_kl - 1) * (no - 1) * no / 2 + code_ij
! Second: decode code_ij and code_kl
!        code_ij = (j - valence - 1) * (j - valence) / 2 + i - valence + 1
!        assume that f = j - valence, g = i - valence + 1, we have
!        code_ij = (f - 1) * f / 2 + g
!        note that i < j, we have g <= f, therefore, for 
!        c = 2 * code_ij = f**2 - f + 2 * g
!        we have (f - 1)**2 < c < (f + 1)**2
!        and thus int(sqrt(c)) = f or int(sqrt(c)) = f - 1
!        if g > int(sqrt(c)), f = int(sqrt(c)) + 1
!        else f = int(sqrt(c)) 
! the decoding procedure for code_kl is much similar

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,code_ij,code_kl,i_tmp
  Integer :: k,l,i,j,i_spin,f,g

! Calculate code_kl and code_ij
  i_tmp = CC_n_occ(i_spin) * (CC_n_occ(i_spin) - 1) / 2
  code_ij = Mod(i_state,i_tmp)
  if ((code_ij.eq.0).and.(i_state.ne.0)) then
    code_ij = i_tmp
  end if

  code_kl = Int((i_state - code_ij) / i_tmp) + 1

! Decode code_kl
!  f = Int(sqrt(2.0D0 * dble(code_kl)))
!  g = code_kl - f * (f - 1) / 2
!  if (g.gt.f) then
!    f = f + 1
!    g = code_kl - f * (f - 1) / 2
!  else if (g.eq.0) then
!    f = f - 1
!    g = code_kl - f * (f - 1) / 2
!  end if

!  l = f + CC_valence
!  k = g + CC_valence - 1

  Call CC_decode_ij(code_kl,k,l,i_spin,i_spin)

! Decode code_ij
!  f = Int(sqrt(2.0D0 * dble(code_ij)))
!  g = code_ij - f * (f - 1) / 2
!  if (g.gt.f) then
!    f = f + 1
!    g = code_ij - f * (f - 1) / 2
!  else if (g.eq.0) then
!    f = f - 1
!    g = code_ij - f * (f - 1) / 2
!  end if

!  j = f + CC_valence
!  i = g + CC_valence - 1

  Call CC_decode_ij(code_ij,i,j,i_spin,i_spin)

  End Subroutine CC_aux_decode_a_xx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_aux_decode_b_xx(i_state,i_spin,c,d,a,b)

! This subroutine calculates four-index (c<d,a<b) for a given index of auxiliary vector b_xx
! First: decode i_state into code_ab and code_cd
!        i_state = (code_cd - 1) * (nv - 1) * nv / 2 + code_ab
! Second: decode code_ab and code_cd
!        code_ab = (b - ne - 2) * (b - ne - 1) / 2 +  a - ne
!        assume that f = b - ne - 1, g = a - ne, we have
!        code_ab = (f - 1) * f / 2 + g
!        note that a < b, we have g <= f, therefore, for 
!        c = 2 * code_ab = f**2 - f + 2 * g
!        we have (f - 1)**2 < c < (f + 1)**2
!        and thus int(sqrt(c)) = f or int(sqrt(c)) = f - 1
!        if g > int(sqrt(c)), f = int(sqrt(c)) + 1
!        else f = int(sqrt(c)) 
! the decoding procedure for code_cd is much similar

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,code_ab,code_cd,i_tmp
  Integer :: a,b,c,d,i_spin,f,g

! Calculate code_cd and code_ab
  i_tmp = CC_n_vir(i_spin) * (CC_n_vir(i_spin) - 1) / 2
  code_ab = Mod(i_state,i_tmp)
  if ((code_ab.eq.0).and.(i_state.ne.0)) then
    code_ab = i_tmp
  end if

  code_cd = Int((i_state - code_ab) / i_tmp) + 1

! Decode code_cd
!  f = Int(sqrt(2.0D0 * dble(code_cd)))
!  g = code_cd - f * (f - 1) / 2
!  if (g.gt.f) then
!    f = f + 1
!    g = code_cd - f * (f - 1) / 2
!  else if (g.eq.0) then
!    f = f - 1
!    g = code_cd - f * (f - 1) / 2
!  end if

!  d = f + CC_n_elec(i_spin) + 1
!  c = g + CC_n_elec(i_spin)

  Call CC_decode_ab(code_cd,c,d,i_spin,i_spin)

! Decode code_ab
!  f = Int(sqrt(2.0D0 * dble(code_ab)))
!  g = code_ab - f * (f - 1) / 2
!  if (g.gt.f) then
!    f = f + 1
!    g = code_ab - f * (f - 1) / 2
!  else if (g.eq.0) then
!    f = f - 1
!    g = code_ab - f * (f - 1) / 2
!  end if

!  b = f + CC_n_elec(i_spin) + 1
!  a = g + CC_n_elec(i_spin)

  Call CC_decode_ab(code_ab,a,b,i_spin,i_spin)

  End Subroutine CC_aux_decode_b_xx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_aux_decode_a_ab(i_state,k,l,i,j)

! This subroutine calculates four-index (k,l,i,j) for a given state of a vector
! in alpha-beta space
! i_state = (((k - valence) * no(2) + l - valence) * no(1) + i - valence) * no(2) 
!           + j - valence + 1

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp,code_kl,code_ij
  Integer :: i,j,k,l

!  j = Mod(i_state,CC_n_occ(2))
!  if ((j.eq.0).and.(i_state.ne.0)) then
!    j = CC_n_occ(2)
!  end if

!  i_tmp = Int((i_state - j) / CC_n_occ(2))

!  j = j + CC_valence - 1

!  i = Mod(i_tmp,CC_n_occ(1)) + CC_valence

!  i_tmp = Int(i_tmp/CC_n_occ(1))

!  l = Mod(i_tmp,CC_n_occ(2)) + CC_valence

!  k = Int(i_tmp/CC_n_occ(2)) + CC_valence

  i_tmp = CC_n_occ(1) * CC_n_occ(2)

  code_ij = Mod(i_state,i_tmp)

  if ((code_ij.eq.0).and.(i_state.ne.0)) then
    code_ij = i_tmp
  end if

  code_kl = Int((i_state - code_ij) / i_tmp) + 1

  Call CC_decode_ij(code_kl,k,l,1,2)

  Call CC_decode_ij(code_ij,i,j,1,2)

  End Subroutine CC_aux_decode_a_ab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_aux_decode_b_ab(i_state,c,d,a,b)

! This subroutine calculates four-index (c,d,a,b) for a given state of b vector
! in alpha-beta space
! i_state = (((c - ne(1) - 1) * nv(2) + d - ne(2) - 1) * nv(1) + a - ne(1) - 1) * nv(2) 
!           + b - ne(2)

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp,code_ab,code_cd
  Integer :: a,b,c,d

!  b = Mod(i_state,CC_n_vir(2))
!  if ((b.eq.0).and.(i_state.ne.0)) then
!    b = CC_n_vir(2)
!  end if

!  i_tmp = Int((i_state - b) / CC_n_vir(2))

!  b = b + CC_n_elec(2)

!  a = Mod(i_tmp,CC_n_vir(1)) + CC_n_elec(1)

!  i_tmp = Int(i_tmp/CC_n_vir(1))

!  d = Mod(i_tmp,CC_n_vir(2)) + CC_n_elec(2)

!  c = Int(i_tmp/CC_n_vir(2)) + CC_n_elec(1)

  i_tmp = CC_n_vir(1) * CC_n_vir(2)

  code_ab = Mod(i_state,i_tmp)
  if ((code_ab.eq.0).and.(i_state.ne.0)) then
    code_ab = i_tmp
  end if

  code_cd = Int((i_state - code_ab) / i_tmp) + 1

  Call CC_decode_ab(code_cd,c,d,1,2)

  Call CC_decode_ab(code_ab,a,b,1,2)

  End Subroutine CC_aux_decode_b_ab

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_aux_decode_h(i_state,i,c,a,k,i_spin,c_spin,a_spin,k_spin)

! (kkk,ccc,iii,aaa)
! i_state = (((k - valence) * nv(c_spin) + c - ne(c_spin) - 1) * no(i_spin) + i - valence) * nv(a_spin) + a - ne(a_spin)

  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state
  Integer :: i,c,a,k,i_spin,c_spin,a_spin,k_spin,i_tmp

  a = Mod(i_state,CC_n_vir(a_spin))
  if ((a.eq.0).and.(i_state.ne.0)) then
    a = CC_n_vir(a_spin)
  end if

  i_tmp = Int((i_state - a) / CC_n_vir(a_spin))

  i = Mod(i_tmp,CC_n_occ(i_spin)) + CC_valence
  i_tmp = Int(i_tmp / CC_n_occ(i_spin))
  c = Mod(i_tmp,CC_n_vir(c_spin)) + 1
  k = Int(i_tmp/CC_n_vir(c_spin)) + CC_valence

  a = a + CC_n_elec(a_spin)
  c = c + CC_n_elec(c_spin)

  End Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_aux_h_pointer(i_aux,iii,ccc,aaa,kkk,i_spin,c_spin,a_spin,k_spin,i_icak)
! i_state = (((k - valence) * nv(c_spin) + c - ne(c_spin) - 1) * no(i_spin) + i - valence) * nv(a_spin) + a - ne(a_spin)

  Use dimensions
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_aux,i_icak,code_ia,code_tp
  Integer :: iii,ccc,aaa,kkk,i_tmp,c_tmp,ainx
  Integer :: i_spin,c_spin,a_spin,k_spin,a_tmp

  Call CC_aux_decode_h(i_aux,i_tmp,ccc,a_tmp,kkk,i_spin,c_spin,a_spin,k_spin)

  ainx = aaa - CC_n_elec(a_spin)

  code_ia = (iii - CC_valence) * CC_n_vir(a_spin) + ainx
  code_tp = (i_tmp - CC_valence) * CC_n_vir(a_spin) + a_tmp - CC_n_elec(a_spin)

  if (code_tp.gt.code_ia) then
    ccc = ccc + 1
    if (ccc.gt.n_states) then
      ccc = CC_n_elec(c_spin) + 1
      kkk = kkk + 1
    end if
  end if

  i_icak = (((kkk - CC_valence) * CC_n_vir(c_spin) + ccc - CC_n_elec(c_spin) - 1) &
            * CC_n_occ(i_spin) + iii - CC_valence) * CC_n_vir(a_spin) + ainx

  End Subroutine CC_aux_h_pointer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_eigenvalue_diff(rrr,sss,r_spin,s_spin,E_diff)

! This subroutine calculates the difference between rrr and sss eigenvalue
! E_diff = KS_eigenvalue(sss) - KS_eigenvalue(rrr)

  Use dimensions
  Use physics
  Use CC_corr

  Implicit None
  Integer :: rrr,sss,r_spin,s_spin
  Double precision :: E_diff

  if (n_spin.eq.1) then
    E_diff = KS_eigenvalue(sss,1,1) - KS_eigenvalue(rrr,1,1)
  else
    E_diff = KS_eigenvalue(sss,s_spin,1) - KS_eigenvalue(rrr,r_spin,1)
  end if

  End Subroutine CC_Calc_eigenvalue_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_integral(rlt,r,s,t,u,r_spin,s_spin,t_spin,u_spin)

! This subroutine calculates the integral <rs||tu>

  Use dimensions
  Use CC_corr

  Implicit None

  Double precision :: rlt,ddot
  Integer :: r,s,t,u,r_spin,s_spin,t_spin,u_spin

  if ((r.eq.s).and.(r_spin.eq.s_spin)) then
    print*,'Calc_int err',r,s,r_spin,s_spin
  end if

  if ((t.eq.u).and.(t_spin.eq.u_spin)) then
    print*,'Calc_int err',t,u,t_spin,u_spin
  end if

  if (r_spin.eq.s_spin) then
    ! In this case, r_spin = s_spin = t_spin = u_spin, in alpha-alpha space or beta-beta space,
    ! both Coulomb and Exchange integrals are calculated.

    if (n_spin.eq.1) then
      ! For close shell molecules, alpha orbitals are identical to beta ones, 
      ! and only alpha orbitals are saved
      rlt = ddot(n_basbas,ovlp_3ks(:,r,t,1),1, &
                          ovlp_3ks(:,s,u,1),1) &
          - ddot(n_basbas,ovlp_3ks(:,r,u,1),1, &
                          ovlp_3ks(:,s,t,1),1)  
    else
      ! For open shell molecules, alpha and beta orbitals are different.
      rlt = ddot(n_basbas,ovlp_3ks(:,r,t,r_spin),1, &
                          ovlp_3ks(:,s,u,r_spin),1) &
          - ddot(n_basbas,ovlp_3ks(:,r,u,r_spin),1, &
                          ovlp_3ks(:,s,t,r_spin),1) 
    end if

  else
! In this case, r_spin /= s_spin, only Coulomb integral is needed.
    if (r_spin.eq.t_spin) then
      ! r_spin = t_spin, s_spin = u_spin
      if (n_spin.eq.1) then
        ! For close shell molecules
        rlt = ddot(n_basbas,ovlp_3ks(:,r,t,1),1, &
                            ovlp_3ks(:,s,u,1),1) 
      else
        ! For open shell calculations
        rlt = ddot(n_basbas,ovlp_3ks(:,r,t,r_spin),1, &
                            ovlp_3ks(:,s,u,s_spin),1) 
      end if
    else
      ! r_spin = u_spin, s_spin = t_spin
      if (n_spin.eq.1) then
        ! For close shell molecules
        rlt = - ddot(n_basbas,ovlp_3ks(:,r,u,1),1, &
                              ovlp_3ks(:,s,t,1),1) 
      else
        ! For open shell calculations
        rlt = - ddot(n_basbas,ovlp_3ks(:,r,u,r_spin),1, &
                              ovlp_3ks(:,s,t,s_spin),1) 
      end if
    end if
  end if

  End Subroutine CC_Calc_integral


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_get_coeff_single(i_spin,iii,aaa,coeff)

  Use dimensions
  Use CC_corr

  Implicit None

  Integer :: i_spin,iii,ainx,aaa,iinx
  Double precision :: coeff

  iinx = iii - CC_valence + 1

  if ((i_spin.eq.1).or.(n_spin.eq.1)) then
    ainx = aaa - CC_n_elec(1)
    coeff = CC_t_1a(iinx,ainx)
  else
    ainx = aaa - CC_n_elec(2)
    coeff = CC_t_1b(iinx,ainx)
  end if

  End Subroutine CC_get_coeff_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_get_coeff_double(iii,jjj,aaa,bbb,i_spin,j_spin,a_spin,b_spin,coeff)

  Use dimensions
  Use CC_corr

  Implicit None

  Integer (kind = CC_ip) :: i_state
  Integer :: iii,jjj,aaa,bbb,i_spin,j_spin,a_spin,b_spin
  Integer :: i_i,i_j,i_a,i_b
  Double precision :: coeff

  if (i_spin+j_spin.ne.a_spin+b_spin) then
    print*,'get_coeff_double err',iii,jjj,aaa,bbb,i_spin,j_spin,a_spin,b_spin
  end if

  coeff = 1.0D0
  if (i_spin.eq.j_spin) then
    if (iii.lt.jjj) then
      i_i = iii
      i_j = jjj
    else
      i_i = jjj
      i_j = iii
      coeff = - coeff
    end if

    if (aaa.lt.bbb) then
      i_a = aaa
      i_b = bbb
    else
      i_a = bbb
      i_b = aaa
      coeff = - coeff
    end if

    if ((i_spin.eq.1).or.(n_spin.eq.1)) then
      Call CC_code_index_xx(i_state,1,i_i,i_j,i_a,i_b)
      coeff = coeff * CC_t_2a(i_state)
    else
      Call CC_code_index_xx(i_state,2,i_i,i_j,i_a,i_b)
      coeff = coeff * CC_t_2b(i_state)
    end if

  else

    if (i_spin.lt.j_spin) then
      i_i = iii
      i_j = jjj
    else
      i_i = jjj
      i_j = iii
      coeff = - coeff
    end if

    if (a_spin.lt.b_spin) then
      i_a = aaa
      i_b = bbb
    else
      i_a = bbb
      i_b = aaa
      coeff = - coeff
    end if

    Call CC_code_index_ab(i_state,i_i,i_j,i_a,i_b)
    if ((i_state.gt.CC_n_ab).or.(i_state.lt.1)) then
      print*,'decode_ab err',iii,jjj,aaa,bbb,i_spin,j_spin,a_spin,b_spin
    end if
    coeff = coeff * CC_t_ab(i_state)

  end if

  End Subroutine CC_get_coeff_double

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Calc_Norm(CC_Norm)

! This subroutine calculates the normalization factor of CC wave function.
! Note in CC energy calculations, the wave function is media-normalized.

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer (kind = CC_ip) :: i_state,i_tmp
  Integer :: i_spin,iii,aaa,jjj,bbb
  Double precision :: CC_Norm,coeff

  CC_Norm = 0.0D0

  do i_spin = 1, n_spin
    i_state = CC_index_w(myid+1,i_spin) - 1
    Call CC_decode_index_s(i_state,i_spin,iii,aaa)
    do i_tmp = 1, CC_mem_w(myid+1,i_spin)
      Call CC_inc_index_s(i_spin,iii,aaa)
      Call CC_get_coeff_single(i_spin,iii,aaa,coeff)
      CC_Norm = CC_Norm + coeff ** 2
    end do

    i_state = CC_index_w(myid+1,3+i_spin) - 1
    do i_tmp = 1, CC_mem_w(myid+1,3+i_spin)
      i_state = i_state + 1
      Call CC_decode_index_xx(i_state,i_spin,iii,jjj,aaa,bbb)
      Call CC_get_coeff_double(iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin,coeff)
      CC_Norm = CC_Norm + coeff ** 2
    end do
  end do

  if (n_spin.eq.1) then
    CC_Norm = CC_Norm * 2.0D0
  end if

  i_state = CC_index_w(myid+1,3) - 1
  Call CC_decode_index_ab(i_state,iii,jjj,aaa,bbb)
  do i_tmp = 1, CC_mem_w(myid+1,3)
    Call CC_inc_index_ab(iii,jjj,aaa,bbb)
    Call CC_get_coeff_double(iii,jjj,aaa,bbb,1,2,1,2,coeff)
    CC_Norm = CC_Norm + coeff ** 2
  end do

  Call sync_real_number(CC_norm)

  CC_norm = sqrt(1.0D0 + CC_norm)

  End Subroutine CC_Calc_Norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_RLE_CalcNorm_w(t_norm)

! This subroutine calculates the normalization factor of the w vector
! (in fact, it is the new trial function)

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None
  Integer :: i_spin,i_local
  Double precision :: t_norm

  ! Calculate the normalized factor of w
  t_norm = 0.0D0

  do i_spin = 1, n_spin
    do i_local = 1, CC_mem_w(myid+1,i_spin)
      if (i_spin.eq.1) then
        t_norm = t_norm + CC_w_1a(i_local) ** 2
      else
        t_norm = t_norm + CC_w_1b(i_local) ** 2
      end if
    end do

    do i_local = 1, CC_mem_w(myid+1,3+i_spin)
      if (i_spin.eq.1) then
        t_norm = t_norm + CC_w_2a(i_local) ** 2
      else
        t_norm = t_norm + CC_w_2b(i_local) ** 2
      end if
    end do
  end do

  if (n_spin.eq.1) then
    t_norm = t_norm * 2.0D0
  end if

  do i_local = 1, CC_mem_w(myid+1,3)
    t_norm = t_norm + CC_w_ab(i_local) ** 2
  end do

  Call sync_real_number(t_norm)

  t_norm = sqrt(t_norm)

  End Subroutine CC_RLE_CalcNorm_w


