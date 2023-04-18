Module CC_3d

!  PURPOSE
!
!  This file may at some point contain all the needed variable declarations
!  for the Coupled-Cluster (CC) calculations (hybrid algorithm)
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

  Use CC_3d_distribution
  !Use CC_3d_hdf5
  Use CC_3d_mem
  Use CC_3d_MPI
  Use dimensions
  Use mpi_tasks
  Use runtime_choices
  Use physics
  Use hartree_fock 

  Implicit None

! Parameters discribe CC equations

  Integer (kind = 8) :: CC_n_config
! Number of configurations in CC calculations

  Integer :: CC_n_s
! Number of singlet excitations in CC calculations

  Integer (kind = 8) :: CC_n_d
! Number of double-excitations in CC calculations

  Integer :: CC_n_elec
! Number of electrons (spin free)

  Integer :: CC_n_occ
! Number of occupied(valence) orbitals (spin free)

  Integer :: CC_n_vir
! Number of virtual orbitals (spin free)

  Integer :: CC_n_fc
! Number of frozen core electrons

  Integer :: CC_valence
  ! The starting index of valence state

  Integer :: CC_n_state
  ! Number of total orbitals included (=n_states if full claculation is taken
  !                                    =n_states - n_frozen if frozen core approximation is used)

  Integer :: CC_n_bas
! Numbers of auxiliary basis functions on each core

  Integer :: CC_sv_strategy
! Saving strategy
! 1 ---- Maximum memory space demand, all vectors needed are saved in memory (replicated) (default)
! 2 ---- Medium memory space demand, no replication of (ka|cd) is saved 
! 3 ---- Minimum memory space demand, only necessary vectors are saved

  Integer :: CC_i_scf
! Number of CC scf iteration

  Double precision :: E_HF, CC_E_corr
! E_HF : HF energy
! CC_E_corr : Correlation energy

! DIIS parameters
  Integer :: CC_DIIS_ndc, CC_DIIS_m_bgn
! CC_DIIS_ndc ---- the dimension of current R matrix in DIIS algorithm
! CC_DIIS_m_bgn ---- the start point of saved basis functions

  Double precision :: CC_Norm
! Normalization coefficient 

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_code(code_fg,f,g,i_step,codetype)

! This subroutine maps an index (f,g) to code_fg
! codetype : 1 ---- code_fg = (f - 1) * i_step + g
!            2 ---- code_fg = (g - 1) * g / 2 + f
!            3 ---- code_fg = (g - 2) * (g - 1) / 2 + f

  Implicit None
  Integer , intent(out) :: code_fg
  Integer , intent(in) :: i_step
  Integer , intent(in) :: f,g
  Integer , intent(in) :: codetype

  if (codetype.eq.1) then
    code_fg = (f - 1) * i_step + g
  else if (codetype.eq.2) then
    code_fg = (g - 1) * g / 2 + f
  else if (codetype.eq.3) then
    code_fg = (g - 2) * (g - 1) / 2 + f
  end if

  End Subroutine CC_3d_code

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_decode(code_fg,f,g,i_step,detype)

! This subroutine maps an index code_fg to (f,g)
! detype : 1 ---- code_fg = (f - 1) * i_step + g
!          2 ---- code_fg = (g - 1) * g / 2 + f
!          3 ---- code_fg = (g - 2) * (g - 1) / 2 + f

! Therefore we have,
! detype : 1 ---- g = Mod(code_fg,i_step) , f = Int(code_fg/i_step) + 1
!          2 ---- Note that f <= g, for c = 2 * code_fg = g**2 - g + 2 * f
!                 we have (g - 1)**2 < c < (g + 1)**2
!                 and thus int(sqrt(c)) = g or int(sqrt(c)) = g - 1
!                 if f > int(sqrt(c)), g = int(sqrt(c)) + 1
!                 else g = int(sqrt(c))
!          3 ---- Similar to 2, the only difference is that f < g, assume that h = g - 1
!                 we get the same expansion by substitute (g - 1) by h.

  Implicit None

  Integer , intent(in) :: code_fg
  Integer , intent(in) :: i_step
  Integer , intent(out) :: f,g
  Integer , intent(in) :: detype

  if (detype.eq.1) then

    g = Mod(code_fg,i_step)
    if ((g.eq.0).and.(code_fg.ne.0)) then
      g = i_step
    end if

    f = Int((code_fg - g) / i_step) + 1

  else if (detype.eq.2) then

    g = Int(sqrt(2.0D0 * dble(code_fg)))
    f = code_fg - g * (g - 1) / 2
    if (f.gt.g) then
      g = g + 1
      f = code_fg - g * (g - 1) / 2
    else if (f.eq.0) then
      g = g - 1
      f = code_fg - g * (g - 1) / 2
    end if

  else if (detype.eq.3) then

    g = Int(sqrt(2.0D0 * dble(code_fg)))
    f = code_fg - g * (g - 1) / 2
    if (f.gt.g) then
      g = g + 1
      f = code_fg - g * (g - 1) / 2
    else if (f.eq.0) then
      g = g - 1
      f = code_fg - g * (g - 1) / 2
    end if

    g = g + 1

  end if

  End Subroutine CC_3d_decode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,k_miss)

! This subroutine gives the k_pattern determined by four k_points: k1,k2,k3 & k4.
! Since momentum conservation law holds, only three of them are free.
! momentum conservation law: k1 + k3 = k2 + k4 (within first Brillioun zone)

  Implicit None
  Integer :: k1,k2,k3,k4
  Integer , optional, intent(in) :: k_miss
  Integer , intent(out) :: k_pattern
  Integer :: k_13,k_24,k_12

  if (present(k_miss)) then

    if (k_miss.eq.1) then
  
      k_24 = kpq_point_list(k2,k4)
      k1 = kq_point_list(k_24,k3)
  
    else if (k_miss.eq.2) then
  
      k_13 = kpq_point_list(k1,k3)
      k2 = kq_point_list(k_13,k4)
  
    else if (k_miss.eq.3) then
  
      k_24 = kpq_point_list(k2,k4)
      k3 = kq_point_list(k_24,k1)
  
    else if (k_miss.eq.4) then
  
      k_13 = kpq_point_list(k1,k3)
      k4 = kq_point_list(k_13,k2)
  
    end if
  end if

  Call CC_check_k(k1,k2,k3,k4)

  Call CC_3d_code(k_12,k1,k2,n_k_points,1)
  Call CC_3d_code(k_pattern,k_12,k3,n_k_points,1)

  End Subroutine CC_3d_determine_k_pattern

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_k_plus(k_13,k1,k3)

  Implicit None
  Integer , intent(in) :: k1,k3
  Integer , optional, intent(out) :: k_13

  k_13 = kpq_point_list(k1,k3)

  End Subroutine CC_3d_k_plus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_k_minus(k_13,k1,k3)

  Implicit None
  Integer , intent(in) :: k_13,k1
  Integer , optional, intent(out) :: k3

  k3 = kq_point_list(k_13,k1)
 
  End Subroutine CC_3d_k_minus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_check_k(k1,k2,k3,k4)

! This subroutine gives the k_points determined by a specific k_pattern.
! momentum conservation law: k1 + k3 = k2 + k4 (within first Brillioun zone)

  Implicit None
  Integer , intent(in) :: k1,k2,k3,k4
  Integer :: k_13,k_24

  k_13 = kpq_point_list(k1,k3)
  k_24 = kpq_point_list(k2,k4)

  if (k_13.ne.k_24) then
    print*,'error in k_pattern',k1,k2,k3,k4
  end if

  End Subroutine CC_check_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_get_t_grp(job,i_grp_task,t_out,local_flag)

  Implicit None
  Character (len=1) , intent(in) :: job
  Integer , intent(in) :: i_grp_task
  Logical , optional :: local_flag
  Double complex , dimension(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                             CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                             CC_mem_k1(i_grp_task)) :: t_out
  Integer (kind=8) :: s_tmp
  Integer :: k_start,k_end,n_k1,source_id,i_tmp,j_tmp
  Integer :: iii,jjj,aaa,bbb,k1,k2,k3,k4

  ! CC_t_d (j,b,kj,i,a,ki,ka)

  if (.not.(present(local_flag))) then
    local_flag = .true.
  end if

  n_k1 = CC_mem_k1(i_grp_task)
  k_start = CC_index_k1(i_grp_task)
  k_end = k_start - 1 + n_k1

  if (CC_sv_strategy.eq.1) then

    if (job.eq.'S') then
      t_out = CC_t_d_A(:,:,:,:,:,:,k_start:k_end)
    else if (job.eq.'D') then

      do k2 = 1, n_k_points
        do k4 = 1, n_k_points
          do iii = 1, CC_n_occ
            do jjj = 1, CC_n_occ
              t_out(iii,:,k2,jjj,:,k4,:) = CC_t_d_A(jjj,:,k4,iii,:,k2,k_start:k_end)
            end do
          end do
        end do
      end do

    end if

  else if (CC_sv_strategy.eq.2) then

    source_id = i_grp_task - 1

    if (source_id.eq.CC_mpi_gid) then
      if (job.eq.'S') then
        t_out = CC_t_d
      else if (job.eq.'D') then

        do k2 = 1, n_k_points
          do k4 = 1, n_k_points
            do iii = 1, CC_n_occ
              do jjj = 1, CC_n_occ
                t_out(iii,:,k2,jjj,:,k4,:) = CC_t_d(jjj,:,k4,iii,:,k2,:)
              end do
            end do
          end do
        end do

      end if
    else
      t_out = 0.0D0
    end if

    if (.not.(local_flag)) then
 
      i_tmp = CC_n_occ * CC_n_vir * n_k_points
      j_tmp = CC_n_occ * CC_mem_aa_D(CC_mpi_did+1) * n_k_points * n_k1
      s_tmp = Int(i_tmp,8) * Int(j_tmp,8)

      Call CC_mpi_complex_bcast(s_tmp, t_out, source_id, CC_mpi_comm_group)
    end if

  end if

  End Subroutine CC_3d_get_t_grp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_get_tau(k1,k2,k3,k4,t_out)

  Implicit None
  Double complex , dimension(CC_mem_aa_D(CC_mpi_did+1),CC_n_vir, &
                             CC_n_occ,CC_n_occ) :: t_out
  Integer :: k_start,k_end,n_k2
  Integer :: iii,jjj,aaa,bbb,a_run,n_a
  Integer :: k1,k2,k3,k4,k_13,k_run

  n_a = CC_mem_aa_D(CC_mpi_did+1)

  ! CC_t_d (j,b,kj,i,a,ki,ka)
  do aaa = 1, CC_mem_aa_D(CC_mpi_did+1)
    do jjj = 1, CC_n_occ
      t_out(aaa,:,:,jjj) = CC_t_d_A(jjj,:,k4,:,aaa,k2,k1)
    end do
  end do

  if (k1.eq.k2) then
    !$OMP PARALLEL Default(Shared) &
    !$OMP Private(aaa,a_run,bbb,iii,jjj)
    !$OMP DO COLLAPSE(4)
    do a_run = 1, n_a
      do bbb = 1, CC_n_vir
        do iii = 1, CC_n_occ
          do jjj = 1, CC_n_occ
            aaa = CC_index_aa_D(CC_mpi_did+1) - 1 + a_run
            t_out(a_run,bbb,iii,jjj) = t_out(a_run,bbb,iii,jjj) &
                                     + CC_t_s(iii,aaa,k1) * CC_t_s(jjj,bbb,k3)
          end do
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end if


  End Subroutine CC_3d_get_tau

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_get_intl_kcld(job,i_grp_task,intl_out,local_flag)

! This subroutine gives the k_points determined by a specific k_pattern.
! momentum conservation law: k1 + k3 = k2 + k4 (within first Brillioun zone)

  Implicit None
  Character (len=1) , intent(in) :: job
  Integer , intent(in) :: i_grp_task
  Logical , optional :: local_flag

  Double complex , dimension(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                             CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                             CC_mem_k1(i_grp_task)) :: intl_out

  Integer (kind=8) :: s_tmp
  Integer :: k_start,k_end,n_k1,source_id,i_tmp,j_tmp
  Integer :: kkk,lll,ccc,ddd,k1,k2,k3,k4

  if (.not.(present(local_flag))) then
    local_flag = .true.
  end if

  ! k1 = kk, k2 = kc, k3 = kl, k4 = kd
  n_k1 = CC_mem_k1(i_grp_task)
  k_start = CC_index_k1(i_grp_task)
  k_end = k_start - 1 + n_k1

  if (CC_sv_strategy.eq.1) then

    do k1 = 1,n_k_points
      do k3 = 1,n_k_points
        do kkk = 1, CC_n_occ
          do lll = 1, CC_n_occ

            if (job.eq.'W') then

              intl_out(lll,:,k3,kkk,:,k1,:) = &
                      2.0D0 * CC_intl_kcld_A(lll,:,k3,kkk,:,k1,k_start:k_end) &
                            - CC_intl_kcld_A(kkk,:,k1,lll,:,k3,k_start:k_end)

            else if (job.eq.'D') then

              intl_out(lll,:,k3,kkk,:,k1,:) = CC_intl_kcld_A(kkk,:,k1,lll,:,k3,k_start:k_end)
 
            end if

          end do
        end do
      end do
    end do

    if (job.eq.'S') then
      intl_out = CC_intl_kcld_A(:,:,:,:,:,:,k_start:k_end)
    end if
  else if (CC_sv_strategy.eq.2) then

    source_id = i_grp_task - 1

    if (source_id.eq.CC_mpi_gid) then

      do k1 = 1,n_k_points
        do k3 = 1,n_k_points
          do kkk = 1, CC_n_occ
            do lll = 1, CC_n_occ

              if (job.eq.'W') then

                intl_out(lll,:,k3,kkk,:,k1,:) = 2.0D0 * CC_intl_kcld(lll,:,k3,kkk,:,k1,:) &
                                                      - CC_intl_kcld(kkk,:,k1,lll,:,k3,:)

              else if (job.eq.'D') then

                intl_out(lll,:,k3,kkk,:,k1,:) = CC_intl_kcld(kkk,:,k1,lll,:,k3,:)
 
              end if

            end do
          end do
        end do
      end do

      if (job.eq.'S') then
        intl_out = CC_intl_kcld
      end if
    else
      intl_out = 0.0D0
    end if

    if (.not.(local_flag)) then
      i_tmp = CC_n_occ * CC_n_vir * n_k_points
      j_tmp = CC_n_occ * CC_mem_aa_D(CC_mpi_did+1) * n_k_points * n_k1
      s_tmp = Int(i_tmp,8) * Int(j_tmp,8)

      Call CC_mpi_complex_bcast(s_tmp, intl_out, source_id, CC_mpi_comm_group)
    end if

  end if

  End Subroutine CC_3d_get_intl_kcld

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_clean_w()

! This subroutine clean w vector

  Implicit None

  CC_w_s = 0.0D0
  CC_w_d_d = 0.0D0
  CC_w_d_nd = 0.0D0

  End Subroutine CC_3d_clean_w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_ev_diff(iii,aaa,ki,ka,E_diff)

  ! This subroutine calculates the difference between rrr and sss eigenvalue
  ! E_diff = KS_eigenvalue(sss,ka) - KS_eigenvalue(rrr,ki)

  Implicit None
  Integer , intent(in) :: iii,aaa,ki,ka
  Double precision , intent(out) :: E_diff
  Integer :: rrr,sss

  sss = aaa + CC_n_elec
  rrr = iii + CC_n_fc
  E_diff = KS_eigenvalue(sss,1,ka) - KS_eigenvalue(rrr,1,ki)

  End Subroutine CC_3d_ev_diff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_w2t()

  Implicit None

  Integer :: k1,k2,k3,k4,k1_start,k1_end,k2_start,k2_end,k_run
  Integer :: k_pattern,n_k1,n_k2,k_13,k_index,k_12,k_34

  Integer (kind = 8) :: s_tmp
  Integer :: i_grp_task,i_dom_task,source_id,i_task_run
  Integer :: iii,jjj,aaa,bbb,i_tmp,j_tmp,code_ab,ab_start,ab_end,a_tmp,n_ab

  Integer :: n_a,a_start,a_end,a_run,b_run,a_index,b_index

  Integer :: i_send,i_recv,req1,req2
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Integer :: errnum

  Logical :: ex_flag1,ex_flag2

  Double complex , dimension(:,:,:,:,:,:) , allocatable :: w_send,w_recv,w_tmp
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: t_tmp

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k1_start = CC_index_k1(CC_mpi_gid+1)
  k1_end = k1_start - 1 + n_k1

  ! For CC_t_s
  CC_t_s = CC_w_s
  
  n_ab = CC_mem_ab_D(CC_mpi_did+1,3)
  ! For CC_t_d_nd
  CC_t_d_nd = CC_w_d_nd
  k_13 = k1_start - 1
  do k_run = 1, n_k1
    k_13 = k_13 + 1
    do k1 = 1, n_k_points
      Call CC_3d_k_minus(k_13,k1,k3)

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(a_run,iii,jjj,code_ab,aaa,bbb)
      !$OMP DO COLLAPSE(3)
      do a_run = 1, n_ab
        do iii = 1, CC_n_occ
          do jjj = 1, CC_n_occ
            code_ab = CC_index_ab_D(CC_mpi_did+1,3) - 1 + a_run
            Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)
            CC_t_d_nd(iii,jjj,a_run,k1,k1,k_run) = CC_t_d_nd(iii,jjj,a_run,k1,k1,k_run) &
                                                 + CC_t_s(iii,aaa,k1) * CC_t_s(jjj,bbb,k3)
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end do
  end do

  ! For CC_t_d_d
  CC_t_d_d = CC_w_d_d

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  k_13 = k1_start - 1
  do k_run = 1, n_k1
    k_13 = k_13 + 1
    do k1 = 1, n_k_points
      Call CC_3d_k_minus(k_13,k1,k3)

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(a_run,iii,jjj,aaa)
      !$OMP DO COLLAPSE(3)
      do a_run = 1, n_a
        do iii = 1, CC_n_occ
          do jjj = 1, CC_n_occ
            aaa = CC_index_aa_D(CC_mpi_did+1) - 1 + a_run
            CC_t_d_d(iii,jjj,a_run,k1,k1,k_run) = CC_t_d_d(iii,jjj,a_run,k1,k1,k_run) &
                                                + CC_t_s(iii,aaa,k1) * CC_t_s(jjj,aaa,k3)
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end do
  end do

  ! For CC_t_d and CC_t_d_A
  ! k1 = ka, k2 = ki, k3 = kb, k4 = kj
  ! For w_nd
  i_grp_task = CC_mpi_gid
  i_send = CC_mpi_gid + 1
  i_recv = CC_mpi_gid + 1

  Allocate(w_tmp(CC_n_occ,CC_n_occ,CC_mem_ab_D(CC_mpi_did+1,3), &
                 n_k_points,n_k_points,n_k1),stat=errnum)
  Call check_allocation(errnum,'w_tmp in CC')

  w_tmp = CC_w_d_nd

  do i_task_run = 1, CC_mpi_group_size

    i_recv = i_recv + 1
    if (i_recv.gt.CC_mpi_group_size) then
      i_recv = 1
    end if

    i_send = i_send - 1
    if (i_send.lt.1) then
      i_send = CC_mpi_group_size
    end if
    
    i_grp_task = i_grp_task + 1
    if (i_grp_task.gt.CC_mpi_group_size) then
      i_grp_task = 1
    end if

    if (i_recv.ne.CC_mpi_gid+1) then

      Allocate(w_recv(CC_n_occ,CC_n_occ,CC_mem_ab_D(CC_mpi_did+1,3), &
                      n_k_points,n_k_points,CC_mem_k1(i_recv)),stat=errnum)
      Call check_allocation(errnum,'w_recv in CC')

      Allocate(w_send(CC_n_occ,CC_n_occ,CC_mem_ab_D(CC_mpi_did+1,3), &
                      n_k_points,n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
      Call check_allocation(errnum,'w_send in CC')

      w_send = CC_w_d_nd

      i_tmp = CC_n_occ * CC_n_occ * CC_mem_ab_D(CC_mpi_did+1,3) * n_k_points ** 2
      s_tmp = Int(i_tmp,8) * Int(CC_mem_k1(CC_mpi_gid+1),8)

      Call CC_mpi_complex_isend(s_tmp,w_send,i_send-1,100,req1,CC_mpi_comm_group)
 
      s_tmp = Int(i_tmp,8) * Int(CC_mem_k1(i_recv),8)
      Call CC_mpi_complex_irecv(s_tmp,w_recv,i_recv-1,100,req2,CC_mpi_comm_group)

    end if

    n_k2 = CC_mem_k1(i_grp_task)
    k2_start = CC_index_k1(i_grp_task)
    k2_end = k2_start - 1 + n_k2

    k_13 = k2_start - 1
    do k_run = 1, n_k2

      k_13 = k_13 + 1

      do k1 = 1, n_k_points

        if (CC_sv_strategy.eq.1) then
          ex_flag1 = .true.
          ex_flag2 = .true.
        else if (CC_sv_strategy.eq.2) then
          ex_flag1 = .false.
          ex_flag2 = .false.
        end if

        Call CC_3d_k_minus(k_13,k1,k3)

        if ((k1.ge.k1_start).and.(k1.le.k1_end)) then
          ex_flag1 = .true.
        end if

        if ((k3.ge.k1_start).and.(k3.le.k1_end)) then
          ex_flag2 = .true.
        end if

        if ((ex_flag1).or.(ex_flag2)) then
          Call CC_3d_w2t_data_exchange(ex_flag1,ex_flag2,k1,k3,w_tmp(:,:,:,:,k1,k_run))
        end if

      end do
    end do

    Deallocate(w_tmp)

    if (i_recv.ne.CC_mpi_gid+1) then

      Call MPI_WAIT(req1,stat1,errnum) 
      Call MPI_WAIT(req2,stat2,errnum)

      Allocate(w_tmp(CC_n_occ,CC_n_occ,CC_mem_ab_D(CC_mpi_did+1,3), &
                     n_k_points,n_k_points,CC_mem_k1(i_recv)),stat=errnum)
      Call check_allocation(errnum,'w_tmp in CC')

      w_tmp = w_recv

      Deallocate(w_recv,w_send)

    end if

  end do

  ! k1 = ka, k2 = ki, k3 = kb, k4 = kj
  ! For w_d
  i_grp_task = CC_mpi_gid
  i_send = CC_mpi_gid + 1
  i_recv = CC_mpi_gid + 1

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k1_start = CC_index_k1(CC_mpi_gid+1)
  k1_end = k1_start - 1 + n_k1

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  Allocate(w_tmp(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1), &
                 n_k_points,n_k_points,n_k1),stat=errnum)
  Call check_allocation(errnum,'w_tmp in CC')

  w_tmp = CC_w_d_d

  do i_task_run = 1, CC_mpi_group_size

    i_recv = i_recv + 1
    if (i_recv.gt.CC_mpi_group_size) then
      i_recv = 1
    end if

    i_send = i_send - 1
    if (i_send.lt.1) then
      i_send = CC_mpi_group_size
    end if
    
    i_grp_task = i_grp_task + 1
    if (i_grp_task.gt.CC_mpi_group_size) then
      i_grp_task = 1
    end if

    if (i_recv.ne.CC_mpi_gid+1) then

      Allocate(w_recv(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1), &
                      n_k_points,n_k_points,CC_mem_k1(i_recv)),stat=errnum)
      Call check_allocation(errnum,'w_recv in CC')

      Allocate(w_send(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1), &
                      n_k_points,n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
      Call check_allocation(errnum,'w_send in CC')

      w_send = CC_w_d_d

      i_tmp = CC_n_occ * CC_n_occ * CC_mem_aa_D(CC_mpi_did+1) * n_k_points ** 2
      s_tmp = Int(i_tmp,8) * Int(CC_mem_k1(CC_mpi_gid+1),8)

      Call CC_mpi_complex_isend(s_tmp,w_send,i_send-1,100,req1,CC_mpi_comm_group)
 
      s_tmp = Int(i_tmp,8) * Int(CC_mem_k1(i_recv),8)
      Call CC_mpi_complex_irecv(s_tmp,w_recv,i_recv-1,100,req2,CC_mpi_comm_group)

    end if

    n_k2 = CC_mem_k1(i_grp_task)
    k2_start = CC_index_k1(i_grp_task)
    k2_end = k2_start - 1 + n_k2

    k_13 = k2_start - 1
    do k_run = 1, n_k2

      k_13 = k_13 + 1

      do k1 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k1,k3)
         
        do k2 = 1, n_k_points

          Call CC_3d_k_minus(k_13,k2,k4)

          if (CC_sv_strategy.eq.1) then

            aaa = a_start - 1
            do a_run = 1, n_a
              aaa = aaa + 1
              do iii = 1, CC_n_occ
                CC_t_d_A(:,aaa,k4,iii,a_run,k2,k1) = w_tmp(iii,:,a_run,k2,k1,k_run)
              end do
            end do

            aaa = a_start - 1
            do a_run = 1, n_a
              aaa = aaa + 1
              CC_t_d_A(:,aaa,k2,:,a_run,k4,k3) = w_tmp(:,:,a_run,k2,k1,k_run)
            end do

          else if (CC_sv_strategy.eq.2) then

            if ((k1.ge.k1_start).and.(k1.le.k1_end)) then
 
              k_index = k1 - CC_index_k1(CC_mpi_gid+1) + 1
              aaa = a_start - 1
              do a_run = 1, n_a
                aaa = aaa + 1
                do iii = 1, CC_n_occ
                  CC_t_d(:,aaa,k4,iii,a_run,k2,k_index) = w_tmp(iii,:,a_run,k2,k1,k_run)
                end do
              end do

            end if

            if ((k3.ge.k1_start).and.(k3.le.k1_end)) then
              k_index = k3 - CC_index_k1(CC_mpi_gid+1) + 1
              aaa = a_start - 1
              do a_run = 1, n_a
                aaa = aaa + 1
                CC_t_d(:,aaa,k2,:,a_run,k4,k_index) = w_tmp(:,:,a_run,k2,k1,k_run)
              end do
            end if
          end if

        end do
      end do
    end do

    Deallocate(w_tmp)

    if (i_recv.ne.CC_mpi_gid+1) then

      Call MPI_WAIT(req1,stat1,errnum) 
      Call MPI_WAIT(req2,stat2,errnum)

      Allocate(w_tmp(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1), &
                     n_k_points,n_k_points,CC_mem_k1(i_recv)),stat=errnum)
      Call check_allocation(errnum,'w_tmp in CC')

      w_tmp = w_recv

      Deallocate(w_recv,w_send)

    end if

  end do

  !write(myid+80,*) 'CC_t_s'
  !do k1 = 1, n_k_points
  !  write(myid+80,*) 'k1',k1
  !  do iii = 1, CC_n_occ
  !    do aaa = 1, CC_n_vir
  !      write(myid+80,*) iii,aaa,CC_t_s(iii,aaa,k1)
  !    end do
  !  end do
  !end do
 
  !write(myid+80,*) 'CC_w_nd'
  !k_13 = CC_index_k1(CC_mpi_gid+1) - 1
  !do k_run = 1, CC_mem_k1(CC_mpi_gid+1)
  !  k_13 = k_13 + 1
  !  do k1 = 1, n_k_points
  !    do k2 = 1, n_k_points
  !      Call CC_3d_k_minus(k_13,k1,k3)
  !      Call CC_3d_k_minus(k_13,k2,k4)
  !      write(myid+80,*) 'k_pattern',k1,k2,k3,k4
  !      do iii = 1, CC_n_occ 
  !        do jjj = 1, CC_n_occ
  !          code_ab = CC_index_ab_D(CC_mpi_did+1,3) - 1
  !          do a_run = 1, CC_mem_ab_D(CC_mpi_did+1,3)
  !            code_ab = code_ab + 1
  !            Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)
  !            write(myid+80,*) iii,jjj,aaa,bbb,CC_w_d_nd(iii,jjj,a_run,k2,k1,k_run)
  !          end do
  !        end do
  !      end do
  !    end do
  !  end do
  !end do

  !write(myid+80,*) 'CC_w_d'
  !k_13 = CC_index_k1(CC_mpi_gid+1) - 1
  !do k_run = 1, CC_mem_k1(CC_mpi_gid+1)
  !  k_13 = k_13 + 1
  !  do k1 = 1, n_k_points
  !    do k2 = 1, n_k_points
  !      Call CC_3d_k_minus(k_13,k1,k3)
  !      Call CC_3d_k_minus(k_13,k2,k4)
  !      write(myid+80,*) 'k_pattern',k1,k2,k3,k4
  !      do iii = 1, CC_n_occ 
  !        do jjj = 1, CC_n_occ
  !          aaa = CC_index_aa_D(CC_mpi_did+1) - 1
  !          do a_run = 1, CC_mem_aa_D(CC_mpi_did+1)
  !            aaa = aaa + 1
  !            bbb = aaa
  !            write(myid+80,*) iii,jjj,aaa,bbb,CC_w_d_d(iii,jjj,a_run,k2,k1,k_run)
  !          end do
  !        end do
  !      end do
  !    end do
  !  end do
  !end do
 
  !write(myid+80,*) 'CC_t_d'
  !do k1 = 1, n_k_points
  !  do k2 = 1, n_k_points
  !    do k3 = 1, n_k_points
  !      Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)
  !      write(myid+80,*) 'k_pattern',k1,k2,k3,k4
  !      do iii = 1, CC_n_occ 
  !        aaa = CC_index_aa_D(CC_mpi_did+1) - 1
  !        do a_run = 1, CC_mem_aa_D(CC_mpi_did+1)
  !          aaa =aaa + 1
  !          do jjj = 1, CC_n_occ
  !            do bbb = 1, CC_n_vir
  !              write(myid+80,*) iii,jjj,aaa,bbb,CC_t_d_A(jjj,bbb,k4,iii,a_run,k2,k1)
  !            end do
  !          end do
  !        end do
  !      end do
  !    end do
  !  end do
  !end do


  End Subroutine CC_3d_w2t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_w2t_data_exchange(ex_flag1,ex_flag2,k1,k3,w_tmp)

  Implicit None

  Integer , intent(in) :: k1,k3
  Logical , intent(in) :: ex_flag1,ex_flag2
  Integer :: k2,k4,k_run,k_13,k_index

  Integer (kind = 8) :: s_tmp
  Integer :: source_id,i_task
  Integer :: iii,jjj,aaa,bbb,code_ab,ab_start,ab_end,n_ab,ab_run,a_start,a_tmp,n_a,a_end
  Integer :: i_tmp,j_tmp,errnum

  Double complex , dimension(CC_n_occ,CC_n_occ, &
                             CC_mem_ab_D(CC_mpi_did+1,3),n_k_points) :: w_tmp

  Double complex , dimension(:,:,:,:) , allocatable :: w_tmp2

  a_start = CC_index_aa_D(CC_mpi_did+1)
  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  Call CC_3d_k_plus(k_13,k1,k3)

  do i_task = 1, CC_mpi_domain_size

    n_ab = CC_mem_ab_D(i_task,3)
    ab_start = CC_index_ab_D(i_task,3)
    ab_end = ab_start - 1 + n_ab

    source_id = i_task - 1
    s_tmp = Int(CC_n_occ,8) * Int(CC_n_occ,8) * Int(n_ab,8) * Int(n_k_points,8)

    Allocate(w_tmp2(CC_n_occ,CC_n_occ,n_ab,n_k_points),stat=errnum)
    Call check_allocation(errnum,'w_tmp2 in CC')

    if (source_id.eq.CC_mpi_did) then
      w_tmp2 = w_tmp
    else
      w_tmp2 = 0.0D0
    end if

    Call CC_mpi_complex_bcast(s_tmp, w_tmp2, source_id, CC_mpi_comm_domain)

    do k2 = 1, n_k_points

      Call CC_3d_k_minus(k_13,k2,k4)

      code_ab = ab_start - 1
      do ab_run = 1, n_ab

        code_ab = code_ab + 1

        Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)

        if (ex_flag1) then
          if ((aaa.ge.a_start).and.(aaa.le.a_end)) then

            a_tmp = aaa - a_start + 1

            if (CC_sv_strategy.eq.1) then

              k_index = k1
              do iii = 1, CC_n_occ
                CC_t_d_A(:,bbb,k4,iii,a_tmp,k2,k_index) = w_tmp2(iii,:,ab_run,k2)
              end do

            else if (CC_sv_strategy.eq.2) then

              k_index = k1 - CC_index_k1(CC_mpi_gid+1) + 1
              do iii = 1, CC_n_occ
                CC_t_d(:,bbb,k4,iii,a_tmp,k2,k_index) = w_tmp2(iii,:,ab_run,k2)
              end do

            end if
          end if
        end if

        if (ex_flag2) then

          if ((bbb.ge.a_start).and.(bbb.le.a_end)) then

            a_tmp = bbb - a_start + 1

            if (CC_sv_strategy.eq.1) then

              k_index = k3
              CC_t_d_A(:,aaa,k2,:,a_tmp,k4,k_index) = w_tmp2(:,:,ab_run,k2)

            else if (CC_sv_strategy.eq.2) then

              k_index = k3 - CC_index_k1(CC_mpi_gid+1) + 1
              CC_t_d(:,aaa,k2,:,a_tmp,k4,k_index) = w_tmp2(:,:,ab_run,k2)

            end if
          end if
        end if

      end do
    end do

    Deallocate(w_tmp2)
  end do

  End Subroutine CC_3d_w2t_data_exchange

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

End Module CC_3d

