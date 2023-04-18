  Subroutine CC_3d_diis()

  Use CC_3d

  Implicit None

  Call CC_3d_Jacob()

  Call CC_3d_DIIS_save('t',CC_DIIS_m_bgn)

  if (CC_i_scf.eq.1) then
    Call CC_3d_DIIS_Bmat()
  end if

  if ((CC_i_scf.ge.2).and.(Mod((CC_i_scf - 1),CC_DIIS_step).eq.0)) then
    Call CC_3d_DIIS_solution()
  end if

  CC_DIIS_ndc = CC_DIIS_ndc + 1

  if (CC_DIIS_ndc.gt.CC_DIIS_n_sv) then
    CC_DIIS_ndc = CC_DIIS_n_sv
  end if

  CC_DIIS_m_bgn = CC_DIIS_m_bgn + 1
  if (CC_DIIS_m_bgn.gt.CC_DIIS_ndc) then
    CC_DIIS_m_bgn = 1
  end if

  Call CC_3d_DIIS_save('r',CC_DIIS_m_bgn)

  End Subroutine CC_3d_diis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_DIIS_solution()

  Use CC_3d

  Implicit None

  Integer :: n_k1,k_start,k_end,n_a,a_start,a_end
  Integer :: iii,jjj,aaa,bbb,kkk,lll
  Integer :: i_t,NNN,INFO,i_bas,errnum
  Integer , dimension(:) , allocatable :: IPIV
  Double complex :: D_element,Rkl,Rlk,ak,t_new,t_best
  Double precision :: et1,et2
  Double complex , dimension(:,:) , allocatable :: Rmat,tau
  Double complex , dimension(:) , allocatable :: avec

  NNN = CC_DIIS_ndc + 1
  Call CC_3d_DIIS_Bmat()

  ! Form R matrix and a vector
  Allocate(Rmat(NNN,NNN),stat=errnum)
  Call check_allocation(errnum, 'Rmat in CC')

  Allocate(tau(NNN,1),stat=errnum)
  Call check_allocation(errnum, 'tau in CC')

  Rmat(1,2:NNN) = -1.0D0
  Rmat(2:NNN,1) = -1.0D0

  tau(2:NNN,1) = 0.0D0

  Rmat(1,1) = 0.0D0
  tau(1,1) = -1.0D0

  do kkk = 2, NNN
    do lll = 2, NNN 
      Rmat(kkk,lll) = CC_DIIS_Bmat(kkk - 1,lll - 1)
    end do
  end do

  Allocate(IPIV(NNN),stat=errnum)
  Call check_allocation(errnum, 'IPIV in CC')

  Call ZGESV(NNN,1,Rmat,NNN,IPIV,tau,NNN,INFO)

  Call CC_3d_clean_w()

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  do kkk = 1, NNN - 1
 
    ! w_s
    CC_w_s = CC_w_s + CC_DIIS_t_s(:,:,:,kkk) * tau(kkk+1,1)

    ! w_d
    CC_w_d_d = CC_w_d_d + CC_DIIS_t_d(:,:,:,:,:,:,kkk) * tau(kkk+1,1)

    ! w_nd
    CC_w_d_nd = CC_w_d_nd + CC_DIIS_t_nd(:,:,:,:,:,:,kkk) * tau(kkk+1,1)

  end do

  Deallocate(Rmat,tau,IPIV)

  End Subroutine CC_3d_DIIS_solution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_DIIS_Bmat()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_run,k_13,k_start,n_k1,k_end
  Integer :: iii,jjj,kkk,lll,nnn,j_tmp,ggg,a_start,a_end
  Integer :: OMP_GET_NUM_THREADS,oid,OMP_GET_THREAD_NUM,errnum
  Integer , dimension(1) :: shp1
  Double complex :: Bij,ZDOTC
  Double complex , dimension(:) , allocatable :: w1,w2,rlt

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + CC_mem_aa_D(CC_mpi_did+1)

  do ggg = 1, CC_DIIS_ndc

    Bij = 0.0D0

    ! w_s
    nnn = CC_n_occ * CC_mem_aa_D(CC_mpi_did+1) * n_k1
    shp1(1) = nnn

    Allocate(w1(nnn),stat=errnum)
    Call check_allocation(errnum,'w1 in CC')

    Allocate(w2(nnn),stat=errnum)
    Call check_allocation(errnum,'w1 in CC')

    w1 = reshape(CC_DIIS_r_s(:,a_start:a_end,k_start:k_end,ggg),shp1)
    w2 = reshape(CC_DIIS_r_s(:,a_start:a_end,k_start:k_end,CC_DIIS_m_bgn),shp1)

    Bij = ZDOTC(nnn,w1,1,w2,1)

    Deallocate(w1,w2)

    ! w_d_nd
    nnn = CC_n_occ * CC_n_occ * CC_mem_ab_D(CC_mpi_did+1,3) * n_k_points ** 2 * n_k1
    shp1(1) = nnn

    Allocate(w1(nnn),stat=errnum)
    Call check_allocation(errnum,'w1 in CC')

    Allocate(w2(nnn),stat=errnum)
    Call check_allocation(errnum,'w1 in CC')

    w1 = reshape(CC_DIIS_r_nd(:,:,:,:,:,:,ggg),shp1)
    w2 = reshape(CC_DIIS_r_nd(:,:,:,:,:,:,CC_DIIS_m_bgn),shp1)

    Bij = Bij + ZDOTC(nnn,w1,1,w2,1)

    Deallocate(w1,w2)

    ! For w_d_d
    !$OMP PARALLEL Default(shared) &
    !$OMP Private(iii,jjj,k2,k1,k_run,k_13,k3,k4,oid)
    !$OMP SINGLE
    nnn = 0
    !$ nnn = OMP_GET_NUM_THREADS()
    if (nnn.eq.0) then
      nnn = 1
    end if

    Allocate(rlt(nnn),stat=errnum)
    Call check_allocation(errnum,'rlt in CC')

    rlt = 0.0D0
    !$OMP END SINGLE

    oid = 0
    !$ oid = OMP_GET_THREAD_NUM()
 
    !$OMP DO COLLAPSE(5)
    do iii = 1, CC_n_occ
      do jjj = 1, CC_n_occ
        do k2 = 1, n_k_points
          do k1 = 1, n_k_points
            do k_run = 1, n_k1

              k_13 = k_start - 1 + k_run
              Call CC_3d_k_minus(k_13,k1,k3)
              Call CC_3d_k_minus(k_13,k2,k4)

              if (k1.lt.k3) then
                rlt(oid+1) = rlt(oid+1) &
                           + dot_product(CC_DIIS_r_d(iii,jjj,:,k2,k1,k_run,ggg), &
                                         CC_DIIS_r_d(iii,jjj,:,k2,k1,k_run,CC_DIIS_m_bgn))
              else if (k1.eq.k3) then
                if (iii.lt.jjj) then
                  rlt(oid+1) = rlt(oid+1) &
                             + dot_product(CC_DIIS_r_d(iii,jjj,:,k2,k1,k_run,ggg), &
                                           CC_DIIS_r_d(iii,jjj,:,k2,k1,k_run,CC_DIIS_m_bgn))
                else if (iii.eq.jjj) then
                  if (k2.le.k4) then
                    rlt(oid+1) = rlt(oid+1) &
                               + dot_product(CC_DIIS_r_d(iii,jjj,:,k2,k1,k_run,ggg), &
                                             CC_DIIS_r_d(iii,jjj,:,k2,k1,k_run,CC_DIIS_m_bgn))
                  end if
                end if
              end if

            end do
          end do
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    do iii = 1, nnn
      Bij = Bij + rlt(iii)
    end do

    Deallocate(rlt)

    Call CC_mpi_complex_number(Bij, MPI_COMM_WORLD)

    CC_DIIS_Bmat(ggg,CC_DIIS_m_bgn) = Bij
    CC_DIIS_Bmat(CC_DIIS_m_bgn,ggg) = Conjg(Bij)

  end do 

  End Subroutine CC_3d_DIIS_Bmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_DIIS_save(jobtype,i_sv)

  Use CC_3d

  Implicit None

  Character (len = 1) , intent(in) :: jobtype
  Integer , intent(in) :: i_sv

  if (jobtype.eq.'t') then

    CC_DIIS_t_s(:,:,:,i_sv) = CC_w_s
    CC_DIIS_t_d(:,:,:,:,:,:,i_sv) = CC_w_d_d
    CC_DIIS_t_nd(:,:,:,:,:,:,i_sv) = CC_w_d_nd

    CC_DIIS_r_s(:,:,:,i_sv) = CC_DIIS_t_s(:,:,:,i_sv) - CC_DIIS_r_s(:,:,:,i_sv)
    CC_DIIS_r_d(:,:,:,:,:,:,i_sv) = CC_DIIS_t_d(:,:,:,:,:,:,i_sv) &
                                  - CC_DIIS_r_d(:,:,:,:,:,:,i_sv)
    CC_DIIS_r_nd(:,:,:,:,:,:,i_sv) = CC_DIIS_t_nd(:,:,:,:,:,:,i_sv) &
                                   - CC_DIIS_r_nd(:,:,:,:,:,:,i_sv)

  else if (jobtype.eq.'r') then

    CC_DIIS_r_s(:,:,:,i_sv) = CC_w_s
    CC_DIIS_r_d(:,:,:,:,:,:,i_sv) = CC_w_d_d
    CC_DIIS_r_nd(:,:,:,:,:,:,i_sv) = CC_w_d_nd

  end if

  End Subroutine CC_3d_DIIS_save

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_E_corr()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_pattern,k_start,k_end,k_run,n_k1
  Integer :: i_tmp,j_tmp,k_tmp,a_run,iii,jjj,aaa,bbb,n_a
  Integer :: errnum
  Integer , dimension(1) :: shp1

  Double complex , dimension(:) , allocatable :: A_vec,B_vec
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_tmp,t_tmp

  Double complex :: rlt,ZDOTU,ZDOTC

  rlt = 0.0D0

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  Allocate(t_tmp(CC_n_occ,CC_n_vir,n_k_points, &
                 CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                 CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'t_tmp in CC')

  Call CC_3d_get_t_grp('S',CC_mpi_gid+1,t_tmp,.true.)

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  k1 = k_start - 1
  do k_run = 1, n_k1

    k1 = k1 + 1

    do k3 = 1, n_k_points

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(aaa,a_run,iii,jjj,bbb)
      !$OMP DO COLLAPSE(4)
      do a_run = 1, n_a
        do iii = 1, CC_n_occ
          do bbb = 1, CC_n_vir
            do jjj = 1, CC_n_occ 
              aaa = CC_index_aa_D(CC_mpi_did+1) - 1 + a_run
              t_tmp(jjj,bbb,k3,iii,a_run,k1,k_run) = t_tmp(jjj,bbb,k3,iii,a_run,k1,k_run) &
                                                   + CC_t_s(jjj,bbb,k3) * CC_t_s(iii,aaa,k1)
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end do
  end do

  i_tmp = CC_n_occ * CC_n_vir * n_k_points * CC_n_occ &
        * CC_mem_aa_D(CC_mpi_did+1) * n_k_points * n_k1

  Allocate(A_vec(i_tmp),stat=errnum)
  Call check_allocation(errnum,'A_vec in CC')

  shp1(1) = i_tmp
  A_vec = reshape(t_tmp,shp1)

  Deallocate(t_tmp)

  Allocate(intl_tmp(CC_n_occ,CC_n_vir,n_k_points, &
                    CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                    CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'t_tmp in CC')

  Call CC_3d_get_intl_kcld('W',CC_mpi_gid+1,intl_tmp,.true.)

  Allocate(B_vec(i_tmp),stat=errnum)
  Call check_allocation(errnum,'B_vec in CC')

  B_vec = reshape(intl_tmp,shp1)

  rlt = ZDOTU(i_tmp,A_vec,1,B_vec,1)

  Call CC_mpi_complex_number(rlt, MPI_COMM_WORLD)
  CC_E_corr = dble(rlt) / dble(n_k_points)

  Deallocate(intl_tmp,A_vec,B_vec)

  End Subroutine CC_3d_E_corr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_calc_Norm()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_pattern,k_start,k_end,k_run
  Integer :: i_tmp,j_tmp,k_tmp,errnum
  Integer , dimension(1) :: shp1

  Double complex , dimension(:) , allocatable :: A_vec,B_vec
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: intl_tmp,t_tmp

  Double complex :: rlt,ZDOTU,ZDOTC

  CC_Norm = 0.0D0

  if (myid.eq.0) then

    i_tmp = CC_n_occ * CC_n_vir * n_k_points

    Allocate(A_vec(i_tmp),stat=errnum)
    Call check_allocation(errnum,'A_tmp in CC')

    shp1(1) = i_tmp

    A_vec = reshape(CC_t_s,shp1)
    B_vec = conjg(A_vec)

    CC_Norm = ZDOTC(i_tmp,A_vec,1,A_vec,1)

    Deallocate(A_vec)

  end if

  CC_Norm = sqrt(1.0D0 + CC_Norm)

  End Subroutine CC_3d_calc_Norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_Jacob()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,n_k1,k_start,k_run,k_13
  Integer :: iii,jjj,aaa,bbb,n_ab,ab_start,code_ab,ab_run,n_a,a_start,a_run
  Double precision :: CC_E_diff_0_s,et1,et2

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1) - 1

  ! For singlet excitation

  !$OMP PARALLEL Default(Shared) &
  !$OMP Private(k1,iii,aaa,et1)
  !$OMP DO COLLAPSE(3)
  do k1 = 1, n_k_points
    do iii = 1, CC_n_occ
      do aaa = 1, CC_n_vir
        Call CC_3d_ev_diff(iii,aaa,k1,k1,et1)
        CC_w_s(iii,aaa,k1) = - CC_w_s(iii,aaa,k1) / et1
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  n_ab = CC_mem_ab_D(CC_mpi_did+1,3)
  ab_start = CC_index_ab_D(CC_mpi_did+1,3)

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)

  ! For double excitation
  ! k1 = ka, k2 = ki, k3 = kb, k4 = kj
  k_13 = k_start
  do k_run = 1, n_k1

    k_13 = k_13 + 1

    do k1 = 1, n_k_points

      Call CC_3d_k_minus(k_13,k1,k3)

      ! For w_d_nd
      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(k2,k4,iii,jjj,ab_run,code_ab,aaa,bbb,et1,et2,CC_E_diff_0_s)
      !$OMP DO COLLAPSE(4)
      do k2 = 1, n_k_points
        do iii = 1, CC_n_occ
          do jjj = 1, CC_n_occ
            do ab_run = 1, n_ab

              Call CC_3d_k_minus(k_13,k2,k4)

              code_ab = ab_start - 1 + ab_run
              Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)
  
              Call CC_3d_ev_diff(iii,aaa,k2,k1,et1)
              Call CC_3d_ev_diff(jjj,bbb,k4,k3,et2)
              CC_E_diff_0_s = - et1 - et2
  
              CC_w_d_nd(iii,jjj,ab_run,k2,k1,k_run) = &
                  (conjg(CC_intl_iajb_nd(iii,jjj,ab_run,k2,k1,k_run)) &
                   + CC_w_d_nd(iii,jjj,ab_run,k2,k1,k_run)) / CC_E_diff_0_s
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL
 
      ! For w_d_d
      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(k2,k4,iii,jjj,a_run,aaa,et1,et2,CC_E_diff_0_s)
      !$OMP DO COLLAPSE(4)
      do k2 = 1, n_k_points
        do iii = 1, CC_n_occ
          do jjj = 1, CC_n_occ
            do a_run = 1, n_a

              Call CC_3d_k_minus(k_13,k2,k4)

              aaa = a_start - 1 + a_run
              Call CC_3d_ev_diff(iii,aaa,k2,k1,et1)
              Call CC_3d_ev_diff(jjj,aaa,k4,k3,et2)
              CC_E_diff_0_s = - et1 - et2
  
              CC_w_d_d(iii,jjj,a_run,k2,k1,k_run) = &
                  (conjg(CC_intl_iajb_d(iii,jjj,a_run,k2,k1,k_run)) &
                   + CC_w_d_d(iii,jjj,a_run,k2,k1,k_run)) / CC_E_diff_0_s
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL
    end do
  end do

  End Subroutine CC_3d_Jacob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

