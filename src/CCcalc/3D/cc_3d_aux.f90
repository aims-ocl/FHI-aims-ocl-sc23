!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_aux_h_g()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_pattern,k_run,k_start,k_end
  Integer :: iii,jjj,aaa,bbb,kkk,lll,ccc,ddd,i_run,i_tmp,j_tmp,c_run
  Integer :: n_dom,a_start,a_end,i_start,i_end,l_start,l_end
  Integer :: errnum
  Integer , dimension(2) :: shp1
  Integer , dimension(3) :: shp2

  Integer (kind=8) :: s_tmp

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: t_tmp,v_tmp
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: tmp,tmp2
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C

  Double complex :: alpha,beta,rlt

  ! Clear tensors
  CC_h_ik = 0.0D0                                                                                 
  CC_h_ca = 0.0D0
  CC_h_ck = 0.0D0
  CC_g_ik = 0.0D0
  CC_g_ca = 0.0D0

  alpha = 1.0D0
  beta = 0.0D0

  Allocate(v_tmp(CC_n_occ,CC_n_vir,n_k_points, &
                 CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                 CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'v_tmp in CC')

  Call CC_3d_get_intl_kcld('W',CC_mpi_gid+1,v_tmp,.true.)

  Allocate(t_tmp(CC_n_occ,CC_n_vir,n_k_points, &
                 CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                 CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'t_tmp in CC')

  Call CC_3d_get_t_grp('S',CC_mpi_gid+1,t_tmp,.true.)

  ! k1 = ka = kc, k2 = kk, k3 = kd, k4 = kl
  k1 = CC_index_k1(CC_mpi_gid+1) - 1
  do k_run = 1, CC_mem_k1(CC_mpi_gid+1)

    k1 = k1 + 1
    ccc = CC_index_aa_D(CC_mpi_did+1) - 1
    do c_run = 1, CC_mem_aa_D(CC_mpi_did+1)

      ccc = ccc + 1

      !$OMP PARALLEL Default(shared) &
      !$OMP Private(k3,iii,lll,ddd)
      !$OMP DO COLLAPSE(4)
      do k3 = 1, n_k_points
        do iii = 1, CC_n_occ
          do ddd = 1, CC_n_vir
            do lll = 1, CC_n_occ
              t_tmp(lll,ddd,k3,iii,c_run,k1,k_run) = t_tmp(lll,ddd,k3,iii,c_run,k1,k_run) &
                                                   + CC_t_s(iii,ccc,k1) * CC_t_s(lll,ddd,k3)
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL
    end do
  end do

  ! Calculate h_ik(iii,kkk,ki=kk)
  Allocate(tmp(CC_n_occ,CC_n_vir,n_k_points,CC_mem_aa_D(CC_mpi_did+1), &
               CC_mem_k1(CC_mpi_gid+1),CC_n_occ,1),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  i_tmp = CC_n_occ * CC_n_vir * CC_mem_aa_D(CC_mpi_did+1) * n_k_points &
        * CC_mem_k1(CC_mpi_gid+1)

  !-------------test-------------
  !rlt = 0.0D0
  !write(80+myid,*) 'h_ik(1,1)'
  !do ccc = 1, CC_mem_aa_D(CC_mpi_did+1)
  !  do lll = 1, CC_n_occ
  !    do ddd = 1, CC_n_vir
  !      write(80+myid,*) 't',1,ccc,lll,ddd,t_tmp(lll,ddd,1,1,ccc,1,1)
  !      write(80+myid,*) 'v',1,ccc,lll,ddd,v_tmp(lll,ddd,1,1,ccc,1,1)
  !      rlt = rlt + v_tmp(lll,ddd,1,1,ccc,1,1) * t_tmp(lll,ddd,1,1,ccc,1,1)
  !    end do
  !  end do
  !end do
  !write(80+myid,*) 'h_ik(1,1)',rlt
  !------------end test-----------

  Allocate(mat_A(i_tmp,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(i_tmp,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(CC_n_occ,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  shp1(1) = i_tmp
  shp1(2) = CC_n_occ

  do k2 = 1, n_k_points

    do iii = 1, CC_n_occ
      tmp(:,:,:,:,:,iii,1) = t_tmp(:,:,:,iii,:,k2,:)
    end do

    mat_A = reshape(tmp,shp1)

    do kkk = 1, CC_n_occ
      tmp(:,:,:,:,:,kkk,1) = v_tmp(:,:,:,kkk,:,k2,:)
    end do

    mat_B = reshape(tmp,shp1)

    Call Zgemm('T','N',CC_n_occ,CC_n_occ,i_tmp,alpha,mat_A,i_tmp, &
               mat_B,i_tmp,beta,mat_C,CC_n_occ)

    CC_h_ik(:,:,k2) = mat_C

  end do

  Deallocate(tmp,mat_A,mat_B,mat_C)

  ! Calculate h_ca(ccc,aaa,kc=ka)
  Allocate(tmp(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
               CC_mem_k1(CC_mpi_gid+1),CC_n_vir,n_k_points),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  Allocate(tmp2(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                CC_mem_k1(CC_mpi_gid+1),CC_n_vir,n_k_points),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  ! k1 = kd, k2 = kl, k3 = ka = kc, k4 = kk
  do k2 = 1, n_k_points

    do k4 = 1, n_k_points

      k1 = CC_index_k1(CC_mpi_gid+1) - 1
      do k_run = 1, CC_mem_k1(CC_mpi_gid+1)

        k1 = k1 + 1

        Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,3)

        do ccc = 1, CC_n_vir
          do kkk = 1, CC_n_occ
            tmp(kkk,:,:,k4,k_run,ccc,k3) = v_tmp(kkk,ccc,k4,:,:,k2,k_run)
          end do
        end do

        do aaa = 1, CC_n_vir
          do kkk = 1, CC_n_occ
            tmp2(kkk,:,:,k4,k_run,aaa,k3) = t_tmp(kkk,aaa,k4,:,:,k2,k_run)
          end do
        end do

      end do
    end do

  end do

  i_tmp = CC_n_occ * CC_n_occ * CC_mem_aa_D(CC_mpi_did+1) * n_k_points &
        * CC_mem_k1(CC_mpi_gid+1)

  Allocate(mat_A(i_tmp,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(i_tmp,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(CC_n_vir,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  shp1(1) = i_tmp
  shp1(2) = CC_n_vir

  do k3 = 1, n_k_points

    mat_A = reshape(tmp(:,:,:,:,:,:,k3),shp1)
    mat_B = reshape(tmp2(:,:,:,:,:,:,k3),shp1)

    Call Zgemm('T','N',CC_n_vir,CC_n_vir,i_tmp,alpha,mat_A,i_tmp, &
               mat_B,i_tmp,beta,mat_C,CC_n_vir)

    CC_h_ca(:,:,k3) = - mat_C

  end do

  Deallocate(tmp,tmp2,mat_A,mat_B,mat_C)

  ! For h_ck(kkk,ccc,kk=kc)
  Allocate(tmp(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1), &
               CC_mem_k1(CC_mpi_gid+1),1),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  i_tmp = CC_n_occ * CC_n_vir * n_k_points 
  j_tmp = CC_n_occ * CC_mem_aa_D(CC_mpi_did+1) * CC_mem_k1(CC_mpi_gid+1)

  Allocate(mat_A(i_tmp,j_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(i_tmp,1),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(j_tmp,1),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + CC_mem_k1(CC_mpi_gid+1)

  k1 = k_start - 1
  do k_run = 1, CC_mem_k1(CC_mpi_gid+1)
    k1 = k1 + 1
    tmp(:,:,:,:,:,k_run,1) = v_tmp(:,:,:,:,:,k1,k_run)
  end do
 
  shp1(1) = i_tmp
  shp1(2) = j_tmp
  mat_A = reshape(tmp,shp1)

  shp1(1) = i_tmp
  shp1(2) = 1
  mat_B = reshape(CC_t_s,shp1)

  Call Zgemm('T','N',j_tmp,1,i_tmp,alpha,mat_A,i_tmp, &
             mat_B,i_tmp,beta,mat_C,j_tmp)

  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + CC_mem_aa_D(CC_mpi_did+1)

  shp2(1) = CC_n_occ
  shp2(2) = CC_mem_aa_D(CC_mpi_did+1)
  shp2(3) = CC_mem_k1(CC_mpi_gid+1)
  CC_h_ck(:,a_start:a_end,k_start:k_end) = reshape(mat_C,shp2)

  Deallocate(tmp,mat_A,mat_B,mat_C)
  Deallocate(t_tmp,v_tmp)

  ! For g_ik(iii,kkk,ki=kk)
  ! k1 = kl = kc, k2 = ki = kk
  Allocate(v_tmp(CC_mem_ii_D(CC_mpi_did+1),CC_n_vir,CC_mem_k1(CC_mpi_gid+1), &
                 CC_n_occ,CC_n_occ,n_k_points,1),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  i_tmp = CC_mem_ii_D(CC_mpi_did+1) * CC_n_vir * CC_mem_k1(CC_mpi_gid+1)
  j_tmp = CC_n_occ * CC_n_occ * n_k_points

  Allocate(mat_A(i_tmp,j_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(i_tmp,1),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(j_tmp,1),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + CC_mem_k1(CC_mpi_gid+1)

  l_start = CC_index_ii_D(CC_mpi_did+1)
  l_end = l_start - 1 + CC_mem_ii_D(CC_mpi_did+1)

  !$OMP PARALLEL Default(Shared) &
  !$OMP Private(k2,iii,kkk,ccc)
  !$OMP DO COLLAPSE(4)
  do k2 = 1, n_k_points
    do iii = 1, CC_n_occ
      do kkk = 1, CC_n_occ
        do ccc = 1, CC_n_vir
          v_tmp(:,ccc,:,iii,kkk,k2,1) = 2.0D0 * CC_intl_kilc(kkk,ccc,k2,iii,:,k2,:) &
                                              - CC_intl_likc(kkk,ccc,k2,iii,:,k2,:)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
 
  shp1(1) = i_tmp
  shp1(2) = j_tmp
  mat_A = reshape(v_tmp,shp1)

  shp1(1) = i_tmp
  shp1(2) = 1
  mat_B = reshape(CC_t_s(l_start:l_end,:,k_start:k_end),shp1)

  Call Zgemm('T','N',j_tmp,1,i_tmp,alpha,mat_A,i_tmp, &
             mat_B,i_tmp,beta,mat_C,j_tmp)

  shp2(1) = CC_n_occ
  shp2(2) = CC_n_occ
  shp2(3) = n_k_points
  CC_g_ik = reshape(mat_C,shp2)

  Deallocate(v_tmp,mat_A,mat_B,mat_C)

  ! For g_ca(ccc,aaa,kc=ka)
  Allocate(v_tmp(CC_n_occ,CC_n_vir,n_k_points,CC_n_vir,CC_mem_aa_D(CC_mpi_did+1), &
               CC_mem_k1(CC_mpi_gid+1),1),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  i_tmp = CC_n_occ * CC_n_vir * n_k_points 
  j_tmp = CC_n_vir * CC_mem_aa_D(CC_mpi_did+1) * CC_mem_k1(CC_mpi_gid+1)

  Allocate(mat_A(i_tmp,j_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(i_tmp,1),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(j_tmp,1),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + CC_mem_k1(CC_mpi_gid+1)

  k1 = k_start - 1
  do k_run = 1, CC_mem_k1(CC_mpi_gid+1)
    k1 = k1 + 1
    !$OMP PARALLEL Default(Shared) &
    !$OMP Private(k3,ccc,ddd)
    !$OMP DO COLLAPSE(3)
    do k3 = 1, n_k_points
      do ccc = 1, CC_n_vir
        do ddd = 1, CC_n_vir
          v_tmp(:,ddd,k3,ccc,:,k_run,1) = 2.0D0 * CC_intl_ackd(:,ddd,k3,ccc,:,k1,k_run) &
                                              - CC_intl_ackd(:,ccc,k3,ddd,:,k3,k_run)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end do
 
  shp1(1) = i_tmp
  shp1(2) = j_tmp
  mat_A = reshape(v_tmp,shp1)

  shp1(1) = i_tmp
  shp1(2) = 1
  mat_B = reshape(CC_t_s,shp1)

  Call Zgemm('T','N',j_tmp,1,i_tmp,alpha,mat_A,i_tmp, &
             mat_B,i_tmp,beta,mat_C,j_tmp)

  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + CC_mem_aa_D(CC_mpi_did+1)

  shp2(1) = CC_n_vir
  shp2(2) = CC_mem_aa_D(CC_mpi_did+1)
  shp2(3) = CC_mem_k1(CC_mpi_gid+1)
  CC_g_ca(:,a_start:a_end,k_start:k_end) = reshape(mat_C,shp2)

  Deallocate(v_tmp,mat_A,mat_B,mat_C)

  ! Synchronize vectors
  s_tmp = Int(CC_n_occ,8) * Int(CC_n_occ,8) * Int(n_k_points,8)
  Call CC_mpi_complex_allreduce(s_tmp, CC_h_ik, MPI_COMM_WORLD) 
  Call CC_mpi_complex_allreduce(s_tmp, CC_g_ik, MPI_COMM_WORLD)

  s_tmp = Int(CC_n_vir,8) * Int(CC_n_vir,8) * Int(n_k_points,8)
  Call CC_mpi_complex_allreduce(s_tmp, CC_h_ca, MPI_COMM_WORLD) 
  Call CC_mpi_complex_allreduce(s_tmp, CC_g_ca, MPI_COMM_WORLD)

  s_tmp = Int(CC_n_occ,8) * Int(CC_n_vir,8) * Int(n_k_points,8)
  Call CC_mpi_complex_allreduce(s_tmp, CC_h_ck, MPI_COMM_WORLD)

  CC_g_ik = CC_g_ik + CC_h_ik
  CC_g_ca = CC_g_ca + CC_h_ca

  !write(80+myid,*) 'h_ik',CC_h_ik
  !write(80+myid,*) 'h_ca',CC_h_ca
  !write(80+myid,*) 'h_ck',CC_h_ck
  !write(80+myid,*) 'g_ik',CC_g_ik
  !write(80+myid,*) 'g_ca',CC_g_ca

  End Subroutine CC_3d_aux_h_g

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_aux_jk()

! Caluclate intermediate j and k 
! k1 = ka, k2 = ki, k3 = kk, k4 = kc, k5 = kl, k6 = kd
 
  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k5,k6,k_pattern,k_run,n_k1,k_start,k_end
  Integer :: i_run,i_tmp,j_tmp,c_run
  Integer :: n_dom,a_start,a_end,i_start,i_end,n_a,a_run
  Integer :: iii,jjj,aaa,bbb,ccc,ddd,kkk,lll
  Integer :: i_task,source_id,i_grp_task,grp_id,errnum
  Integer , dimension(2) :: shp1
  Integer , dimension(6) :: shp1i
  Integer , dimension(4) :: shp2i

  Integer (kind=8) :: s_tmp

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: t_tmp,v_tmp
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: tmp,tmp2
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C

  Double complex :: alpha,beta

!  write(100,*) 'j_aux'

  ! First term (ai|kc) & (ki|ac)
  !$OMP PARALLEL Default(Shared)
  !$OMP WORKSHARE
  CC_j_aux = 2.0D0 * CC_intl_aikc - CC_intl_kiac
  !$OMP END WORKSHARE
  !$OMP END PARALLEL
 
  CC_k_aux = CC_intl_kiac

  !print*,'intl_k'
  !do iii = 1, CC_n_occ
  !  do aaa = 1, CC_n_vir
  !    do kkk = 1, CC_n_occ
  !      do ccc = 1, CC_n_vir
  !        print*,iii,aaa,kkk,ccc,CC_k_aux(kkk,ccc,1,iii,aaa,1,1)
  !      end do
  !    end do
  !  end do
  !end do

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  alpha = 1.0D0
  beta = 0.0D0

  !print*,'k before 2'
  !do iii = 1, CC_n_occ
  !  do aaa = 1, CC_n_vir
  !    do kkk = 1, CC_n_occ
  !      do ccc = 1, CC_n_vir
  !        print*,iii,aaa,kkk,ccc,CC_k_aux(kkk,ccc,1,iii,aaa,1,1)
  !      end do
  !    end do
  !  end do
  !end do


  ! Second term (ki|lc)
  ! k1 = ka = kl, k2 = ki, k3 = kk, k4 = kc
  do i_task = 1, CC_mpi_domain_size
    
    n_dom = CC_mem_ii_D(i_task)
    i_start = CC_index_ii_D(i_task)
    i_end = i_start - 1 + n_dom

    Allocate(v_tmp(n_dom,CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                   n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
    Call check_allocation(errnum,'v_tmp in CC')

    source_id = i_task - 1
    j_tmp = CC_n_occ * CC_n_vir * n_k_points * CC_n_occ * n_k_points
    s_tmp = Int(j_tmp,8) * Int(n_dom,8) * Int(n_k1,8)

    n_a = CC_mem_aa_D(CC_mpi_did+1)
    a_start = CC_index_aa_D(CC_mpi_did+1) 
    a_end = a_start - 1 + n_a

    Allocate(mat_A(n_dom,j_tmp),stat=errnum)
    Call check_allocation(errnum,'mat_A in CC')

    Allocate(mat_B(n_dom,n_a),stat=errnum)
    Call check_allocation(errnum,'mat_B in CC')

    Allocate(mat_C(j_tmp,n_a),stat=errnum)
    Call check_allocation(errnum,'mat_C in CC')

    Allocate(tmp(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_k_points,n_a,1),stat=errnum)
    Call check_allocation(errnum,'tmp in CC')

    ! For j_aux
    if (source_id.eq.CC_mpi_did) then
      do lll = 1, n_dom
        v_tmp(lll,:,:,:,:,:,:) = 2.0D0 * CC_intl_likc(:,:,:,:,lll,:,:) &
                                       - CC_intl_kilc(:,:,:,:,lll,:,:)
      end do
    else
      v_tmp = 0.0D0
    end if

    Call CC_mpi_complex_bcast(s_tmp, v_tmp, source_id, CC_mpi_comm_domain)

    k1 = k_start - 1
    do k_run = 1, n_k1

      k1 = k1 + 1

      shp1(1) = n_dom
      shp1(2) = j_tmp
      mat_A = reshape(v_tmp(:,:,:,:,:,:,k_run),shp1)

      mat_B = CC_t_s(i_start:i_end,a_start:a_end,k1)

      Call Zgemm('T','N',j_tmp,n_a,n_dom,alpha,mat_A,n_dom, &
                 mat_B,n_dom,beta,mat_C,j_tmp)

      shp1i(1) = CC_n_occ
      shp1i(2) = CC_n_vir
      shp1i(3) = n_k_points
      shp1i(4) = CC_n_occ
      shp1i(5) = n_k_points
      shp1i(6) = n_a

      tmp(:,:,:,:,:,:,1) = reshape(mat_C,shp1i)

      do aaa = 1, n_a
        CC_j_aux(:,:,:,:,aaa,:,k_run) = &
               CC_j_aux(:,:,:,:,aaa,:,k_run) - tmp(:,:,:,:,:,aaa,1)
      end do
    end do

    ! for k_aux
    if (source_id.eq.CC_mpi_did) then
      do lll = 1, n_dom
        v_tmp(lll,:,:,:,:,:,:) = CC_intl_kilc(:,:,:,:,lll,:,:)
      end do
    else
      v_tmp = 0.0D0
    end if

    Call CC_mpi_complex_bcast(s_tmp, v_tmp, source_id, CC_mpi_comm_domain)

    k1 = k_start - 1
    do k_run = 1, n_k1

      k1 = k1 + 1

      shp1(1) = n_dom
      shp1(2) = j_tmp
      mat_A = reshape(v_tmp(:,:,:,:,:,:,k_run),shp1)

      mat_B = CC_t_s(i_start:i_end,a_start:a_end,k1)

      Call Zgemm('T','N',j_tmp,n_a,n_dom,alpha,mat_A,n_dom, &
                 mat_B,n_dom,beta,mat_C,j_tmp)

      shp1i(1) = CC_n_occ
      shp1i(2) = CC_n_vir
      shp1i(3) = n_k_points
      shp1i(4) = CC_n_occ
      shp1i(5) = n_k_points
      shp1i(6) = n_a

      tmp(:,:,:,:,:,:,1) = reshape(mat_C,shp1i)

      do aaa = 1, n_a
        CC_k_aux(:,:,:,:,aaa,:,k_run) = &
               CC_k_aux(:,:,:,:,aaa,:,k_run) - tmp(:,:,:,:,:,aaa,1)
      end do
    end do

    Deallocate(v_tmp,mat_A,mat_B,mat_C,tmp)

  end do

  !print*,'second intl_j'
  !do aaa = 1, CC_n_vir
  !  do ccc = 1, CC_n_vir
  !    print*,aaa,1,1,ccc,CC_j_aux(1,ccc,1,1,aaa,1,1)
  !  end do
  !end do
  !print*,'k before 3'
  !do iii = 1, CC_n_occ
  !  do aaa = 1, CC_n_vir
  !    do kkk = 1, CC_n_occ
  !      do ccc = 1, CC_n_vir
  !        print*,iii,aaa,kkk,ccc,CC_k_aux(kkk,ccc,1,iii,aaa,1,1)
  !      end do
  !    end do
  !  end do
  !end do


  ! Third term (ac|kd)
  ! k1 = ka, k2 = ki = kd, k3 = kk, k4 = kc
  Allocate(v_tmp(CC_n_vir,CC_n_occ,CC_n_vir,CC_mem_aa_D(CC_mpi_did+1), &
                 1,1,1),stat=errnum)
  Call check_allocation(errnum,'v_tmp in CC')

  Allocate(tmp(CC_n_occ,CC_n_occ,CC_n_vir,CC_mem_aa_D(CC_mpi_did+1), &
               1,1,1),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  i_tmp = CC_n_occ * CC_n_vir * CC_mem_aa_D(CC_mpi_did+1)

  Allocate(mat_A(CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(CC_n_vir,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(CC_n_occ,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + CC_mem_k1(CC_mpi_gid+1)

  k1 = k_start - 1
  do k_run = 1, CC_mem_k1(CC_mpi_gid+1)
    k1 = k1 + 1
    do k2 = 1, n_k_points
      do k3 = 1, n_k_points

        Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

        ! For j_aux

        !$OMP PARALLEL Default(shared) &
        !$OMP Private(ccc,ddd)
        !$OMP DO COLLAPSE(2)
        do ccc = 1, CC_n_vir
          do ddd = 1, CC_n_vir
            v_tmp(ddd,:,ccc,:,1,1,1) = 2.0D0 * CC_intl_ackd(:,ccc,k3,ddd,:,k2,k_run) &
                                             - CC_intl_ackd(:,ddd,k3,ccc,:,k4,k_run)
          end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        mat_A = CC_t_s(:,:,k2)

        shp1(1) = CC_n_vir
        shp1(2) = i_tmp
        mat_B = reshape(v_tmp,shp1)

        Call Zgemm('N','N',CC_n_occ,i_tmp,CC_n_vir,alpha,mat_A,CC_n_occ, &
                   mat_B,CC_n_vir,beta,mat_C,CC_n_occ)

        shp2i(1) = CC_n_occ
        shp2i(2) = CC_n_occ
        shp2i(3) = CC_n_vir
        shp2i(4) = CC_mem_aa_D(CC_mpi_did+1)

        tmp(:,:,:,:,1,1,1) = reshape(mat_C,shp2i)

        do iii = 1, CC_n_occ
          CC_j_aux(:,:,k3,iii,:,k2,k_run) = &
                 CC_j_aux(:,:,k3,iii,:,k2,k_run) + tmp(iii,:,:,:,1,1,1)
        end do

        ! For k_aux
        do ddd = 1, CC_n_vir
          v_tmp(ddd,:,:,:,1,1,1) = CC_intl_ackd(:,ddd,k3,:,:,k4,k_run)
        end do

        mat_A = CC_t_s(:,:,k2)

        shp1(1) = CC_n_vir
        shp1(2) = i_tmp
        mat_B = reshape(v_tmp,shp1)

        Call Zgemm('N','N',CC_n_occ,i_tmp,CC_n_vir,alpha,mat_A,CC_n_occ, &
                   mat_B,CC_n_vir,beta,mat_C,CC_n_occ)

        shp2i(1) = CC_n_occ
        shp2i(2) = CC_n_occ
        shp2i(3) = CC_n_vir
        shp2i(4) = CC_mem_aa_D(CC_mpi_did+1)

        tmp(:,:,:,:,1,1,1) = reshape(mat_C,shp2i)

        do iii = 1, CC_n_occ
          CC_k_aux(:,:,k3,iii,:,k2,k_run) = &
                 CC_k_aux(:,:,k3,iii,:,k2,k_run) + tmp(iii,:,:,:,1,1,1)
        end do

      end do
    end do
  end do
 
  Deallocate(v_tmp,mat_A,mat_B,mat_C,tmp)

  !print*,'after intl j'
  !do aaa = 1, CC_n_vir
  !  do ccc = 1, CC_n_vir
  !    print*,aaa,1,1,ccc,CC_j_aux(1,ccc,1,1,aaa,1,1)
  !  end do
  !end do
 
  End Subroutine CC_3d_aux_jk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_aux_a()

! Caluclate intermediate a
! k1 = kl, k2 = ki, k3 = kk, k4 = kj, k5 = kc, k6 = kd
 
  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k5,k6,k_pattern,k_run,n_k1,k_start,k_end,k_tmp,k_13
  Integer :: ij_run,c_start,c_end,n_c,c_run
  Integer :: n_cd,n_kl,kl_start,kl_end,kl_run,cd_run,code_kl
  Integer :: iii,jjj,aaa,bbb,kkk,lll,ccc,ddd,i_run,i_tmp,j_tmp,l_tmp,ij_tmp
  Integer :: i_task,source_id,i_grp_task,grp_id,target_id
  Integer :: errnum
  Integer , dimension(2) :: shp1
  Integer , dimension(3) :: shp1i
  Integer , dimension(4) :: shp2i

  Integer (kind=8) :: s_tmp

  Double complex , dimension(:,:,:,:,:,:) , allocatable :: aux_tmp
  Double complex , dimension(:,:,:) , allocatable :: kilc_tmp
  Double complex , dimension(:,:,:) , allocatable :: tmp
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C
  Double complex , dimension(:,:) , allocatable :: mat_D,mat_E,mat_F
  Double complex , dimension(:,:,:,:) , allocatable :: t_plus_nd,t_minus_nd
  Double complex , dimension(:,:,:,:) , allocatable :: t_plus_d,t_minus_d
  Double complex , dimension(:,:,:) , allocatable :: v_plus_nd,v_minus_nd
  Double complex , dimension(:,:,:) , allocatable :: v_plus_d,v_minus_d
  Double complex , dimension(:,:,:) , allocatable :: rlt_S,rlt_A


  Double complex :: alpha,beta,ZDOTU

!  write(100,*) 'a_aux'

  ! First term (ki|lj) 
  CC_a_aux = CC_intl_kilj

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  n_cd = CC_mem_ab_D(CC_mpi_did+1,3)
  n_c = CC_mem_aa_D(CC_mpi_did+1)

  alpha = 1.0D0
  beta = 0.0D0

  n_kl = CC_mem_ij_D(CC_mpi_did+1,2)
  kl_start = CC_index_ij_D(CC_mpi_did+1,2)
  kl_end = kl_start - 1 + n_kl

  ! Second term (ki|lc) & (li|kc)
  Allocate(kilc_tmp(n_kl,CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'kilc_tmp in CC')

  i_tmp = n_kl * CC_n_occ

  Allocate(mat_A(i_tmp,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(CC_n_vir,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(i_tmp,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  Allocate(tmp(n_kl,CC_n_occ,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  ! k1 = kl, k2 = ki, k3 = kk, k4 = kj, k5 = kc, k6 = kd
  k_13 = k_start - 1
  do k_run = 1, n_k1

    k_13 = k_13 + 1
    do k1 = 1, n_k_points
      Call CC_3d_k_minus(k_13,k1,k3)
      
      do k2 = 1, n_k_points
        Call CC_3d_k_minus(k_13,k2,k4)
     
        ! For (ki|lc)
        kilc_tmp = CC_intl_kilc_T(:,k1,:,:,k2,k_run,1)

        shp1(1) = i_tmp
        shp1(2) = CC_n_vir

        mat_A = reshape(kilc_tmp,shp1)
        do iii = 1, CC_n_occ
          mat_B(:,iii) = CC_t_s(iii,:,k4)
        end do

        Call Zgemm('N','N',i_tmp,CC_n_occ,CC_n_vir,alpha,mat_A,i_tmp, &
                   mat_B,CC_n_vir,beta,mat_C,i_tmp)

        shp1i(1) = n_kl
        shp1i(2) = CC_n_occ
        shp1i(3) = CC_n_occ

        tmp = reshape(mat_C,shp1i)

        CC_a_aux(:,k1,:,:,k2,k_run) = CC_a_aux(:,k1,:,:,k2,k_run) + tmp

        ! For (li|kc)
        kilc_tmp = CC_intl_kilc_T(:,k1,:,:,k2,k_run,2)

        shp1(1) = i_tmp
        shp1(2) = CC_n_vir

        mat_A = reshape(kilc_tmp,shp1)
        do iii = 1, CC_n_occ
          mat_B(:,iii) = CC_t_s(iii,:,k2)
        end do

        Call Zgemm('N','N',i_tmp,CC_n_occ,CC_n_vir,alpha,mat_A,i_tmp, &
                   mat_B,CC_n_vir,beta,mat_C,i_tmp)

        shp1i(1) = n_kl
        shp1i(2) = CC_n_occ
        shp1i(3) = CC_n_occ

        tmp = reshape(mat_C,shp1i)

        do iii = 1, CC_n_occ
          CC_a_aux(:,k1,iii,:,k2,k_run) = CC_a_aux(:,k1,iii,:,k2,k_run) + tmp(:,:,iii)
        end do

      end do
    end do
  end do

  Deallocate(mat_A,mat_B,mat_C,tmp,kilc_tmp)

  ! Third term (kc|ld)
  ! k1 = kl, k2 = ki, k3 = kk, k4 = kj, k5 = kc, k6 = kd
 
  do i_task = 1, CC_mpi_domain_size
  
    n_kl = CC_mem_ij_D(i_task,2)
    kl_start = CC_index_ij_D(i_task,2)
    kl_end = kl_start - 1 + n_kl

    Allocate(aux_tmp(n_kl,n_k_points,CC_n_occ,CC_n_occ,n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'aux_tmp in CC')

    aux_tmp = 0.0D0

    ij_tmp = CC_n_occ * (CC_n_occ + 1) / 2
    i_tmp = ij_tmp * n_k_points

    Allocate(t_plus_nd(n_cd,n_k_points,ij_tmp,n_k_points),stat=errnum)
    Call check_allocation(errnum,'t_plus_nd in CC')

    Allocate(t_minus_nd(n_cd,n_k_points,ij_tmp,n_k_points),stat=errnum)
    Call check_allocation(errnum,'t_minus_nd in CC')

    Allocate(v_plus_nd(n_cd,n_k_points,n_kl),stat=errnum)
    Call check_allocation(errnum,'v_plus_nd in CC')

    Allocate(v_minus_nd(n_cd,n_k_points,n_kl),stat=errnum)
    Call check_allocation(errnum,'v_minus_nd in CC')

    Allocate(t_plus_d(n_c,n_k_points,ij_tmp,n_k_points),stat=errnum)
    Call check_allocation(errnum,'t_plus_d in CC')

    Allocate(t_minus_d(n_c,n_k_points,ij_tmp,n_k_points),stat=errnum)
    Call check_allocation(errnum,'t_minus_d in CC')

    Allocate(v_plus_d(n_c,n_k_points,n_kl),stat=errnum)
    Call check_allocation(errnum,'v_plus_d in CC')

    Allocate(v_minus_d(n_c,n_k_points,n_kl),stat=errnum)
    Call check_allocation(errnum,'v_minus_d in CC')

    Allocate(rlt_S(n_kl,ij_tmp,n_k_points),stat=errnum)
    Call check_allocation(errnum,'v_minus_d in CC')

    Allocate(rlt_A(n_kl,ij_tmp,n_k_points),stat=errnum)
    Call check_allocation(errnum,'v_minus_d in CC')

    j_tmp = n_cd * n_k_points
    Allocate(mat_A(j_tmp,n_kl),stat=errnum)
    Call check_allocation(errnum,'mat_A in CC')

    Allocate(mat_B(j_tmp,i_tmp),stat=errnum)
    Call check_allocation(errnum,'mat_B in CC')

    Allocate(mat_C(n_kl,i_tmp),stat=errnum)
    Call check_allocation(errnum,'mat_C in CC')

    l_tmp = n_c * n_k_points
    Allocate(mat_D(l_tmp,n_kl),stat=errnum)
    Call check_allocation(errnum,'mat_D in CC')

    Allocate(mat_E(l_tmp,i_tmp),stat=errnum)
    Call check_allocation(errnum,'mat_E in CC')

    Allocate(mat_F(n_kl,i_tmp),stat=errnum)
    Call check_allocation(errnum,'mat_F in CC')

    k_13 = k_start - 1
    do k_run = 1, n_k1

      k_13 = k_13 + 1

      !$OMP PARALLEL Default(shared) &
      !$OMP Private(k2,k4,ij_run,iii,jjj)
      !$OMP DO COLLAPSE(2)
      do k2 = 1, n_k_points
        do ij_run = 1, ij_tmp

          Call CC_3d_k_minus(k_13,k2,k4)
          Call CC_3d_decode(ij_run,iii,jjj,CC_n_occ,2)

          t_plus_nd(:,:,ij_run,k2) = 0.5D0 * (CC_t_d_nd(iii,jjj,:,k2,:,k_run) &
                                            + CC_t_d_nd(jjj,iii,:,k4,:,k_run))

          t_minus_nd(:,:,ij_run,k2) = 0.5D0 * (CC_t_d_nd(iii,jjj,:,k2,:,k_run) &
                                             - CC_t_d_nd(jjj,iii,:,k4,:,k_run))

          t_plus_d(:,:,ij_run,k2) = 0.5D0 * (CC_t_d_d(iii,jjj,:,k2,:,k_run) &
                                           + CC_t_d_d(jjj,iii,:,k4,:,k_run))

          t_minus_d(:,:,ij_run,k2) = 0.5D0 * (CC_t_d_d(iii,jjj,:,k2,:,k_run) &
                                            - CC_t_d_d(jjj,iii,:,k4,:,k_run))

        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL 

      do k1 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k1,k3)

        code_kl = kl_start - 1
        do kl_run = 1, n_kl

          code_kl = code_kl + 1
          Call CC_3d_decode(code_kl,kkk,lll,CC_n_occ,2)

          v_plus_nd(:,:,kl_run) = 0.5D0 * (CC_intl_iajb_nd(kkk,lll,:,k1,:,k_run) &
                                         + CC_intl_iajb_nd(lll,kkk,:,k3,:,k_run))

          v_minus_nd(:,:,kl_run) = 0.5D0 * (CC_intl_iajb_nd(kkk,lll,:,k1,:,k_run) &
                                          - CC_intl_iajb_nd(lll,kkk,:,k3,:,k_run))

          v_plus_d(:,:,kl_run) = 0.5D0 * (CC_intl_iajb_d(kkk,lll,:,k1,:,k_run) &
                                        + CC_intl_iajb_d(lll,kkk,:,k3,:,k_run))

          v_minus_d(:,:,kl_run) = 0.5D0 * (CC_intl_iajb_d(kkk,lll,:,k1,:,k_run) &
                                         - CC_intl_iajb_d(lll,kkk,:,k3,:,k_run))

        end do

        shp1(1) = j_tmp
        shp1(2) = n_kl
        mat_A = reshape(v_plus_nd,shp1)

        shp1(1) = j_tmp
        shp1(2) = i_tmp
        mat_B = reshape(t_plus_nd,shp1)

        alpha = 2.0D0
        beta = 0.0D0

        Call Zgemm('T','N',n_kl,i_tmp,j_tmp,alpha,mat_A,j_tmp, &
                   mat_B,j_tmp,beta,mat_C,n_kl)

        shp1i(1) = n_kl
        shp1i(2) = ij_tmp
        shp1i(3) = n_k_points

        rlt_S = reshape(mat_C,shp1i)

        shp1(1) = j_tmp
        shp1(2) = n_kl
        mat_A = reshape(v_minus_nd,shp1)

        shp1(1) = j_tmp
        shp1(2) = i_tmp
        mat_B = reshape(t_minus_nd,shp1)

        alpha = 2.0D0
        beta = 0.0D0

        Call Zgemm('T','N',n_kl,i_tmp,j_tmp,alpha,mat_A,j_tmp, &
                   mat_B,j_tmp,beta,mat_C,n_kl)

        shp1i(1) = n_kl
        shp1i(2) = ij_tmp
        shp1i(3) = n_k_points

        rlt_A = reshape(mat_C,shp1i)

        shp1(1) = l_tmp
        shp1(2) = n_kl
        mat_D = reshape(v_plus_d,shp1)

        shp1(1) = l_tmp
        shp1(2) = i_tmp
        mat_E = reshape(t_plus_d,shp1)

        alpha = 1.0D0
        beta = 0.0D0

        Call Zgemm('T','N',n_kl,i_tmp,l_tmp,alpha,mat_D,l_tmp, &
                   mat_E,l_tmp,beta,mat_F,n_kl)

        shp1i(1) = n_kl
        shp1i(2) = ij_tmp
        shp1i(3) = n_k_points

        rlt_S = rlt_S + reshape(mat_F,shp1i)

        shp1(1) = l_tmp
        shp1(2) = n_kl
        mat_D = reshape(v_minus_d,shp1)

        shp1(1) = l_tmp
        shp1(2) = i_tmp
        mat_E = reshape(t_minus_d,shp1)

        alpha = 1.0D0
        beta = 0.0D0

        Call Zgemm('T','N',n_kl,i_tmp,l_tmp,alpha,mat_D,l_tmp, &
                   mat_E,l_tmp,beta,mat_F,n_kl)

        shp1i(1) = n_kl
        shp1i(2) = ij_tmp
        shp1i(3) = n_k_points

        rlt_A = rlt_A + reshape(mat_F,shp1i)

        !$OMP PARALLEL Default(Shared) &
        !$OMP Private(k2,k4,ij_run,iii,jjj)
        !$OMP DO COLLAPSE(2)
        do k2 = 1, n_k_points
          do ij_run = 1, ij_tmp
            Call CC_3d_k_minus(k_13,k2,k4)
            Call CC_3d_decode(ij_run,iii,jjj,CC_n_occ,2)
            aux_tmp(:,k1,iii,jjj,k2,k_run) = rlt_S(:,ij_run,k2) + rlt_A(:,ij_run,k2)
            aux_tmp(:,k1,jjj,iii,k4,k_run) = rlt_S(:,ij_run,k2) - rlt_A(:,ij_run,k2)
          end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
      end do
    end do

    s_tmp = Int(CC_n_occ**2,8) * Int(n_kl,8) * Int(n_k1,8) * Int(n_k_points**2,8)
    target_id = i_task - 1

    Call CC_mpi_complex_reduce(s_tmp,aux_tmp,target_id,CC_mpi_comm_domain)

    if (target_id.eq.CC_mpi_did) then
      !$OMP PARALLEL Default(Shared)
      !$OMP WORKSHARE
      CC_a_aux = CC_a_aux + aux_tmp
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    end if

    Deallocate(t_plus_nd,t_minus_nd,t_plus_d,t_minus_d)
    Deallocate(v_plus_nd,v_minus_nd,v_plus_d,v_minus_d)
    Deallocate(mat_A,mat_B,mat_C,mat_D,mat_E,mat_F,aux_tmp,rlt_S,rlt_A)
  end do

  End Subroutine CC_3d_aux_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_add_intl_klcd_to_jk()

! Add integrals (kc|ld) to aux j and k
 
  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k5,k6,k_pattern,k1_run,n_k1,n_k1_task,k1_start,k1_end
  Integer :: k4_start,k4_end,k4_run
  Integer :: iii,aaa,jjj,bbb,kkk,lll,ccc,ddd,i_run,i_tmp,j_tmp,c_run,l_tmp
  Integer :: n_dom,n_dom_task,a_start,a_end,i_start,i_end,a_run,c_start,c_end
  Integer :: i_task,source_id,i_grp_task,grp_id,errnum
  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp1i
  Integer , dimension(3) :: shp2
  Integer , dimension(4) :: shp2i

  Integer (kind=8) :: s_tmp

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: t_tmp,v_tmp,v_tmp2
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: tmp,v_work,t_work
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C

  Double complex :: alpha,beta


  !print*,'k before kcld'
  !do iii = 1, CC_n_occ
  !  do aaa = 1, CC_n_vir
  !    do kkk = 1, CC_n_occ
  !      do ccc = 1, CC_n_vir
  !        print*,iii,aaa,kkk,ccc,CC_k_aux(kkk,ccc,1,iii,aaa,1,1)
  !      end do
  !    end do
  !  end do
  !end do

! k1 = ka, k2 = ki, k3 = kk, k4 = kc, k5 = kl, k6 = kd
  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k1_start = CC_index_k1(CC_mpi_gid+1)
  k1_end = k1_start - 1 + n_k1

  n_dom = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_dom

  Allocate(t_tmp(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_dom, &
                 n_k_points,n_k1),stat=errnum)
  Call check_allocation(errnum,'t_tmp in CC')

  Call CC_3d_get_t_grp('D',CC_mpi_gid+1,t_tmp,.true.)

  do i_grp_task = 1, CC_mpi_group_size

    n_k1_task = CC_mem_k1(i_grp_task)
    k4_start = CC_index_k1(i_grp_task)
    k4_end = k4_start - 1 + n_k1_task

    Allocate(v_tmp(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_dom, &
                   n_k_points,n_k1_task),stat=errnum)
    Call check_allocation(errnum,'v_tmp in CC')

    Call CC_3d_get_intl_kcld('D',i_grp_task,v_tmp,.false.)

    do i_task = 1, CC_mpi_domain_size

      n_dom_task = CC_mem_aa_D(i_task)
      c_start = CC_index_aa_D(i_task)
      c_end = c_start - 1 + n_dom_task

      Allocate(v_tmp2(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_dom_task, &
                      n_k_points,n_k1_task),stat=errnum)
      Call check_allocation(errnum,'v_tmp2 in CC')

      source_id = i_task - 1

      i_tmp = CC_n_occ * CC_n_vir * n_k_points
      j_tmp = CC_n_occ * n_dom_task * n_k_points * n_k1_task

      if (CC_mpi_did.eq.source_id) then
        v_tmp2 = v_tmp
      else
        v_tmp2 = 0.0D0
      end if

      s_tmp = Int(i_tmp,8) * Int(j_tmp,8)
      Call CC_mpi_complex_bcast(s_tmp, v_tmp2, source_id, CC_mpi_comm_domain)

      ! add to k
      ! k1 = ka, k2 = ki, k3 = kk, k4 = kc, k5 = kl, k6 = kd
      alpha = 0.5D0
      beta = 0.0D0

      ! v_work
      Allocate(v_work(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_dom_task,1,1),stat=errnum)
      Call check_allocation(errnum,'v_work in CC')

      ! t_work
      Allocate(t_work(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_dom,1,1),stat=errnum)
      Call check_allocation(errnum,'t_work in CC')

      i_tmp = CC_n_occ * CC_n_vir * n_k_points
      j_tmp = CC_n_occ * n_dom_task
      l_tmp = CC_n_occ * n_dom

      Allocate(mat_A(i_tmp,j_tmp),stat=errnum)
      Call check_allocation(errnum,'mat_A in CC')

      Allocate(mat_B(i_tmp,l_tmp),stat=errnum)
      Call check_allocation(errnum,'mat_B in CC')

      Allocate(mat_C(j_tmp,l_tmp),stat=errnum)
      Call check_allocation(errnum,'mat_C in CC')

      Allocate(tmp(CC_n_occ,n_dom_task,CC_n_occ,n_dom,1,1,1),stat=errnum)
      Call check_allocation(errnum,'tmp in CC')

      k1 = k1_start - 1
      do k1_run = 1, n_k1
        k1 = k1 + 1

        do k2 = 1, n_k_points

          k4 = k4_start - 1
          do k4_run = 1, n_k1_task

            k4 = k4 + 1

            Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,3)

            v_work(:,:,:,:,:,1,1) = v_tmp2(:,:,:,:,:,k3,k4_run)
            t_work(:,:,:,:,:,1,1) = t_tmp(:,:,:,:,:,k2,k1_run)

            ! For kl = ka
            !$OMP PARALLEL Default(shared) &
            !$OMP Private(iii,lll,ddd,aaa,a_run)
            !$OMP DO COLLAPSE(4)
            do iii = 1, CC_n_occ
              do ddd = 1, CC_n_vir
                do lll = 1, CC_n_occ
                  do a_run = 1, n_dom
                    aaa = a_start - 1 + a_run
                    t_work(lll,ddd,k1,iii,a_run,1,1) = t_work(lll,ddd,k1,iii,a_run,1,1) &
                                      + 2.0D0 * CC_t_s(iii,ddd,k2) * CC_t_s(lll,aaa,k1)
                  end do
                end do
              end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL
          
            !i_tmp = CC_n_occ * CC_n_vir * n_k_points
            !j_tmp = CC_n_occ * n_dom_task
            !l_tmp = CC_n_occ * n_dom

            shp1(1) = i_tmp
            shp1(2) = j_tmp
            mat_A = reshape(v_work,shp1)

            shp1(1) = i_tmp
            shp1(2) = l_tmp
            mat_B = reshape(t_work,shp1)

            Call Zgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,mat_A,i_tmp, &
                       mat_B,i_tmp,beta,mat_C,j_tmp)

            shp1i(1) = CC_n_occ
            shp1i(2) = n_dom_task
            shp1i(3) = CC_n_occ
            shp1i(4) = n_dom

            tmp(:,:,:,:,1,1,1) = reshape(mat_C,shp1i)

            CC_k_aux(:,c_start:c_end,k3,:,:,k2,k1_run) = &
                   CC_k_aux(:,c_start:c_end,k3,:,:,k2,k1_run) - tmp(:,:,:,:,1,1,1)

          end do

        end do
      end do

      ! add to j
      ! k1 = ka, k2 = ki, k3 = kk, k4 = kc, k5 = kl, k6 = kd
      alpha = 0.5D0
      beta = 0.0D0

      k1 = k1_start - 1
      do k1_run = 1, n_k1
        k1 = k1 + 1

        do k2 = 1, n_k_points

          k4 = k4_start - 1
          do k4_run = 1, n_k1_task

            k4 = k4 + 1

            Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,3)

            !$OMP PARALLEL Default(shared) &
            !$OMP Private(kkk,lll,iii,k5)
            do k5 = 1, n_k_points

              !$OMP DO COLLAPSE(2)
              do kkk = 1, CC_n_occ
                do lll = 1, CC_n_occ
                  v_work(lll,:,k5,kkk,:,1,1) = 2.0D0 * v_tmp2(kkk,:,k3,lll,:,k5,k4_run) &
                                                     - v_tmp2(lll,:,k5,kkk,:,k3,k4_run)
                end do
              end do
              !$OMP END DO

              !$OMP DO COLLAPSE(2)
              do iii = 1, CC_n_occ
                do lll = 1, CC_n_occ
                  t_work(lll,:,k5,iii,:,1,1) = 2.0D0 * t_tmp(iii,:,k2,lll,:,k5,k1_run) &
                                                     - t_tmp(lll,:,k5,iii,:,k2,k1_run)
                end do
              end do
              !$OMP END DO

            end do
            !$OMP END PARALLEL

            ! For kl = ka
            !$OMP PARALLEL Default(shared) &
            !$OMP Private(iii,lll,ddd,aaa,a_run)
            !$OMP DO COLLAPSE(4)
            do iii = 1, CC_n_occ
              do ddd = 1, CC_n_vir
                do lll = 1, CC_n_occ
                  do a_run = 1, n_dom
                    aaa = a_start - 1 + a_run
                    t_work(lll,ddd,k1,iii,a_run,1,1) = t_work(lll,ddd,k1,iii,a_run,1,1) &
                                     - 2.0D0 * CC_t_s(iii,ddd,k2) * CC_t_s(lll,aaa,k1)
                  end do
                end do
              end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL
          
            !i_tmp = CC_n_occ * CC_n_vir * n_k_points
            !j_tmp = CC_n_occ * n_dom_task
            !l_tmp = CC_n_occ * n_dom

            shp1(1) = i_tmp
            shp1(2) = j_tmp
            mat_A = reshape(v_work(:,:,:,:,:,1,1),shp1)

            shp1(1) = i_tmp
            shp1(2) = l_tmp
            mat_B = reshape(t_work(:,:,:,:,:,1,1),shp1)

            Call Zgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,mat_A,i_tmp, &
                       mat_B,i_tmp,beta,mat_C,j_tmp)
 
            shp1i(1) = CC_n_occ
            shp1i(2) = n_dom_task
            shp1i(3) = CC_n_occ
            shp1i(4) = n_dom

            tmp(:,:,:,:,1,1,1) = reshape(mat_C,shp1i)

            CC_j_aux(:,c_start:c_end,k3,:,:,k2,k1_run) = &
                   CC_j_aux(:,c_start:c_end,k3,:,:,k2,k1_run) + tmp(:,:,:,:,1,1,1)

          end do

        end do
      end do

      Deallocate(v_tmp2,t_work,v_work,mat_A,mat_B,mat_C,tmp)

    end do

    Deallocate(v_tmp)

  end do

  !print*,'k after kcld'
  !do iii = 1, CC_n_occ
  !  do aaa = 1, CC_n_vir
  !    do kkk = 1, CC_n_occ
  !      do ccc = 1, CC_n_vir
  !        print*,aaa,1,1,ccc,CC_j_aux(kkk,ccc,1,iii,aaa,1,1)
  !      end do
  !    end do
  !  end do
  !end do



  End Subroutine CC_3d_add_intl_klcd_to_jk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_aux_b_nd(i_task,k1,k3,k_13,CC_b_plus_nd,CC_b_minus_nd, &
                                              CC_b_plus_d,CC_b_minus_d)

  Use CC_3d

  Implicit None

  Integer , intent(in) :: i_task,k1,k3,k_13

  Integer :: k2,k4,k5,k6,k_run,n_k1,k_start,k_end,k_pattern

  Integer :: aaa,bbb,ccc,ddd,kkk,a_run,c_run,code_ab,code_cd
  Integer :: n_ab,ab_start,ab_end,n_cd,cd_start,cd_end,n_a,a_start,a_end
  Integer :: nnn,i_tmp,j_tmp
  Integer :: i_task_run,target_id,errnum

  Integer (kind=8) :: s_tmp

  Double complex , dimension(CC_mem_ab_D(CC_mpi_did+1,3),n_k_points, &
                             CC_mem_ab_D(i_task,3)) :: CC_b_plus_nd,CC_b_minus_nd

  Double complex , dimension(CC_mem_aa_D(CC_mpi_did+1),n_k_points,   &
                             CC_mem_ab_D(i_task,3)) :: CC_b_plus_d,CC_b_minus_d

  Double complex , dimension(:,:,:,:) , allocatable :: t_work,v_work
  Double complex , dimension(:,:,:,:) , allocatable :: tmp,tmp2
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C
  Double complex , dimension(:,:,:,:) , allocatable :: intl_rlt

  Double complex :: alpha,beta,ZDOTU,rlt1,rlt2

  Integer :: OMP_GET_NUM_THREADS,o_num

  ! k1 = ka, k2 = kc, k3 = kb, k4 = kd
  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)

  n_ab = CC_mem_ab_D(i_task,3)
  ab_start = CC_index_ab_D(i_task,3)
  ab_end = ab_start - 1 + n_ab

  ! First term (ac|bd)
  if (CC_abcd_sv_flag) then

    code_ab = ab_start - 1
    do c_run = 1, n_ab

      code_ab = code_ab + 1
      Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)

      CC_b_plus_nd(:,:,c_run) = 0.5D0 * (CC_intl_acbd_nd(:,:,aaa,bbb,k1,k_13) &
                                       + CC_intl_acbd_nd(:,:,bbb,aaa,k3,k_13))

      CC_b_minus_nd(:,:,c_run) = 0.5D0 * (CC_intl_acbd_nd(:,:,aaa,bbb,k1,k_13) &
                                        - CC_intl_acbd_nd(:,:,bbb,aaa,k3,k_13))

      CC_b_plus_d(:,:,c_run) = 0.5D0 * (CC_intl_acbd_d(:,:,aaa,bbb,k1,k_13) &
                                      + CC_intl_acbd_d(:,:,bbb,aaa,k3,k_13))

      CC_b_minus_d(:,:,c_run) = 0.5D0 * (CC_intl_acbd_d(:,:,aaa,bbb,k1,k_13) &
                                       - CC_intl_acbd_d(:,:,bbb,aaa,k3,k_13))

    end do

  else

    nnn = CC_mem_bas(CC_mpi_did+1)

    do i_task_run = 1, CC_mpi_domain_size

      n_cd = CC_mem_ab_D(i_task_run,3)
      cd_start = CC_index_ab_D(i_task_run,3)
      cd_end = cd_start - 1 + n_cd

      n_a = CC_mem_aa_D(i_task_run)
      a_start = CC_index_aa_D(i_task_run)
      a_end = a_start - 1 + n_a

      Allocate(intl_rlt(n_cd,n_k_points,n_ab,2),stat=errnum)
      Call check_allocation(errnum,'intl_rlt in CC')

      !$ o_num = OMP_GET_NUM_THREADS()
      !!$ Call mkl_set_num_threads(1)

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(k2,c_run,a_run,k4,code_cd,code_ab,ccc,ddd,aaa,bbb)
      !$OMP DO COLLAPSE(3)
      do k2 = 1, n_k_points
        do c_run = 1, n_cd
          do a_run = 1, n_ab

            Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

            code_cd = cd_start - 1 + c_run
            code_ab = ab_start - 1 + a_run
            Call CC_3d_decode(code_cd,ccc,ddd,CC_n_vir,3)
            Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)

            aaa = aaa + CC_n_occ
            bbb = bbb + CC_n_occ
            ccc = ccc + CC_n_occ
            ddd = ddd + CC_n_occ

            intl_rlt(c_run,k2,a_run,1) = ZDOTU(nnn,CC_RI_L(:,aaa,ccc,k1,k2),1, &
                                                   CC_RI_R(:,bbb,ddd,k3,k4),1)

            intl_rlt(c_run,k2,a_run,2) = ZDOTU(nnn,CC_RI_L(:,bbb,ccc,k3,k2),1, &
                                                   CC_RI_R(:,aaa,ddd,k1,k4),1)

          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      !!$ Call mkl_set_num_threads(o_num)

      s_tmp = Int(n_cd,8) * Int(n_ab*2,8) * Int(n_k_points,8)
      target_id = i_task_run - 1

      Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

      if (CC_mpi_did.eq.target_id) then

        CC_b_plus_nd = 0.5D0 * (intl_rlt(:,:,:,1) + intl_rlt(:,:,:,2))
        CC_b_minus_nd = 0.5D0 * (intl_rlt(:,:,:,1) - intl_rlt(:,:,:,2))

      end if

      Deallocate(intl_rlt)

      Allocate(intl_rlt(n_a,n_k_points,n_ab,2),stat=errnum)
      Call check_allocation(errnum,'intl_rlt in CC')

      !!$ Call mkl_set_num_threads(1)

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(k2,c_run,a_run,k4,ccc,code_ab,aaa,bbb)
      !$OMP DO COLLAPSE(3)
      do k2 = 1, n_k_points
        do c_run = 1, n_a
          do a_run = 1, n_ab

            Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

            ccc = a_start - 1 + c_run
            code_ab = ab_start - 1 + a_run
            Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)

            aaa = aaa + CC_n_occ
            bbb = bbb + CC_n_occ
            ccc = ccc + CC_n_occ

            intl_rlt(c_run,k2,a_run,1) = ZDOTU(nnn,CC_RI_L(:,aaa,ccc,k1,k2),1, &
                                                   CC_RI_R(:,bbb,ccc,k3,k4),1)

            intl_rlt(c_run,k2,a_run,2) = ZDOTU(nnn,CC_RI_L(:,bbb,ccc,k3,k2),1, &
                                                   CC_RI_R(:,aaa,ccc,k1,k4),1)

          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      !!$ Call mkl_set_num_threads(o_num)

      s_tmp = Int(n_a,8) * Int(n_ab*2,8) * Int(n_k_points,8)
      target_id = i_task_run - 1

      Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

      if (CC_mpi_did.eq.target_id) then

        CC_b_plus_d = 0.5D0 * (intl_rlt(:,:,:,1) + intl_rlt(:,:,:,2))
        CC_b_minus_d = 0.5D0 * (intl_rlt(:,:,:,1) - intl_rlt(:,:,:,2))

      end if

      Deallocate(intl_rlt)

    end do
   
  end if

  alpha = 1.0D0
  beta = 0.0D0

  ! Second term (ac|kd) & (kc|bd)
  ! k1 = ka, k2 = kc, k3 = kb, k4 = kd
  n_cd = CC_mem_ab_D(CC_mpi_did+1,3)
  cd_start = CC_index_ab_D(CC_mpi_did+1,3)
  cd_end = cd_start - 1 + n_cd

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  !!$ Call mkl_set_num_threads(1)

  !$OMP PARALLEL Default(Shared) &
  !$OMP Private(k2,c_run,a_run,k4,code_ab,aaa,bbb,rlt1,rlt2)

  !$OMP DO COLLAPSE(3)
  do k2 = 1, n_k_points
    do c_run = 1, n_cd
      do a_run = 1, n_ab

        code_ab = ab_start - 1 + a_run
        Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)

        rlt1 = 0.0D0
        rlt2 = 0.0D0

        ! (ac|kd) 
        rlt1 = rlt1 - ZDOTU(CC_n_occ,CC_intl_ackd_nd(:,c_run,aaa,k2,k1,k_13,1),1, &
                                                              CC_t_s(:,bbb,k3),1)

        rlt2 = rlt2 - ZDOTU(CC_n_occ,CC_intl_ackd_nd(:,c_run,bbb,k2,k3,k_13,1),1, &
                                                              CC_t_s(:,aaa,k1),1)
 
        ! (kc|bd)
        rlt1 = rlt1 - ZDOTU(CC_n_occ,CC_intl_ackd_nd(:,c_run,bbb,k2,k1,k_13,2),1, &
                                                              CC_t_s(:,aaa,k1),1)

        rlt2 = rlt2 - ZDOTU(CC_n_occ,CC_intl_ackd_nd(:,c_run,aaa,k2,k3,k_13,2),1, &
                                                              CC_t_s(:,bbb,k3),1)

        CC_b_plus_nd(c_run,k2,a_run) = CC_b_plus_nd(c_run,k2,a_run) &
                                       + 0.5D0 * (rlt1 + rlt2)
        CC_b_minus_nd(c_run,k2,a_run) = CC_b_minus_nd(c_run,k2,a_run) &
                                        + 0.5D0 * (rlt1 - rlt2)
 
      end do
    end do
  end do
  !$OMP END DO

  !$OMP DO COLLAPSE(3)
  do k2 = 1, n_k_points
    do c_run = 1, n_a
      do a_run = 1, n_ab

        Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

        code_ab = ab_start - 1 + a_run
        Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)

        rlt1 = 0.0D0
        rlt2 = 0.0D0

        ! (ac|kd) 
        rlt1 = rlt1 - ZDOTU(CC_n_occ,CC_intl_ackd_d(:,c_run,aaa,k2,k1,k_13),1, &
                                                           CC_t_s(:,bbb,k3),1)

        rlt2 = rlt2 - ZDOTU(CC_n_occ,CC_intl_ackd_d(:,c_run,bbb,k2,k3,k_13),1, &
                                                           CC_t_s(:,aaa,k1),1)
 
        ! (kc|bd)
        rlt1 = rlt1 - ZDOTU(CC_n_occ,CC_intl_ackd_d(:,c_run,bbb,k4,k3,k_13),1, &
                                                          CC_t_s(:,aaa,k1),1)

        rlt2 = rlt2 - ZDOTU(CC_n_occ,CC_intl_ackd_d(:,c_run,aaa,k4,k1,k_13),1, &
                                                           CC_t_s(:,bbb,k3),1)

        CC_b_plus_d(c_run,k2,a_run) = CC_b_plus_d(c_run,k2,a_run) &
                                      + 0.5D0 * (rlt1 + rlt2)
        CC_b_minus_d(c_run,k2,a_run) = CC_b_minus_d(c_run,k2,a_run) &
                                       + 0.5D0 * (rlt1 - rlt2)
 
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !!$ Call mkl_set_num_threads(o_num)

  !write(80+myid,*) 'b_aux'
  !do k2 = 1, n_k_points
  !  do c_run = 1, n_cd
  !    do a_run = 1, n_ab

  !      code_ab = ab_start - 1 + a_run
  !      Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)

  !      Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

  !      Call CC_3d_decode(c_run,ccc,ddd,CC_n_vir,3)

  !      write(80+myid,*) aaa,bbb,ccc,ddd
  !      write(80+myid,*) CC_b_plus_nd(c_run,k2,a_run)
  !      write(80+myid,*) CC_b_minus_nd(c_run,k2,a_run)
  !    end do
  !  end do
  !end do
 
  !do k2 = 1, n_k_points
  !  do c_run = 1, n_a
  !    do a_run = 1, n_ab

  !      code_ab = ab_start - 1 + a_run
  !      Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)

  !      Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

  !      ccc = a_start - 1 + c_run
  !      ddd = ccc

  !      write(80+myid,*) aaa,bbb,ccc,ddd
  !      write(80+myid,*) CC_b_plus_d(c_run,k2,a_run)
  !      write(80+myid,*) CC_b_minus_d(c_run,k2,a_run)
  !    end do
  !  end do
  !end do
 
  End Subroutine CC_3d_aux_b_nd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_aux_b_d(i_task,k1,k3,k_13,CC_b_plus_nd,CC_b_minus_nd, &
                                             CC_b_plus_d,CC_b_minus_d)

  Use CC_3d

  Implicit None

  Integer , intent(in) :: i_task,k1,k3,k_13

  Integer :: k2,k4,k5,k6,k_run,n_k1,k_start,k_end,k_pattern

  Integer :: aaa,bbb,ccc,ddd,kkk,a_run,c_run,code_ab,code_cd
  Integer :: n_cd,cd_start,cd_end
  Integer :: n_a,a_start,a_end,n_c,c_start,c_end
  Integer :: nnn,i_tmp,j_tmp
  Integer :: i_task_run,target_id,errnum

  Integer (kind=8) :: s_tmp

  Double complex , dimension(CC_mem_ab_D(CC_mpi_did+1,3),n_k_points, &
                             CC_mem_aa_D(i_task)) :: CC_b_plus_nd,CC_b_minus_nd

  Double complex , dimension(CC_mem_aa_D(CC_mpi_did+1),n_k_points,   &
                             CC_mem_aa_D(i_task)) :: CC_b_plus_d,CC_b_minus_d

  Double complex , dimension(:,:,:,:) , allocatable :: t_work,v_work
  Double complex , dimension(:,:,:,:) , allocatable :: tmp,tmp2
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C
  Double complex , dimension(:,:,:,:) , allocatable :: intl_rlt

  Double complex :: alpha,beta,ZDOTU,rlt1,rlt2

  Integer :: OMP_GET_NUM_THREADS,o_num

  ! k1 = ka, k2 = kc, k3 = kb, k4 = kd
  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)

  n_c = CC_mem_aa_D(i_task)
  c_start = CC_index_aa_D(i_task)
  c_end = c_start - 1 + n_c

  ! First term (ac|bd)
  if (CC_abcd_sv_flag) then

    aaa = c_start - 1
    do c_run = 1, n_c

      aaa = aaa + 1
      bbb = aaa

      CC_b_plus_nd(:,:,c_run) = 0.5D0 * (CC_intl_acbd_nd(:,:,aaa,bbb,k1,k_13) &
                                         + CC_intl_acbd_nd(:,:,bbb,aaa,k3,k_13))

      CC_b_minus_nd(:,:,c_run) = 0.5D0 * (CC_intl_acbd_nd(:,:,aaa,bbb,k1,k_13) &
                                          - CC_intl_acbd_nd(:,:,bbb,aaa,k3,k_13))

      CC_b_plus_d(:,:,c_run) = 0.5D0 * (CC_intl_acbd_d(:,:,aaa,bbb,k1,k_13) &
                                        + CC_intl_acbd_d(:,:,bbb,aaa,k3,k_13))

      CC_b_minus_d(:,:,c_run) = 0.5D0 * (CC_intl_acbd_d(:,:,aaa,bbb,k1,k_13) &
                                         - CC_intl_acbd_d(:,:,bbb,aaa,k3,k_13))

    end do

  else

    nnn = CC_mem_bas(CC_mpi_did+1)

    do i_task_run = 1, CC_mpi_domain_size

      n_cd = CC_mem_ab_D(i_task_run,3)
      cd_start = CC_index_ab_D(i_task_run,3)
      cd_end = cd_start - 1 + n_cd

      n_a = CC_mem_aa_D(i_task_run)
      a_start = CC_index_aa_D(i_task_run)
      a_end = a_start - 1 + n_a

      Allocate(intl_rlt(n_cd,n_k_points,n_c,2),stat=errnum)
      Call check_allocation(errnum,'intl_rlt in CC')

      !$ o_num = OMP_GET_NUM_THREADS()
      !!$ Call mkl_set_num_threads(1)

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(k2,c_run,a_run,k4,code_cd,ccc,ddd,aaa,bbb)
      !$OMP DO COLLAPSE(3)
      do k2 = 1, n_k_points
        do c_run = 1, n_cd
          do a_run = 1, n_c

            Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

            code_cd = cd_start - 1 + c_run
            aaa = c_start - 1 + a_run
            Call CC_3d_decode(code_cd,ccc,ddd,CC_n_vir,3)

            aaa = aaa + CC_n_occ
            bbb = aaa
            ccc = ccc + CC_n_occ
            ddd = ddd + CC_n_occ

            intl_rlt(c_run,k2,a_run,1) = ZDOTU(nnn,CC_RI_L(:,aaa,ccc,k1,k2),1, &
                                                   CC_RI_R(:,bbb,ddd,k3,k4),1)

            intl_rlt(c_run,k2,a_run,2) = ZDOTU(nnn,CC_RI_L(:,bbb,ccc,k3,k2),1, &
                                                   CC_RI_R(:,aaa,ddd,k1,k4),1)

          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      !!$ Call mkl_set_num_threads(o_num)

      s_tmp = Int(n_cd,8) * Int(n_c*2,8) * Int(n_k_points,8)
      target_id = i_task_run - 1

      Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

      if (CC_mpi_did.eq.target_id) then

        CC_b_plus_nd = 0.5D0 * (intl_rlt(:,:,:,1) + intl_rlt(:,:,:,2))
        CC_b_minus_nd = 0.5D0 * (intl_rlt(:,:,:,1) - intl_rlt(:,:,:,2))

      end if

      Deallocate(intl_rlt)

      Allocate(intl_rlt(n_a,n_k_points,n_c,2),stat=errnum)
      Call check_allocation(errnum,'intl_rlt in CC')

      !!$ Call mkl_set_num_threads(1)

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(k2,c_run,a_run,k4,ccc,aaa,bbb)
      !$OMP DO COLLAPSE(3)
      do k2 = 1, n_k_points
        do c_run = 1, n_a
          do a_run = 1, n_c

            Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

            ccc = a_start - 1 + c_run + CC_n_occ
            aaa = c_start - 1 + a_run + CC_n_occ
            bbb = aaa

            intl_rlt(c_run,k2,a_run,1) = ZDOTU(nnn,CC_RI_L(:,aaa,ccc,k1,k2),1, &
                                                   CC_RI_R(:,bbb,ccc,k3,k4),1)

            intl_rlt(c_run,k2,a_run,2) = ZDOTU(nnn,CC_RI_L(:,bbb,ccc,k3,k2),1, &
                                                   CC_RI_R(:,aaa,ccc,k1,k4),1)

          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      !!$ Call mkl_set_num_threads(o_num)

      s_tmp = Int(n_a,8) * Int(n_c*2,8) * Int(n_k_points,8)
      target_id = i_task_run - 1

      Call CC_mpi_complex_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

      if (CC_mpi_did.eq.target_id) then

        CC_b_plus_d = 0.5D0 * (intl_rlt(:,:,:,1) + intl_rlt(:,:,:,2))
        CC_b_minus_d = 0.5D0 * (intl_rlt(:,:,:,1) - intl_rlt(:,:,:,2))

      end if

      Deallocate(intl_rlt)

    end do
   
  end if

  alpha = 1.0D0
  beta = 0.0D0

  ! Second term (ac|kd) & (kc|bd)
  ! k1 = ka, k2 = kc, k3 = kb, k4 = kd
  n_cd = CC_mem_ab_D(CC_mpi_did+1,3)
  cd_start = CC_index_ab_D(CC_mpi_did+1,3)
  cd_end = cd_start - 1 + n_cd

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  !!$ Call mkl_set_num_threads(1)

  !$OMP PARALLEL Default(Shared) &
  !$OMP Private(k2,k4,c_run,a_run,aaa,bbb,rlt1,rlt2)

  !$OMP DO COLLAPSE(3)
  do k2 = 1, n_k_points
    do c_run = 1, n_cd
      do a_run = 1, n_c

        aaa = c_start - 1 + a_run
        bbb = aaa

        rlt1 = 0.0D0
        rlt2 = 0.0D0

        ! (ac|kd) 
        rlt1 = rlt1 - ZDOTU(CC_n_occ,CC_intl_ackd_nd(:,c_run,aaa,k2,k1,k_13,1),1, &
                                                              CC_t_s(:,bbb,k3),1)

        rlt2 = rlt2 - ZDOTU(CC_n_occ,CC_intl_ackd_nd(:,c_run,bbb,k2,k3,k_13,1),1, &
                                                              CC_t_s(:,aaa,k1),1)
 
        ! (kc|bd)
        rlt1 = rlt1 - ZDOTU(CC_n_occ,CC_intl_ackd_nd(:,c_run,bbb,k2,k1,k_13,2),1, &
                                                              CC_t_s(:,aaa,k1),1)

        rlt2 = rlt2 - ZDOTU(CC_n_occ,CC_intl_ackd_nd(:,c_run,aaa,k2,k3,k_13,2),1, &
                                                              CC_t_s(:,bbb,k3),1)

        CC_b_plus_nd(c_run,k2,a_run) = CC_b_plus_nd(c_run,k2,a_run) &
                                       + 0.5D0 * (rlt1 + rlt2)
        CC_b_minus_nd(c_run,k2,a_run) = CC_b_minus_nd(c_run,k2,a_run) &
                                        + 0.5D0 * (rlt1 - rlt2)
 
      end do
    end do
  end do
  !$OMP END DO

  !$OMP DO COLLAPSE(3)
  do k2 = 1, n_k_points
    do c_run = 1, n_a
      do a_run = 1, n_c

        Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

        aaa = c_start - 1 + a_run
        bbb = aaa

        rlt1 = 0.0D0
        rlt2 = 0.0D0

        ! (ac|kd) 
        rlt1 = rlt1 - ZDOTU(CC_n_occ,CC_intl_ackd_d(:,c_run,aaa,k2,k1,k_13),1, &
                                                           CC_t_s(:,bbb,k3),1)

        rlt2 = rlt2 - ZDOTU(CC_n_occ,CC_intl_ackd_d(:,c_run,bbb,k2,k3,k_13),1, &
                                                           CC_t_s(:,aaa,k1),1)
 
        ! (kc|bd)
        rlt1 = rlt1 - ZDOTU(CC_n_occ,CC_intl_ackd_d(:,c_run,bbb,k4,k3,k_13),1, &
                                                          CC_t_s(:,aaa,k1),1)

        rlt2 = rlt2 - ZDOTU(CC_n_occ,CC_intl_ackd_d(:,c_run,aaa,k4,k1,k_13),1, &
                                                           CC_t_s(:,bbb,k3),1)

        CC_b_plus_d(c_run,k2,a_run) = CC_b_plus_d(c_run,k2,a_run) &
                                      + 0.5D0 * (rlt1 + rlt2)
        CC_b_minus_d(c_run,k2,a_run) = CC_b_minus_d(c_run,k2,a_run) &
                                       + 0.5D0 * (rlt1 - rlt2)
 
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !!$ Call mkl_set_num_threads(o_num)

  !do k2 = 1, n_k_points
  !  do c_run = 1, n_cd
  !    do a_run = 1, n_c

  !      Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

  !      Call CC_3d_decode(c_run,ccc,ddd,CC_n_vir,3)
  !      aaa = c_start - 1 + a_run
  !      bbb = aaa

  !      print*,aaa,bbb,ccc,ddd,CC_b_plus_nd(c_run,k2,a_run),CC_b_minus_nd(c_run,k2,a_run)
  !    end do
  !  end do
  !end do
 
  !do k2 = 1, n_k_points
  !  do c_run = 1, n_a
  !    do a_run = 1, n_c

  !      Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

  !      aaa = c_start - 1 + a_run
  !      bbb = aaa
  !      ccc = a_start - 1 + c_run
  !      ddd = ccc

  !      print*,aaa,bbb,ccc,ddd,CC_b_plus_d(c_run,k2,a_run),CC_b_minus_d(c_run,k2,a_run)
  !    end do
  !  end do
  !end do
 
  End Subroutine CC_3d_aux_b_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

