!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_calc_w()

  Use CC_3d

  Implicit None

  Integer (kind = 8) :: s_tmp
  Integer :: k1,k2,k3,k4,k_pattern,k_run,k_13
  Integer :: errnum
  Integer :: i_tmp,j_tmp,i_grp,i_dom,d_start,g_start
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,code_kl,kl_run,a_run,code_ab

  Call CC_3d_clean_w()

  ! Calculate intermediate tensors h_ik, h_ca, g_ik, g_ca, and h_ck
  Call CC_3d_aux_h_g()

  ! Calculate intermediate tensors j and k
  Call CC_3d_aux_jk()
  Call CC_3d_add_intl_klcd_to_jk()

  !write(90+myid,*) 'j_aux'
  !k1 = CC_index_k1(CC_mpi_gid+1) - 1
  !do k_run = 1, CC_mem_k1(CC_mpi_gid+1)
  !  k1 = k1 + 1
  !  do k2 = 1, n_k_points
  !    do k3 = 1, n_k_points
  !      Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)
  !      write(90+myid,*) 'k_pattern',k_pattern,k1,k2,k3,k4
  !      do iii = 1, CC_n_occ
  !        do aaa = 1, CC_mem_aa_D(CC_mpi_did+1)
  !          do kkk = 1, CC_n_occ
  !            do ccc = 1, CC_n_vir
  !              write(90+myid,*) iii,aaa+CC_index_aa_D(CC_mpi_did+1)-1,&
  !                   kkk,ccc,CC_j_aux(kkk,ccc,k3,iii,aaa,k2,k_run)
  !            end do
  !          end do
  !        end do
  !      end do
  !    end do
  !  end do
  !end do

  !write(90+myid,*) 'k_aux'
  !k1 = CC_index_k1(CC_mpi_gid+1) - 1
  !do k_run = 1, CC_mem_k1(CC_mpi_gid+1)
  !  k1 = k1 + 1
  !  do k2 = 1, n_k_points
  !    do k3 = 1, n_k_points
  !      Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)
  !      write(90+myid,*) 'k_pattern',k_pattern,k1,k2,k3,k4
  !      do iii = 1, CC_n_occ
  !        do aaa = 1, CC_mem_aa_D(CC_mpi_did+1)
  !          do kkk = 1, CC_n_occ
  !            do ccc = 1, CC_n_vir
  !              write(90+myid,*) iii,aaa+CC_index_aa_D(CC_mpi_did+1)-1,&
  !                   kkk,ccc,CC_k_aux(kkk,ccc,k3,iii,aaa,k2,k_run)
  !            end do
  !          end do
  !        end do
  !      end do
  !    end do
  !  end do
  !end do

  ! Calculate intermediate array a
  Call CC_3d_aux_a()

  !write(90+myid,*) 'a_aux'
  !k_13 = CC_index_k1(CC_mpi_gid+1) - 1
  !do k_run = 1, CC_mem_k1(CC_mpi_gid+1)
  !  k_13 = k_13 + 1
  !  do k1 = 1, n_k_points
  !    Call CC_3d_k_minus(k_13,k1,k3)
  !    do k2 = 1, n_k_points

  !      Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)
  !      write(90+myid,*) 'k_pattern',k_pattern,k1,k2,k3,k4
  !      do iii = 1, CC_n_occ
  !        do jjj = 1, CC_n_occ
  !          code_kl = CC_index_ij_D(CC_mpi_did+1,2) - 1
  !          do kl_run = 1, CC_mem_ij_D(CC_mpi_did+1,2)
  !            code_kl = code_kl + 1
  !            Call CC_3d_decode(code_kl,kkk,lll,CC_n_occ,2)
  !            write(90+myid,*) iii,jjj,kkk,lll,CC_a_aux(kl_run,k1,iii,jjj,k2,k_run)
  !          end do
  !        end do
  !      end do
  !    end do
  !  end do
  !end do

  ! add aux_a to w_d
  Call CC_3d_add_a2w()

  !write(myid+80,*) 'a2w'
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
  !            write(myid+80,*) iii,jjj,aaa,aaa,CC_w_d_d(iii,jjj,a_run,k2,k1,k_run)
  !          end do
  !        end do
  !      end do
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

  ! add aux_b to w_d
  Call CC_3d_add_b2w()

  !write(myid+80,*) 'b2w'
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
  !            write(myid+80,*) iii,jjj,aaa,aaa,CC_w_d_d(iii,jjj,a_run,k2,k1,k_run)
  !          end do
  !        end do
  !      end do
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

  ! Tensor contraction
  Call CC_3d_tensor_contraction()

  ! Calculate w_s
  Call CC_3d_calc_w_s()

  ! add contraction tensor to w
  Call CC_3d_add_ct2w()

  !write(myid+90,*) 'w_d before Jacob'
  !k_13 = CC_index_k1(CC_mpi_gid+1) - 1
  !do k_run = 1, CC_mem_k1(CC_mpi_gid+1)
  !  k_13 = k_13 + 1
  !  do k1 = 1, n_k_points
  !    do k2 = 1, n_k_points
  !      Call CC_3d_k_minus(k_13,k1,k3)
  !      Call CC_3d_k_minus(k_13,k2,k4)
  !      write(myid+90,*) 'k_pattern',k1,k2,k3,k4
  !      do iii = 1, CC_n_occ 
  !        do jjj = 1, CC_n_occ
  !          aaa = CC_index_aa_D(CC_mpi_did+1) - 1
  !          do a_run = 1, CC_mem_aa_D(CC_mpi_did+1)
  !            aaa = aaa + 1
  !            write(myid+90,*) iii,jjj,aaa,aaa,CC_w_d_d(iii,jjj,a_run,k2,k1,k_run)
  !          end do
  !        end do
  !      end do
  !      do iii = 1, CC_n_occ 
  !        do jjj = 1, CC_n_occ
  !          code_ab = CC_index_ab_D(CC_mpi_did+1,3) - 1
  !          do a_run = 1, CC_mem_ab_D(CC_mpi_did+1,3)
  !            code_ab = code_ab + 1
  !            Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)
  !            write(myid+90,*) iii,jjj,aaa,bbb,CC_w_d_nd(iii,jjj,a_run,k2,k1,k_run)
  !          end do
  !        end do
  !      end do
  !    end do
  !  end do
  !end do

  End Subroutine CC_3d_calc_w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_tensor_contraction()

! Tensor contraction
! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kk, k6 = kc
 
  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_run,k_13,iii,jjj,aaa,bbb,a_run

  CC_ten_con = 0.0D0

  Call CC_3d_tc_g_ca()
  Call CC_3d_tc_g_ik()
  Call CC_3d_tc_aibc()
  Call CC_3d_tc_kjac()
  Call CC_3d_tc_aikj()
  Call CC_3d_tc_aikc()
  Call CC_3d_tc_jks()
  
  End Subroutine CC_3d_tensor_contraction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_tc_g_ca()

! Tensor contraction
! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kk, k6 = kc
 
  Use CC_3d
      
  Implicit None

  Integer :: n_k1,k_start,n_a,a_start,a_end,i_tmp
  Integer :: k_run,k1,k2,k3,k4,k_pattern
  Integer :: iii,jjj,aaa,bbb
  Integer :: errnum

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp1i

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: t_tmp,tmp,t_work
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C

  Double complex :: alpha,beta

  alpha = 1.0D0
  beta = 0.0D0

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  Allocate(t_tmp(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                 n_a,n_k_points,n_k1),stat=errnum)
  Call check_allocation(errnum,'t_tmp in CC')

  Call CC_3d_get_t_grp('S',CC_mpi_gid+1,t_tmp,.true.)

  ! g_ca 
  ! Tensor contraction of g_ca(ccc,aaa,kc=ka)
  Allocate(t_work(CC_n_vir,CC_n_occ,CC_n_occ,n_a,1,1,1),stat=errnum)
  Call check_allocation(errnum,'t_work in CC')

  i_tmp = CC_n_occ * CC_n_occ * n_a

  Allocate(mat_A(CC_n_vir,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(CC_n_vir,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(CC_n_vir,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  Allocate(tmp(CC_n_vir,CC_n_occ,CC_n_occ,n_a,1,1,1),stat=errnum)
  Call check_allocation(errnum,'t_work in CC')

  k1 = k_start - 1
  do k_run = 1, n_k1
    k1 = k1 + 1

    do k2 = 1, n_k_points
      do k4 = 1, n_k_points

        Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,3)

        mat_A = CC_g_ca(:,:,k3)
 
        do jjj = 1, CC_n_occ
          t_work(:,jjj,:,:,1,1,1) = t_tmp(jjj,:,k4,:,:,k2,k_run)
        end do

        shp1(1) = CC_n_vir
        shp1(2) = i_tmp
        mat_B = reshape(t_work,shp1)

        Call Zgemm('T','N',CC_n_vir,i_tmp,CC_n_vir,alpha,mat_A,CC_n_vir, &
                   mat_B,CC_n_vir,beta,mat_C,CC_n_vir)

        shp1i(1) = CC_n_vir
        shp1i(2) = CC_n_occ
        shp1i(3) = CC_n_occ
        shp1i(4) = n_a
        tmp(:,:,:,:,1,1,1) = reshape(mat_C,shp1i)

        !$OMP PARALLEL Default(Shared) &
        !$OMP Private(jjj,bbb)
        !$OMP DO COLLAPSE(2)
        do jjj = 1, CC_n_occ
          do bbb = 1, CC_n_vir
            CC_ten_con(jjj,bbb,k4,:,:,k2,k_run) = &
                      CC_ten_con(jjj,bbb,k4,:,:,k2,k_run) + tmp(bbb,jjj,:,:,1,1,1)
          end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
      end do
    end do
  end do

  Deallocate(t_work,mat_A,mat_B,mat_C,tmp,t_tmp)

  End Subroutine CC_3d_tc_g_ca

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_tc_g_ik()

! Tensor contraction
! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kk, k6 = kc
 
  Use CC_3d

  Implicit None

  Integer :: n_k1,k_start,n_a,a_start,a_end,i_tmp
  Integer :: k_run,k1,k2,k3,k4,k_pattern
  Integer :: iii,jjj,aaa,bbb,kkk
  Integer :: errnum

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp1i
  Integer , dimension(6) :: shp2i

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: t_tmp,tmp,t_work
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C

  Double complex :: alpha,beta

  alpha = 1.0D0
  beta = 0.0D0

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  Allocate(t_tmp(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                 n_a,n_k_points,n_k1),stat=errnum)
  Call check_allocation(errnum,'t_tmp in CC')

  Call CC_3d_get_t_grp('S',CC_mpi_gid+1,t_tmp,.true.)

  ! g_ik
  ! Tensor contraction of g_ik(iii,kkk,ki=kk)
  Allocate(t_work(CC_n_occ,CC_n_occ,CC_n_vir,n_k_points, &
                  n_a,n_k1,1),stat=errnum)
  Call check_allocation(errnum,'t_work in CC')

  i_tmp = CC_n_occ * CC_n_vir * n_k_points * n_a * n_k1

  Allocate(mat_A(CC_n_occ,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(CC_n_occ,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(CC_n_occ,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  Allocate(tmp(CC_n_occ,CC_n_occ,CC_n_vir,n_k_points, &
               n_a,n_k1,1),stat=errnum)
  Call check_allocation(errnum,'t_work in CC')

  do k2 = 1, n_k_points

    mat_A = CC_g_ik(:,:,k2)
 
    do kkk = 1, CC_n_occ
      t_work(kkk,:,:,:,:,:,1) = t_tmp(:,:,:,kkk,:,k2,:)
    end do

    shp1(1) = CC_n_occ
    shp1(2) = i_tmp
    mat_B = reshape(t_work,shp1)

    Call Zgemm('N','N',CC_n_occ,i_tmp,CC_n_occ,alpha,mat_A,CC_n_occ, &
               mat_B,CC_n_occ,beta,mat_C,CC_n_occ)

    shp2i(1) = CC_n_occ
    shp2i(2) = CC_n_occ
    shp2i(3) = CC_n_vir
    shp2i(4) = n_k_points
    shp2i(5) = n_a
    shp2i(6) = n_k1
    tmp(:,:,:,:,:,:,1) = reshape(mat_C,shp2i)

    do iii = 1, CC_n_occ
      CC_ten_con(:,:,:,iii,:,k2,:) = CC_ten_con(:,:,:,iii,:,k2,:) &
                                   - tmp(iii,:,:,:,:,:,1)
    end do
  end do

  Deallocate(t_work,mat_A,mat_B,mat_C,tmp,t_tmp)

  End Subroutine CC_3d_tc_g_ik

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_tc_aibc()

! Tensor contraction
! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kk, k6 = kc
 
  Use CC_3d

  Implicit None

  Integer :: n_k1,k_start,n_a,a_start,a_end,i_tmp
  Integer :: k_run,k1,k2,k3,k4
  Integer :: iii,jjj,aaa,bbb
  Integer :: errnum

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp1i
  Integer , dimension(6) :: shp2i

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: tmp,v_work
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C

  Double complex :: alpha,beta

  alpha = 1.0D0
  beta = 0.0D0

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  ! (ai|bc)
  ! Tensor contraction of (ai|bc) (c,b,kc,i,a,ki,ka)
  Allocate(v_work(CC_n_vir,CC_n_vir,CC_n_occ,n_a, &
                  n_k_points,n_k1,1),stat=errnum)
  Call check_allocation(errnum,'v_work in CC')

  i_tmp = CC_n_vir * CC_n_occ * n_a * n_k_points * n_k1

  Allocate(mat_A(CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(CC_n_vir,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(CC_n_occ,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  Allocate(tmp(CC_n_occ,CC_n_vir,CC_n_occ,n_a, &
               n_k_points,n_k1,1),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  do k4 = 1, n_k_points

    mat_A = CC_t_s(:,:,k4)
 
    v_work(:,:,:,:,:,:,1) = CC_intl_aibc(:,:,k4,:,:,:,:)

    shp1(1) = CC_n_vir
    shp1(2) = i_tmp
    mat_B = reshape(v_work,shp1)

    Call Zgemm('N','N',CC_n_occ,i_tmp,CC_n_vir,alpha,mat_A,CC_n_occ, &
               mat_B,CC_n_vir,beta,mat_C,CC_n_occ)

    shp2i(1) = CC_n_occ
    shp2i(2) = CC_n_vir
    shp2i(3) = CC_n_occ
    shp2i(4) = n_a
    shp2i(5) = n_k_points
    shp2i(6) = n_k1
    tmp(:,:,:,:,:,:,1) = reshape(mat_C,shp2i)

    CC_ten_con(:,:,k4,:,:,:,:) = CC_ten_con(:,:,k4,:,:,:,:) &
                               + tmp(:,:,:,:,:,:,1)
  end do

  Deallocate(v_work,mat_A,mat_B,mat_C,tmp)

  End Subroutine CC_3d_tc_aibc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_tc_kjac()

! Tensor contraction
! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kk, k6 = kc
 
  Use CC_3d

  Implicit None

  Integer :: n_k1,k_start,n_a,a_start,a_end,i_tmp
  Integer :: k_run,k1,k2,k3,k4,k_pattern
  Integer :: iii,jjj,aaa,bbb,ccc
  Integer :: errnum

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp1i
  Integer , dimension(6) :: shp2i

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: tmp,tmp2,v_work,t_work
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C,mat_D,mat_E,mat_F

  Double complex :: alpha,beta

  alpha = 1.0D0
  beta = 0.0D0

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  ! (kj|ac)
  ! Tensor contraction of (kj|ac) (k,c,kk,i,a,ki,ka)
  Allocate(v_work(CC_n_vir,CC_n_occ,CC_n_occ,n_a, &
                  1,1,1),stat=errnum)
  Call check_allocation(errnum,'v_work in CC')

  i_tmp = CC_n_occ * CC_n_occ * n_a

  Allocate(mat_A(CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(CC_n_vir,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(CC_n_occ,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  Allocate(tmp(CC_n_occ,CC_n_occ,CC_n_occ,n_a,1,1,1),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  Allocate(t_work(CC_n_occ,CC_n_occ,CC_n_occ,n_a, &
                  1,1,1),stat=errnum)
  Call check_allocation(errnum,'v_work in CC')

  Allocate(mat_D(CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_E(CC_n_occ,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_F(CC_n_vir,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  Allocate(tmp2(CC_n_vir,CC_n_occ,CC_n_occ,n_a,1,1,1),stat=errnum)
  Call check_allocation(errnum,'tmp2 in CC')

  k1 = k_start - 1

  do k_run = 1, n_k1
    k1 = k1 + 1

    do k2 = 1, n_k_points
      do k4 = 1, n_k_points

        Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,3)

        mat_A = CC_t_s(:,:,k2)
 
        do ccc = 1, CC_n_vir
          v_work(ccc,:,:,:,1,1,1) = CC_intl_kiac(:,ccc,k3,:,:,k4,k_run)
        end do

        shp1(1) = CC_n_vir
        shp1(2) = i_tmp
        mat_B = reshape(v_work,shp1)

        Call Zgemm('N','N',CC_n_occ,i_tmp,CC_n_vir,alpha,mat_A,CC_n_occ, &
                   mat_B,CC_n_vir,beta,mat_C,CC_n_occ)

        shp1i(1) = CC_n_occ ! i
        shp1i(2) = CC_n_occ ! k
        shp1i(3) = CC_n_occ ! j
        shp1i(4) = n_a      ! a
        tmp(:,:,:,:,1,1,1) = reshape(mat_C,shp1i)

        mat_D = CC_t_s(:,:,k3)
 
        do iii = 1, CC_n_occ
          t_work(:,:,iii,:,1,1,1) = tmp(iii,:,:,:,1,1,1)
        end do

        shp1(1) = CC_n_occ
        shp1(2) = i_tmp
        mat_E = reshape(t_work,shp1)

        Call Zgemm('T','N',CC_n_vir,i_tmp,CC_n_occ,alpha,mat_D,CC_n_occ, &
                   mat_E,CC_n_occ,beta,mat_F,CC_n_vir)

        shp1i(1) = CC_n_vir ! b
        shp1i(2) = CC_n_occ ! j
        shp1i(3) = CC_n_occ ! i
        shp1i(4) = n_a      ! a
        tmp2(:,:,:,:,1,1,1) = reshape(mat_F,shp1i)

        do jjj = 1, CC_n_occ
          CC_ten_con(jjj,:,k4,:,:,k2,k_run) = CC_ten_con(jjj,:,k4,:,:,k2,k_run) &
                                            - tmp2(:,jjj,:,:,1,1,1)
        end do
      end do
    end do
  end do

  Deallocate(t_work,v_work,mat_A,mat_B,mat_C,mat_D,mat_E,mat_F,tmp,tmp2)

  End Subroutine CC_3d_tc_kjac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_tc_aikj()

! Tensor contraction
! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kk, k6 = kc
 
  Use CC_3d

  Implicit None

  Integer :: n_k1,k_start,n_a,a_start,a_end,i_tmp
  Integer :: k_run,k1,k2,k3,k4,k_pattern
  Integer :: iii,jjj,aaa,bbb,kkk
  Integer :: errnum

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp1i
  Integer , dimension(6) :: shp2i

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: tmp,v_work
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C

  Double complex :: alpha,beta

  alpha = 1.0D0
  beta = 0.0D0

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  ! (ai|kj)
  ! Tensor contraction of (ai|kj) (j,k,kj,i,a,ki,ka)
  Allocate(v_work(CC_n_occ,CC_n_occ,CC_n_occ,n_a,1,1,1),stat=errnum)
  Call check_allocation(errnum,'v_work in CC')

  i_tmp = CC_n_occ * CC_n_occ * n_a

  Allocate(mat_A(CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(CC_n_occ,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(CC_n_vir,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  Allocate(tmp(CC_n_vir,CC_n_occ,CC_n_occ,n_a,1,1,1),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  k1 = k_start - 1
  do k_run = 1, n_k1

    k1 = k1 + 1
    do k2 = 1, n_k_points
      do k4 = 1, n_k_points

        Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,3)
        mat_A = CC_t_s(:,:,k3)
 
        do kkk = 1, CC_n_occ
          v_work(kkk,:,:,:,1,1,1) = CC_intl_aikj(:,kkk,k4,:,:,k2,k_run)
        end do

        shp1(1) = CC_n_occ
        shp1(2) = i_tmp
        mat_B = reshape(v_work,shp1)

        Call Zgemm('T','N',CC_n_vir,i_tmp,CC_n_occ,alpha,mat_A,CC_n_occ, &
                   mat_B,CC_n_occ,beta,mat_C,CC_n_vir)

        shp1i(1) = CC_n_vir
        shp1i(2) = CC_n_occ
        shp1i(3) = CC_n_occ
        shp1i(4) = n_a
        tmp(:,:,:,:,1,1,1) = reshape(mat_C,shp1i)

        do jjj = 1, CC_n_occ
          CC_ten_con(jjj,:,k4,:,:,k2,k_run) = CC_ten_con(jjj,:,k4,:,:,k2,k_run) &
                                            - tmp(:,jjj,:,:,1,1,1)
        end do
      end do
    end do
  end do

  Deallocate(v_work,mat_A,mat_B,mat_C,tmp)

  End Subroutine CC_3d_tc_aikj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_tc_aikc()

! Tensor contraction
! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kk, k6 = kc
 
  Use CC_3d

  Implicit None

  Integer :: n_k1,k_start,n_a,a_start,a_end,i_tmp
  Integer :: k_run,k1,k2,k3,k4,k_pattern
  Integer :: iii,jjj,aaa,bbb,kkk,ccc
  Integer :: errnum

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp1i
  Integer , dimension(6) :: shp2i

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: tmp,tmp2,v_work,t_work,t_mp2
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C,mat_D,mat_E,mat_F

  Double complex :: alpha,beta

  alpha = 1.0D0
  beta = 0.0D0

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  ! (ai|kc)
  ! Tensor contraction of (ai|kc) (k,c,kk,i,a,ki,ka)
  Allocate(v_work(CC_n_vir,CC_n_occ,CC_n_occ,n_a, &
                  1,1,1),stat=errnum)
  Call check_allocation(errnum,'v_work in CC')

  i_tmp = CC_n_occ * CC_n_occ * n_a

  Allocate(mat_A(CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(CC_n_vir,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(CC_n_occ,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  Allocate(tmp(CC_n_occ,CC_n_occ,CC_n_occ,n_a,1,1,1),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  Allocate(t_work(CC_n_occ,CC_n_occ,CC_n_occ,n_a, &
                  1,1,1),stat=errnum)
  Call check_allocation(errnum,'v_work in CC')

  Allocate(mat_D(CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_E(CC_n_occ,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_F(CC_n_vir,i_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  Allocate(tmp2(CC_n_vir,CC_n_occ,CC_n_occ,n_a,1,1,1),stat=errnum)
  Call check_allocation(errnum,'tmp2 in CC')

  k1 = k_start - 1

  do k_run = 1, n_k1
    k1 = k1 + 1

    do k2 = 1, n_k_points
      do k4 = 1, n_k_points

        Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,3)

        mat_A = CC_t_s(:,:,k4)
 
        do ccc = 1, CC_n_vir
          v_work(ccc,:,:,:,1,1,1) = CC_intl_aikc(:,ccc,k3,:,:,k2,k_run)
        end do

        shp1(1) = CC_n_vir
        shp1(2) = i_tmp
        mat_B = reshape(v_work,shp1)

        Call Zgemm('N','N',CC_n_occ,i_tmp,CC_n_vir,alpha,mat_A,CC_n_occ, &
                   mat_B,CC_n_vir,beta,mat_C,CC_n_occ)

        shp1i(1) = CC_n_occ ! j
        shp1i(2) = CC_n_occ ! k
        shp1i(3) = CC_n_occ ! i
        shp1i(4) = n_a      ! a
        tmp(:,:,:,:,1,1,1) = reshape(mat_C,shp1i)

        mat_D = CC_t_s(:,:,k3)
 
        do kkk = 1, CC_n_occ
          t_work(kkk,:,:,:,1,1,1) = tmp(:,kkk,:,:,1,1,1)
        end do

        shp1(1) = CC_n_occ
        shp1(2) = i_tmp
        mat_E = reshape(t_work,shp1)

        Call Zgemm('T','N',CC_n_vir,i_tmp,CC_n_occ,alpha,mat_D,CC_n_occ, &
                   mat_E,CC_n_occ,beta,mat_F,CC_n_vir)

        shp1i(1) = CC_n_vir ! b
        shp1i(2) = CC_n_occ ! j
        shp1i(3) = CC_n_occ ! i
        shp1i(4) = n_a      ! a
        tmp2(:,:,:,:,1,1,1) = reshape(mat_F,shp1i)

        do jjj = 1, CC_n_occ
          CC_ten_con(jjj,:,k4,:,:,k2,k_run) = CC_ten_con(jjj,:,k4,:,:,k2,k_run) &
                                            - tmp2(:,jjj,:,:,1,1,1)
        end do
      end do
    end do
  end do

  Deallocate(t_work,v_work,mat_A,mat_B,mat_C,mat_D,mat_E,mat_F,tmp,tmp2)

  End Subroutine CC_3d_tc_aikc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_tc_jks()

! Tensor contraction
! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kk, k6 = kc
 
  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k5,k6,k1_run,k2_run,k3_run,k_pattern
  Integer :: n_k1,k_start,n_a,a_start,a_end,n_c,c_start,c_end,a_run,c_run
  Integer :: i_start,i_end,l_run,n_l
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,i_tmp,j_tmp,l_tmp

  Integer :: i_grp_task,n_k1_task,k_start_task,k_end_task,i_task
  Integer :: errnum,source_id
 
  Integer (kind=8) :: s_tmp

  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp1i
  Integer , dimension(6) :: shp2i

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: tmp,v_work,t_work,j_work
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: t_tmp,t_tmp2
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C,mat_D,mat_E,mat_F

  Double complex :: alpha,beta

  alpha = 1.0D0
  beta = 0.0D0

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  i_start = CC_index_ii_D(CC_mpi_did+1)
  i_end = i_start - 1 + CC_mem_ii_D(CC_mpi_did+1)

  ! Tensor contraction of aux j, k and w_s
  do i_grp_task = 1, CC_mpi_group_size

    n_k1_task = CC_mem_k1(i_grp_task)
    k_start_task = CC_index_k1(i_grp_task)
    k_end_task = k_start_task - 1 + n_k1_task

    Allocate(t_tmp(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_a, &
                   n_k_points,n_k1_task),stat=errnum)
    Call check_allocation(errnum,'t_tmp in CC')

    Call CC_3d_get_t_grp('S',i_grp_task,t_tmp,.false.)

    do i_task = 1, CC_mpi_domain_size
      
      n_c = CC_mem_aa_D(i_task)
      c_start = CC_index_aa_D(i_task)
      c_end = c_start - 1 + n_c

      Allocate(t_tmp2(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_c, &
                      n_k_points,n_k1_task),stat=errnum)
      Call check_allocation(errnum,'t_tmp in CC')

      i_tmp = CC_n_occ * CC_n_vir * n_k_points * CC_n_occ * n_c
      j_tmp = n_k_points * n_k1_task
      s_tmp = Int(i_tmp,8) * Int(j_tmp,8)
      source_id = i_task - 1

      if (source_id.eq.CC_mpi_did) then
        t_tmp2 = t_tmp
      else
        t_tmp2 = 0.0D0
      end if

      Call CC_mpi_complex_bcast(s_tmp, t_tmp2, source_id, CC_mpi_comm_domain)

      Allocate(j_work(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_a,1,1),stat=errnum)
      Call check_allocation(errnum,'j_work in CC')

      Allocate(t_work(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ,n_c,1,1),stat=errnum)
      Call check_allocation(errnum,'t_work in CC')

      i_tmp = CC_n_occ * CC_n_vir * n_k_points
      j_tmp = CC_n_occ * n_c
      l_tmp = CC_n_occ * n_a

      Allocate(mat_A(i_tmp,j_tmp),stat=errnum)
      Call check_allocation(errnum,'mat_A in CC')

      Allocate(mat_B(i_tmp,l_tmp),stat=errnum)
      Call check_allocation(errnum,'mat_B in CC')

      Allocate(mat_C(j_tmp,l_tmp),stat=errnum)
      Call check_allocation(errnum,'mat_C in CC')

      Allocate(tmp(CC_n_occ,n_c,CC_n_occ,n_a,1,1,1),stat=errnum)
      Call check_allocation(errnum,'tmp in CC')

      ! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kk, k6 = kc
      k1 = k_start - 1
      do k1_run = 1, n_k1

        k1 = k1 + 1

        do k2 = 1, n_k_points
          
          k3 = k_start_task - 1
          do k3_run = 1, n_k1_task
            k3 = k3 + 1

            Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

            ! Tensor contraction of j

            !$OMP PARALLEL Default(Shared) &
            !$OMP Private(k5,jjj,kkk)
            !$OMP DO COLLAPSE(3)
            do k5 = 1, n_k_points
              do jjj = 1, CC_n_occ
                do kkk = 1, CC_n_occ
                  t_work(kkk,:,k5,jjj,:,1,1) = t_tmp2(kkk,:,k5,jjj,:,k4,k3_run) &
                                     - 0.5D0 * t_tmp2(jjj,:,k4,kkk,:,k5,k3_run)
                end do
              end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL

            j_work(:,:,:,:,:,1,1) = CC_j_aux(:,:,:,:,:,k2,k1_run)

            shp1(1) = i_tmp
            shp1(2) = j_tmp
            mat_A = reshape(t_work,shp1)

            shp1(1) = i_tmp
            shp1(2) = l_tmp
            mat_B = reshape(j_work,shp1)

            Call Zgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,mat_A,i_tmp, &
                       mat_B,i_tmp,beta,mat_C,j_tmp)

            shp1i(1) = CC_n_occ ! j
            shp1i(2) = n_c      ! b
            shp1i(3) = CC_n_occ ! i
            shp1i(4) = n_a      ! a
            tmp(:,:,:,:,1,1,1) = reshape(mat_C,shp1i)

            CC_ten_con(:,c_start:c_end,k4,:,:,k2,k1_run) = &
                    CC_ten_con(:,c_start:c_end,k4,:,:,k2,k1_run) + tmp(:,:,:,:,1,1,1)

            ! Tensor contraction of k
            do k5 = 1, n_k_points
              do jjj = 1, CC_n_occ
                do kkk = 1, CC_n_occ
                  t_work(kkk,:,k5,jjj,:,1,1) = t_tmp2(jjj,:,k4,kkk,:,k5,k3_run)
                end do
              end do
            end do

            j_work(:,:,:,:,:,1,1) = CC_k_aux(:,:,:,:,:,k2,k1_run)

            shp1(1) = i_tmp
            shp1(2) = j_tmp
            mat_A = reshape(t_work,shp1)

            shp1(1) = i_tmp
            shp1(2) = l_tmp
            mat_B = reshape(j_work,shp1)

            Call Zgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,mat_A,i_tmp, &
                       mat_B,i_tmp,beta,mat_C,j_tmp)

            shp1i(1) = CC_n_occ ! j
            shp1i(2) = n_c      ! b
            shp1i(3) = CC_n_occ ! i
            shp1i(4) = n_a      ! a
            tmp(:,:,:,:,1,1,1) = reshape(mat_C,shp1i)

            CC_ten_con(:,c_start:c_end,k4,:,:,k2,k1_run) = &
              CC_ten_con(:,c_start:c_end,k4,:,:,k2,k1_run) - 0.5D0 * tmp(:,:,:,:,1,1,1)

            !$OMP PARALLEL Default(Shared) &
            !$OMP Private(iii,jjj)
            !$OMP DO COLLAPSE(2)
            do iii = 1, CC_n_occ
              do jjj = 1, CC_n_occ
                CC_ten_con(iii,c_start:c_end,k2,jjj,:,k4,k1_run) = &
                    CC_ten_con(iii,c_start:c_end,k2,jjj,:,k4,k1_run) - tmp(jjj,:,iii,:,1,1,1)
              end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL
          end do
        end do
      end do

      Deallocate(t_work,j_work,mat_A,mat_B,mat_C,tmp)

      ! Tensor contraction in w_s
      ! (ac|kd)
      Allocate(v_work(CC_n_occ,CC_n_vir,n_k_points, & 
                      n_c,n_k1_task,n_a,1),stat=errnum)
      Call check_allocation(errnum,'v_work in CC')

      Allocate(t_work(CC_n_occ,CC_n_vir,n_k_points, &
                      n_c,n_k1_task,CC_n_occ,1),stat=errnum)
      Call check_allocation(errnum,'t_work in CC')

      i_tmp = CC_n_occ * CC_n_vir * n_k_points * n_c * n_k1_task
      j_tmp = CC_n_occ 
      l_tmp = n_a

      Allocate(mat_A(i_tmp,j_tmp),stat=errnum)
      Call check_allocation(errnum,'mat_A in CC')

      Allocate(mat_B(i_tmp,l_tmp),stat=errnum)
      Call check_allocation(errnum,'mat_B in CC')

      Allocate(mat_C(j_tmp,l_tmp),stat=errnum)
      Call check_allocation(errnum,'mat_C in CC')

      ! k1 = ka = ki, k2 = kc, k3 = kk, k4 = kd
      k1 = k_start - 1
      do k1_run = 1, n_k1
        k1 = k1 + 1

        !$OMP PARALLEL Default(Shared) &
        !$OMP Private(k3,k2,k2_run,k4,ccc,c_run,ddd)
        do k3 = 1, n_k_points

          k2 = k_start_task - 1
          do k2_run = 1, n_k1_task
            k2 = k2 + 1

            Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)
            !$OMP DO COLLAPSE(2)
            do c_run = 1, n_c
              do ddd = 1, CC_n_vir
                ccc = c_start - 1 + c_run
                v_work(:,ddd,k3,c_run,k2_run,:,1) = &
                          2.0D0 * CC_intl_ackd(:,ddd,k3,ccc,:,k2,k1_run) &
                                - CC_intl_ackd(:,ccc,k3,ddd,:,k4,k1_run)
              end do
            end do
            !$OMP END DO
          end do
        end do
        !$OMP END PARALLEL 

        do iii = 1, CC_n_occ
          t_work(:,:,:,:,:,iii,1) = t_tmp2(:,:,:,iii,:,k1,:)
        end do

        if (i_grp_task.eq.CC_mpi_gid+1) then
          !$OMP PARALLEL Default(Shared) &
          !$OMP Private(k3,iii,kkk,c_run,ddd,ccc)
          !$OMP DO COLLAPSE(5)
          do k3 = 1, n_k_points
            do iii = 1, CC_n_occ
              do kkk = 1, CC_n_occ
                do c_run = 1, n_c
                  do ddd = 1, CC_n_vir
                    ccc = c_start - 1 + c_run
                    t_work(kkk,ddd,k3,c_run,k1_run,iii,1) = &
                              t_work(kkk,ddd,k3,c_run,k1_run,iii,1) &
                            + CC_t_s(kkk,ddd,k3) * CC_t_s(iii,ccc,k1)
                  end do
                end do
              end do
            end do
          end do
          !$OMP END DO
          !$OMP END PARALLEL
        end if

        shp1(1) = i_tmp
        shp1(2) = j_tmp
        mat_A = reshape(t_work,shp1)

        shp1(1) = i_tmp
        shp1(2) = l_tmp
        mat_B = reshape(v_work,shp1)

        Call Zgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,mat_A,i_tmp, &
                   mat_B,i_tmp,beta,mat_C,j_tmp)

        CC_w_s(:,a_start:a_end,k1) = CC_w_s(:,a_start:a_end,k1) + mat_C

      end do

      Deallocate(v_work,t_work,mat_A,mat_B,mat_C)

      ! (ki|lc)
      ! v_work(k,c,kk,l,kl,i,1)
      Allocate(v_work(CC_n_occ,CC_n_vir,n_k_points, & 
                      CC_mem_ii_D(CC_mpi_did+1),n_k1_task,CC_n_occ,1),stat=errnum)
      Call check_allocation(errnum,'v_work in CC')

      ! t_work(k,c,kk,l,kl,a,1)
      Allocate(t_work(CC_n_occ,CC_n_vir,n_k_points, &
                      CC_mem_ii_D(CC_mpi_did+1),n_k1_task,n_c,1),stat=errnum)
      Call check_allocation(errnum,'t_work in CC')

      i_tmp = CC_n_occ * CC_n_vir * n_k_points * CC_mem_ii_D(CC_mpi_did+1) * n_k1_task
      j_tmp = CC_n_occ
      l_tmp = n_c

      Allocate(mat_A(i_tmp,j_tmp),stat=errnum)
      Call check_allocation(errnum,'mat_A in CC')

      Allocate(mat_B(i_tmp,l_tmp),stat=errnum)
      Call check_allocation(errnum,'mat_B in CC')

      Allocate(mat_C(j_tmp,l_tmp),stat=errnum)
      Call check_allocation(errnum,'mat_C in CC')

      ! k1 = ka = ki, k2 = kl, k3 = kk, k4 = kc
      k1 = k_start_task - 1
      do k1_run = 1, n_k1_task
        k1 = k1 + 1

        do iii = 1, CC_n_occ
          v_work(:,:,:,:,:,iii,1) = 2.0D0 * CC_intl_kilc(:,:,:,iii,:,k1,:) &
                                          - CC_intl_likc(:,:,:,iii,:,k1,:)
        end do

        k2 = k_start - 1
        do k2_run = 1, n_k1
          k2 = k2 + 1
          do k3 = 1, n_k_points
            lll = CC_index_ii_D(CC_mpi_did+1) - 1
            do l_run = 1, CC_mem_ii_D(CC_mpi_did+1)
              lll = lll + 1
              do kkk = 1, CC_n_occ
                t_work(kkk,:,k3,l_run,k2_run,:,1) = t_tmp2(lll,:,k2,kkk,:,k3,k1_run)
              end do
            end do
          end do

          n_l = CC_mem_ii_D(CC_mpi_did+1)
          !$OMP PARALLEL Default(Shared) &
          !$OMP Private(aaa,a_run,kkk,l_run,lll,ccc)
          !$OMP DO COLLAPSE(4)
          !!$OMP DO
          do a_run = 1, n_c
            do kkk = 1, CC_n_occ
              do l_run = 1, n_l
                do ccc = 1, CC_n_vir

                  lll = CC_index_ii_D(CC_mpi_did+1) - 1 + l_run
                  aaa = c_start - 1 + a_run
                  t_work(kkk,ccc,k1,l_run,k2_run,a_run,1) = &
                         t_work(kkk,ccc,k1,l_run,k2_run,a_run,1) &
                       + CC_t_s(kkk,aaa,k1) * CC_t_s(lll,ccc,k2)
                end do
              end do
            end do
          end do
          !$OMP END DO
          !$OMP END PARALLEL
        end do

        shp1(1) = i_tmp
        shp1(2) = j_tmp
        mat_A = reshape(v_work,shp1)

        shp1(1) = i_tmp
        shp1(2) = l_tmp
        mat_B = reshape(t_work,shp1)

        Call Zgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,mat_A,i_tmp, &
                   mat_B,i_tmp,beta,mat_C,j_tmp)

        CC_w_s(:,c_start:c_end,k1) = CC_w_s(:,c_start:c_end,k1) - mat_C

      end do

      Deallocate(t_work,v_work,mat_A,mat_B,mat_C)

      Deallocate(t_tmp2)
    end do

    Deallocate(t_tmp)

  end do

  End Subroutine CC_3d_tc_jks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_calc_w_s()

! (iii,aaa,k1)

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_pattern,k_run,n_k1,k_start,k_end
  Integer :: iii,jjj,aaa,bbb,kkk,lll,ccc,ddd,i_run,i_tmp,j_tmp,a_run,c_run,n_a
  Integer :: n_dom,a_start,a_end,i_start,i_end
  Integer :: i_task,source_id,errnum
  Integer , dimension(2) :: shp1
  Integer , dimension(3) :: shp1i

  Integer (kind=8) :: s_tmp

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: t_tmp,v_tmp,t_work,v_work,t_tmp2
  Double complex , dimension(:,:,:) , allocatable :: tmp,w_tmp
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C

  Double complex :: alpha,beta

  alpha = 1.0D0
  beta = 0.0D0

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + CC_mem_aa_D(CC_mpi_did+1)

  ! h_ca & h_ik
  if (CC_mpi_did.eq.0) then

    k1 = k_start - 1
    do k_run = 1, n_k1
      k1 = k1 + 1

      ! Tensor contraction of h_ca(ccc,aaa,kc=ka)
      Allocate(mat_A(CC_n_occ,CC_n_vir),stat=errnum)
      Call check_allocation(errnum,'mat_A in CC')

      Allocate(mat_B(CC_n_vir,CC_n_vir),stat=errnum)
      Call check_allocation(errnum,'mat_B in CC')

      Allocate(mat_C(CC_n_occ,CC_n_vir),stat=errnum)
      Call check_allocation(errnum,'mat_C in CC')

      mat_A = CC_t_s(:,:,k1)
      mat_B = CC_h_ca(:,:,k1)

      Call Zgemm('N','N',CC_n_occ,CC_n_vir,CC_n_vir,alpha,mat_A,CC_n_occ, &
                 mat_B,CC_n_vir,beta,mat_C,CC_n_occ)

      CC_w_s(:,:,k1) = CC_w_s(:,:,k1) + mat_C

      Deallocate(mat_A,mat_B)

      ! Tensor contraction of h_ik(iii,kkk,ki=kk)
      Allocate(mat_A(CC_n_occ,CC_n_occ),stat=errnum)
      Call check_allocation(errnum,'mat_A in CC')

      Allocate(mat_B(CC_n_occ,CC_n_vir),stat=errnum)
      Call check_allocation(errnum,'mat_B in CC')

      mat_A = CC_h_ik(:,:,k1)
      mat_B = CC_t_s(:,:,k1)

      Call Zgemm('N','N',CC_n_occ,CC_n_vir,CC_n_occ,alpha,mat_A,CC_n_occ, &
                 mat_B,CC_n_occ,beta,mat_C,CC_n_occ)

      CC_w_s(:,:,k1) = CC_w_s(:,:,k1) - mat_C

      Deallocate(mat_A,mat_B,mat_C)

    end do

  end if

  ! h_ck
  Allocate(t_tmp(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                 CC_mem_aa_D(CC_mpi_did+1),n_k_points,n_k1),stat=errnum)
  Call check_allocation(errnum,'t_tmp in CC')

  Call CC_3d_get_t_grp('S',CC_mpi_gid+1,t_tmp,.true.)

  Allocate(t_work(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                  CC_mem_aa_D(CC_mpi_did+1),1,1),stat=errnum)
  Call check_allocation(errnum,'t_work in CC')

  i_tmp = CC_n_occ * CC_n_vir * n_k_points
  j_tmp = CC_n_occ * CC_mem_aa_D(CC_mpi_did+1)

  Allocate(mat_A(i_tmp,j_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(i_tmp,1),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(j_tmp,1),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  Allocate(tmp(CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),1),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  k1 = k_start - 1
  do k_run = 1, n_k1
    k1 = k1 + 1
    do k3 = 1, n_k_points

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(iii,kkk)
      !$OMP DO COLLAPSE(2)
      do iii = 1, CC_n_occ
        do kkk  = 1, CC_n_occ
          t_work(kkk,:,k3,iii,:,1,1) = 2.0D0 * t_tmp(kkk,:,k3,iii,:,k1,k_run) &
                                             - t_tmp(iii,:,k1,kkk,:,k3,k_run)
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end do

    n_a = CC_mem_aa_D(CC_mpi_did+1)
    !$OMP PARALLEL Default(Shared) &
    !$OMP Private(a_run,aaa,iii,kkk,ccc)
    !$OMP DO COLLAPSE(4)
    do a_run = 1, n_a
      do iii = 1, CC_n_occ
        do kkk  = 1, CC_n_occ
          do ccc = 1, CC_n_vir
            aaa = a_start - 1 + a_run
            t_work(kkk,ccc,k1,iii,a_run,1,1) = t_work(kkk,ccc,k1,iii,a_run,1,1) &
                                             + CC_t_s(iii,ccc,k1) * CC_t_s(kkk,aaa,k1)
          end do
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    shp1(1) = i_tmp
    shp1(2) = j_tmp
    mat_A = reshape(t_work,shp1)

    shp1(1) = i_tmp
    shp1(2) = 1
    mat_B = reshape(CC_h_ck,shp1)

    Call Zgemm('T','N',j_tmp,1,i_tmp,alpha,mat_A,i_tmp, &
               mat_B,i_tmp,beta,mat_C,j_tmp)

    shp1(1) = CC_n_occ
    shp1(2) = CC_mem_aa_D(CC_mpi_did+1)
    tmp(:,:,1) = reshape(mat_C,shp1)

    CC_w_s(:,a_start:a_end,k1) = CC_w_s(:,a_start:a_end,k1) + tmp(:,:,1)

  end do

  Deallocate(t_work,mat_A,mat_B,mat_C,tmp)


  ! (ai|kc) & (ki|ac)
  Allocate(v_work(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                  CC_mem_aa_D(CC_mpi_did+1),1,1),stat=errnum)
  Call check_allocation(errnum,'v_work in CC')

  i_tmp = CC_n_occ * CC_n_vir * n_k_points
  j_tmp = CC_n_occ * CC_mem_aa_D(CC_mpi_did+1)

  Allocate(mat_A(i_tmp,j_tmp),stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(mat_B(i_tmp,1),stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(mat_C(j_tmp,1),stat=errnum)
  Call check_allocation(errnum,'mat_C in CC')

  Allocate(tmp(CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),1),stat=errnum)
  Call check_allocation(errnum,'tmp in CC')

  k1 = k_start - 1
  do k_run = 1, n_k1
    k1 = k1 + 1

    v_work(:,:,:,:,:,1,1) = 2.0D0 * CC_intl_aikc(:,:,:,:,:,k1,k_run) &
                                  - CC_intl_kiac(:,:,:,:,:,k1,k_run)

    shp1(1) = i_tmp
    shp1(2) = j_tmp
    mat_A = reshape(v_work,shp1)

    shp1(1) = i_tmp
    shp1(2) = 1
    mat_B = reshape(CC_t_s,shp1)

    Call Zgemm('T','N',j_tmp,1,i_tmp,alpha,mat_A,i_tmp, &
               mat_B,i_tmp,beta,mat_C,j_tmp)

    shp1(1) = CC_n_occ
    shp1(2) = CC_mem_aa_D(CC_mpi_did+1)
    tmp(:,:,1) = reshape(mat_C,shp1)

    CC_w_s(:,a_start:a_end,k1) = CC_w_s(:,a_start:a_end,k1) + tmp(:,:,1)

  end do

  Deallocate(v_work,mat_A,mat_B,mat_C,tmp)

  s_tmp = Int(CC_n_occ * CC_n_vir * n_k_points,8)
  Call CC_mpi_complex_allreduce(s_tmp, CC_w_s, MPI_COMM_WORLD)
 
  End Subroutine CC_3d_calc_w_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_add_a2w()

! Tensor contraction
! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kk, k6 = kl
 
  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k5,k6,k_pattern,k_run,n_k1,k_start,k_end,k_13
  Integer :: n_ab,n_a,code_kl,kl_run,n_kl,kl_start,code_ij
  Integer :: iii,jjj,kkk,lll,aaa,bbb,i_run,i_tmp,j_tmp,c_run,l_tmp,ij_tmp,ij_run
  Integer :: i_task,target_id,errnum,source_id
  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp1i
  Integer , dimension(3) :: shp2i

  Integer (kind=8) :: s_tmp

  Double complex , dimension(:,:,:,:,:,:) , allocatable :: aux_tmp,w_tmp
  Double complex , dimension(:,:,:,:) , allocatable :: t_plus,t_minus
  Double complex , dimension(:,:,:,:) , allocatable :: a_plus,a_minus
  Double complex , dimension(:,:,:,:) , allocatable :: rlt_S,rlt_A
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C,mat_D,mat_E,mat_F

  Double complex :: alpha,beta,coeff

  alpha = 1.0D0
  beta = 0.0D0

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)

  n_ab = CC_mem_ab_D(CC_mpi_did+1,3)
  n_a = CC_mem_aa_D(CC_mpi_did+1)

  do i_task = 1, CC_mpi_domain_size

    n_kl = CC_mem_ij_D(i_task,2)
    kl_start = CC_index_ij_D(i_task,2)

    Allocate(aux_tmp(n_kl,n_k_points,CC_n_occ,CC_n_occ,n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'t_plus in CC')

    source_id = i_task - 1

    if (CC_mpi_did.eq.source_id) then
      aux_tmp = CC_a_aux
    else
      aux_tmp = 0.0D0
    end if

    s_tmp = Int(n_kl,8) * Int(n_k1,8) * Int(n_k_points**2,8) * Int(CC_n_occ**2,8)

    Call CC_mpi_complex_bcast(s_tmp, aux_tmp, source_id, CC_mpi_comm_domain)

    ! For w_d_k2
    ! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kk, k6 = kl
    Allocate(w_tmp(CC_n_occ,CC_n_occ,n_ab,n_k_points,n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'t_plus in CC')

    ij_tmp = CC_n_occ * (CC_n_occ + 1) / 2
    i_tmp = ij_tmp * n_k_points

    Allocate(t_plus(n_kl,n_k_points,n_ab,n_k_points),stat=errnum)
    Call check_allocation(errnum,'t_plus in CC')

    Allocate(t_minus(n_kl,n_k_points,n_ab,n_k_points),stat=errnum)
    Call check_allocation(errnum,'t_minus in CC')

    Allocate(a_plus(n_kl,n_k_points,ij_tmp,n_k_points),stat=errnum)
    Call check_allocation(errnum,'a_plus in CC')

    Allocate(a_minus(n_kl,n_k_points,ij_tmp,n_k_points),stat=errnum)
    Call check_allocation(errnum,'a_minus in CC')

    Allocate(rlt_S(ij_tmp,n_k_points,n_ab,n_k_points),stat=errnum)
    Call check_allocation(errnum,'a_minus in CC')

    Allocate(rlt_A(ij_tmp,n_k_points,n_ab,n_k_points),stat=errnum)
    Call check_allocation(errnum,'a_minus in CC')

    j_tmp = n_kl * n_k_points
    l_tmp = n_ab * n_k_points
    Allocate(mat_A(j_tmp,i_tmp),stat=errnum)
    Call check_allocation(errnum,'mat_A in CC')

    Allocate(mat_B(j_tmp,l_tmp),stat=errnum)
    Call check_allocation(errnum,'mat_B in CC')

    Allocate(mat_C(i_tmp,l_tmp),stat=errnum)
    Call check_allocation(errnum,'mat_C in CC')

    k_13 = k_start - 1
    do k_run = 1, n_k1

      k_13 = k_13 + 1
 
      do k5 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k5,k6)

        code_kl = kl_start - 1
        do kl_run = 1, n_kl

          code_kl = code_kl + 1
          Call CC_3d_decode(code_kl,kkk,lll,CC_n_occ,2)

          if (kkk.ne.lll) then
            coeff = 1.0D0
          else
            coeff = 0.5D0
          end if

          t_plus(kl_run,k5,:,:) = coeff * (CC_t_d_nd(kkk,lll,:,k5,:,k_run) &
                                         + CC_t_d_nd(lll,kkk,:,k6,:,k_run))

          t_minus(kl_run,k5,:,:) = coeff * (CC_t_d_nd(kkk,lll,:,k5,:,k_run) &
                                          - CC_t_d_nd(lll,kkk,:,k6,:,k_run))

        end do
      end do

      do k2 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k2,k4)

        do ij_run = 1, ij_tmp

          Call CC_3d_decode(ij_run,iii,jjj,CC_n_occ,2)

          a_plus(:,:,ij_run,k2) = 0.5D0 * (aux_tmp(:,:,iii,jjj,k2,k_run) &
                                         + aux_tmp(:,:,jjj,iii,k4,k_run))

          a_minus(:,:,ij_run,k2) = 0.5D0 * (aux_tmp(:,:,iii,jjj,k2,k_run) &
                                          - aux_tmp(:,:,jjj,iii,k4,k_run))

        end do
      end do

      shp1(1) = j_tmp
      shp1(2) = i_tmp
      mat_A = reshape(a_plus,shp1)

      shp1(1) = j_tmp
      shp1(2) = l_tmp
      mat_B = reshape(t_plus,shp1)

      Call Zgemm('T','N',i_tmp,l_tmp,j_tmp,alpha,mat_A,j_tmp, &
                 mat_B,j_tmp,beta,mat_C,i_tmp)

      shp1i(1) = ij_tmp
      shp1i(2) = n_k_points
      shp1i(3) = n_ab
      shp1i(4) = n_k_points
      rlt_S = reshape(mat_C,shp1i)

      shp1(1) = j_tmp
      shp1(2) = i_tmp
      mat_A = reshape(a_minus,shp1)

      shp1(1) = j_tmp
      shp1(2) = l_tmp
      mat_B = reshape(t_minus,shp1)

      Call Zgemm('T','N',i_tmp,l_tmp,j_tmp,alpha,mat_A,j_tmp, &
                 mat_B,j_tmp,beta,mat_C,i_tmp)

      shp1i(1) = ij_tmp
      shp1i(2) = n_k_points
      shp1i(3) = n_ab
      shp1i(4) = n_k_points
      rlt_A = reshape(mat_C,shp1i)

      do k2 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k2,k4)

        do ij_run = 1, ij_tmp

          Call CC_3d_decode(ij_run,iii,jjj,CC_n_occ,2)

          w_tmp(iii,jjj,:,k2,:,k_run) = rlt_S(ij_run,k2,:,:) + rlt_A(ij_run,k2,:,:)
          w_tmp(jjj,iii,:,k4,:,k_run) = rlt_S(ij_run,k2,:,:) - rlt_A(ij_run,k2,:,:)

        end do
      end do

    end do

    CC_w_d_nd = CC_w_d_nd + w_tmp

    Deallocate(w_tmp,t_plus,t_minus,a_plus,a_minus,rlt_S,rlt_A,mat_A,mat_B,mat_C)
 
    ! For w_d
    ! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kk, k6 = kl
    Allocate(w_tmp(CC_n_occ,CC_n_occ,n_a,n_k_points,n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'t_plus in CC')

    ij_tmp = CC_n_occ * (CC_n_occ + 1) / 2
    i_tmp = ij_tmp * n_k_points

    Allocate(t_plus(n_kl,n_k_points,n_a,n_k_points),stat=errnum)
    Call check_allocation(errnum,'t_plus in CC')

    Allocate(t_minus(n_kl,n_k_points,n_a,n_k_points),stat=errnum)
    Call check_allocation(errnum,'t_minus in CC')

    Allocate(a_plus(n_kl,n_k_points,ij_tmp,n_k_points),stat=errnum)
    Call check_allocation(errnum,'a_plus in CC')

    Allocate(a_minus(n_kl,n_k_points,ij_tmp,n_k_points),stat=errnum)
    Call check_allocation(errnum,'a_minus in CC')

    Allocate(rlt_S(ij_tmp,n_k_points,n_a,n_k_points),stat=errnum)
    Call check_allocation(errnum,'a_minus in CC')

    Allocate(rlt_A(ij_tmp,n_k_points,n_a,n_k_points),stat=errnum)
    Call check_allocation(errnum,'a_minus in CC')

    j_tmp = n_kl * n_k_points
    l_tmp = n_a * n_k_points
    Allocate(mat_A(j_tmp,i_tmp),stat=errnum)
    Call check_allocation(errnum,'mat_A in CC')

    Allocate(mat_B(j_tmp,l_tmp),stat=errnum)
    Call check_allocation(errnum,'mat_B in CC')

    Allocate(mat_C(i_tmp,l_tmp),stat=errnum)
    Call check_allocation(errnum,'mat_C in CC')

    k_13 = k_start - 1
    do k_run = 1, n_k1

      k_13 = k_13 + 1
 
      do k5 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k5,k6)

        code_kl = kl_start - 1
        do kl_run = 1, n_kl

          code_kl = code_kl + 1
          Call CC_3d_decode(code_kl,kkk,lll,CC_n_occ,2)

          if (kkk.ne.lll) then
            coeff = 1.0D0
          else
            coeff = 0.5D0
          end if

          t_plus(kl_run,k5,:,:) = coeff * (CC_t_d_d(kkk,lll,:,k5,:,k_run) &
                                         + CC_t_d_d(lll,kkk,:,k6,:,k_run))

          t_minus(kl_run,k5,:,:) = coeff * (CC_t_d_d(kkk,lll,:,k5,:,k_run) &
                                          - CC_t_d_d(lll,kkk,:,k6,:,k_run))

        end do
      end do

      do k2 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k2,k4)

        do ij_run = 1, ij_tmp

          Call CC_3d_decode(ij_run,iii,jjj,CC_n_occ,2)

          a_plus(:,:,ij_run,k2) = 0.5D0 * (aux_tmp(:,:,iii,jjj,k2,k_run) &
                                         + aux_tmp(:,:,jjj,iii,k4,k_run))

          a_minus(:,:,ij_run,k2) = 0.5D0 * (aux_tmp(:,:,iii,jjj,k2,k_run) &
                                          - aux_tmp(:,:,jjj,iii,k4,k_run))

        end do
      end do

      shp1(1) = j_tmp
      shp1(2) = i_tmp
      mat_A = reshape(a_plus,shp1)

      shp1(1) = j_tmp
      shp1(2) = l_tmp
      mat_B = reshape(t_plus,shp1)

      Call Zgemm('T','N',i_tmp,l_tmp,j_tmp,alpha,mat_A,j_tmp, &
                 mat_B,j_tmp,beta,mat_C,i_tmp)

      shp1i(1) = ij_tmp
      shp1i(2) = n_k_points
      shp1i(3) = n_a
      shp1i(4) = n_k_points
      rlt_S = reshape(mat_C,shp1i)

      shp1(1) = j_tmp
      shp1(2) = i_tmp
      mat_A = reshape(a_minus,shp1)

      shp1(1) = j_tmp
      shp1(2) = l_tmp
      mat_B = reshape(t_minus,shp1)

      Call Zgemm('T','N',i_tmp,l_tmp,j_tmp,alpha,mat_A,j_tmp, &
                 mat_B,j_tmp,beta,mat_C,i_tmp)

      shp1i(1) = ij_tmp
      shp1i(2) = n_k_points
      shp1i(3) = n_a
      shp1i(4) = n_k_points
      rlt_A = reshape(mat_C,shp1i)

      do k2 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k2,k4)

        do ij_run = 1, ij_tmp

          Call CC_3d_decode(ij_run,iii,jjj,CC_n_occ,2)

          w_tmp(iii,jjj,:,k2,:,k_run) = rlt_S(ij_run,k2,:,:) + rlt_A(ij_run,k2,:,:)
          w_tmp(jjj,iii,:,k4,:,k_run) = rlt_S(ij_run,k2,:,:) - rlt_A(ij_run,k2,:,:)

        end do
      end do

    end do

    CC_w_d_d = CC_w_d_d + w_tmp
 
    Deallocate(w_tmp,t_plus,t_minus,a_plus,a_minus,rlt_S,rlt_A,mat_A,mat_B,mat_C,aux_tmp)
 
  end do

  End Subroutine CC_3d_add_a2w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_add_b2w()

! k1 = ka, k2 = ki, k3 = kb, k4 = kj, k5 = kc, k6 = kd

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k5,k6,k_run,n_k1,k_start,k_end,k_13
  Integer :: n_ab,ab_start,ab_end,n_cd,cd_start,cd_end
  Integer :: n_a,a_start,a_end,n_c,c_start,c_end
  Integer :: iii,jjj,aaa,bbb,i_tmp,j_tmp
  Integer :: c1,c2,c3
  Integer :: i_task,target_id,errnum
  Integer , dimension(2) :: shp1
  Integer , dimension(4) :: shp1i
  Integer , dimension(3) :: shp2i

  Integer (kind=8) :: s_tmp

  Integer , dimension(CC_mem_k1(CC_mpi_gid+1)) :: n_calc

  Double complex , dimension(:,:,:) , allocatable :: CC_b_minus_nd,CC_b_plus_nd
  Double complex , dimension(:,:,:) , allocatable :: CC_b_minus_d,CC_b_plus_d
  Double complex , dimension(:,:,:) , allocatable :: t_plus_nd,t_minus_nd
  Double complex , dimension(:,:,:) , allocatable :: t_plus_d,t_minus_d
  Double complex , dimension(:,:) , allocatable :: mat_A,mat_B,mat_C,mat_D,mat_E,mat_F
  Double complex , dimension(:,:) , allocatable :: rlt_S,rlt_A
  Double complex , dimension(:,:,:,:,:,:) , allocatable :: wnd_tmp,wd_tmp

  Double complex :: alpha,beta

  !print*,'b_aux'
  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)

  n_calc = 0
  k_13 = k_start - 1
  do k_run = 1, n_k1
    k_13 = k_13 + 1
    do k2 = 1, n_k_points
      Call CC_3d_k_minus(k_13,k2,k4)
      do iii = 1, CC_n_occ
        do jjj = 1, CC_n_occ
          Call CC_3d_code(c1,k2,iii,CC_n_occ,1)
          Call CC_3d_code(c2,k4,jjj,CC_n_occ,1)
          if (c1.le.c2) then
            n_calc(k_run) = n_calc(k_run) + 1
          end if
        end do
      end do
    end do
  end do

  ! For w_nd
  n_cd = CC_mem_ab_D(CC_mpi_did+1,3)
  cd_start = CC_index_ab_D(CC_mpi_did+1,3)
  cd_end = cd_start - 1 + n_cd

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  do i_task = 1, CC_mpi_domain_size

    n_ab = CC_mem_ab_D(i_task,3)
    ab_start = CC_index_ab_D(i_task,3)
    ab_end = ab_start - 1 + n_ab

    n_c = CC_mem_aa_D(i_task)
    c_start = CC_index_aa_D(i_task)
    c_end = c_start - 1 + n_c

    Allocate(wnd_tmp(CC_n_occ,CC_n_occ,n_ab,n_k_points,n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'wnd_tmp in CC')

    wnd_tmp = 0.0D0

    Allocate(wd_tmp(CC_n_occ,CC_n_occ,n_c,n_k_points,n_k_points,n_k1),stat=errnum)
    Call check_allocation(errnum,'wnd_tmp in CC')

    wd_tmp = 0.0D0

    k_13 = k_start - 1
    do k_run = 1, n_k1
      k_13 = k_13 + 1

      ! get tau(+) and tau(-)
      Allocate(t_plus_nd(n_cd,n_k_points,n_calc(k_run)),stat=errnum)
      Call check_allocation(errnum,'t_plus_nd in CC')

      Allocate(t_minus_nd(n_cd,n_k_points,n_calc(k_run)),stat=errnum)
      Call check_allocation(errnum,'t_minus_nd in CC')

      c3 = 0
      do k2 = 1, n_k_points
        Call CC_3d_k_minus(k_13,k2,k4)
        do iii = 1, CC_n_occ
          do jjj = 1, CC_n_occ
            Call CC_3d_code(c1,k2,iii,CC_n_occ,1)
            Call CC_3d_code(c2,k4,jjj,CC_n_occ,1)
            if (c1.le.c2) then
              c3 = c3 + 1
              t_plus_nd(:,:,c3) = 0.5D0 * (CC_t_d_nd(iii,jjj,:,k2,:,k_run) &
                                           + CC_t_d_nd(jjj,iii,:,k4,:,k_run))
              t_minus_nd(:,:,c3) = 0.5D0 * (CC_t_d_nd(iii,jjj,:,k2,:,k_run) &
                                            - CC_t_d_nd(jjj,iii,:,k4,:,k_run))
            end if
          end do
        end do
      end do

      Allocate(t_plus_d(n_a,n_k_points,n_calc(k_run)),stat=errnum)
      Call check_allocation(errnum,'t_plus_d in CC')

      Allocate(t_minus_d(n_a,n_k_points,n_calc(k_run)),stat=errnum)
      Call check_allocation(errnum,'t_minus_d in CC')

      c3 = 0
      do k2 = 1, n_k_points
        Call CC_3d_k_minus(k_13,k2,k4)
        do iii = 1, CC_n_occ
          do jjj = 1, CC_n_occ
            Call CC_3d_code(c1,k2,iii,CC_n_occ,1)
            Call CC_3d_code(c2,k4,jjj,CC_n_occ,1)
            if (c1.le.c2) then
              c3 = c3 + 1
              t_plus_d(:,:,c3) = 0.5D0 * (CC_t_d_d(iii,jjj,:,k2,:,k_run) &
                                          + CC_t_d_d(jjj,iii,:,k4,:,k_run))
              t_minus_d(:,:,c3) = 0.5D0 * (CC_t_d_d(iii,jjj,:,k2,:,k_run) &
                                           - CC_t_d_d(jjj,iii,:,k4,:,k_run))
            end if
          end do
        end do
      end do

      ! For w_nd
      Allocate(CC_b_plus_nd(n_cd,n_k_points,n_ab),stat=errnum)
      Call check_allocation(errnum,'b_minus in CC')

      Allocate(CC_b_minus_nd(n_cd,n_k_points,n_ab),stat=errnum)
      Call check_allocation(errnum,'b_minus in CC')

      Allocate(CC_b_plus_d(n_a,n_k_points,n_ab),stat=errnum)
      Call check_allocation(errnum,'b_minus in CC')

      Allocate(CC_b_minus_d(n_a,n_k_points,n_ab),stat=errnum)
      Call check_allocation(errnum,'b_minus in CC')

      i_tmp = n_cd * n_k_points

      Allocate(mat_A(i_tmp,n_calc(k_run)),stat=errnum)
      Call check_allocation(errnum,'mat_A in CC')

      Allocate(mat_B(i_tmp,n_ab),stat=errnum)
      Call check_allocation(errnum,'mat_B in CC')

      Allocate(mat_C(n_calc(k_run),n_ab),stat=errnum)
      Call check_allocation(errnum,'mat_C in CC')

      Allocate(rlt_S(n_calc(k_run),n_ab),stat=errnum)
      Call check_allocation(errnum,'rlt_S in CC')

      Allocate(rlt_A(n_calc(k_run),n_ab),stat=errnum)
      Call check_allocation(errnum,'rlt_A in CC')

      j_tmp = n_a * n_k_points

      Allocate(mat_D(j_tmp,n_calc(k_run)),stat=errnum)
      Call check_allocation(errnum,'mat_D in CC')

      Allocate(mat_E(j_tmp,n_ab),stat=errnum)
      Call check_allocation(errnum,'mat_E in CC')

      Allocate(mat_F(n_calc(k_run),n_ab),stat=errnum)
      Call check_allocation(errnum,'mat_F in CC')

      do k1 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k1,k3)
        Call CC_3d_aux_b_nd(i_task,k1,k3,k_run,CC_b_plus_nd,CC_b_minus_nd, &
                                               CC_b_plus_d,CC_b_minus_d)

        shp1(1) = i_tmp
        shp1(2) = n_calc(k_run)
        mat_A = reshape(t_plus_nd,shp1)

        shp1(1) = i_tmp
        shp1(2) = n_ab
        mat_B = reshape(CC_b_plus_nd,shp1)

        alpha = 2.0D0
        beta = 0.0D0

        Call Zgemm('T','N',n_calc(k_run),n_ab,i_tmp,alpha,mat_A,i_tmp, &
                   mat_B,i_tmp,beta,mat_C,n_calc(k_run))

        rlt_S = mat_C

        shp1(1) = i_tmp
        shp1(2) = n_calc(k_run)
        mat_A = reshape(t_minus_nd,shp1)

        shp1(1) = i_tmp
        shp1(2) = n_ab
        mat_B = reshape(CC_b_minus_nd,shp1)

        alpha = 2.0D0
        beta = 0.0D0

        Call Zgemm('T','N',n_calc(k_run),n_ab,i_tmp,alpha,mat_A,i_tmp, &
                   mat_B,i_tmp,beta,mat_C,n_calc(k_run))

        rlt_A = mat_C

        shp1(1) = j_tmp
        shp1(2) = n_calc(k_run)
        mat_D = reshape(t_plus_d,shp1)

        shp1(1) = j_tmp
        shp1(2) = n_ab
        mat_E = reshape(CC_b_plus_d,shp1)

        alpha = 1.0D0
        beta = 0.0D0

        Call Zgemm('T','N',n_calc(k_run),n_ab,j_tmp,alpha,mat_D,j_tmp, &
                   mat_E,j_tmp,beta,mat_F,n_calc(k_run))

        rlt_S = rlt_S + mat_F

        shp1(1) = j_tmp
        shp1(2) = n_calc(k_run)
        mat_D = reshape(t_minus_d,shp1)

        shp1(1) = j_tmp
        shp1(2) = n_ab
        mat_E = reshape(CC_b_minus_d,shp1)

        alpha = 1.0D0
        beta = 0.0D0

        Call Zgemm('T','N',n_calc(k_run),n_ab,j_tmp,alpha,mat_D,j_tmp, &
                   mat_E,j_tmp,beta,mat_F,n_calc(k_run))

        rlt_A = rlt_A + mat_F

        c3 = 0
        do k2 = 1, n_k_points
          Call CC_3d_k_minus(k_13,k2,k4)
          do iii = 1, CC_n_occ
            do jjj = 1, CC_n_occ
              Call CC_3d_code(c1,k2,iii,CC_n_occ,1)
              Call CC_3d_code(c2,k4,jjj,CC_n_occ,1)
              if (c1.le.c2) then
                c3 = c3 + 1
                wnd_tmp(iii,jjj,:,k2,k1,k_run) = wnd_tmp(iii,jjj,:,k2,k1,k_run) &
                                               + rlt_S(c3,:) + rlt_A(c3,:)

                if (c1.ne.c2) then
                  wnd_tmp(jjj,iii,:,k4,k1,k_run) = wnd_tmp(jjj,iii,:,k4,k1,k_run) &
                                                 + rlt_S(c3,:) - rlt_A(c3,:)
                end if
              end if
            end do
          end do
        end do
        
      end do

      Deallocate(mat_A,mat_B,mat_C,mat_D,mat_E,mat_F,rlt_S,rlt_A)
      Deallocate(CC_b_plus_nd,CC_b_minus_nd,CC_b_plus_d,CC_b_minus_d)

      ! For w_d
      Allocate(CC_b_plus_nd(n_cd,n_k_points,n_c),stat=errnum)
      Call check_allocation(errnum,'b_minus in CC')

      Allocate(CC_b_minus_nd(n_cd,n_k_points,n_c),stat=errnum)
      Call check_allocation(errnum,'b_minus in CC')

      Allocate(CC_b_plus_d(n_a,n_k_points,n_c),stat=errnum)
      Call check_allocation(errnum,'b_minus in CC')

      Allocate(CC_b_minus_d(n_a,n_k_points,n_c),stat=errnum)
      Call check_allocation(errnum,'b_minus in CC')

      i_tmp = n_cd * n_k_points

      Allocate(mat_A(i_tmp,n_calc(k_run)),stat=errnum)
      Call check_allocation(errnum,'mat_A in CC')

      Allocate(mat_B(i_tmp,n_c),stat=errnum)
      Call check_allocation(errnum,'mat_B in CC')

      Allocate(mat_C(n_calc(k_run),n_c),stat=errnum)
      Call check_allocation(errnum,'mat_C in CC')

      Allocate(rlt_S(n_calc(k_run),n_c),stat=errnum)
      Call check_allocation(errnum,'rlt_S in CC')

      Allocate(rlt_A(n_calc(k_run),n_c),stat=errnum)
      Call check_allocation(errnum,'rlt_A in CC')

      j_tmp = n_a * n_k_points

      Allocate(mat_D(j_tmp,n_calc(k_run)),stat=errnum)
      Call check_allocation(errnum,'mat_D in CC')

      Allocate(mat_E(j_tmp,n_c),stat=errnum)
      Call check_allocation(errnum,'mat_E in CC')

      Allocate(mat_F(n_calc(k_run),n_c),stat=errnum)
      Call check_allocation(errnum,'mat_F in CC')

      do k1 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k1,k3)
        Call CC_3d_aux_b_d(i_task,k1,k3,k_run,CC_b_plus_nd,CC_b_minus_nd, &
                                              CC_b_plus_d,CC_b_minus_d)

        shp1(1) = i_tmp
        shp1(2) = n_calc(k_run)
        mat_A = reshape(t_plus_nd,shp1)

        shp1(1) = i_tmp
        shp1(2) = n_c
        mat_B = reshape(CC_b_plus_nd,shp1)

        alpha = 2.0D0
        beta = 0.0D0

        Call Zgemm('T','N',n_calc(k_run),n_c,i_tmp,alpha,mat_A,i_tmp, &
                   mat_B,i_tmp,beta,mat_C,n_calc(k_run))

        rlt_S = mat_C

        shp1(1) = i_tmp
        shp1(2) = n_calc(k_run)
        mat_A = reshape(t_minus_nd,shp1)

        shp1(1) = i_tmp
        shp1(2) = n_c
        mat_B = reshape(CC_b_minus_nd,shp1)

        alpha = 2.0D0
        beta = 0.0D0

        Call Zgemm('T','N',n_calc(k_run),n_c,i_tmp,alpha,mat_A,i_tmp, &
                   mat_B,i_tmp,beta,mat_C,n_calc(k_run))

        rlt_A = mat_C

        shp1(1) = j_tmp
        shp1(2) = n_calc(k_run)
        mat_D = reshape(t_plus_d,shp1)

        shp1(1) = j_tmp
        shp1(2) = n_c
        mat_E = reshape(CC_b_plus_d,shp1)

        alpha = 1.0D0
        beta = 0.0D0

        Call Zgemm('T','N',n_calc(k_run),n_c,j_tmp,alpha,mat_D,j_tmp, &
                   mat_E,j_tmp,beta,mat_F,n_calc(k_run))

        rlt_S = rlt_S + mat_F

        shp1(1) = j_tmp
        shp1(2) = n_calc(k_run)
        mat_D = reshape(t_minus_d,shp1)

        shp1(1) = j_tmp
        shp1(2) = n_c
        mat_E = reshape(CC_b_minus_d,shp1)

        alpha = 1.0D0
        beta = 0.0D0

        Call Zgemm('T','N',n_calc(k_run),n_c,j_tmp,alpha,mat_D,j_tmp, &
                   mat_E,j_tmp,beta,mat_F,n_calc(k_run))

        rlt_A = rlt_A + mat_F

        c3 = 0
        do k2 = 1, n_k_points
          Call CC_3d_k_minus(k_13,k2,k4)
          do iii = 1, CC_n_occ
            do jjj = 1, CC_n_occ
              Call CC_3d_code(c1,k2,iii,CC_n_occ,1)
              Call CC_3d_code(c2,k4,jjj,CC_n_occ,1)
              if (c1.le.c2) then
                c3 = c3 + 1
                wd_tmp(iii,jjj,:,k2,k1,k_run) = wd_tmp(iii,jjj,:,k2,k1,k_run) &
                                              + rlt_S(c3,:) + rlt_A(c3,:)

                if (c1.ne.c2) then
                  wd_tmp(jjj,iii,:,k4,k1,k_run) = wd_tmp(jjj,iii,:,k4,k1,k_run) &
                                                + rlt_S(c3,:) - rlt_A(c3,:)
                end if
              end if
            end do
          end do
        end do
        
      end do

      Deallocate(mat_A,mat_B,mat_C,mat_D,mat_E,mat_F,rlt_S,rlt_A)
      Deallocate(CC_b_plus_nd,CC_b_minus_nd,CC_b_plus_d,CC_b_minus_d)
      Deallocate(t_plus_nd,t_minus_nd,t_plus_d,t_minus_d)
 
    end do

    target_id = i_task - 1

    i_tmp = CC_n_occ * CC_n_occ * n_k_points * n_k_points * n_k1

    s_tmp = Int(i_tmp,8) * Int(n_ab,8)
    Call CC_mpi_complex_reduce(s_tmp,wnd_tmp,target_id,CC_mpi_comm_domain)

    s_tmp = Int(i_tmp,8) * Int(n_c,8)
    Call CC_mpi_complex_reduce(s_tmp,wd_tmp,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      CC_w_d_nd = CC_w_d_nd + wnd_tmp
      CC_w_d_d = CC_w_d_d + wd_tmp
    end if

    Deallocate(wnd_tmp,wd_tmp)

  end do

  End Subroutine CC_3d_add_b2w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_add_ct2w()

  Use CC_3d

  Implicit None

  Integer :: k1,k2,k3,k4,k_13,k_run,n_k1,k_start,k_end,k_tmp
  Integer :: k_pattern,k_recv_start,n_k_recv
  Integer :: n_ab,ab_start,ab_end,n_cd,cd_start,cd_end
  Integer :: n_a,a_start,a_end,n_c,c_start,c_end
  Integer :: iii,jjj,aaa,bbb,a_run,code_ab,code_ij
  Integer :: c1,c2,c3,i_tmp
  Integer :: i_task,source_id,i_grp_task,i_task_run,errnum

  Integer :: req1,req2,i_send,i_recv
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Integer (kind=8) :: s_tmp

  Integer , dimension(CC_mem_k1(CC_mpi_gid+1)) :: n_calc

  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: tc_tmp,tc_send,tc_recv
  Double complex , dimension(:,:,:,:,:,:,:) , allocatable :: tc_tmp2

  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)
  k_end = k_start - 1 + n_k1

  n_ab = CC_mem_ab_D(CC_mpi_did+1,3)
  ab_start = CC_index_ab_D(CC_mpi_did+1,3)
  ab_end = ab_start - 1 + n_ab

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  i_grp_task = CC_mpi_gid
  i_send = CC_mpi_gid + 1
  i_recv = CC_mpi_gid + 1

  Allocate(tc_tmp(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                  CC_mem_aa_D(CC_mpi_did+1),n_k_points,n_k1),stat=errnum)
  Call check_allocation(errnum,'tc_tmp in CC')

  tc_tmp = CC_ten_con

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

    n_k_recv = CC_mem_k1(i_grp_task)
    k_recv_start = CC_index_k1(i_grp_task)

    if (i_recv.ne.CC_mpi_gid+1) then

      Allocate(tc_send(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                       CC_mem_aa_D(CC_mpi_did+1),n_k_points,n_k1),stat=errnum)
      Call check_allocation(errnum,'tc_send in CC')

      Allocate(tc_recv(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                       CC_mem_aa_D(CC_mpi_did+1),n_k_points,CC_mem_k1(i_recv)),stat=errnum)
      Call check_allocation(errnum,'tc_recv in CC')

      tc_send = CC_ten_con

      i_tmp = CC_n_occ * CC_n_occ * CC_mem_aa_D(CC_mpi_did+1) * CC_n_vir
      s_tmp = Int(i_tmp,8) * Int(n_k1,8) * Int(n_k_points**2,8)

      Call CC_mpi_complex_isend(s_tmp,tc_send,i_send-1,100,req1,CC_mpi_comm_group)
 
      s_tmp = Int(i_tmp,8) * Int(CC_mem_k1(i_recv),8) * Int(n_k_points**2,8)
      Call CC_mpi_complex_irecv(s_tmp,tc_recv,i_recv-1,100,req2,CC_mpi_comm_group)

    end if

    n_a = CC_mem_aa_D(CC_mpi_did+1)
    a_start = CC_index_aa_D(CC_mpi_did+1)
    a_end = a_start - 1 + n_a

    ! For w_d
    k1 = k_recv_start - 1
    do k_run = 1, n_k_recv

      k1 = k1 + 1
      do k3 = 1, n_k_points

        Call CC_3d_k_plus(k_tmp,k1,k3)
        k_13 = k_tmp - k_start + 1

        if ((k_13.ge.1).and.(k_13.le.n_k1)) then
          do k2 = 1, n_k_points
            do a_run = 1, n_a
              Call CC_3d_k_minus(k_tmp,k2,k4)
              do jjj = 1, CC_n_occ
                aaa = a_start - 1 + a_run
                CC_w_d_d(:,jjj,a_run,k2,k1,k_13) = CC_w_d_d(:,jjj,a_run,k2,k1,k_13) &
                                                 + tc_tmp(jjj,aaa,k4,:,a_run,k2,k_run)
              end do
              CC_w_d_d(:,:,a_run,k4,k3,k_13) = CC_w_d_d(:,:,a_run,k4,k3,k_13) &
                                             + tc_tmp(:,aaa,k4,:,a_run,k2,k_run)
            end do
          end do
        end if
      end do
    end do

    ! For w_nd
    do i_task = 1, CC_mpi_domain_size

      n_a = CC_mem_aa_D(i_task)
      a_start = CC_index_aa_D(i_task)
      a_end = a_start - 1 + n_a

      Allocate(tc_tmp2(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                       n_a,n_k_points,n_k_recv),stat=errnum)
      Call check_allocation(errnum,'tc_tmp2 in CC')

      source_id = i_task - 1
      i_tmp = CC_n_occ * CC_n_occ * n_a * CC_n_vir
      s_tmp = Int(i_tmp,8) * Int(n_k_recv,8) * Int(n_k_points**2,8)

      if (source_id.eq.CC_mpi_did) then
        tc_tmp2 = tc_tmp
      else
        tc_tmp2 = 0.0D0
      end if

      Call  CC_mpi_complex_bcast(s_tmp, tc_tmp2, source_id, CC_mpi_comm_domain)

      k1 = k_recv_start - 1
      do k_run = 1, n_k_recv

        k1 = k1 + 1
        do k3 = 1, n_k_points

          Call CC_3d_k_plus(k_13,k1,k3)

          k_13 = k_13 - k_start + 1

          if ((k_13.ge.1).and.(k_13.le.n_k1)) then

            do k2 = 1, n_k_points

              Call CC_3d_determine_k_pattern(k_pattern,k1,k2,k3,k4,4)

              aaa = a_start - 1
              do a_run = 1, n_a
                aaa = aaa + 1

                do bbb = aaa+1, CC_n_vir
                  Call CC_3d_code(code_ab,aaa,bbb,CC_n_vir,3)

                  code_ab = code_ab - ab_start + 1
                  if ((code_ab.ge.1).and.(code_ab.le.n_ab)) then
                    do jjj = 1, CC_n_occ
                      CC_w_d_nd(:,jjj,code_ab,k2,k1,k_13) = CC_w_d_nd(:,jjj,code_ab,k2,k1,k_13) &
                                                          + tc_tmp2(jjj,bbb,k4,:,a_run,k2,k_run)
                    end do
                  end if
                end do

                do bbb = 1, aaa - 1
                  Call CC_3d_code(code_ab,bbb,aaa,CC_n_vir,3)

                  code_ab = code_ab - ab_start + 1
                  if ((code_ab.ge.1).and.(code_ab.le.n_ab)) then
                    CC_w_d_nd(:,:,code_ab,k4,k3,k_13) = CC_w_d_nd(:,:,code_ab,k4,k3,k_13) &
                                                      + tc_tmp2(:,bbb,k4,:,a_run,k2,k_run)
                  end if
                end do
              end do
            end do
          end if
        end do
      end do

      Deallocate(tc_tmp2)
    end do

    Deallocate(tc_tmp)

    if (i_recv.ne.CC_mpi_gid+1) then

      Call MPI_WAIT(req1,stat1,errnum) 
      Call MPI_WAIT(req2,stat2,errnum)

      Allocate(tc_tmp(CC_n_occ,CC_n_vir,n_k_points,CC_n_occ, &
                      CC_mem_aa_D(CC_mpi_did+1),n_k_points,CC_mem_k1(i_recv)),stat=errnum)

      Call check_allocation(errnum,'tc_tmp in CC')

      tc_tmp = tc_recv

      Deallocate(tc_recv,tc_send)

    end if
  end do

  End Subroutine CC_3d_add_ct2w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

