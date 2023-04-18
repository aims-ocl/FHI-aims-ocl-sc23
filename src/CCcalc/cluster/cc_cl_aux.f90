!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_aux_h_g()

! (iii,kkk)

  Use CC_cl

  Implicit None

  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,i_tmp,j_tmp
  Integer :: n_c,c_start,c_end,n_d,d_start,d_end,n_l,l_start,l_end
  Integer :: c_run,d_run,l_run

  Integer (kind=8) :: s_tmp
  Integer :: errnum

  Double precision :: alpha,beta
  Double precision , dimension(:,:,:,:) , allocatable :: tau,intl
  Double precision , dimension(:,:,:,:) , allocatable :: MatA,matB
  Double precision , dimension(:,:) , allocatable :: MatC

  ! Clear axiliary vectors
  CC_h_ik = 0.0D0
  CC_h_ca = 0.0D0
  CC_h_ck = 0.0D0
  CC_g_ik = 0.0D0
  CC_g_ca = 0.0D0

  alpha = 1.0D0
  beta = 0.0D0

  n_l = CC_mem_ii_G(CC_mpi_gid+1)
  l_start = CC_index_ii_G(CC_mpi_gid+1)
  l_end = l_start - 1 + n_l

  n_d = CC_mem_aa_D(CC_mpi_did+1)
  d_start = CC_index_aa_D(CC_mpi_did+1)
  d_end = d_start - 1 + n_d

  n_c = CC_mem_aa_G(CC_mpi_gid+1)
  c_start = CC_index_aa_G(CC_mpi_gid+1)
  c_end = c_start - 1 + n_c

  if (n_l.gt.0) then
    ! Calculate tau and intl_iajb (l,d,k,c)
    Allocate(tau(n_l,n_d,CC_n_occ,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'tau in CC')

    Allocate(intl(n_l,n_d,CC_n_occ,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'intl in CC')

    !$OMP PARALLEL Default(Shared) &
    !$OMP Private(l_run,d_run,kkk,lll,ccc,ddd)
    !$OMP DO
    do ccc = 1, CC_n_vir
      do kkk = 1, CC_n_occ
        do d_run = 1, n_d
          do l_run = 1, n_l
            lll = l_start - 1 + l_run
            ddd = d_start - 1 + d_run
            tau(l_run,d_run,kkk,ccc) = CC_t_d(lll,kkk,d_run,ccc) &
                                     + CC_t_s(kkk,ccc) * CC_t_s(lll,ddd)

            intl(l_run,d_run,kkk,ccc) = 2.0D0 * CC_intl_iajb(lll,kkk,d_run,ccc) &
                                              - CC_intl_iajb(kkk,lll,d_run,ccc)
          end do
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! h_ik(iii,kkk)
    Allocate(MatA(n_l,n_d,CC_n_vir,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'matA in CC')

    Allocate(MatB(n_l,n_d,CC_n_vir,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'matB in CC')

    Allocate(MatC(CC_n_occ,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'matC in CC')

    do iii = 1, CC_n_occ
      MatA(:,:,:,iii) = tau(:,:,iii,:)
    end do

    do kkk = 1, CC_n_occ
      MatB(:,:,:,kkk) = intl(:,:,kkk,:)
    end do

    i_tmp = n_l * n_d * CC_n_vir

    Call Dgemm('T','N',CC_n_occ,CC_n_occ,i_tmp,alpha,MatA,i_tmp, &
               MatB,i_tmp,beta,MatC,CC_n_occ)

    CC_h_ik = MatC

    Deallocate(MatA,MatB,MatC)

    ! h_ca(ccc,aaa)
    Allocate(MatC(CC_n_vir,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'matC in CC')

    i_tmp = n_d * n_l * CC_n_occ

    Call Dgemm('T','N',CC_n_vir,CC_n_vir,i_tmp,alpha,intl,i_tmp, &
               tau,i_tmp,beta,MatC,CC_n_vir)

    CC_h_ca = - MatC

    Deallocate(MatC)

    ! h_ck(kkk,ccc)
    i_tmp = n_l * n_d
    j_tmp = CC_n_occ * CC_n_vir

    Allocate(MatB(n_l,n_d,1,1),stat=errnum)
    Call check_allocation(errnum,'matB in CC')

    Allocate(MatC(CC_n_occ,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'matC in CC')

    MatB(:,:,1,1) = CC_t_s(l_start:l_end,d_start:d_end)

    Call Dgemm('T','N',j_tmp,1,i_tmp,alpha,intl,i_tmp, &
               MatB,i_tmp,beta,MatC,j_tmp)

    CC_h_ck = MatC

    Deallocate(MatB,MatC)
    Deallocate(intl,tau)

    ! g_ik(iii,kkk)
    i_tmp = n_l * n_d
    j_tmp = CC_n_occ * CC_n_occ

    ! (l,c,i,k)
    Allocate(MatA(n_l,n_d,CC_n_occ,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'matA in CC')

    Allocate(MatB(n_l,n_d,1,1),stat=errnum)
    Call check_allocation(errnum,'matB in CC')

    Allocate(MatC(CC_n_occ,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'matC in CC')

    !$OMP PARALLEL Default(Shared) &
    !$OMP Private(l_run,kkk,iii,lll,ccc)
    !$OMP DO
    do kkk = 1, CC_n_occ
      do iii = 1, CC_n_occ
        do ccc = 1, n_d
          do l_run = 1, n_l

            lll = l_start - 1 + l_run
            MatA(l_run,ccc,iii,kkk) = 2.0D0 * CC_intl_kilc(kkk,ccc,iii,lll) &
                                            - CC_intl_kilc(lll,ccc,iii,kkk)
          end do
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    MatB(:,:,1,1) = CC_t_s(l_start:l_end,d_start:d_end)

    Call Dgemm('T','N',j_tmp,1,i_tmp,alpha,MatA,i_tmp, &
               MatB,i_tmp,beta,MatC,j_tmp)

    CC_g_ik = MatC

    Deallocate(MatA,MatB,MatC)

  end if

  ! g_ca(ccc,aaa)
  ! first part
  i_tmp = CC_n_vir * CC_n_occ

  ! (k,d,c,a)
  Allocate(MatA(CC_n_occ,CC_n_vir,n_d,1),stat=errnum)
  Call check_allocation(errnum,'matA in CC')

  Allocate(MatC(n_d,1),stat=errnum)
  Call check_allocation(errnum,'matC in CC')

  do aaa = c_start, c_end

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(ccc,ddd)
    !$OMP DO
    do ccc = 1, n_d
      do ddd = 1, CC_n_vir
        MatA(:,ddd,ccc,1) = CC_intl_ackd(:,ccc,ddd,aaa)
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL 

    alpha = 2.0D0
    beta = 0.0D0

    Call Dgemm('T','N',n_d,1,i_tmp,alpha,MatA,i_tmp, &
               CC_t_s,i_tmp,beta,MatC,n_d)

    CC_g_ca(d_start:d_end,aaa) = CC_g_ca(d_start:d_end,aaa) + MatC(:,1)

  end do

  Deallocate(MatA,MatC)

  ! second part
  i_tmp = n_d * CC_n_occ
  j_tmp = n_c * CC_n_vir

  Allocate(MatA(CC_n_occ,n_d,1,1),stat=errnum)
  Call check_allocation(errnum,'matA in CC')

  Allocate(MatC(CC_n_vir,n_c),stat=errnum)
  Call check_allocation(errnum,'matC in CC')

  MatA(:,:,1,1) = CC_t_s(:,d_start:d_end)

  alpha = -1.0D0
  beta = 0.0D0

  Call Dgemm('T','N',j_tmp,1,i_tmp,alpha,CC_intl_ackd(:,:,:,c_start:c_end),i_tmp, &
             MatA,i_tmp,beta,MatC,j_tmp)

  CC_g_ca(:,c_start:c_end) = CC_g_ca(:,c_start:c_end) + MatC

  Deallocate(MatA,MatC)

  ! Synchronize vectors
  s_tmp = Int(CC_n_occ,8) * Int(CC_n_occ,8)
  Call CC_mpi_allreduce(s_tmp, CC_h_ik, MPI_COMM_WORLD) 
  Call CC_mpi_allreduce(s_tmp, CC_g_ik, MPI_COMM_WORLD)

  s_tmp = Int(CC_n_vir,8) * Int(CC_n_vir,8)
  Call CC_mpi_allreduce(s_tmp, CC_h_ca, MPI_COMM_WORLD) 
  Call CC_mpi_allreduce(s_tmp, CC_g_ca, MPI_COMM_WORLD)

  s_tmp = Int(CC_n_occ,8) * Int(CC_n_vir,8)
  Call CC_mpi_allreduce(s_tmp, CC_h_ck, MPI_COMM_WORLD)

  CC_g_ik = CC_g_ik + CC_h_ik
  CC_g_ca = CC_g_ca + CC_h_ca

  End Subroutine CC_cl_aux_h_g

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_aux_j()

! Caluclate intermediate (2 * j - k) (kkk,ccc,iii,aaa) 

  Use CC_cl

  Implicit None

  Integer :: n_c,c_start,c_end,n_a,a_start,a_end,n_d,d_start,d_end
  Integer :: i_run,k_run,a_run,c_run,d_run
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,i_tmp,j_tmp,k_tmp,l_tmp
  Integer (kind = 8) :: s_tmp
  Integer :: errnum

  Integer :: i_task,i_task_run
  Integer :: req1,req2,i_send,i_recv
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Double precision :: alpha,beta
  Double precision , dimension(:,:,:,:) , allocatable :: matA,matB,matC
  Double precision , dimension(:,:,:,:) , allocatable :: t_tmp,v_tmp,j_send,j_recv,aux_tmp

  n_c = CC_mem_aa_D(CC_mpi_did+1)
  c_start = CC_index_aa_D(CC_mpi_did+1)
  c_end = c_start - 1 + n_c

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  alpha = 1.0D0
  beta = 0.0D0

  if (.not.(allocated(CC_j_aux))) then
    ! (k,c,i,a)
    Allocate(CC_j_aux(CC_n_occ,n_c,CC_n_occ,n_a),stat=errnum)
  end if

  ! First term 2(ai|kc) - (ki|ac)
  !$OMP PARALLEL Default(shared) &
  !$OMP Private(a_run,iii,ccc,aaa,kkk)
  !$OMP DO
  do a_run = 1, n_a
    do iii = 1, CC_n_occ
      do ccc = 1, n_c
        do kkk = 1, CC_n_occ

          aaa = a_start - 1 + a_run
          CC_j_aux(kkk,ccc,iii,a_run) = 2.0D0 * CC_intl_iajb(kkk,iii,ccc,aaa) &
                                              - CC_intl_kiac(kkk,ccc,iii,aaa)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Second term
  !(l,k,c,i)
  Allocate(matA(CC_n_occ,CC_n_occ,n_c,CC_n_occ), stat=errnum)
  Call check_allocation(errnum,'matA in CC')

  Allocate(matB(CC_n_occ,n_a,1,1), stat=errnum)
  Call check_allocation(errnum,'matB in CC')

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(lll,iii,ccc,kkk)
  !$OMP DO
  do iii = 1, CC_n_occ
    do ccc = 1, n_c
      do kkk = 1, CC_n_occ
        do lll = 1, CC_n_occ
          matA(lll,kkk,ccc,iii) = 2.0D0 * CC_intl_kilc(lll,ccc,iii,kkk) &
                                        - CC_intl_kilc(kkk,ccc,iii,lll)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  matB(:,:,1,1) = CC_t_s(:,a_start:a_end)

  alpha = -1.0D0
  beta = 1.0D0

  j_tmp = CC_n_occ * n_c * CC_n_occ

  Call Dgemm('T','N',j_tmp,n_a,CC_n_occ,alpha,matA,CC_n_occ, &
             matB,CC_n_occ,beta,CC_j_aux,j_tmp)

  Deallocate(matA,matB)

  ! Third term second part
  i_tmp = CC_n_occ * n_c

  alpha = -1.0D0
  beta = 1.0D0

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(a_run,aaa)
  !$OMP DO Schedule(Dynamic)
  do a_run = 1, n_a

    aaa = a_start - 1 + a_run

    Call Dgemm('N','T',i_tmp,CC_n_occ,CC_n_vir,alpha,CC_intl_ackd(:,:,:,aaa), &
               i_tmp,CC_t_s,CC_n_occ,beta,CC_j_aux(:,:,:,a_run),i_tmp)

  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Third term first part & Forth term
  ! Get T
  ! t_tmp(l,d,i,a)
  Allocate(t_tmp(CC_n_occ,n_c,CC_n_occ,n_a),stat=errnum)
  Call check_allocation(errnum,'t_tmp in CC')

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(a_run,d_run,aaa,iii,ddd,lll)
  !$OMP DO
  do a_run = 1, n_a
    do iii = 1, CC_n_occ
      do d_run = 1, n_c
        do lll = 1, CC_n_occ
          aaa = a_start - 1 + a_run
          ddd = c_start - 1 + d_run
          t_tmp(lll,d_run,iii,a_run) = CC_t_d(lll,iii,d_run,aaa)         &
                                     - 0.5D0 * CC_t_d(iii,lll,d_run,aaa) &
                                     - CC_t_s(iii,ddd) * CC_t_s(lll,aaa)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Get V
  ! v_tmp(l,d,k,c)
  Allocate(v_tmp(CC_n_occ,n_c,CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'v_tmp in CC')

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(ccc,kkk,ddd,lll)
  !$OMP DO
  do ccc = 1, CC_n_vir
    do kkk = 1, CC_n_occ
      do ddd = 1, n_c
        do lll = 1, CC_n_occ
          v_tmp(lll,ddd,kkk,ccc) = 2.0D0 * CC_intl_iajb(lll,kkk,ddd,ccc) &
                                         - CC_intl_iajb(kkk,lll,ddd,ccc)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  n_d = CC_mem_aa_D(CC_mpi_did+1)
  d_start = CC_index_aa_D(CC_mpi_did+1)
  d_end = d_start - 1 + n_d

  i_task = CC_mpi_did + 1
  i_send = CC_mpi_did + 1
  i_recv = CC_mpi_did + 1

  do i_task_run = 1, CC_mpi_domain_size

    i_recv = i_recv - 1
    if (i_recv.lt.1) then
      i_recv = CC_mpi_domain_size
    end if

    i_send = i_send + 1
    if (i_send.gt.CC_mpi_domain_size) then
      i_send = 1
    end if

    i_task = i_task + 1
    if (i_task.gt.CC_mpi_domain_size) then
      i_task = 1
    end if

    n_c = CC_mem_aa_D(i_task)
    c_start = CC_index_aa_D(i_task)
    c_end = c_start - 1 + n_c

    !(k,c,i,a)
    Allocate(aux_tmp(CC_n_occ,n_c,CC_n_occ,n_a),stat=errnum)
    Call check_allocation(errnum,'aux_tmp in CC')

    ! Third term first part (k,c,d)
    Allocate(matA(CC_n_occ,n_c,n_d,1),stat=errnum)
    Call check_allocation(errnum,'matA in CC')

    i_tmp = CC_n_occ * n_c

    alpha = 2.0D0
    beta = 0.0D0

    do a_run = 1, n_a

      aaa = a_start - 1 + a_run

      !$OMP PARALLEL Default(shared) &
      !$OMP Private(c_run,ccc,ddd)
      !$OMP DO
      do c_run = 1, n_c
        do ddd = 1, n_d
          ccc = c_start - 1 + c_run
          matA(:,c_run,ddd,1) = CC_intl_ackd(:,ddd,ccc,aaa)
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      Call Dgemm('N','T',i_tmp,CC_n_occ,n_d,alpha,matA,i_tmp, &
                 CC_t_s(:,d_start:d_end),CC_n_occ,beta,aux_tmp(:,:,:,a_run),i_tmp)

    end do

    Deallocate(matA)

    ! Forth term
    alpha = 1.0D0
    beta = 1.0D0

    i_tmp = CC_n_occ * n_d
    j_tmp = CC_n_occ * n_c
    l_tmp = CC_n_occ * n_a

    Call Dgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,v_tmp(:,:,:,c_start:c_end),i_tmp, &
               t_tmp,i_tmp,beta,aux_tmp,j_tmp)

    if (i_task_run.ne.1) then

      Call MPI_WAIT(req1,stat1,errnum)
      Call MPI_WAIT(req2,stat2,errnum)

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(aaa,iii,ccc,kkk)
      !$OMP DO
      do aaa = 1, n_a
        do iii = 1, CC_n_occ
          do ccc = 1, n_d
            do kkk = 1, CC_n_occ
              CC_j_aux(kkk,ccc,iii,aaa) = CC_j_aux(kkk,ccc,iii,aaa) &
                                        + j_recv(kkk,ccc,iii,aaa)
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      Deallocate(j_send,j_recv)

    end if

    if (i_task_run.ne.CC_mpi_domain_size) then

      Allocate(j_send(CC_n_occ,n_c,CC_n_occ,n_a),stat=errnum)
      Call check_allocation(errnum,'j_send in CC')

      j_send = aux_tmp

      Allocate(j_recv(CC_n_occ,n_d,CC_n_occ,n_a),stat=errnum)
      Call check_allocation(errnum,'j_recv in CC')

      i_tmp = CC_n_occ * CC_n_occ * n_a
      s_tmp = Int(i_tmp,8) * Int(n_c,8)

      Call CC_mpi_real_isend(s_tmp,j_send,i_send-1,100,req1,CC_mpi_comm_domain)

      s_tmp = Int(i_tmp,8) * Int(n_d,8)
      Call CC_mpi_real_irecv(s_tmp,j_recv,i_recv-1,100,req2,CC_mpi_comm_domain)

    else

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(aaa,iii,ccc,kkk)
      !$OMP DO
      do aaa = 1, n_a
        do iii = 1, CC_n_occ
          do ccc = 1, n_d
            do kkk = 1, CC_n_occ
              CC_j_aux(kkk,ccc,iii,aaa) = CC_j_aux(kkk,ccc,iii,aaa) &
                                        + aux_tmp(kkk,ccc,iii,aaa)
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end if

    Deallocate(aux_tmp)

  end do

  Deallocate(t_tmp,v_tmp)

  End Subroutine CC_cl_aux_j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_aux_k()

! Caluclate intermediate k (kkk,ccc,iii,aaa) 

  Use CC_cl

  Implicit None

  Integer :: n_c,c_start,c_end,n_a,a_start,a_end,n_d,d_start,d_end
  Integer :: i_run,k_run,a_run,c_run,d_run
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,i_tmp,j_tmp,k_tmp,l_tmp
  Integer (kind = 8) :: s_tmp
  Integer :: errnum

  Integer :: i_task,i_task_run
  Integer :: req1,req2,i_send,i_recv
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Double precision :: alpha,beta
  Double precision , dimension(:,:,:,:) , allocatable :: matA,matB,matC
  Double precision , dimension(:,:,:,:) , allocatable :: t_tmp,v_tmp,k_send,k_recv,aux_tmp

  n_c = CC_mem_aa_D(CC_mpi_did+1)
  c_start = CC_index_aa_D(CC_mpi_did+1)
  c_end = c_start - 1 + n_c

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  alpha = 1.0D0
  beta = 0.0D0

  if (.not.(allocated(CC_k_aux))) then
    ! (k,c,i,a)
    Allocate(CC_k_aux(CC_n_occ,n_c,CC_n_occ,n_a),stat=errnum)
  end if

  ! First term (ki|ac)
  CC_k_aux = CC_intl_kiac(:,:,:,a_start:a_end)

  ! Second term
  j_tmp = CC_n_occ * n_c * CC_n_occ
  !(l,k,c,i)
  alpha = -1.0D0
  beta = 1.0D0

  Call Dgemm('N','N',j_tmp,n_a,CC_n_occ,alpha,CC_intl_kilc,j_tmp, &
             CC_t_s(:,a_start:a_end),CC_n_occ,beta,CC_k_aux,j_tmp)

  ! Third term 
  i_tmp = CC_n_occ * n_c

  alpha = 1.0D0
  beta = 1.0D0

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(a_run,aaa)
  !$OMP DO
  do a_run = 1, n_a

    aaa = a_start - 1 + a_run

    Call Dgemm('N','T',i_tmp,CC_n_occ,CC_n_vir,alpha,CC_intl_ackd(:,:,:,aaa), &
               i_tmp,CC_t_s,CC_n_occ,beta,CC_k_aux(:,:,:,a_run),i_tmp)

  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Third term first part & Forth term
  ! Get T
  ! t_tmp(l,d,i,a)
  Allocate(t_tmp(CC_n_occ,n_c,CC_n_occ,n_a),stat=errnum)
  Call check_allocation(errnum,'t_tmp in CC')

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(a_run,d_run,aaa,iii,ddd,lll)
  !$OMP DO
  do a_run = 1, n_a
    do iii = 1, CC_n_occ
      do d_run = 1, n_c
        do lll = 1, CC_n_occ
          aaa = a_start - 1 + a_run
          ddd = c_start - 1 + d_run
          t_tmp(lll,d_run,iii,a_run) = 0.5D0 * CC_t_d(iii,lll,d_run,aaa) &
                                     + CC_t_s(iii,ddd) * CC_t_s(lll,aaa)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Get V
  ! v_tmp(l,d,k,c)
  Allocate(v_tmp(CC_n_occ,n_c,CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'v_tmp in CC')

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(ccc,kkk,ddd,lll)
  !$OMP DO
  do ccc = 1, CC_n_vir
    do kkk = 1, CC_n_occ
      do ddd = 1, n_c
        do lll = 1, CC_n_occ
          v_tmp(lll,ddd,kkk,ccc) = CC_intl_iajb(kkk,lll,ddd,ccc)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  n_d = CC_mem_aa_D(CC_mpi_did+1)
  d_start = CC_index_aa_D(CC_mpi_did+1)
  d_end = d_start - 1 + n_d

  i_task = CC_mpi_did + 1
  i_send = CC_mpi_did + 1
  i_recv = CC_mpi_did + 1

  do i_task_run = 1, CC_mpi_domain_size

    i_recv = i_recv - 1
    if (i_recv.lt.1) then
      i_recv = CC_mpi_domain_size
    end if

    i_send = i_send + 1
    if (i_send.gt.CC_mpi_domain_size) then
      i_send = 1
    end if

    i_task = i_task + 1
    if (i_task.gt.CC_mpi_domain_size) then
      i_task = 1
    end if

    n_c = CC_mem_aa_D(i_task)
    c_start = CC_index_aa_D(i_task)
    c_end = c_start - 1 + n_c

    !(k,c,i,a)
    Allocate(aux_tmp(CC_n_occ,n_c,CC_n_occ,n_a),stat=errnum)
    Call check_allocation(errnum,'matC in CC')

    ! Forth term
    alpha = -1.0D0
    beta = 0.0D0

    i_tmp = CC_n_occ * n_d
    j_tmp = CC_n_occ * n_c
    l_tmp = CC_n_occ * n_a

    Call Dgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,v_tmp(:,:,:,c_start:c_end),i_tmp, &
               t_tmp,i_tmp,beta,aux_tmp,j_tmp)

    !if (myid.eq.0) then
    !  print*,'matC',matC(n_k,CC_n_vir,n_i,n_a)
    !end if

    if (i_task_run.ne.1) then

      Call MPI_WAIT(req1,stat1,errnum)
      Call MPI_WAIT(req2,stat2,errnum)

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(aaa,iii,ccc,kkk)
      !$OMP DO
      do aaa = 1, n_a
        do iii = 1, CC_n_occ
          do ccc = 1, n_d
            do kkk = 1, CC_n_occ
              CC_k_aux(kkk,ccc,iii,aaa) = CC_k_aux(kkk,ccc,iii,aaa) &
                                        + k_recv(kkk,ccc,iii,aaa)
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      Deallocate(k_send,k_recv)

    end if

    if (i_task_run.ne.CC_mpi_domain_size) then

      Allocate(k_send(CC_n_occ,n_c,CC_n_occ,n_a),stat=errnum)
      Call check_allocation(errnum,'k_send in CC')

      k_send = aux_tmp

      Allocate(k_recv(CC_n_occ,n_d,CC_n_occ,n_a),stat=errnum)
      Call check_allocation(errnum,'k_recv in CC')

      i_tmp = CC_n_occ * CC_n_occ * n_a
      s_tmp = Int(i_tmp,8) * Int(n_c,8)

      Call CC_mpi_real_isend(s_tmp,k_send,i_send-1,100,req1,CC_mpi_comm_domain)

      s_tmp = Int(i_tmp,8) * Int(n_d,8)
      Call CC_mpi_real_irecv(s_tmp,k_recv,i_recv-1,100,req2,CC_mpi_comm_domain)

    else

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(aaa,iii,ccc,kkk)
      !$OMP DO
      do aaa = 1, n_a
        do iii = 1, CC_n_occ
          do ccc = 1, n_d
            do kkk = 1, CC_n_occ
              CC_k_aux(kkk,ccc,iii,aaa) = CC_k_aux(kkk,ccc,iii,aaa) &
                                        + aux_tmp(kkk,ccc,iii,aaa)
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end if

    Deallocate(aux_tmp)

  end do

  Deallocate(t_tmp,v_tmp)

  !if (myid.eq.0) then
  !  print*,'k',CC_k_aux(n_k,CC_n_vir,CC_n_occ,n_a)
  !end if

  End Subroutine CC_cl_aux_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_aux_a()

! Calculate aux_a(k,l,ij,1)

  Use CC_cl

  Implicit None

  Integer :: n_i,i_start,i_end,n_kl,kl_start,kl_end,n_c,c_start,c_end
  Integer :: i_run,k_run,a_run,c_run,d_run,kl,i_kl
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,i_tmp,j_tmp,k_tmp,l_tmp
  Integer (kind = 8) :: s_tmp
  Integer :: errnum

  Integer :: i_task,i_task_run
  Integer :: req1,req2,i_send,i_recv
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Double precision :: alpha,beta
  Double precision , dimension(:,:,:,:) , allocatable :: matA,matB,matC
  Double precision , dimension(:,:,:,:) , allocatable :: t_tmp,a_send,a_recv,aux_tmp
  Double precision , dimension(:,:,:,:) , allocatable :: v_tmp

  n_c = CC_mem_aa_D(CC_mpi_did+1)
  c_start = CC_index_aa_D(CC_mpi_did+1)
  c_end = c_start - 1 + n_c

  n_i = CC_mem_ii_G(CC_mpi_gid+1)
  i_start = CC_index_ii_G(CC_mpi_gid+1)
  i_end = i_start - 1 + n_i

  n_kl = CC_mem_kl_D(CC_mpi_did+1)
  kl_start = CC_index_kl_D(CC_mpi_did+1)
  kl_end = kl_start - 1 + n_kl

  if (.not.(allocated(CC_a_aux))) then
    ! (kl,i,j)
    Allocate(CC_a_aux(n_kl,CC_n_occ,CC_n_occ,1),stat=errnum)
    Call check_allocation(errnum,'CC_a_aux in CC')

  end if

  ! First term 
  CC_a_aux = CC_intl_kilj

  if (n_i.gt.0) then

    alpha = 1.0D0
    beta = 0.0D0

    ! Get T
    ! t_tmp(c,d,i,j)
    Allocate(t_tmp(n_c,CC_n_vir,n_i,CC_n_occ),stat=errnum)
    Call check_allocation(errnum,'t_tmp in CC')

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(i_run,c_run,ccc,ddd,iii,jjj)
    !$OMP DO
    do jjj = 1, CC_n_occ
      do i_run = 1, n_i
        do ddd = 1, CC_n_vir
          do c_run = 1, n_c
            ccc = c_start - 1 + c_run
            iii = i_start - 1 + i_run
            t_tmp(c_run,ddd,i_run,jjj) = CC_t_d(iii,jjj,c_run,ddd) &
                                       + CC_t_s(iii,ccc) * CC_t_s(jjj,ddd)
          end do
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! Get V
    ! v_tmp(c,d,kl)
    Allocate(v_tmp(n_c,CC_n_vir,CC_n_occ*CC_n_occ,1),stat=errnum)
    Call check_allocation(errnum,'v_tmp in CC')

    do kkk = 1, CC_n_occ
      do lll = 1, CC_n_occ
        Call CC_cl_code(i_kl,kkk,lll,CC_n_occ,1)
        v_tmp(:,:,i_kl,1) = CC_intl_iajb(kkk,lll,:,:)
      end do
    end do

    i_task = CC_mpi_did + 1
    i_send = CC_mpi_did + 1
    i_recv = CC_mpi_did + 1

    do i_task_run = 1, CC_mpi_domain_size

      i_recv = i_recv - 1
      if (i_recv.lt.1) then
        i_recv = CC_mpi_domain_size
      end if

      i_send = i_send + 1
      if (i_send.gt.CC_mpi_domain_size) then
        i_send = 1
      end if

      i_task = i_task + 1
      if (i_task.gt.CC_mpi_domain_size) then
        i_task = 1
      end if

      n_kl = CC_mem_kl_D(i_task)
      kl_start = CC_index_kl_D(i_task)
      kl_end = kl_start - 1 + n_kl

      !(kl,i,j,1)
      Allocate(aux_tmp(n_kl,n_i,CC_n_occ,1),stat=errnum)
      Call check_allocation(errnum,'aux_tmp in CC')

      ! Second term
      ! matA(c,kl,i)
      Allocate(matA(n_c,n_kl,n_i,1), stat=errnum)
      Call check_allocation(errnum,'matA in CC')

      Allocate(matB(n_c,CC_n_occ,1,1), stat=errnum)
      Call check_allocation(errnum,'matB in CC')

      ! matC(k,l,i,j)
      Allocate(matC(n_kl,n_i,CC_n_occ,1), stat=errnum)
      Call check_allocation(errnum,'matC in CC')
      
      do k_run = 1, n_kl

        kl = kl_start - 1 + k_run

        Call CC_cl_decode(kl,kkk,lll,CC_n_occ,1)

        matA(:,k_run,:,1) = CC_intl_kilc(kkk,:,i_start:i_end,lll)

      end do

      do jjj = 1, CC_n_occ
        matB(:,jjj,1,1) = CC_t_s(jjj,c_start:c_end)
      end do

      j_tmp = n_kl * n_i

      alpha = 1.0D0
      beta = 0.0D0

      Call Dgemm('T','N',j_tmp,CC_n_occ,n_c,alpha,matA,n_c, &
                 matB,n_c,beta,aux_tmp,j_tmp)

      Deallocate(matA,matB,matC)

      ! Third term
      ! matA(c,kl,j)
      Allocate(matA(n_c,n_kl,CC_n_occ,1), stat=errnum)
      Call check_allocation(errnum,'matA in CC')

      Allocate(matB(n_c,n_i,1,1), stat=errnum)
      Call check_allocation(errnum,'matB in CC')

      ! matC(k,l,i,j)
      Allocate(matC(n_kl,CC_n_occ,n_i,1), stat=errnum)
      Call check_allocation(errnum,'matC in CC')
 
      do k_run = 1, n_kl

        kl = kl_start - 1 + k_run

        Call CC_cl_decode(kl,kkk,lll,CC_n_occ,1)

        matA(:,k_run,:,1) = CC_intl_kilc(lll,:,:,kkk)

      end do

      do i_run = 1, n_i
        iii = i_start - 1 + i_run
        matB(:,i_run,1,1) = CC_t_s(iii,c_start:c_end)
      end do

      j_tmp = n_kl * CC_n_occ

      Call Dgemm('T','N',j_tmp,n_i,n_c,alpha,matA,n_c, &
                 matB,n_c,beta,matC,j_tmp)

      do jjj = 1, CC_n_occ
        aux_tmp(:,:,jjj,1) = aux_tmp(:,:,jjj,1) + matC(:,jjj,:,1)
      end do

      Deallocate(matA,matB,matC)

      ! Forth term
      alpha = 1.0D0
      beta = 1.0D0

      i_tmp = CC_n_vir * n_c
      j_tmp = n_kl
      l_tmp = n_i * CC_n_occ

      Call Dgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,v_tmp(:,:,kl_start:kl_end,1),i_tmp, &
                 t_tmp,i_tmp,beta,aux_tmp,j_tmp)

      !if (myid.eq.0) then
      !  print*,'matC',matC(n_k,CC_n_vir,n_i,n_a)
      !end if

      if (i_task_run.ne.1) then

        Call MPI_WAIT(req1,stat1,errnum)
        Call MPI_WAIT(req2,stat2,errnum)

        CC_a_aux(:,i_start:i_end,:,1) = CC_a_aux(:,i_start:i_end,:,1) &
                                      + a_recv(:,:,:,1)

        Deallocate(a_send,a_recv)

      end if

      if (i_task_run.ne.CC_mpi_domain_size) then

        Allocate(a_send(n_kl,n_i,CC_n_occ,1),stat=errnum)
        Call check_allocation(errnum,'j_send in CC')

        a_send = aux_tmp

        Allocate(a_recv(CC_mem_kl_D(CC_mpi_did+1),n_i,CC_n_occ,1),stat=errnum)
        Call check_allocation(errnum,'j_recv in CC')

        i_tmp = n_i * CC_n_occ 
        s_tmp = Int(i_tmp,8) * Int(n_kl,8)

        Call CC_mpi_real_isend(s_tmp,a_send,i_send-1,100,req1,CC_mpi_comm_domain)

        s_tmp = Int(i_tmp,8) * Int(CC_mem_kl_D(CC_mpi_did+1),8)
        Call CC_mpi_real_irecv(s_tmp,a_recv,i_recv-1,100,req2,CC_mpi_comm_domain)

      else

        CC_a_aux(:,i_start:i_end,:,1) = CC_a_aux(:,i_start:i_end,:,1)&
                                      + aux_tmp(:,:,:,1)

      end if

      Deallocate(aux_tmp)

    end do

    Deallocate(t_tmp,v_tmp)

  end if

  i_tmp = CC_mem_kl_D(CC_mpi_did+1)
  j_tmp = CC_n_occ * CC_n_occ
  s_tmp = Int(i_tmp,8) * Int(j_tmp,8)

  Call CC_mpi_allreduce(s_tmp,CC_a_aux,CC_mpi_comm_group)

  End Subroutine CC_cl_aux_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_aux_b(a_finish,n_a,CC_b_aux,RI_B,ctime_calc,wtime_calc,&
                                                    ctime_comm,wtime_comm)

  Use timing
  Use CC_cl

  Implicit None

  Integer , intent(in) :: a_finish,n_a
  Double precision , dimension(CC_n_vir,CC_mem_aa_D(CC_mpi_did+1),n_a,CC_n_vir) :: CC_b_aux
  Double precision , dimension(CC_n_bas,CC_mem_aa_D(CC_mpi_did+1),CC_n_vir) :: RI_B

  Integer :: i_finish

  Integer :: a_start,a_end,i_a,a_run
  Integer :: c_start,c_end,n_c,d_start,d_end,n_d,c_run

  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,i_tmp,j_tmp
  Integer :: errnum

  Integer :: oid,nth,OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS,calc_nth

  Integer (kind = 8) :: s_tmp

  Integer :: i_task,i_task_run
  Integer :: req1,req2,i_send,i_recv
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Double precision :: c_pre,w_pre,c_now,w_now,ctime_calc,wtime_calc,ctime_comm,wtime_comm
  Double precision :: ctime_1,wtime_1,ctime_2,wtime_2,ctime_3,wtime_3

  Double precision :: alpha,beta
  Double precision , dimension(:,:,:) , allocatable :: RI_tmp,RI_send,RI_recv
  Double precision , dimension(:,:,:,:) , allocatable :: b_tmp

  ctime_2 = 0.0D0
  wtime_2 = 0.0D0
  ctime_3 = 0.0D0
  wtime_3 = 0.0D0

  n_d = CC_mem_aa_D(CC_mpi_did+1)
  d_start = CC_index_aa_D(CC_mpi_did+1)
  d_end = d_start - 1 + n_d

  a_start = a_finish +1
  a_end = a_finish + n_a

  i_tmp = n_d * CC_n_vir * CC_n_vir

  alpha = -1.0D0
  beta = 0.0D0

  ! b_tmp(d,c,b,a)
  Allocate(b_tmp(n_d,CC_n_vir,CC_n_vir,n_a),stat=errnum)
  Call check_allocation(errnum,'t_tmp in CC')

  b_tmp = 0.0D0

  Call get_timestamps(c_pre, w_pre)

  Call Dgemm('T','N',i_tmp,n_a,CC_n_occ,alpha,CC_intl_ackd,CC_n_occ, &
             CC_t_s(:,a_start:a_end),CC_n_occ,beta,b_tmp,i_tmp)

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(aaa,bbb,ccc,ddd)
  !$OMP DO
  do aaa = 1, n_a
    do bbb = 1, CC_n_vir
      do ccc = 1, CC_n_vir
        do ddd = 1, n_d
          CC_b_aux(ccc,ddd,aaa,bbb) = b_tmp(ddd,ccc,bbb,aaa)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  Deallocate(b_tmp)

  !if (CC_abcd_sv_flag) then
  !  CC_b_aux = CC_b_aux + CC_intl_acbd(:,:,ab_start:ab_end,1)
  !end if

  Call get_timestamps(c_now, w_now)
  ctime_calc = ctime_calc + c_now - c_pre
  wtime_calc = wtime_calc + w_now - w_pre
  ctime_1 = c_now - c_pre
  wtime_1 = w_now - w_pre

  !RI_tmp(bas,ccc,aaa)
  Allocate(RI_tmp(CC_n_bas,n_d,n_a),stat=errnum)
  Call check_allocation(errnum,'RI_tmp in CC1')

  !RI_send(bas,ccc,aaa)
  Allocate(RI_send(CC_n_bas,n_d,n_a),stat=errnum)
  Call check_allocation(errnum,'RI_send in CC1')

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(a_run,aaa,ccc)
  !$OMP DO
  do a_run = 1, n_a
    do ccc = 1, n_d
      aaa = a_finish + a_run + CC_n_occ
      RI_tmp(:,ccc,a_run) = CC_RI_B(:,ccc,aaa)
      RI_send(:,ccc,a_run) = CC_RI_B(:,ccc,aaa)
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  i_task = CC_mpi_did
  i_send = CC_mpi_did + 1
  i_recv = CC_mpi_did + 1

  do i_task_run = 1, CC_mpi_domain_size

    i_recv = i_recv + 1
    if (i_recv.gt.CC_mpi_domain_size) then
      i_recv = 1
    end if

    i_send = i_send - 1
    if (i_send.lt.1) then
      i_send = CC_mpi_domain_size
    end if

    i_task = i_task + 1
    if (i_task.gt.CC_mpi_domain_size) then
      i_task = 1
    end if

    n_c = CC_mem_aa_D(i_task)
    c_start = CC_index_aa_D(i_task)
    c_end = c_start - 1 + n_c

    if (i_recv.ne.CC_mpi_did+1) then

      Allocate(RI_recv(CC_n_bas,CC_mem_aa_D(i_recv),n_a),stat=errnum)
      Call check_allocation(errnum,'RI_recv in CC')

      i_tmp = CC_n_bas * n_a
      s_tmp = Int(CC_mem_aa_D(CC_mpi_did+1),8) * Int(i_tmp,8)

      Call CC_mpi_real_isend(s_tmp,RI_send,i_send-1,100,req1,CC_mpi_comm_domain)

      s_tmp = Int(CC_mem_aa_D(i_recv),8) * Int(i_tmp,8)
      Call CC_mpi_real_irecv(s_tmp,RI_recv,i_recv-1,100,req2,CC_mpi_comm_domain)

    end if

    ! b_tmp(c,a,d,b)
    Allocate(b_tmp(n_c,n_a,n_d,CC_n_vir), stat=errnum)
    Call check_allocation(errnum,'b_tmp in CC')

    i_tmp = n_c * n_a
    j_tmp = n_d * CC_n_vir

    Call get_timestamps(c_pre, w_pre)

    alpha = 1.0D0
    beta = 0.0D0

    Call Dgemm('T','N',i_tmp,j_tmp,CC_n_bas,alpha, &
               RI_tmp,CC_n_bas,RI_B,CC_n_bas,beta,b_tmp,i_tmp)

    !$OMP PARALLEL Default(shared) &
    !$OMP private(aaa,bbb,ccc,ddd,c_run)
    !$OMP DO
    do aaa = 1, n_a
      do bbb = 1, CC_n_vir
        do c_run = 1, n_c
          do ddd = 1, n_d
            ccc = c_start - 1 + c_run
            CC_b_aux(ccc,ddd,aaa,bbb) = CC_b_aux(ccc,ddd,aaa,bbb) &
                                      + b_tmp(c_run,aaa,ddd,bbb)
          end do
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    Call get_timestamps(c_now, w_now)
    ctime_calc = ctime_calc + c_now - c_pre
    wtime_calc = wtime_calc + w_now - w_pre

    Deallocate(RI_tmp,b_tmp)

    if (i_recv.ne.CC_mpi_did+1) then

      Call get_timestamps(c_pre, w_pre)

      Call MPI_WAIT(req1,stat1,errnum)
      Call MPI_WAIT(req2,stat2,errnum)

      Call get_timestamps(c_now, w_now)
      ctime_comm = ctime_comm + c_now - c_pre
      wtime_comm = wtime_comm + w_now - w_pre

      Allocate(RI_tmp(CC_n_bas,CC_mem_aa_D(i_recv),n_a),stat=errnum)
      Call check_allocation(errnum,'RI_tmp in CC2')

      !$OMP PARALLEL Default(shared) &
      !$OMP private(aaa,ccc)
      !$OMP DO
      do aaa = 1, n_a
        do ccc = 1, CC_mem_aa_D(i_recv)
          RI_tmp(:,ccc,aaa) = RI_recv(:,ccc,aaa)
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      Deallocate(RI_recv)

    end if

  end do

  Deallocate(RI_send)

  End Subroutine CC_cl_aux_b

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

