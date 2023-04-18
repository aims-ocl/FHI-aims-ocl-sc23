  Subroutine CC_cl_PT()

  Use CC_cl

  Implicit None

  Call CC_cl_PT_initial()
  Call CC_cl_PT_calc()

  End Subroutine CC_cl_PT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_PT_initial()

  Use CC_cl

  Implicit None

  Integer :: i_tmp

  Logical :: resex
  Integer :: errnum

  ! Clean up unnecessary tensors
  Deallocate(CC_h_ik,CC_g_ik,CC_h_ca,CC_g_ca,CC_h_ck)

  Deallocate(CC_DIIS_Bmat,CC_DIIS_tau)
  Deallocate(CC_DIIS_r_s,CC_DIIS_t_s)
  Deallocate(CC_DIIS_r_d,CC_DIIS_t_d)

  if (allocated(CC_a_aux)) then
    Deallocate(CC_a_aux)
  end if

  if (allocated(CC_j_aux)) then
    Deallocate(CC_j_aux)
  end if

  if (allocated(CC_k_aux)) then
    Deallocate(CC_k_aux)
  end if

  if (allocated(CC_tau_m)) then
    Deallocate(CC_tau_m)
  end if

  if (allocated(CC_tau_p)) then
    Deallocate(CC_tau_p)
  end if

  Deallocate(CC_intl_kiac)
  Deallocate(CC_intl_kilj,CC_intl_ackd)

  if (allocated(CC_intl_acbd)) then
    Deallocate(CC_intl_acbd)
  end if

  ! redistribution of tensors
  Call CC_cl_PT_redistribution()

  ! Clean up unnecessary tensors
  if (allocated(CC_intl_iajb)) then
    Deallocate(CC_intl_iajb)
  end if

  if (allocated(CC_intl_kilc)) then
    Deallocate(CC_intl_kilc)
  end if

  if (allocated(CC_t_d)) then
    Deallocate(CC_t_d)
  end if

  if (allocated(CC_RI)) then
    Deallocate(CC_RI)
  end if

  ! Create distribution vectors
  Allocate(CC_mem_ijk(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_ijk in CC')

  Allocate(CC_index_ijk(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_ijk in CC')

  Call CC_cl_pyramidal((CC_n_occ-2),i_tmp)
  if (i_tmp.ge.1) then
    Call CC_cl_vec_distr(1,i_tmp,CC_mpi_group_size,CC_mem_ijk,CC_index_ijk)
  else
    CC_mem_ijk = 0
    CC_index_ijk = 1
  end if

  Allocate(CC_mem_abc(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_abc in CC')

  Allocate(CC_index_abc(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_abc in CC')

  Call CC_cl_pyramidal(CC_n_vir,i_tmp)
  Call CC_cl_vec_distr(1,i_tmp,CC_mpi_domain_size,CC_mem_abc,CC_index_abc)

  ! Create vectors
  ! 1 ---- (a,b,c)
  ! 2 ---- (b,c,a)
  ! 3 ---- (c,a,b)
  ! 4 ---- (a,c,b)
  ! 5 ---- (b,a,c)
  ! 6 ---- (c,b,a)
  Allocate(CC_PT_W_abc(CC_mem_abc(CC_mpi_did+1),6,1), stat=errnum)
  Call check_allocation(errnum,'CC_PT_W_abc in CC')

  ! 1 ---- U(b,c,a) ---- (a,b,c)
  ! 2 ---- U(c,a,b) ---- (b,c,a)
  ! 3 ---- U(a,b,c) ---- (c,a,b)
  Allocate(CC_PT_U(CC_n_vir,CC_n_vir,CC_mem_aa_D(CC_mpi_did+1),3), stat=errnum)
  Call check_allocation(errnum,'CC_PT_U in CC')

  ! 1 ---- t(d,c,j,k)
  ! 2 ---- t(d,b,k,j)
  ! 3 ---- t(d,b,i,j)
  ! 4 ---- t(d,a,j,i)
  ! 5 ---- t(d,a,k,i)
  ! 6 ---- t(d,c,i,k)
  Allocate(CC_PT_T_abc(CC_n_vir,CC_n_vir,6,1), stat=errnum)
  Call check_allocation(errnum,'CC_PT_T_abc in CC')

  Allocate(CC_PT_T_send(CC_n_vir,CC_mem_aa_D(1),6,CC_mpi_domain_size,1), stat=errnum)
  Call check_allocation(errnum,'CC_PT_T_send in CC')

  Allocate(CC_PT_T_recv(CC_n_vir,CC_mem_aa_D(1),6,CC_mpi_domain_size,1), stat=errnum)
  Call check_allocation(errnum,'CC_PT_T_recv in CC')

  ! 1 ---- intl_jlkc(l,c,j,k)
  ! 2 ---- intl_jlkc(l,b,k,j)
  ! 3 ---- intl_jlkc(l,b,i,j)
  ! 4 ---- intl_jlkc(l,a,j,i)
  ! 5 ---- intl_jlkc(l,a,k,i)
  ! 6 ---- intl_jlkc(l,c,i,k)
  Allocate(CC_PT_jlkc_abc(CC_n_occ,CC_n_vir,6), stat=errnum)
  Call check_allocation(errnum,'CC_PT_jlkc_abc in CC')

  Allocate(CC_PT_jlkc_send(CC_n_occ,CC_mem_aa_D(1),6,CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_PT_jlkc_send in CC')

  Allocate(CC_PT_jlkc_recv(CC_n_occ,CC_mem_aa_D(1),6,CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_PT_jlkc_recv in CC')

  ! 1 ---- intl_iajb(j,k,b,c)
  ! 2 ---- intl_iajb(i,k,a,c)
  ! 3 ---- intl_iajb(i,j,a,b)
  Allocate(CC_PT_iajb_abc(CC_n_vir,CC_n_vir,3), stat=errnum)
  Call check_allocation(errnum,'CC_PT_iajb_abc in CC')

  Allocate(CC_PT_iajb_send(CC_n_vir,CC_mem_aa_D(1),3,CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_PT_iajb_send in CC')

  Allocate(CC_PT_iajb_recv(CC_n_vir,CC_mem_aa_D(1),3,CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_PT_iajb_recv in CC')

  End Subroutine CC_cl_PT_initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_PT_calc()

  Use timing
  Use CC_cl

  Implicit None

  Integer :: errnum,i_step,i_tmp

  Integer , dimension(CC_mpi_domain_size,2) :: req_t,req_jlkc,req_iajb

  Double precision :: CC_clock_stamp,CC_time_stamp,CC_clock_pre,CC_time_elapse
  Double precision :: ctime_calc,wtime_calc,ctime_comm,wtime_comm,ctime_calcu,wtime_calcu
  Double precision :: ctime_calcw,wtime_calcw,ctime_commw,wtime_commw
  Double precision :: ctime_pre,wtime_pre,ctime_now,wtime_now

  Double precision :: CC_PT_percent
  Integer , dimension(CC_PT_n_out) :: CC_PT_outstep
  Logical :: CC_PT_finish_flag
  Logical , dimension(CC_PT_n_out) :: out_flag

  Integer :: iii,jjj,kkk,i_next,j_next,k_next
  Integer :: n_ijk,n_calc,n_ij,i_calc
  Double precision :: rlt

  ctime_calc = 0.0D0
  wtime_calc = 0.0D0
  ctime_comm = 0.0D0
  wtime_comm = 0.0D0
  ctime_calcu = 0.0D0
  wtime_calcu = 0.0D0
  ctime_calcw = 0.0D0
  wtime_calcw = 0.0D0
  ctime_commw = 0.0D0
  wtime_commw = 0.0D0

  !print*,'CCSD(T) start'

  n_ijk = CC_mem_ijk(CC_mpi_gid+1)
  n_ij = CC_mem_ij_G(CC_mpi_gid+1)
  n_calc = n_ijk + n_ij

  !print*,'n_ijk',n_ijk,'n_ij',n_ij,'n_calc',n_calc

  i_step = Int(n_calc / CC_PT_n_out)

  do iii = 1, CC_PT_n_out
    CC_PT_outstep(iii) = i_step * iii
    CC_PT_outstep(iii) = CC_PT_outstep(iii) + min(mod(n_calc,CC_PT_n_out),iii)
    !print*,'outstep',iii,CC_PT_outstep(iii),i_step
  end do

  CC_PT_finish_flag = .false.

  i_calc = 0
  i_calc = CC_res_inf(CC_mpi_gid+2)

  out_flag = .true.
  do iii = 1, CC_PT_n_out - 1
    if (CC_PT_outstep(iii).eq.CC_PT_outstep(iii+1)) then
      out_flag(iii) = .false.
    end if
  end do

  do iii = 1, CC_PT_n_out
    if (i_calc.ge.CC_PT_outstep(iii)) then
      out_flag(iii) = .false.
    end if
  end do

  if (i_calc.lt.n_calc) then

    Call get_timestamps(CC_time_stamp, CC_clock_stamp)
    CC_time_elapse = 0.0D0

    Call CC_cl_PT_get_ijk(i_calc+1,iii,jjj,kkk)

    ! First iteration
    Call CC_cl_PT_send(iii,jjj,kkk,req_t,req_jlkc,req_iajb)

    do while(.not.(CC_PT_finish_flag))

      i_calc = i_calc + 1

      Call CC_cl_PT_get_ijk(i_calc,iii,jjj,kkk)

      !print*,'iii,jjj,kkk',iii,jjj,kkk

      Call get_timestamps(ctime_pre, wtime_pre)

      Call CC_cl_PT_recv(req_t,req_jlkc,req_iajb)

      Call get_timestamps(ctime_now, wtime_now)
      ctime_comm = ctime_comm + ctime_now - ctime_pre
      wtime_comm = wtime_comm + wtime_now - wtime_pre

      if (kkk.eq.jjj) then

        i_next = iii
        j_next = iii
        k_next = kkk
        Call CC_cl_PT_send(i_next,j_next,k_next,req_t,req_jlkc,req_iajb)

      else if (i_calc.lt.n_calc) then

        Call CC_cl_PT_get_ijk(i_calc+1,i_next,j_next,k_next)
        Call CC_cl_PT_send(i_next,j_next,k_next,req_t,req_jlkc,req_iajb)

      end if

      Call CC_cl_PT_calc_W(iii,jjj,kkk,ctime_calcu,wtime_calcu)
      Call CC_cl_PT_sync_W(ctime_commw,wtime_commw,ctime_calcw,wtime_calcw)
      Call CC_cl_PT_calc_corr(iii,jjj,kkk,rlt)

      CC_E_PT = CC_E_PT + rlt

      if (kkk.eq.jjj) then

        Call get_timestamps(ctime_pre, wtime_pre)

        Call CC_cl_PT_recv(req_t,req_jlkc,req_iajb)

        Call get_timestamps(ctime_now, wtime_now)
        ctime_comm = ctime_comm + ctime_now - ctime_pre
        wtime_comm = wtime_comm + wtime_now - wtime_pre

        if (i_calc.lt.n_calc) then

          Call CC_cl_PT_get_ijk(i_calc+1,i_next,j_next,k_next)
          Call CC_cl_PT_send(i_next,j_next,k_next,req_t,req_jlkc,req_iajb)

        end if

        jjj = iii
        Call CC_cl_PT_calc_W(iii,jjj,kkk,ctime_calcu,wtime_calcu)
        Call CC_cl_PT_sync_W(ctime_commw,wtime_commw,ctime_calcw,wtime_calcw)
        Call CC_cl_PT_calc_corr(iii,jjj,kkk,rlt)

        CC_E_PT = CC_E_PT + rlt

      end if

      if (i_calc.eq.n_calc) then

        CC_PT_finish_flag = .true.

      end if

      CC_clock_pre = CC_clock_stamp
      Call get_timestamps(CC_time_stamp, CC_clock_stamp)

      CC_time_elapse = CC_time_elapse + CC_clock_stamp - CC_clock_pre

      Call CC_cl_PT_output(i_calc,CC_PT_outstep,out_flag,CC_time_elapse)

    end do

  end if

  Call CC_mpi_real_number(CC_E_PT,MPI_COMM_WORLD)

  if (myid.eq.0) then
    write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of calc U', &
                                               ctime_calcu, 's (cpu)', &
                                               wtime_calcu, 's (wall clock)'
  end if

  if (myid.eq.0) then
    write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of calc W', &
                                               ctime_calcw, 's (cpu)', &
                                               wtime_calcw, 's (wall clock)'
  end if

  if (myid.eq.0) then
    write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of comm W', &
                                               ctime_commw, 's (cpu)', &
                                               wtime_commw, 's (wall clock)'
  end if

  if (myid.eq.0) then
    write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of comm', &
                                               ctime_comm, 's (cpu)', &
                                               wtime_comm, 's (wall clock)'
  end if

  End Subroutine CC_cl_PT_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_PT_get_ijk(i_calc,iii,jjj,kkk)

  Use CC_cl

  Implicit None

  Integer , intent(in) :: i_calc
  Integer , intent(out) :: iii,jjj,kkk
  Integer :: n_ijk,i_state,ijk_start,ij_start

  n_ijk = CC_mem_ijk(CC_mpi_gid+1)
  ijk_start = CC_index_ijk(CC_mpi_gid+1)

  ij_start = CC_index_ij_G(CC_mpi_gid+1)

  if (i_calc.le.n_ijk) then
    i_state = ijk_start - 1 + i_calc
    Call CC_cl_decode_ijk_ne(i_state,iii,jjj,kkk)
  else
    i_state = ij_start - 1 + i_calc - n_ijk
    Call CC_cl_decode(i_state,jjj,iii,CC_n_occ,3)
    kkk = jjj
  end if

  End Subroutine CC_cl_PT_get_ijk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_PT_calc_W(iii,jjj,kkk,ctime_calc,wtime_calc)

  Use timing
  Use CC_cl

  Implicit None

  Integer , intent(in) :: iii,jjj,kkk

  Double precision :: ctime_calc,wtime_calc
  Double precision :: ctime_pre,wtime_pre,ctime_now,wtime_now

  Integer :: a_start,a_end,n_a

  Integer :: errnum

  Integer :: n_abc,i_abc,code_abc,abc_start
  Integer :: aaa,bbb,ccc,ddd,i_tmp,j_tmp,aa1,bb1,cc1

  Integer (kind=8) :: s_tmp

  Double precision :: alpha,beta
  Double precision , dimension(:,:,:) , allocatable :: matC

  Call get_timestamps(ctime_pre, wtime_pre)

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  n_abc = CC_mem_abc(CC_mpi_did+1)

  CC_PT_U = 0.0D0

  i_tmp = CC_n_vir * n_a

  Allocate(matC(CC_n_vir,CC_n_vir,n_a), stat=errnum)
  Call check_allocation(errnum,'matC in CC')

  !aa1 = 1
  !bb1 = 1
  !cc1 = 1
  !print*,'U',iii,jjj,kkk,aa1,bb1,cc1

  ! Calculate U
  ! 1 --- t(k,j,c,d) # (db|ia) - (jl|kc) # t(l,i,b,a) = U(c,b,a)
  alpha = 1.0D0
  beta = 0.0D0
  Call Dgemm('T','N',CC_n_vir,i_tmp,CC_n_vir,alpha,CC_PT_T_abc(:,:,1,1),CC_n_vir, &
             CC_intl_dbia(:,:,:,iii),CC_n_vir,beta,matC,CC_n_vir)

  alpha = - 1.0D0
  beta = 1.0D0
  Call Dgemm('T','N',CC_n_vir,i_tmp,CC_n_occ,alpha,CC_PT_jlkc_abc(:,:,1),CC_n_occ, &
             CC_t_d_I(:,:,:,iii),CC_n_occ,beta,matC,CC_n_vir)

  !write(60+myid,*) '111',matC(cc1,bb1,aa1)
  !write(60+myid,*) CC_PT_T_abc(:,cc1,1,1)
  !write(60+myid,*) CC_intl_dbia(:,bb1,aa1,iii)

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(aaa,bbb,ccc)
  !$OMP DO
  do aaa = 1, n_a
    do ccc = 1, CC_n_vir
      do bbb = 1, CC_n_vir
        CC_PT_U(bbb,ccc,aaa,1) = CC_PT_U(bbb,ccc,aaa,1) + matC(ccc,bbb,aaa)
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  if (iii.eq.jjj) then

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(aaa,bbb,ccc)
    !$OMP DO
    do aaa = 1, CC_n_vir
      do ccc = 1, CC_n_vir
        do bbb = 1, n_a
          CC_PT_U(ccc,aaa,bbb,2) = CC_PT_U(ccc,aaa,bbb,2) + matC(ccc,aaa,bbb)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  else if (jjj.eq.kkk) then

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(aaa,bbb,ccc)
    !$OMP DO
    do aaa = 1, n_a
      do ccc = 1, CC_n_vir
        do bbb = 1, CC_n_vir
          CC_PT_U(bbb,ccc,aaa,1) = CC_PT_U(bbb,ccc,aaa,1) + matC(bbb,ccc,aaa)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end if

  ! 3 --- t(j,i,b,d) # (da|kc) - (il|jb) # t(l,k,a,c) = U(b,a,c)
  alpha = 1.0D0
  beta = 0.0D0
  Call Dgemm('T','N',CC_n_vir,i_tmp,CC_n_vir,alpha,CC_PT_T_abc(:,:,3,1),CC_n_vir, &
             CC_intl_dbia(:,:,:,kkk),CC_n_vir,beta,matC,CC_n_vir)

  alpha = - 1.0D0
  beta = 1.0D0
  Call Dgemm('T','N',CC_n_vir,i_tmp,CC_n_occ,alpha,CC_PT_jlkc_abc(:,:,3),CC_n_occ, &
             CC_t_d_I(:,:,:,kkk),CC_n_occ,beta,matC,CC_n_vir)

  !write(60+myid,*) '111',matC(bb1,aa1,cc1)

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(aaa,bbb,ccc)
  !$OMP DO
  do aaa = 1, CC_n_vir
    do ccc = 1, n_a
      do bbb = 1, CC_n_vir
        CC_PT_U(aaa,bbb,ccc,3) = CC_PT_U(aaa,bbb,ccc,3) + matC(bbb,aaa,ccc)
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  if (iii.eq.jjj) then

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(aaa,bbb,ccc)
    !$OMP DO
    do aaa = 1, CC_n_vir
      do ccc = 1, n_a
        do bbb = 1, CC_n_vir
          CC_PT_U(aaa,bbb,ccc,3) = CC_PT_U(aaa,bbb,ccc,3) + matC(aaa,bbb,ccc)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  else if (jjj.eq.kkk) then

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(aaa,bbb,ccc)
    !$OMP DO
    do aaa = 1, CC_n_vir
      do ccc = 1, CC_n_vir
        do bbb = 1, n_a
          CC_PT_U(ccc,aaa,bbb,2) = CC_PT_U(ccc,aaa,bbb,2) + matC(ccc,aaa,bbb)
        end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
      
  end if

  ! 5 --- t(i,k,a,d) # (dc|jb) - (kl|ia) # t(l,j,c,b) = U(a,c,b)
  alpha = 1.0D0
  beta = 0.0D0
  Call Dgemm('T','N',CC_n_vir,i_tmp,CC_n_vir,alpha,CC_PT_T_abc(:,:,5,1),CC_n_vir, &
             CC_intl_dbia(:,:,:,jjj),CC_n_vir,beta,matC,CC_n_vir)

  alpha = - 1.0D0
  beta = 1.0D0
  Call Dgemm('T','N',CC_n_vir,i_tmp,CC_n_occ,alpha,CC_PT_jlkc_abc(:,:,5),CC_n_occ, &
             CC_t_d_I(:,:,:,jjj),CC_n_occ,beta,matC,CC_n_vir)

  !write(60+myid,*) '111',matC(aa1,cc1,bb1)

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(aaa,bbb,ccc)
  !$OMP DO
  do aaa = 1, CC_n_vir
    do ccc = 1, CC_n_vir
      do bbb = 1, n_a
        CC_PT_U(ccc,aaa,bbb,2) = CC_PT_U(ccc,aaa,bbb,2) + matC(aaa,ccc,bbb)
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  if (iii.eq.jjj) then

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(aaa,bbb,ccc)
    !$OMP DO
    do aaa = 1, n_a
      do ccc = 1, CC_n_vir
        do bbb = 1, CC_n_vir
          CC_PT_U(bbb,ccc,aaa,1) = CC_PT_U(bbb,ccc,aaa,1) + matC(bbb,ccc,aaa)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  else if (jjj.eq.kkk) then

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(aaa,bbb,ccc)
    !$OMP DO
    do aaa = 1, CC_n_vir
      do ccc = 1, n_a
        do bbb = 1, CC_n_vir
          CC_PT_U(aaa,bbb,ccc,3) = CC_PT_U(aaa,bbb,ccc,3) + matC(aaa,bbb,ccc)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end if

  if ((iii.ne.jjj).and.(jjj.ne.kkk)) then

    ! 2 --- t(j,k,b,d) # (dc|ia) - (kl|jb) # t(l,i,c,a) = U(b,c,a)
    alpha = 1.0D0
    beta = 0.0D0
    Call Dgemm('T','N',CC_n_vir,i_tmp,CC_n_vir,alpha,CC_PT_T_abc(:,:,2,1),CC_n_vir, &
               CC_intl_dbia(:,:,:,iii),CC_n_vir,beta,matC,CC_n_vir)

    alpha = - 1.0D0
    beta = 1.0D0
    Call Dgemm('T','N',CC_n_vir,i_tmp,CC_n_occ,alpha,CC_PT_jlkc_abc(:,:,2),CC_n_occ, &
               CC_t_d_I(:,:,:,iii),CC_n_occ,beta,matC,CC_n_vir)

    !write(60+myid,*) '111',matC(bb1,cc1,aa1)

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(aaa,bbb,ccc)
    !$OMP DO
    do aaa = 1, n_a
      do ccc = 1, CC_n_vir
        do bbb = 1, CC_n_vir
          CC_PT_U(bbb,ccc,aaa,1) = CC_PT_U(bbb,ccc,aaa,1) + matC(bbb,ccc,aaa)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! 4 --- t(i,j,a,d) # (db|kc) - (jl|ia) # t(l,k,b,c) = U(a,b,c)
    alpha = 1.0D0
    beta = 0.0D0
    Call Dgemm('T','N',CC_n_vir,i_tmp,CC_n_vir,alpha,CC_PT_T_abc(:,:,4,1),CC_n_vir, &
               CC_intl_dbia(:,:,:,kkk),CC_n_vir,beta,matC,CC_n_vir)

    alpha = - 1.0D0
    beta = 1.0D0
    Call Dgemm('T','N',CC_n_vir,i_tmp,CC_n_occ,alpha,CC_PT_jlkc_abc(:,:,4),CC_n_occ, &
               CC_t_d_I(:,:,:,kkk),CC_n_occ,beta,matC,CC_n_vir)

    !write(60+myid,*) '111',matC(aa1,bb1,cc1)

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(aaa,bbb,ccc)
    !$OMP DO
    do aaa = 1, CC_n_vir
      do ccc = 1, n_a
        do bbb = 1, CC_n_vir
          CC_PT_U(aaa,bbb,ccc,3) = CC_PT_U(aaa,bbb,ccc,3) + matC(aaa,bbb,ccc)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! 6 --- t(k,i,c,d) # (da|jb) - (il|kc) # t(l,j,a,b) = U(c,a,b)
    alpha = 1.0D0
    beta = 0.0D0
    Call Dgemm('T','N',CC_n_vir,i_tmp,CC_n_vir,alpha,CC_PT_T_abc(:,:,6,1),CC_n_vir, &
               CC_intl_dbia(:,:,:,jjj),CC_n_vir,beta,matC,CC_n_vir)

    alpha = - 1.0D0
    beta = 1.0D0
    Call Dgemm('T','N',CC_n_vir,i_tmp,CC_n_occ,alpha,CC_PT_jlkc_abc(:,:,6),CC_n_occ, &
               CC_t_d_I(:,:,:,jjj),CC_n_occ,beta,matC,CC_n_vir)

    !write(60+myid,*) '111',matC(cc1,aa1,bb1)

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(aaa,bbb,ccc)
    !$OMP DO
    do aaa = 1, CC_n_vir
      do ccc = 1, CC_n_vir
        do bbb = 1, n_a
          CC_PT_U(ccc,aaa,bbb,2) = CC_PT_U(ccc,aaa,bbb,2) + matC(ccc,aaa,bbb)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end if

  Deallocate(matC)

  Call get_timestamps(ctime_now, wtime_now)
  ctime_calc = ctime_calc + ctime_now - ctime_pre
  wtime_calc = wtime_calc + wtime_now - wtime_pre

  End Subroutine CC_cl_PT_calc_W

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_PT_sync_W(ctime_comm,wtime_comm,ctime_calc,wtime_calc)

  Use timing
  Use CC_cl

  Implicit None

  Double precision :: ctime_comm,wtime_comm,ctime_calc,wtime_calc
  Double precision :: ctime_pre,wtime_pre,ctime_now,wtime_now

  Integer :: a_start,a_end,n_a

  Integer :: errnum

  Integer :: i_task,i_task_run
  Integer :: req1,req2,i_send,i_recv
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Integer :: n_abc,i_abc,code_abc,abc_start
  Integer :: aaa,bbb,ccc,ddd,i_tmp,j_tmp

  Integer (kind=8) :: s_tmp

  Double precision :: alpha,beta
  Double precision , dimension(:,:) , allocatable :: w_tmp,w_send,w_recv

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  n_abc = CC_mem_abc(CC_mpi_did+1)

  CC_PT_W_abc = 0.0D0

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

    Allocate(w_tmp(CC_mem_abc(i_task),6), stat=errnum)
    Call check_allocation(errnum,'w_tmp in CCT')

    w_tmp = 0.0D0

    Call get_timestamps(ctime_pre, wtime_pre)

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(i_abc,aaa,bbb,ccc,ddd,code_abc)
    !$OMP DO Schedule(Dynamic)
    do i_abc = 1, CC_mem_abc(i_task)

      code_abc = CC_index_abc(i_task) - 1 + i_abc
      Call CC_cl_decode_ijk(code_abc,aaa,bbb,ccc)

      ! add U to W
      ! W(1) = W(a,b,c) is the reference change W(1) to W(n) by swaping a,b,c,
      ! one gets the formula for W(n)
      ! for example, W(2) = W(b,c,a), it can be obtained from W(1) by a permutation:
      ! a --> b, b --> c, c --> a, then we have:
      ! W(2) = U(c,a,b,1) + U(a,b,c,2) + U(b,c,a,3)
      ! 
      ! for W(1) = W(a,b,c) = U(b,c,a,1) + U(c,a,b,2) + U(a,b,c,3) --- ref
      !     W(2) = W(b,c,a) = U(c,a,b,1) + U(a,b,c,2) + U(b,c,a,3)
      !     W(3) = W(c,a,b) = U(a,b,c,1) + U(b,c,a,2) + U(c,a,b,3)
      !     W(4) = W(a,c,b) = U(c,b,a,1) + U(b,a,c,2) + U(a,c,b,3)
      !     W(5) = W(b,a,c) = U(a,c,b,1) + U(c,b,a,2) + U(b,a,c,3)
      !     W(6) = W(c,b,a) = U(b,a,c,1) + U(a,c,b,2) + U(c,b,a,3)
      ddd = aaa - a_start + 1
      if ((ddd.ge.1).and.(ddd.le.n_a)) then
        ! 1 ---- W(a,b,c) <---- U(b,c,a,1)
        w_tmp(i_abc,1) = w_tmp(i_abc,1) + CC_PT_U(bbb,ccc,ddd,1)
        ! 2 ---- W(b,c,a) <---- U(b,c,a,3)
        w_tmp(i_abc,2) = w_tmp(i_abc,2) + CC_PT_U(bbb,ccc,ddd,3)
        ! 3 ---- W(c,a,b) <---- U(b,c,a,2)
        w_tmp(i_abc,3) = w_tmp(i_abc,3) + CC_PT_U(bbb,ccc,ddd,2)
        ! 4 ---- W(a,c,b) <---- U(c,b,a,1)
        w_tmp(i_abc,4) = w_tmp(i_abc,4) + CC_PT_U(ccc,bbb,ddd,1)
        ! 5 ---- W(b,a,c) <---- U(c,b,a,2)
        w_tmp(i_abc,5) = w_tmp(i_abc,5) + CC_PT_U(ccc,bbb,ddd,2)
        ! 6 ---- W(c,b,a) <---- U(c,b,a,3)
        w_tmp(i_abc,6) = w_tmp(i_abc,6) + CC_PT_U(ccc,bbb,ddd,3)
      end if

      ddd = bbb - a_start + 1
      if ((ddd.ge.1).and.(ddd.le.n_a)) then
        ! 1 ---- W(a,b,c) <---- U(c,a,b,2)
        w_tmp(i_abc,1) = w_tmp(i_abc,1) + CC_PT_U(ccc,aaa,ddd,2)
        ! 2 ---- W(b,c,a) <---- U(c,a,b,1)
        w_tmp(i_abc,2) = w_tmp(i_abc,2) + CC_PT_U(ccc,aaa,ddd,1)
        ! 3 ---- W(c,a,b) <---- U(c,a,b,3)
        w_tmp(i_abc,3) = w_tmp(i_abc,3) + CC_PT_U(ccc,aaa,ddd,3)
        ! 4 ---- W(a,c,b) <---- U(a,c,b,3)
        w_tmp(i_abc,4) = w_tmp(i_abc,4) + CC_PT_U(aaa,ccc,ddd,3)
        ! 5 ---- W(b,a,c) <---- U(a,c,b,1)
        w_tmp(i_abc,5) = w_tmp(i_abc,5) + CC_PT_U(aaa,ccc,ddd,1)
        ! 6 ---- W(c,b,a) <---- U(a,c,b,2)
        w_tmp(i_abc,6) = w_tmp(i_abc,6) + CC_PT_U(aaa,ccc,ddd,2)
      end if

      ddd = ccc - a_start + 1
      if ((ddd.ge.1).and.(ddd.le.n_a)) then
        ! 1 ---- W(a,b,c) <---- U(a,b,c,3)
        w_tmp(i_abc,1) = w_tmp(i_abc,1) + CC_PT_U(aaa,bbb,ddd,3)
        ! 2 ---- W(b,c,a) <---- U(a,b,c,2)
        w_tmp(i_abc,2) = w_tmp(i_abc,2) + CC_PT_U(aaa,bbb,ddd,2)
        ! 3 ---- W(c,a,b) <---- U(a,b,c,1)
        w_tmp(i_abc,3) = w_tmp(i_abc,3) + CC_PT_U(aaa,bbb,ddd,1)
        ! 4 ---- W(a,c,b) <---- U(b,a,c,2)
        w_tmp(i_abc,4) = w_tmp(i_abc,4) + CC_PT_U(bbb,aaa,ddd,2)
        ! 5 ---- W(b,a,c) <---- U(b,a,c,3)
        w_tmp(i_abc,5) = w_tmp(i_abc,5) + CC_PT_U(bbb,aaa,ddd,3)
        ! 6 ---- W(c,b,a) <---- U(b,a,c,1)
        w_tmp(i_abc,6) = w_tmp(i_abc,6) + CC_PT_U(bbb,aaa,ddd,1)
      end if

      !write(60+myid,*) 'W' 
      !write(60+myid,*) iii,jjj,kkk,aaa,bbb,ccc,W_tmp(i_abc,1,i_use)
      !write(60+myid,*) iii,jjj,kkk,bbb,ccc,aaa,W_tmp(i_abc,2,i_use) 
      !write(60+myid,*) iii,jjj,kkk,ccc,aaa,bbb,W_tmp(i_abc,3,i_use)
      !write(60+myid,*) iii,jjj,kkk,aaa,ccc,bbb,W_tmp(i_abc,4,i_use)
      !write(60+myid,*) iii,jjj,kkk,bbb,aaa,ccc,W_tmp(i_abc,5,i_use)
      !write(60+myid,*) iii,jjj,kkk,ccc,bbb,aaa,W_tmp(i_abc,6,i_use)

    end do
    !$OMP END DO
    !$OMP END PARALLEL

    Call get_timestamps(ctime_now, wtime_now)
    ctime_calc = ctime_calc + ctime_now - ctime_pre
    wtime_calc = wtime_calc + wtime_now - wtime_pre

    if (i_task_run.ne.1) then

      Call get_timestamps(ctime_pre, wtime_pre)

      Call MPI_WAIT(req1,stat1,errnum)
      Call MPI_WAIT(req2,stat2,errnum)

      Call get_timestamps(ctime_now, wtime_now)
      ctime_comm = ctime_comm + ctime_now - ctime_pre
      wtime_comm = wtime_comm + wtime_now - wtime_pre

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(i_abc)
      !$OMP DO
      do  i_abc = 1, n_abc
        CC_PT_W_abc(i_abc,:,1) = CC_PT_W_abc(i_abc,:,1) &
                               + w_recv(i_abc,:)
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      Deallocate(w_send,w_recv)

    end if

    if (i_task_run.ne.CC_mpi_domain_size) then

      Allocate(w_send(CC_mem_abc(i_task),6),stat=errnum)
      Call check_allocation(errnum,'w_send in CC')

      w_send = w_tmp

      Allocate(w_recv(CC_mem_abc(CC_mpi_did+1),6),stat=errnum)
      Call check_allocation(errnum,'w_recv in CC')

      s_tmp = Int(CC_mem_abc(i_task),8) * 6

      Call CC_mpi_real_isend(s_tmp,w_send,i_send-1,200,req1,CC_mpi_comm_domain)

      s_tmp = Int(CC_mem_abc(CC_mpi_did+1),8) * 6
      Call CC_mpi_real_irecv(s_tmp,w_recv,i_recv-1,200,req2,CC_mpi_comm_domain)

    else

      !$OMP PARALLEL Default(Shared) &
      !$OMP Private(i_abc)
      !$OMP DO
      do  i_abc = 1, n_abc
        CC_PT_W_abc(i_abc,:,1) = CC_PT_W_abc(i_abc,:,1) &
                               + w_tmp(i_abc,:)
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end if

    Deallocate(w_tmp)

  end do

  !code_abc = CC_index_abc(CC_mpi_did+1) - 1
  !do i_abc = 1, CC_mem_abc(CC_mpi_did+1)
  !  code_abc = code_abc + 1
  !  Call CC_cl_decode_ijk(code_abc,aaa,bbb,ccc)
  !  write(60+myid,*) iii,jjj,kkk,aaa,bbb,ccc,CC_PT_W_abc(i_abc,1,1)
  !  write(60+myid,*) iii,jjj,kkk,bbb,ccc,aaa,CC_PT_W_abc(i_abc,2,1) 
  !  write(60+myid,*) iii,jjj,kkk,ccc,aaa,bbb,CC_PT_W_abc(i_abc,3,1)
  !  write(60+myid,*) iii,jjj,kkk,aaa,ccc,bbb,CC_PT_W_abc(i_abc,4,1)
  !  write(60+myid,*) iii,jjj,kkk,bbb,aaa,ccc,CC_PT_W_abc(i_abc,5,1)
  !  write(60+myid,*) iii,jjj,kkk,ccc,bbb,aaa,CC_PT_W_abc(i_abc,6,1)
  !end do

  End Subroutine CC_cl_PT_sync_W

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_PT_calc_corr(iii,jjj,kkk,rlt)

  Use CC_cl

  Implicit None

  Integer , intent(in) :: iii,jjj,kkk
  Double precision :: delta_ijk,D_ijk,rlt

  Integer :: aaa,bbb,ccc,i_abc,n_abc,i_run,code_abc,abc_start,outline
  Double precision :: D_abc,delta_abc,xxx,yyy,zzz,www1,www2,rlt1,rlt2
  Double precision , dimension(6) :: V_abc

  n_abc = CC_mem_abc(CC_mpi_did+1)
  outline = Int(n_abc/2) + mod(n_abc,2)
  abc_start = CC_index_abc(CC_mpi_did+1)

  Call CC_cl_D_ijk(iii,jjj,kkk,delta_ijk,D_ijk)

  rlt = 0.0D0
  !$OMP PARALLEL Default(shared) Reduction(+:rlt) &
  !$OMP Private(i_abc,code_abc,aaa,bbb,ccc,i_run,delta_abc,D_abc,V_abc, &
  !$OMP         xxx,yyy,zzz,www1,www2,rlt1)
  !$OMP DO Schedule(Dynamic)
  do i_abc = 1, n_abc

    code_abc = abc_start - 1 + i_abc

    Call CC_cl_decode_ijk(code_abc,aaa,bbb,ccc)
    Call CC_cl_D_abc(aaa,bbb,ccc,D_abc)

    delta_abc = 1.0D0
    if (aaa.eq.bbb) then
      delta_abc = delta_abc + 1.0D0
    end if

    if (bbb.eq.ccc) then
      delta_abc = delta_abc + 1.0D0
    end if

    !write(myid+70,*) 'term',iii,jjj,kkk,aaa,bbb,ccc
    !print*,'W'
    !print*,iii,jjj,kkk,aaa,bbb,ccc,CC_PT_W_abc(iii,jjj,kkk) 
    !print*,iii,kkk,jjj,aaa,bbb,ccc,CC_PT_W_abc(iii,kkk,jjj)
    !print*,jjj,iii,kkk,aaa,bbb,ccc,CC_PT_W_abc(jjj,iii,kkk)
    !print*,kkk,iii,jjj,aaa,bbb,ccc,CC_PT_W_abc(kkk,iii,jjj)
    !print*,jjj,kkk,iii,aaa,bbb,ccc,CC_PT_W_abc(jjj,kkk,iii)
    !print*,kkk,jjj,iii,aaa,bbb,ccc,CC_PT_W_abc(kkk,jjj,iii)

    ! get V_ijk
    Call CC_cl_PT_calc_V(V_abc(1),iii,jjj,kkk,aaa,bbb,ccc)
    Call CC_cl_PT_calc_V(V_abc(2),iii,jjj,kkk,bbb,ccc,aaa)
    Call CC_cl_PT_calc_V(V_abc(3),iii,jjj,kkk,ccc,aaa,bbb)
    Call CC_cl_PT_calc_V(V_abc(4),iii,jjj,kkk,aaa,ccc,bbb)
    Call CC_cl_PT_calc_V(V_abc(5),iii,jjj,kkk,bbb,aaa,ccc)
    Call CC_cl_PT_calc_V(V_abc(6),iii,jjj,kkk,ccc,bbb,aaa)

    xxx = 0.0D0
    do i_run = 1, 6
      V_abc(i_run) = (V_abc(i_run) + CC_PT_W_abc(i_abc,i_run,1)) / delta_abc
      xxx = xxx + V_abc(i_run) * CC_PT_W_abc(i_abc,i_run,1)
    end do

    www1 = CC_PT_W_abc(i_abc,1,1) + CC_PT_W_abc(i_abc,2,1) + CC_PT_W_abc(i_abc,3,1)
    www2 = CC_PT_W_abc(i_abc,4,1) + CC_PT_W_abc(i_abc,5,1) + CC_PT_W_abc(i_abc,6,1)

    yyy = V_abc(1) + V_abc(2) + V_abc(3)
    zzz = V_abc(4) + V_abc(5) + V_abc(6)

    !print*,'V'
    !print*,iii,jjj,kkk,aaa,bbb,ccc,V_ijk(1)
    !print*,iii,kkk,jjj,aaa,bbb,ccc,V_ijk(2)
    !print*,jjj,iii,kkk,aaa,bbb,ccc,V_ijk(3)
    !print*,kkk,iii,jjj,aaa,bbb,ccc,V_ijk(4)
    !print*,jjj,kkk,iii,aaa,bbb,ccc,V_ijk(5)
    !print*,kkk,jjj,iii,aaa,bbb,ccc,V_ijk(6)

    !if ((iii.eq.4).and.(jjj.eq.3).and.(kkk.eq.1)) then
    !  if ((aaa.eq.13).and.(bbb.eq.6).and.(ccc.eq.1)) then
    !    print*,'V'
    !    print*,iii,jjj,kkk,aaa,bbb,ccc,V_abc(1) 
    !    print*,iii,jjj,kkk,bbb,ccc,aaa,V_abc(2)
    !    print*,iii,jjj,kkk,ccc,aaa,bbb,V_abc(3)
    !    print*,iii,jjj,kkk,aaa,ccc,bbb,V_abc(4)
    !    print*,iii,jjj,kkk,bbb,aaa,ccc,V_abc(5)
    !    print*,iii,jjj,kkk,ccc,bbb,aaa,V_abc(6)
   

    !    print*, 'XYZ',xxx,yyy,zzz

    !    print*, 'W1',www1,'W2',www2
    !    print*, 'delta_ijk',delta_ijk
    !    print*, 'Dijkabc',D_ijk - D_abc

    !    print*, 'XYZ',xxx,yyy,zzz

    !    print*, 'W1',www1,'W2',www2
    !    print*, 'delta_ijk',delta_ijk
    !    print*, 'Dijkabc',D_ijk - D_abc
    !  end if
    !end if
 
    rlt1 =  ((yyy - 2.0D0 * zzz) * www1 + (zzz - 2.0D0 * yyy) * www2 &
             + 3.0D0 * xxx) * delta_ijk / (D_ijk - D_abc)
    !write(myid+70,*) 'rlt',rlt1

    !write(80+myid,*) iii,jjj,kkk,aaa,bbb,ccc,rlt1
    rlt = rlt + rlt1
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !write(50+myid,*) 'E_c',iii,jjj,kkk,rlt

  End Subroutine CC_cl_PT_calc_corr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_PT_calc_V(V_abc,iii,jjj,kkk,aaa,bbb,ccc)

  Use CC_cl

  Implicit None

  Double precision , intent(out) :: V_abc
  Integer , intent(in) :: iii,jjj,kkk,aaa,bbb,ccc

  V_abc = CC_PT_iajb_abc(bbb,ccc,1) * CC_t_s(iii,aaa) &
        + CC_PT_iajb_abc(aaa,ccc,2) * CC_t_s(jjj,bbb) &
        + CC_PT_iajb_abc(aaa,bbb,3) * CC_t_s(kkk,ccc)

  !if ((iii.eq.4).and.(jjj.eq.3).and.(kkk.eq.1)) then
  !  if ((aaa.eq.1).and.(bbb.eq.6).and.(ccc.eq.13)) then
 
  !    print*,'iajb',iii,jjj,kkk,aaa,bbb,ccc
  !    print*,jjj,kkk,bbb,ccc,CC_PT_iajb_abc(bbb,ccc,1),CC_t_s(iii,aaa)
  !    print*,iii,kkk,aaa,ccc,CC_PT_iajb_abc(aaa,ccc,2),CC_t_s(jjj,bbb)
  !    print*,iii,jjj,aaa,bbb,CC_PT_iajb_abc(aaa,bbb,3),CC_t_s(kkk,ccc)
  !  end if
  !end if

  End Subroutine CC_cl_PT_calc_V

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_PT_output(i_calc,CC_PT_outstep,out_flag,CC_elapsed_time)

  Use CC_cl

  Implicit None

  Integer , intent(in) :: i_calc

  Double precision :: CC_elapsed_time

  Integer , parameter :: ures = 100

  Integer :: iii,i_tmp
  Double precision :: CC_PT_percent
  Integer , dimension(CC_PT_n_out) :: CC_PT_outstep
  Logical , dimension(CC_PT_n_out) :: out_flag

  CC_PT_percent = dble(i_calc) / dble(CC_PT_outstep(CC_PT_n_out)) * 100.0D0

  do iii = 1, CC_PT_n_out
    if ((i_calc.ge.CC_PT_outstep(iii)).and.(out_flag(iii))) then
      out_flag(iii) = .false.
      if (myid.eq.0) then
        write(use_unit,"(2X,A,F6.2,A,F12.3,A)") 'Perturbative triples calculation: ', &
                                        CC_PT_percent,'% completed, elapsed time: ', &
                                        CC_elapsed_time,' s'
      end if

      CC_elapsed_time = 0.0D0

      if (CC_restart_flag) then
        CC_res_inf = 0
        CC_res_inf(CC_mpi_gid+2) = i_calc
        Call CC_mpi_int_allreduce(CC_n_domain+1,CC_res_inf,CC_mpi_comm_group)
        Call CC_mpi_real_number(CC_E_PT,MPI_COMM_WORLD)
        CC_res_inf(1) = -1
        if (myid.eq.0) then
          open(unit=ures,file='CCrestart')
          do i_tmp = 1, CC_n_domain+1
            write(ures,*) CC_res_inf(i_tmp)
          end do
          write(ures,*) CC_E_PT
          close(ures)
        else
          CC_E_PT = 0.0D0
        end if
        if (iii.eq.CC_max_cycle) then
          Call aims_stop
        end if
      end if
    end if
  end do

  End Subroutine CC_cl_PT_output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_PT_send(iii,jjj,kkk,req_t,req_jlkc,req_iajb)

  Use CC_cl

  Implicit None

  Integer , intent(in) :: iii,jjj,kkk
  Integer , dimension(CC_mpi_domain_size,2) :: req_t,req_jlkc,req_iajb

  Integer , parameter :: tag_t = 100
  Integer , parameter :: tag_jlkc = 101
  Integer , parameter :: tag_iajb = 102

  Integer :: i_task,i_task_run,i_run,i_send,i_recv
  Integer :: i_abc,i_tmp,abc_end,n_a
  Integer :: aaa,bbb,ccc

  Integer :: errnum
  Integer (kind=8) :: s_tmp

  n_a = CC_mem_aa_D(CC_mpi_did+1)

  req_t = MPI_REQUEST_NULL
  req_jlkc = MPI_REQUEST_NULL
  req_iajb = MPI_REQUEST_NULL

  CC_PT_T_send = 0.0D0
  CC_PT_jlkc_send = 0.0D0
  CC_PT_iajb_send = 0.0D0
  CC_PT_T_recv = 0.0D0
  CC_PT_jlkc_recv = 0.0D0
  CC_PT_iajb_recv = 0.0D0

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

    if (CC_mpi_did+1.ne.i_send) then

      ! get (ia|jb)
      CC_PT_iajb_send(:,1:n_a,1,i_send) = CC_intl_iajb_A(jjj,kkk,:,:)
      CC_PT_iajb_send(:,1:n_a,2,i_send) = CC_intl_iajb_A(iii,kkk,:,:)
      CC_PT_iajb_send(:,1:n_a,3,i_send) = CC_intl_iajb_A(iii,jjj,:,:)

      i_tmp = CC_n_vir * CC_mem_aa_D(1) * 3
      s_tmp = Int(i_tmp,8)

      Call CC_mpi_real_isend(s_tmp,CC_PT_iajb_send(:,:,:,i_send),i_send-1, &
                             tag_iajb,req_iajb(i_send,1),CC_mpi_comm_domain)

      Call CC_mpi_real_irecv(s_tmp,CC_PT_iajb_recv(:,:,:,i_recv),i_recv-1, &
                             tag_iajb,req_iajb(i_recv,2),CC_mpi_comm_domain)

      ! get T
      CC_PT_T_send(:,1:n_a,1,i_send,1) = CC_t_d_I(jjj,:,:,kkk)
      CC_PT_T_send(:,1:n_a,2,i_send,1) = CC_t_d_I(kkk,:,:,jjj)
      CC_PT_T_send(:,1:n_a,3,i_send,1) = CC_t_d_I(iii,:,:,jjj)
      CC_PT_T_send(:,1:n_a,4,i_send,1) = CC_t_d_I(jjj,:,:,iii)
      CC_PT_T_send(:,1:n_a,5,i_send,1) = CC_t_d_I(kkk,:,:,iii)
      CC_PT_T_send(:,1:n_a,6,i_send,1) = CC_t_d_I(iii,:,:,kkk)

      i_tmp = CC_n_vir * CC_mem_aa_D(1) * 6
      s_tmp = Int(i_tmp,8)

      Call CC_mpi_real_isend(s_tmp,CC_PT_T_send(:,:,:,i_send,1),i_send-1, &
                             tag_t,req_t(i_send,1),CC_mpi_comm_domain)

      Call CC_mpi_real_irecv(s_tmp,CC_PT_T_recv(:,:,:,i_recv,1),i_recv-1, &
                             tag_t,req_t(i_recv,2),CC_mpi_comm_domain)

      ! get (jl|kc)
      CC_PT_jlkc_send(:,1:n_a,1,i_send) =  CC_intl_jlkc(:,:,jjj,kkk)
      CC_PT_jlkc_send(:,1:n_a,2,i_send) =  CC_intl_jlkc(:,:,kkk,jjj)
      CC_PT_jlkc_send(:,1:n_a,3,i_send) =  CC_intl_jlkc(:,:,iii,jjj)
      CC_PT_jlkc_send(:,1:n_a,4,i_send) =  CC_intl_jlkc(:,:,jjj,iii)
      CC_PT_jlkc_send(:,1:n_a,5,i_send) =  CC_intl_jlkc(:,:,kkk,iii)
      CC_PT_jlkc_send(:,1:n_a,6,i_send) =  CC_intl_jlkc(:,:,iii,kkk)

      i_tmp = CC_n_occ * CC_mem_aa_D(1) * 6 
      s_tmp = Int(i_tmp,8)

      Call CC_mpi_real_isend(s_tmp,CC_PT_jlkc_send(:,:,:,i_send),i_send-1, &
                             tag_jlkc,req_jlkc(i_send,1),CC_mpi_comm_domain)

      Call CC_mpi_real_irecv(s_tmp,CC_PT_jlkc_recv(:,:,:,i_recv),i_recv-1, &
                             tag_jlkc,req_jlkc(i_recv,2),CC_mpi_comm_domain)

    else

      ! get (ia|jb)
      CC_PT_iajb_recv(:,1:n_a,1,i_send) = CC_intl_iajb_A(jjj,kkk,:,:)
      CC_PT_iajb_recv(:,1:n_a,2,i_send) = CC_intl_iajb_A(iii,kkk,:,:)
      CC_PT_iajb_recv(:,1:n_a,3,i_send) = CC_intl_iajb_A(iii,jjj,:,:)

      ! get T
      CC_PT_T_recv(:,1:n_a,1,i_send,1) = CC_t_d_I(jjj,:,:,kkk)
      CC_PT_T_recv(:,1:n_a,2,i_send,1) = CC_t_d_I(kkk,:,:,jjj)
      CC_PT_T_recv(:,1:n_a,3,i_send,1) = CC_t_d_I(iii,:,:,jjj)
      CC_PT_T_recv(:,1:n_a,4,i_send,1) = CC_t_d_I(jjj,:,:,iii)
      CC_PT_T_recv(:,1:n_a,5,i_send,1) = CC_t_d_I(kkk,:,:,iii)
      CC_PT_T_recv(:,1:n_a,6,i_send,1) = CC_t_d_I(iii,:,:,kkk)

      ! get (jl|kc)
      CC_PT_jlkc_recv(:,1:n_a,1,i_send) =  CC_intl_jlkc(:,:,jjj,kkk)
      CC_PT_jlkc_recv(:,1:n_a,2,i_send) =  CC_intl_jlkc(:,:,kkk,jjj)
      CC_PT_jlkc_recv(:,1:n_a,3,i_send) =  CC_intl_jlkc(:,:,iii,jjj)
      CC_PT_jlkc_recv(:,1:n_a,4,i_send) =  CC_intl_jlkc(:,:,jjj,iii)
      CC_PT_jlkc_recv(:,1:n_a,5,i_send) =  CC_intl_jlkc(:,:,kkk,iii)
      CC_PT_jlkc_recv(:,1:n_a,6,i_send) =  CC_intl_jlkc(:,:,iii,kkk)

    end if
  end do

  End Subroutine CC_cl_PT_send

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_PT_recv(req_t,req_jlkc,req_iajb)

  Use CC_cl

  Implicit None

  Integer , dimension(CC_mpi_domain_size,2) :: req_t,req_jlkc,req_iajb
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Integer :: i_task,n_a,a_start,a_end

  Integer :: errnum

  do i_task = 1, CC_mpi_domain_size
    if (CC_mpi_did+1.ne.i_task) then
      Call MPI_WAIT(req_t(i_task,1),stat1,errnum)
      Call MPI_WAIT(req_t(i_task,2),stat2,errnum)
      Call MPI_WAIT(req_jlkc(i_task,1),stat1,errnum)
      Call MPI_WAIT(req_jlkc(i_task,2),stat2,errnum)
      Call MPI_WAIT(req_iajb(i_task,1),stat1,errnum)
      Call MPI_WAIT(req_iajb(i_task,2),stat2,errnum)
    end if
  end do

  do i_task = 1, CC_mpi_domain_size

    n_a = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task)
    a_end = a_start - 1 + n_a

    CC_PT_T_abc(:,a_start:a_end,:,1) = CC_PT_T_recv(:,1:n_a,:,i_task,1)
    CC_PT_jlkc_abc(:,a_start:a_end,:) = CC_PT_jlkc_recv(:,1:n_a,:,i_task)
    CC_PT_iajb_abc(:,a_start:a_end,:) = CC_PT_iajb_recv(:,1:n_a,:,i_task)

  end do

  End Subroutine CC_cl_PT_recv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_PT_redistribution()

  Use CC_cl

  Implicit None

  Integer :: n_a,a_start,a_end,n_b,b_start,b_end
  Integer :: i_task,code_ia,i_ind,code_kl,i_tmp,j_tmp,nnn,a_run
  Integer :: n_i,i_start,i_end,i_run
  Integer :: kkk,lll,ccc,ddd,aaa,bbb,iii,jjj
  Integer (kind=8) :: s_tmp
  Integer :: source_id,errnum,target_id

  Double precision :: alpha,beta
  Double precision , dimension(:,:,:) , allocatable :: mat_A,mat_B
  Double precision , dimension(:,:,:,:) , allocatable :: intl_rlt

  Double precision , dimension(:,:,:,:) , allocatable :: CC_comm_tmp

  ! CC_intl_jlkc(l,c,j,k) ---- (jl|kc)
  Allocate(CC_intl_jlkc(CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),CC_n_occ, & 
                        CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'CC_intl_jlkc in CC')

  ! CC_intl_dbia(d,b,a,i) ---- (db|ia)
  Allocate(CC_intl_dbia(CC_n_vir,CC_n_vir,CC_mem_aa_D(CC_mpi_did+1), &
                        CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'CC_intl_dbia in CC')

  ! CC_t_d_I(j,b,a,i)
  Allocate(CC_t_d_I(CC_n_occ,CC_n_vir,CC_mem_aa_D(CC_mpi_did+1), &
                    CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'CC_t_d_I in CC')

  ! CC_intl_iajb(i,j,a,b) ---- (ia|jb)
  Allocate(CC_intl_iajb_A(CC_n_occ,CC_n_occ,CC_n_vir, &
                          CC_mem_aa_D(CC_mpi_did+1)),stat=errnum)
  Call check_allocation(errnum,'CC_intl_iajb in CC')

  n_a = CC_mem_aa_D(CC_mpi_did+1)
  a_start = CC_index_aa_D(CC_mpi_did+1)
  a_end = a_start - 1 + n_a

  ! redistribution of t(j,b,a,i)
  do iii = 1, CC_n_occ
    do bbb = 1, CC_n_vir
      CC_t_d_I(:,bbb,:,iii) = CC_t_d(iii,:,:,bbb)
    end do
  end do

  ! redistibution of (ia|jb)
  do iii = 1, CC_n_occ
    do aaa = 1, CC_n_vir
      CC_intl_iajb_A(iii,:,aaa,:) = CC_intl_iajb(:,iii,:,aaa)
    end do
  end do

  ! redistribution of (ki|lc) (k,c,i,l) ----> (jl|kc) (l,c,j,k)
  do ccc = 1, n_a
    do lll = 1, CC_n_occ
      CC_intl_jlkc(lll,ccc,:,:) = CC_intl_kilc(:,ccc,lll,:)
    end do
  end do

  ! Calculate (db|ia)
  nnn = CC_mem_bas(CC_mpi_did+1) 
  n_b = CC_mem_aa_G(CC_mpi_gid+1)
  b_start = CC_index_aa_G(CC_mpi_gid+1)
  b_end = b_start - 1 + n_b

  CC_intl_dbia = 0.0D0

  !print*,'n_b',n_b
  !print*,'b_start',b_start
  !print*,'b_end',b_end

  alpha = 1.0D0
  beta = 0.0D0

  do i_task = 1, CC_mpi_domain_size

    source_id = i_task - 1

    n_a = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task)
    a_end = a_start - 1 + n_a

    i_tmp = CC_n_vir * n_b
    j_tmp = CC_n_occ * n_a

    Allocate(mat_A(nnn,CC_n_vir,n_b), stat=errnum)
    Call check_allocation(errnum,'mat_A in CC')

    Allocate(mat_B(nnn,CC_n_occ,n_a), stat=errnum)
    Call check_allocation(errnum,'mat_B in CC')

    Allocate(intl_rlt(CC_n_vir,n_b,CC_n_occ,n_a), stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')

    mat_A = CC_RI(:,CC_n_occ+1:CC_n_occ+CC_n_vir,CC_n_occ+b_start:CC_n_occ+b_end,1)
    mat_B = CC_RI(:,1:CC_n_occ,CC_n_occ+a_start:CC_n_occ+a_end,1)

    Call Dgemm('T','N',i_tmp,j_tmp,nnn,alpha,mat_A,nnn, &
               mat_B,nnn,beta,intl_rlt,i_tmp)

    target_id = i_task - 1
    s_tmp = Int(i_tmp,8) * Int(j_tmp,8)
    Call CC_mpi_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      ! (d,b,i,a) ----> (d,b,a,i)
      do iii = 1, CC_n_occ
        CC_intl_dbia(:,b_start:b_end,:,iii) = intl_rlt(:,:,iii,:)
      end do
    end if

    Deallocate(mat_A,mat_B,intl_rlt)
  end do

  i_tmp = CC_n_vir * CC_n_vir
  j_tmp = CC_n_occ * CC_mem_aa_D(CC_mpi_did+1)
  s_tmp = Int(i_tmp,8) * Int(j_tmp,8)
  Call CC_mpi_allreduce(s_tmp, CC_intl_dbia, CC_mpi_comm_group)

  !! CC_t_d_I(j,b,a,i)
  !write(myid+70,*) 't'
  !do iii = 1, CC_n_occ
  !  aaa = CC_index_aa_D(CC_mpi_did+1) - 1
  !  do i_run = 1, CC_mem_aa_D(CC_mpi_did+1)
  !    aaa = aaa + 1
  !    do bbb = 1, CC_n_vir
  !      do jjj = 1, CC_n_occ
  !        write(myid+70,*) jjj,bbb,aaa,iii,CC_t_d_I(jjj,bbb,i_run,iii) 
  !      end do
  !    end do
  !  end do
  !end do

  !! CC_intl_dbia(d,b,a,i) ---- (db|ia)
  !write(myid+70,*) 'dbia'
  !do iii = 1, CC_n_occ
  !  aaa = CC_index_aa_D(CC_mpi_did+1) - 1
  !  do i_run = 1, CC_mem_aa_D(CC_mpi_did+1)
  !    aaa = aaa + 1
  !    do bbb = 1, CC_n_vir
  !      do ddd = 1, CC_n_vir
  !        write(myid+70,*) ddd,bbb,aaa,iii,CC_intl_dbia(ddd,bbb,i_run,iii)
  !      end do
  !    end do
  !  end do
  !end do

  !! CC_intl_jlkc(l,c,j,k) ---- (jl|kc)
  !write(myid+70,*) 'jlkc'
  !do kkk = 1, CC_n_occ
  !  do jjj = 1, CC_n_occ
  !    ccc = CC_index_aa_D(CC_mpi_did+1) - 1
  !    do i_run = 1, CC_mem_aa_D(CC_mpi_did+1)
  !      ccc = ccc + 1
  !      do lll = 1, CC_n_occ
  !        write(myid+70,*) lll,ccc,jjj,kkk,CC_intl_jlkc(lll,i_run,jjj,kkk)
  !      end do
  !    end do
  !  end do
  !end do

  !! CC_intl_iajb(i,j,a,b) ---- (ia|jb)
  !write(myid+70,*) 'iajb'
  !do iii = 1, CC_n_occ
  !  do jjj = 1, CC_n_occ
  !    do aaa = 1, CC_n_vir
  !      bbb = CC_index_aa_D(CC_mpi_did+1) - 1
  !      do i_run = 1, CC_mem_aa_D(CC_mpi_did+1)
  !        bbb = bbb + 1
  !        write(myid+70,*) iii,jjj,aaa,bbb,CC_intl_iajb(iii,jjj,aaa,i_run)
  !      end do
  !    end do
  !  end do
  !end do

  !! CC_t_d_I(j,b,a,i)
  !aaa = CC_index_aa_D(CC_mpi_did+1) - 1
  !do i_run = 1, CC_mem_aa_D(CC_mpi_did+1)
  !  aaa = aaa + 1
  !  do iii = 1, CC_n_occ
  !    do bbb = 1, CC_n_vir
  !      do jjj = 1, CC_n_occ
  !        write(myid+70,*) jjj,bbb,aaa,iii,CC_t_d_I(jjj,bbb,i_run,iii) 
  !      end do
  !    end do
  !  end do
  !end do

  ! CC_intl_dbia(d,b,a,i) ---- (db|ia)
  !aaa = CC_index_aa_D(CC_mpi_did+1) - 1
  !do i_run = 1, CC_mem_aa_D(CC_mpi_did+1)
  !  aaa = aaa + 1
  !  do iii = 1, CC_n_occ
  !    do bbb = 1, CC_n_vir
  !      do ddd = 1, CC_n_vir
  !        write(myid+80,*) ddd,bbb,aaa,iii,CC_intl_dbia(ddd,bbb,i_run,iii)
  !      end do
  !    end do
  !  end do
  !end do

  !! CC_intl_jlkc(l,c,j,k) ---- (jl|kc)
  !ccc = CC_index_aa_D(CC_mpi_did+1) - 1
  !do i_run = 1, CC_mem_aa_D(CC_mpi_did+1)
  !  ccc = ccc + 1
  !  do kkk = 1, CC_n_occ
  !    do jjj = 1, CC_n_occ
  !      do lll = 1, CC_n_occ
  !        write(myid+90,*) lll,ccc,jjj,kkk,CC_intl_jlkc(lll,i_run,jjj,kkk)
  !      end do
  !    end do
  !  end do
  !end do

  !! CC_intl_iajb(i,j,a,b) ---- (ia|jb)
  !bbb = CC_index_aa_D(CC_mpi_did+1) - 1
  !do i_run = 1, CC_mem_aa_D(CC_mpi_did+1)
  !  bbb = bbb + 1
  !  do iii = 1, CC_n_occ
  !    do jjj = 1, CC_n_occ
  !      do aaa = 1, CC_n_vir
  !        write(myid+100,*) iii,jjj,aaa,bbb,CC_intl_iajb(iii,jjj,aaa,i_run)
  !      end do
  !    end do
  !  end do
  !end do

  !a_start = CC_index_aa_D(CC_mpi_did+1)
  !a_end = CC_index_aa_D(CC_mpi_did+1) - 1 + CC_mem_aa_D(CC_mpi_did+1)
  !if ((a_start.le.13).and.(a_end.ge.13)) then
  !  print*,'did',CC_mpi_did
  !  print*,'iajb',4,1,1,13,CC_intl_iajb(4,1,1,13-a_start+1)
  !end if

  End Subroutine CC_cl_PT_redistribution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

