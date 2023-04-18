!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_calc_w()

    Use CC_cl
    Use timing

    Implicit None

    Integer :: errnum
    Integer :: j_tmp,s_tmp
    Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd
    Integer :: n_b,b_start,b_end,i_run,n_a,a_start,a_end,a_run,b_run,c_run
    Integer :: n_ab,ab_start,ab_end,ab_run,ab,ij_run,n_ij,ij_start,ij_end,i_ij
    Double precision :: ctime, wtime, c_pre, w_pre

    n_b = CC_mem_aa_D(CC_mpi_did+1)
    b_start = CC_index_aa_D(CC_mpi_did+1)
    b_end = b_start - 1 + n_b

    n_a = CC_mem_aa_G(CC_mpi_gid+1)
    a_start = CC_index_aa_G(CC_mpi_gid+1)
    a_end = a_start - 1 + n_a

    n_ab = CC_mem_ab_G(CC_mpi_gid+1)
    ab_start = CC_index_ab_G(CC_mpi_gid+1)
    ab_end = ab_start - 1 + n_ab

    n_ij = CC_mem_kl_D(CC_mpi_did+1)
    ij_start = CC_index_kl_D(CC_mpi_did+1)
    ij_end = ij_start - 1 + n_ij

    !write(70+myid,*) 'start'

    Call CC_cl_clean_w()

    ! Calculate intermediate arrays h_ik, h_ca, g_ik, g_ca, and h_ck
    Call get_timestamps(c_pre, w_pre)
    Call CC_cl_aux_h_g()
    Call get_timestamps(ctime, wtime)
    if (myid.eq.0) then
      write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of evaluating h & g', &
                                               ctime - c_pre, 's (cpu)', &
                                               wtime - w_pre, 's (wall clock)'
    end if

    !if (myid.eq.0) then
    !  print*,'h_ik'
    !  do iii = 1, CC_n_occ
    !    do kkk = 1, CC_n_occ
    !      print*,iii,kkk,CC_h_ik(iii,kkk)
    !    end do
    !  end do

    !  print*,'h_ca'
    !  do ccc = 1, CC_n_vir
    !    do aaa = 1, CC_n_vir
    !      print*,ccc,aaa,CC_h_ca(ccc,aaa)
    !    end do
    !  end do

    !  print*,'h_ck'
    !  do kkk = 1, CC_n_occ
    !    do ccc = 1, CC_n_vir
    !      print*,kkk,ccc,CC_h_ck(kkk,ccc)
    !    end do
    !  end do

    !  print*,'g_ik'
    !  do iii = 1, CC_n_occ
    !    do kkk = 1, CC_n_occ
    !      print*,iii,kkk,CC_g_ik(iii,kkk)
    !    end do
    !  end do

    !  print*,'g_ca'
    !  do ccc = 1, CC_n_vir
    !    do aaa = 1, CC_n_vir
    !      print*,ccc,aaa,CC_g_ca(ccc,aaa)
    !   end do
    !  end do
    !end if


    ! Calculate w_s
    Call get_timestamps(c_pre, w_pre)

    Call CC_cl_calc_w_s()

    Call get_timestamps(ctime, wtime)
    if (myid.eq.0) then
      write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of evaluating singlet excitaions', &
                                               ctime - c_pre, 's (cpu)', &
                                               wtime - w_pre, 's (wall clock)'
    end if

    !if (myid.eq.0) then
    !  print*,'w_s'
    !  do iii = 1, CC_n_occ
    !    do aaa = 1, CC_n_vir
    !      print*,CC_w_s(iii,aaa)
    !    end do
    !  end do
    !end if

    ! Calculate w_d
    Call get_timestamps(c_pre, w_pre)

    Call CC_cl_calc_w_d()

    Call get_timestamps(ctime, wtime)
    if (myid.eq.0) then
      write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of evaluating mosaics', &
                                               ctime - c_pre, 's (cpu)', &
                                               wtime - w_pre, 's (wall clock)'
    end if

    !write(80+myid,*) 'w_d coeff add i'
    !do b_run = 1, n_b
    !  do a_run = 1, n_a
    !    do jjj = 1, CC_n_occ
    !      do iii = 1, CC_n_occ
    !        bbb = b_start - 1 + b_run 
    !        aaa = a_start - 1 + a_run
    !        write(80+myid,*) iii,jjj,aaa,bbb,CC_ten_con(iii,jjj,a_run,b_run)
    !      end do
    !    end do
    !  end do
    !end do

    ! Calculate intermediate array j
    Call get_timestamps(c_pre, w_pre)

    Call CC_cl_aux_j()

    Call get_timestamps(ctime, wtime)
    if (myid.eq.0) then
      write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of evaluating intermediate J', &
                                               ctime - c_pre, 's (cpu)', &
                                               wtime - w_pre, 's (wall clock)'
    end if

    !write(80+myid,*) 'j_aux'
    !do kkk = 1, CC_n_occ
    !  do c_run = 1, n_b
    !    do iii = 1, CC_n_occ
    !      do a_run = 1, n_a
    !        ccc = b_start - 1 + c_run
    !        aaa = a_start - 1 + a_run
    !        write(80+myid,*) kkk,ccc,iii,aaa,CC_j_aux(kkk,c_run,iii,a_run)
    !      end do
    !    end do
    !  end do
    !end do

    ! Add j to w_d
    Call get_timestamps(c_pre, w_pre)

    Call CC_cl_add_j2w()

    Call get_timestamps(ctime, wtime)
    if (myid.eq.0) then
      write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of evaluating ring term', &
                                               ctime - c_pre, 's (cpu)', &
                                               wtime - w_pre, 's (wall clock)'
    end if

    !write(80+myid,*) 'w_d coeff add j'
    !do b_run = 1, n_b
    !  do a_run = 1, n_a
    !    do jjj = 1, CC_n_occ
    !      do iii = 1, CC_n_occ
    !        bbb = b_start - 1 + b_run 
    !        aaa = a_start - 1 + a_run
    !        write(80+myid,*) iii,jjj,aaa,bbb,CC_ten_con(iii,jjj,a_run,b_run)
    !      end do
    !    end do
    !  end do
    !end do

    ! Calculate intermediate array k
    Call get_timestamps(c_pre, w_pre)

    Call CC_cl_aux_k()

    Call get_timestamps(ctime, wtime)
    if (myid.eq.0) then
      write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of evaluating intermediate K', &
                                               ctime - c_pre, 's (cpu)', &
                                               wtime - w_pre, 's (wall clock)'
    end if

    !write(80+myid,*) 'k_aux'
    !do kkk = 1, CC_n_occ
    !  do c_run = 1, n_b
    !    do iii = 1, CC_n_occ
    !      do a_run = 1, n_a
    !        ccc = b_start - 1 + c_run
    !        aaa = a_start - 1 + a_run
    !        write(80+myid,*) kkk,ccc,iii,aaa,CC_k_aux(kkk,c_run,iii,a_run)
    !      end do
    !    end do
    !  end do
    !end do

    ! Add k to w_d
    Call get_timestamps(c_pre, w_pre)

    Call CC_cl_add_k2w()

    Call get_timestamps(ctime, wtime)
    if (myid.eq.0) then
      write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of evaluating cross-ring term', &
                                               ctime - c_pre, 's (cpu)', &
                                               wtime - w_pre, 's (wall clock)'
    end if

    !write(80+myid,*) 'w_d coeff add k'
    !do b_run = 1, n_b
    !  do a_run = 1, n_a
    !    do jjj = 1, CC_n_occ
    !      do iii = 1, CC_n_occ
    !        bbb = b_start - 1 + b_run 
    !        aaa = a_start - 1 + a_run
    !        write(80+myid,*) iii,jjj,aaa,bbb,CC_ten_con(iii,jjj,a_run,b_run)
    !      end do
    !    end do
    !  end do
    !end do

    ! Calculate intermediate a
    Call get_timestamps(c_pre, w_pre)

    Call CC_cl_aux_a()

    Call get_timestamps(ctime, wtime)
    if (myid.eq.0) then
      write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of evaluating intermediate A', &
                                               ctime - c_pre, 's (cpu)', &
                                               wtime - w_pre, 's (wall clock)'
    end if

    !write(80+myid,*) 'a_aux'
    !do iii = 1, CC_n_occ
    !  do jjj = 1, CC_n_occ
    !    do ij_run = 1, n_ij
    !      i_ij = ij_start - 1 + ij_run
    !      Call CC_cl_decode(i_ij,kkk,lll,CC_n_occ,1)
    !      write(80+myid,*) iii,jjj,kkk,lll,CC_a_aux(ij_run,iii,jjj,1)
    !    end do
    !  end do
    !end do

    Call get_timestamps(c_pre, w_pre)

    Call CC_cl_add_a2w()

    Call get_timestamps(ctime, wtime)
    if (myid.eq.0) then
      write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of evaluating HH ladder term', &
                                               ctime - c_pre, 's (cpu)', &
                                               wtime - w_pre, 's (wall clock)'
    end if

    !write(80+myid,*) 'w_d coeff add a'
    !do b_run = 1, n_b
    !  do a_run = 1, n_a
    !    do jjj = 1, CC_n_occ
    !      do iii = 1, CC_n_occ
    !        bbb = b_start - 1 + b_run 
    !        aaa = a_start - 1 + a_run
    !        write(80+myid,*) iii,jjj,aaa,bbb,CC_ten_con(iii,jjj,a_run,b_run)
    !      end do
    !    end do
    !  end do
    !end do

    Call get_timestamps(c_pre, w_pre)

    Call CC_cl_add_b2w()

    Call get_timestamps(ctime, wtime)
    if (myid.eq.0) then
      write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of evaluating PP ladder term', &
                                               ctime - c_pre, 's (cpu)', &
                                               wtime - w_pre, 's (wall clock)'
    end if

    Call get_timestamps(c_pre, w_pre)

    Call CC_cl_add_tc2w()

    Call get_timestamps(ctime, wtime)
    if (myid.eq.0) then
      write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of communicating asymmetric term', &
                                               ctime - c_pre, 's (cpu)', &
                                               wtime - w_pre, 's (wall clock)'
    end if

  !write(90+myid,*) 'w_d tc'
  !do a_run = 1, CC_mem_aa_G(CC_mpi_gid+1)
  !  do bbb = 1, CC_n_vir
  !    do ij_run = 1, CC_mem_ij_D(CC_mpi_did+1)
  !      aaa = CC_index_aa_G(CC_mpi_gid+1) - 1 + a_run
  !      Call CC_cl_decode(ij_run,iii,jjj,CC_n_occ,2)
  !      write(90+myid,*) iii,jjj,aaa,bbb,CC_w_d(ij_run,aaa,bbb,1)
  !    end do
  !  end do
  !end do
 

!  write(70+myid,*) 'w_s coeff'
!  do iii = 1, CC_n_occ
!    do aaa = 1, CC_n_vir
!      write(70+myid,*) iii,aaa,CC_w_s(iii,aaa)
!    end do
!  end do
!
!  write(70+myid,*) 'w_d coeff add all'
!  do iii = 1, CC_n_occ
!    do jjj = iii, CC_n_occ
!      do ab_run = 1, n_ab
!        ab = ab_start - 1 + ab_run
!        Call CC_cl_decode(ab,aaa,bbb,CC_n_vir,2)
!        if (aaa.eq.bbb) then
!          write(70+myid,*) iii,jjj,aaa,aaa,CC_w_d(iii,jjj,ab_run,1)
!        end if
!      end do
!    end do
!  end do
!
!  write(70+myid,*) 'w_nd coeff add all'
!  do iii = 1, CC_n_occ
!    do jjj = 1, CC_n_occ
!      do ab_run = 1, n_ab
!        ab = ab_start - 1 + ab_run
!        Call CC_cl_decode(ab,aaa,bbb,CC_n_vir,2)
!        if (aaa.ne.bbb) then
!          write(70+myid,*) iii,jjj,aaa,bbb,CC_w_d(iii,jjj,ab_run,1)
!        end if
!      end do
!    end do
!  end do

  End Subroutine CC_cl_calc_w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_calc_w_s()

! (iii,aaa)

  Use CC_cl

  Implicit None

  Integer (kind = 8) :: s_tmp

  Integer :: n_i,i_start,i_end,n_a,a_start,a_end,n_c,c_start,c_end
  Integer :: i_tmp,j_tmp,l_tmp,i_run,j_run,k_run,a_run,b_run,c_run
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,errnum
  Double precision :: alpha,beta

  Double precision , dimension(:,:,:,:) , allocatable :: matA,matB,matD
  Double precision , dimension(:,:) , allocatable :: matC

!  print*,'w_s'

  n_i = CC_mem_ii_G(CC_mpi_gid+1)
  i_start = CC_index_ii_G(CC_mpi_gid+1)
  i_end = i_start - 1 + n_i

  n_c = CC_mem_aa_D(CC_mpi_did+1)
  c_start = CC_index_aa_D(CC_mpi_did+1)
  c_end = c_start - 1 + n_c

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  alpha = 1.0D0
  beta = 0.0D0

  ! First term h_ca(c,a) # t(i,c)
  Allocate(matA(CC_n_occ,n_c,1,1),stat=errnum)
  Call check_allocation(errnum,'matA in CC')

  Allocate(matB(n_c,n_a,1,1),stat=errnum)
  Call check_allocation(errnum,'matB in CC')

  Allocate(matC(CC_n_occ,n_a),stat=errnum)
  Call check_allocation(errnum,'matC in CC')

  matA(:,:,1,1) = CC_t_s(:,c_start:c_end)

  matB(:,:,1,1) = CC_h_ca(c_start:c_end,a_start:a_end)

  Call Dgemm('N','N',CC_n_occ,n_a,n_c,alpha,MatA,CC_n_occ, &
             MatB,n_c,beta,MatC,CC_n_occ)

  CC_w_s(:,a_start:a_end) = CC_w_s(:,a_start:a_end) + matC

  Deallocate(matA,matB,matC)

  if (n_i.gt.0) then
    ! Second term h_ik(i,k) # t(k,a)
    Allocate(matA(CC_n_occ,n_i,1,1),stat=errnum)
    Call check_allocation(errnum,'matA in CC')

    Allocate(matB(n_i,n_c,1,1),stat=errnum)
    Call check_allocation(errnum,'matB in CC')

    Allocate(matC(CC_n_occ,n_c),stat=errnum)
    Call check_allocation(errnum,'matC in CC')

    matA(:,:,1,1) = CC_h_ik(:,i_start:i_end)
    matB(:,:,1,1) = CC_t_s(i_start:i_end,c_start:c_end)

    Call Dgemm('N','N',CC_n_occ,n_c,n_i,alpha,MatA,CC_n_occ, &
               MatB,n_i,beta,MatC,CC_n_occ)

    CC_w_s(:,c_start:c_end) = CC_w_s(:,c_start:c_end) - matC

    Deallocate(matA,matB,matC)
  end if

  ! Third term h_ck(k,c) # (2 * t(k,i,c,a) - t(i,k,c,a) + t(i,c) * t(k,a))
  ! matA t(k,c,i,a)
  Allocate(matA(CC_n_occ,n_c,CC_n_occ,n_a),stat=errnum)
  Call check_allocation(errnum,'matA in CC')

  Allocate(matC(CC_n_occ,n_a),stat=errnum)
  Call check_allocation(errnum,'matC in CC')

  !$OMP PARALLEL Default(Shared) & 
  !$OMP Private(a_run,c_run,kkk,ccc,iii,aaa)
  !$OMP DO
  do a_run = 1, n_a
    do iii = 1, CC_n_occ
      do c_run = 1, n_c
        do kkk = 1, CC_n_occ
          ccc = c_start - 1 + c_run
          aaa = a_start - 1 + a_run
          matA(kkk,c_run,iii,a_run) = 2.0D0 * CC_t_d(kkk,iii,c_run,aaa) &
                                            - CC_t_d(iii,kkk,c_run,aaa) &
                                            + CC_t_s(iii,ccc) * CC_t_s(kkk,aaa)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  i_tmp = CC_n_occ * n_c
  j_tmp = CC_n_occ * n_a

  Call Dgemm('T','N',j_tmp,1,i_tmp,alpha,MatA,i_tmp, &
             CC_h_ck(:,c_start:c_end),i_tmp,beta,MatC,j_tmp)

  CC_w_s(:,a_start:a_end) = CC_w_s(:,a_start:a_end) + matC

  Deallocate(matA,matC)

  ! Forth term t(k,c) # (2 * (kc|ia) - (ki|ac))
  ! matA (k,c,i,a)
  Allocate(matA(CC_n_occ,n_c,CC_n_occ,n_a),stat=errnum)
  Call check_allocation(errnum,'matA in CC')

  Allocate(matB(CC_n_occ,n_c,1,1),stat=errnum)
  Call check_allocation(errnum,'matB in CC')

  Allocate(matC(CC_n_occ,n_a),stat=errnum)
  Call check_allocation(errnum,'matC in CC')

  !$OMP PARALLEL Default(Shared) & 
  !$OMP Private(a_run,kkk,ccc,iii,aaa)
  !$OMP DO
  do a_run = 1, n_a
    do iii = 1, CC_n_occ
      do ccc = 1, n_c
        do kkk = 1, CC_n_occ
          aaa = a_start - 1 + a_run
          matA(kkk,ccc,iii,a_run) = 2.0D0 * CC_intl_iajb(kkk,iii,ccc,aaa) &
                                          - CC_intl_kiac(kkk,ccc,iii,aaa)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  matB(:,:,1,1) = CC_t_s(:,c_start:c_end)

  i_tmp = CC_n_occ * n_c
  j_tmp = CC_n_occ * n_a

  Call Dgemm('T','N',j_tmp,1,i_tmp,alpha,MatA,i_tmp, &
             matB,i_tmp,beta,MatC,j_tmp)

  CC_w_s(:,a_start:a_end) = CC_w_s(:,a_start:a_end) + matC

  Deallocate(matA,matB,matC)


  ! Fifth term (t(k,i,c,d) + t(k,c) * t(i,d)) # (2 * (kc|ad) - (kd|ac))
  ! equal to (ac|kd) # (2 * tau(k,i,d,c) - tau(k,i,c,d))
  ! matA (k,c,d,i)
  Allocate(matA(CC_n_occ,n_c,CC_n_vir,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'matA in CC')

  Allocate(matC(CC_n_occ,n_a),stat=errnum)
  Call check_allocation(errnum,'matC in CC')

  ! t(k,c,i,d) + t(k,c) + t(i,d)
  !$OMP PARALLEL Default(Shared) & 
  !$OMP Private(c_run,kkk,ccc,iii,ddd)
  !$OMP DO
  do iii = 1, CC_n_occ
    do ddd = 1, CC_n_vir
      do c_run = 1, n_c
        do kkk = 1, CC_n_occ
          ccc = c_start - 1 + c_run
          matA(kkk,c_run,ddd,iii) = 2.0 * (CC_t_d(iii,kkk,c_run,ddd) &
                                           + CC_t_s(kkk,ddd) * CC_t_s(iii,ccc)) &
                                        - CC_t_d(kkk,iii,c_run,ddd)  &
                                        - CC_t_s(kkk,ccc) * CC_t_s(iii,ddd)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  i_tmp = CC_n_occ * n_c * CC_n_vir

  alpha = 1.0D0
  beta = 0.0D0

  Call Dgemm('T','N',CC_n_occ,n_a,i_tmp,alpha,MatA,i_tmp,  &
             CC_intl_ackd(:,:,:,a_start:a_end),i_tmp,beta,matC,CC_n_occ)

  CC_w_s(:,a_start:a_end) = CC_w_s(:,a_start:a_end) + matC

  Deallocate(matA,matC)

  ! Sixth term (t(k,c,l,a) + t(k,c) * t(l,a)) # (2 * (li|kc) - (ki|lc))
  ! matA (k,l,c,i)
  Allocate(matA(CC_n_occ,CC_n_occ,n_c,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'matA in CC')

  ! matB (k,l,c,a)
  Allocate(matB(CC_n_occ,CC_n_occ,n_c,n_a),stat=errnum)
  Call check_allocation(errnum,'matB in CC')

  Allocate(matC(CC_n_occ,n_a),stat=errnum)
  Call check_allocation(errnum,'matC in CC')

  ! 2 * (li|kc) - (ki|lc)
  !$OMP PARALLEL Default(Shared) & 
  !$OMP Private(kkk,iii,lll,ccc)
  !$OMP DO
  do iii = 1, CC_n_occ
    do ccc = 1, n_c
      do lll = 1, CC_n_occ
        do kkk = 1, CC_n_occ
          matA(kkk,lll,ccc,iii) = 2.0D0 * CC_intl_kilc(lll,ccc,iii,kkk) &
                                        - CC_intl_kilc(kkk,ccc,iii,lll)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! t(k,c,l,a) + t(k,c) * t(l,a) 
  !$OMP PARALLEL Default(Shared) & 
  !$OMP Private(c_run,a_run,kkk,ccc,aaa,lll)
  !$OMP DO
  do a_run = 1, n_a
    do c_run = 1, n_c
      do lll = 1, CC_n_occ
        do kkk = 1, CC_n_occ
          aaa = a_start - 1 + a_run
          ccc = c_start - 1 + c_run
          matB(kkk,lll,c_run,a_run) = CC_t_d(kkk,lll,c_run,aaa) &
                                    + CC_t_s(kkk,ccc) * CC_t_s(lll,aaa)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  i_tmp = n_c * CC_n_occ * CC_n_occ

  Call Dgemm('T','N',CC_n_occ,n_a,i_tmp,alpha,MatA,i_tmp, &
             matB,i_tmp,beta,MatC,CC_n_occ)

  CC_w_s(:,a_start:a_end) = CC_w_s(:,a_start:a_end) - matC

  Deallocate(matA,matB,matC)

  s_tmp = Int(CC_n_occ * CC_n_vir,8)

  Call CC_mpi_allreduce(s_tmp,CC_w_s,MPI_COMM_WORLD) 

  End Subroutine CC_cl_calc_w_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_calc_w_d()

  Use CC_cl

  Implicit None

  Integer :: n_c,c_start,c_end,n_a,a_start,a_end
  Integer :: i_tmp,j_tmp,l_tmp,i_run,j_run,k_run,a_run,b_run,c_run
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,errnum
  Double precision :: alpha,beta

  Double precision , dimension(:,:,:,:) , allocatable :: matA,matB,matC

!  print*,'w_d'

  n_c = CC_mem_aa_D(CC_mpi_did+1)
  c_start = CC_index_aa_D(CC_mpi_did+1)
  c_end = c_start - 1 + n_c

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  !$OMP PARALLEL Default(Shared) & 
  !$OMP Private(iii,jjj,aaa,bbb,a_run)
  !$OMP DO
  do bbb = 1, n_c
    do a_run = 1, n_a
      do iii = 1, CC_n_occ
        do jjj = 1, CC_n_occ
          aaa = a_start - 1 + a_run
          CC_ten_con(iii,jjj,a_run,bbb) = 0.5D0 * CC_intl_iajb(jjj,iii,bbb,aaa)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  alpha = 1.0D0
  beta = 0.0D0

  ! First term g_ca(c,a) # t(j,i,b,c)
  ! CC_t_d (j,i,b,c)
  ! matC(j,i,b,a)
  Allocate(matC(CC_n_occ,CC_n_occ,n_c,n_a),stat=errnum)
  Call check_allocation(errnum,'matC in CC')

  i_tmp = CC_n_occ * CC_n_occ * n_c

  Call Dgemm('N','N',i_tmp,n_a,CC_n_vir,alpha,CC_t_d,i_tmp, &
             CC_g_ca(:,a_start:a_end),CC_n_vir,beta,MatC,i_tmp)

  !$OMP PARALLEL Default(Shared) & 
  !$OMP Private(iii,jjj,aaa,bbb)
  !$OMP DO
  do bbb = 1, n_c
    do aaa = 1, n_a
      do iii = 1, CC_n_occ
        do jjj = 1, CC_n_occ
          CC_ten_con(iii,jjj,aaa,bbb) = CC_ten_con(iii,jjj,aaa,bbb) &
                                      + matC(jjj,iii,bbb,aaa)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  Deallocate(matC)

  ! Second term g_ik(j,k) # t(k,i,b,a)
  ! matC(j,i,b,a)
  Allocate(matC(CC_n_occ,CC_n_occ,n_c,n_a),stat=errnum)
  Call check_allocation(errnum,'matC in CC')

  i_tmp = CC_n_occ * n_a * n_c

  Call Dgemm('N','N',CC_n_occ,i_tmp,CC_n_occ,alpha,CC_g_ik,CC_n_occ, &
             CC_t_d(:,:,:,a_start:a_end),CC_n_occ,beta,MatC,CC_n_occ)

  !$OMP PARALLEL Default(Shared) & 
  !$OMP Private(iii,jjj,aaa,bbb)
  !$OMP DO
  do bbb = 1, n_c
    do aaa = 1, n_a
      do jjj = 1, CC_n_occ
        do iii = 1, CC_n_occ
          CC_ten_con(iii,jjj,aaa,bbb) = CC_ten_con(iii,jjj,aaa,bbb) &
                                      - matC(jjj,iii,bbb,aaa)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  Deallocate(matC)

  ! Third term t(j,c) # (cb|ia) - t(j,c) # (cb|ik) # t(k,a)
  ! Calculate t(k,a) # (cb|ik)
  ! intl(k,b,i,c)
  ! matA (b,i,c,a)
  Allocate(matA(n_c,CC_n_occ,CC_n_vir,1),stat=errnum)
  Call check_allocation(errnum,'matA in CC')

  ! matB(i,b,c,a)
  Allocate(matB(CC_n_occ,n_c,CC_n_vir,1),stat=errnum)
  Call check_allocation(errnum,'matB in CC')

  ! matC(i,b,j,a)
  Allocate(matC(CC_n_occ,n_c,CC_n_occ,1),stat=errnum)
  Call check_allocation(errnum,'matC in CC')

  aaa = a_start - 1
  do a_run = 1, n_a
    aaa = aaa + 1

    i_tmp = CC_n_occ * n_c * CC_n_vir

    Call Dgemm('T','N',i_tmp,1,CC_n_occ,alpha,CC_intl_kiac,CC_n_occ, &
               CC_t_s(:,aaa),CC_n_occ,beta,MatA,i_tmp)

    ! Calculate (cb|ia) - t(k,a) # (cb|ik)
    ! intl_ackd(b,a,i,c)
    !$OMP PARALLEL Default(Shared) & 
    !$OMP Private(bbb,iii,ccc)
    !$OMP DO
    do bbb = 1, n_c
      do iii = 1, CC_n_occ
        do ccc = 1, CC_n_vir
          matB(iii,bbb,ccc,1) = CC_intl_ackd(iii,bbb,aaa,ccc) - matA(bbb,iii,ccc,1)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    i_tmp = CC_n_occ * n_c
    Call Dgemm('N','T',i_tmp,CC_n_occ,CC_n_vir,alpha,matB,i_tmp, &
               CC_t_s,CC_n_occ,beta,MatC,i_tmp)

    !$OMP PARALLEL Default(Shared) & 
    !$OMP Private(bbb,iii,jjj)
    !$OMP DO
    do bbb = 1, n_c
      do iii = 1, CC_n_occ
        do jjj = 1, CC_n_occ
          CC_ten_con(iii,jjj,a_run,bbb) = CC_ten_con(iii,jjj,a_run,bbb) &
                                        + matC(iii,bbb,jjj,1)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end do

  Deallocate(matA,matB,matC)

  ! Forth term - t(k,a) # (ik|jb) - t(i,c) # (jb|kc)) # t(k,a)
  ! Calculate t(i,c) # (jb|kc)
  ! intl(j,k,b,c)
  ! matA (j,k,b,i)
  Allocate(matA(CC_n_occ,CC_n_occ,n_c,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'matA in CC')

  i_tmp = CC_n_occ * n_c * CC_n_occ

  Call Dgemm('N','T',i_tmp,CC_n_occ,CC_n_vir,alpha,CC_intl_iajb,i_tmp, &
             CC_t_s,CC_n_occ,beta,MatA,i_tmp)

  ! matB(i,j,b,k)
  Allocate(matB(CC_n_occ,CC_n_occ,n_c,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'matB in CC')

  ! Calculate (ik|jb) + t(i,c) # (jb|kc)
  ! intl_likc(k,b,i,j)
  !$OMP PARALLEL Default(Shared) & 
  !$OMP Private(kkk,bbb,iii,jjj)
  !$OMP DO
  do kkk = 1, CC_n_occ
    do bbb = 1, n_c
      do jjj = 1, CC_n_occ
        do iii = 1, CC_n_occ
          matB(iii,jjj,bbb,kkk) = CC_intl_kilc(kkk,bbb,iii,jjj) &
                                + matA(jjj,kkk,bbb,iii)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  Deallocate(matA)

  ! matC(i,j,b,a)
  Allocate(matC(CC_n_occ,CC_n_occ,n_c,n_a),stat=errnum)
  Call check_allocation(errnum,'matC in CC')

  i_tmp = CC_n_occ * CC_n_occ * n_c
  Call Dgemm('N','N',i_tmp,n_a,CC_n_occ,alpha,matB,i_tmp, &
             CC_t_s(:,a_start:a_end),CC_n_occ,beta,MatC,i_tmp)

  !$OMP PARALLEL Default(Shared) & 
  !$OMP Private(aaa,bbb,iii,jjj)
  !$OMP DO
  do bbb = 1, n_c
    do aaa = 1, n_a
      do jjj = 1, CC_n_occ
        do iii = 1, CC_n_occ
          CC_ten_con(iii,jjj,aaa,bbb) = CC_ten_con(iii,jjj,aaa,bbb) &
                                      - matC(iii,jjj,bbb,aaa)
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  Deallocate(matB,matC)

  End Subroutine CC_cl_calc_w_d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_add_j2w()

  Use CC_cl

  Implicit None

  Integer (kind = 8) :: s_tmp

  Integer :: n_b,b_start,b_end,n_a,a_start,a_end,n_c,c_start,c_end
  Integer :: i_tmp,j_tmp,l_tmp
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,b_run
  Integer :: errnum
  Double precision :: alpha,beta

  Integer :: i_task,i_task_run
  Integer :: req1,req2,i_send,i_recv
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Double precision , dimension(:,:,:,:) , allocatable :: matA,matB,matC
  Double precision , dimension(:,:,:,:) , allocatable :: t_tmp,t_send,t_recv

  n_c = CC_mem_aa_D(CC_mpi_did+1)
  c_start = CC_index_aa_D(CC_mpi_did+1)
  c_end = c_start - 1 + n_c

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  alpha = 1.0D0
  beta = 0.0D0

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

    n_b = CC_mem_aa_D(i_task)
    b_start = CC_index_aa_D(i_task)
    b_end = b_start - 1 + n_b

    ! Get T
    ! matB(k,c,j,b)
    Allocate(matB(CC_n_occ,n_c,CC_n_occ,n_b),stat=errnum)
    Call check_allocation(errnum,'matB in CC')

    ! t(k,j,c,b) - 0.5 * t(j,k,c,b)
    !$OMP PARALLEL Default(Shared) & 
    !$OMP Private(b_run,bbb,jjj,ccc,kkk)
    !$OMP DO
    do b_run = 1, n_b
      do jjj = 1, CC_n_occ
        do ccc = 1, n_c
          do kkk = 1, CC_n_occ
            bbb = b_start - 1 + b_run
            matB(kkk,ccc,jjj,b_run) = CC_t_d(kkk,jjj,ccc,bbb) &
                          - 0.5D0 * CC_t_d(jjj,kkk,ccc,bbb)
          end do
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! matC(i,a,j,b)
    Allocate(matC(CC_n_occ,n_a,CC_n_occ,n_b),stat=errnum)
    Call check_allocation(errnum,'matC in CC')

    i_tmp = CC_n_occ * n_c
    j_tmp = CC_n_occ * n_a
    l_tmp = CC_n_occ * n_b

    Call Dgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,CC_j_aux,i_tmp, &
               MatB,i_tmp,beta,MatC,j_tmp)

    if (i_task_run.ne.1) then

      Call MPI_WAIT(req1,stat1,errnum)
      Call MPI_WAIT(req2,stat2,errnum)

      !$OMP PARALLEL Default(shared) &
      !$OMP Private(iii,jjj,aaa,bbb)
      !$OMP DO
      do bbb = 1, n_c
        do aaa = 1, n_a
          do jjj = 1, CC_n_occ
            do iii = 1, CC_n_occ
              CC_ten_con(iii,jjj,aaa,bbb) = CC_ten_con(iii,jjj,aaa,bbb) &
                                          + t_recv(iii,aaa,jjj,bbb)
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      Deallocate(t_send,t_recv)

    end if

    if (i_task_run.ne.CC_mpi_domain_size) then

      Allocate(t_send(CC_n_occ,n_a,CC_n_occ,n_b),stat=errnum)
      Call check_allocation(errnum,'t_send in CC')

      t_send = matC

      Allocate(t_recv(CC_n_occ,n_a,CC_n_occ,n_c),stat=errnum)
      Call check_allocation(errnum,'t_recv in CC')

      i_tmp = CC_n_occ * CC_n_occ * n_a
      s_tmp = Int(i_tmp,8) * Int(n_b,8)

      Call CC_mpi_real_isend(s_tmp,t_send,i_send-1,100,req1,CC_mpi_comm_domain)

      s_tmp = Int(i_tmp,8) * Int(n_c,8)
      Call CC_mpi_real_irecv(s_tmp,t_recv,i_recv-1,100,req2,CC_mpi_comm_domain)

    else

      !$OMP PARALLEL Default(shared) &
      !$OMP Private(iii,jjj,aaa,bbb)
      !$OMP DO
      do bbb = 1, n_c
        do aaa = 1, n_a
          do jjj = 1, CC_n_occ
            do iii = 1, CC_n_occ
              CC_ten_con(iii,jjj,aaa,bbb) = CC_ten_con(iii,jjj,aaa,bbb) &
                                          + matC(iii,aaa,jjj,bbb)
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end if

    Deallocate(matB,matC)

  end do

  Deallocate(CC_j_aux)

  End Subroutine CC_cl_add_j2w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_add_k2w()

  Use CC_cl

  Implicit None

  Integer (kind = 8) :: s_tmp

  Integer :: n_b,b_start,b_end,n_a,a_start,a_end,n_c,c_start,c_end
  Integer :: i_tmp,j_tmp,l_tmp
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,b_run
  Integer :: errnum
  Double precision :: alpha,beta

  Integer :: i_task,i_task_run
  Integer :: req1,req2,i_send,i_recv
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Double precision , dimension(:,:,:,:) , allocatable :: matA,matB,matC
  Double precision , dimension(:,:,:,:) , allocatable :: t_tmp,t_send,t_recv

  n_c = CC_mem_aa_D(CC_mpi_did+1)
  c_start = CC_index_aa_D(CC_mpi_did+1)
  c_end = c_start - 1 + n_c

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  alpha = 1.0D0
  beta = 0.0D0

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

    n_b = CC_mem_aa_D(i_task)
    b_start = CC_index_aa_D(i_task)
    b_end = b_start - 1 + n_b

    ! Get T
    ! matB(k,c,j,b)
    Allocate(matB(CC_n_occ,n_c,CC_n_occ,n_b),stat=errnum)
    Call check_allocation(errnum,'matB in CC')

    ! t(k,j,c,b) - 0.5 * t(j,k,c,b)
    !$OMP PARALLEL Default(Shared) & 
    !$OMP Private(b_run,bbb,jjj,ccc,kkk)
    !$OMP DO
    do b_run = 1, n_b
      do jjj = 1, CC_n_occ
        do ccc = 1, n_c
          do kkk = 1, CC_n_occ
            bbb = b_start - 1 + b_run
            matB(kkk,ccc,jjj,b_run) = CC_t_d(jjj,kkk,ccc,bbb)
          end do
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! matC(i,a,j,b)
    Allocate(matC(CC_n_occ,n_a,CC_n_occ,n_b),stat=errnum)
    Call check_allocation(errnum,'matC in CC')

    i_tmp = CC_n_occ * n_c
    j_tmp = CC_n_occ * n_a
    l_tmp = CC_n_occ * n_b

    Call Dgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,CC_k_aux,i_tmp, &
               MatB,i_tmp,beta,MatC,j_tmp)

    if (i_task_run.ne.1) then

      Call MPI_WAIT(req1,stat1,errnum)
      Call MPI_WAIT(req2,stat2,errnum)

      !$OMP PARALLEL Default(shared) &
      !$OMP Private(iii,jjj,aaa,bbb)
      !$OMP DO
      do bbb = 1, n_c
        do aaa = 1, n_a
          do jjj = 1, CC_n_occ
            do iii = 1, CC_n_occ
              CC_ten_con(iii,jjj,aaa,bbb) = CC_ten_con(iii,jjj,aaa,bbb) &
                                          - 0.5D0 * t_recv(iii,aaa,jjj,bbb) &
                                          - t_recv(jjj,aaa,iii,bbb)
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      Deallocate(t_send,t_recv)

    end if

    if (i_task_run.ne.CC_mpi_domain_size) then

      Allocate(t_send(CC_n_occ,n_a,CC_n_occ,n_b),stat=errnum)
      Call check_allocation(errnum,'t_send in CC')

      t_send = matC

      Allocate(t_recv(CC_n_occ,n_a,CC_n_occ,n_c),stat=errnum)
      Call check_allocation(errnum,'t_recv in CC')

      i_tmp = CC_n_occ * CC_n_occ * n_a
      s_tmp = Int(i_tmp,8) * Int(n_b,8)

      Call CC_mpi_real_isend(s_tmp,t_send,i_send-1,100,req1,CC_mpi_comm_domain)

      s_tmp = Int(i_tmp,8) * Int(n_c,8)
      Call CC_mpi_real_irecv(s_tmp,t_recv,i_recv-1,100,req2,CC_mpi_comm_domain)

    else

      !$OMP PARALLEL Default(shared) &
      !$OMP Private(iii,jjj,aaa,bbb)
      !$OMP DO
      do bbb = 1, n_c
        do aaa = 1, n_a
          do jjj = 1, CC_n_occ
            do iii = 1, CC_n_occ
              CC_ten_con(iii,jjj,aaa,bbb) = CC_ten_con(iii,jjj,aaa,bbb) &
                                          - 0.5D0 * matC(iii,aaa,jjj,bbb) &
                                          - matC(jjj,aaa,iii,bbb)
            end do
          end do
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL

    end if

    Deallocate(matB,matC)

  end do

  Deallocate(CC_k_aux)

  End Subroutine CC_cl_add_k2w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_add_a2w()

  Use CC_cl

  Implicit None

  Integer (kind = 8) :: s_tmp

  Integer :: n_a,a_start,a_end,n_b,b_start,b_end,n_kl,kl_start,kl_end,code_kl,kl_run
  Integer :: i_tmp,j_tmp,l_tmp,i_run,j_run,k_run,a_run,b_run,c_run
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd
  Integer :: errnum
  Double precision :: alpha,beta

  Integer :: i_task,i_task_run
  Integer :: req1,req2,i_send,i_recv
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Double precision , dimension(:,:,:,:) , allocatable :: matA,matB,matC
  Double precision , dimension(:,:,:,:) , allocatable :: t_tmp,a_recv,a_tmp

  n_b = CC_mem_aa_D(CC_mpi_did+1)
  b_start = CC_index_aa_D(CC_mpi_did+1)
  b_end = b_start - 1 + n_b

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  alpha = 0.5D0
  beta = 1.0D0

  !print*,'a2w'

  Allocate(a_tmp(CC_mem_kl_D(CC_mpi_did+1),CC_n_occ,CC_n_occ,1),stat=errnum)
  Call check_allocation(errnum,'a_tmp in CC')

  a_tmp = CC_a_aux

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

    n_kl = CC_mem_kl_D(i_task)
    kl_start = CC_index_kl_D(i_task)
    kl_end = kl_start - 1 + n_kl

    if (i_recv.ne.CC_mpi_did+1) then

      Allocate(a_recv(CC_mem_kl_D(i_recv),CC_n_occ,CC_n_occ,1),stat=errnum)
      Call check_allocation(errnum,'t_recv in CC')

      i_tmp = CC_n_occ * CC_n_occ
      s_tmp = Int(i_tmp,8) * Int(CC_mem_kl_D(CC_mpi_did+1),8)

      Call CC_mpi_real_isend(s_tmp,CC_a_aux,i_send-1,100,req1,CC_mpi_comm_domain)

      s_tmp = Int(i_tmp,8) * Int(CC_mem_kl_D(i_recv),8)
      Call CC_mpi_real_irecv(s_tmp,a_recv,i_recv-1,100,req2,CC_mpi_comm_domain)

    end if

    Allocate(t_tmp(n_kl,n_a,n_b,1),stat=errnum)
    Call check_allocation(errnum,'t_tmp in CC')

    !$OMP PARALLEL Default(shared) &
    !$OMP Private(b_run,a_run,kl_run,aaa,bbb,code_kl,kkk,lll)
    !$OMP DO
    do b_run = 1, n_b
      do a_run = 1, n_a
        do kl_run = 1, n_kl

          aaa = a_start - 1 + a_run
          bbb = b_start - 1 + b_run
          code_kl = kl_start - 1 + kl_run
          Call CC_cl_decode(code_kl,kkk,lll,CC_n_occ,1)

          t_tmp(kl_run,a_run,b_run,1) = CC_t_d(lll,kkk,b_run,aaa) &
                                      + CC_t_s(lll,bbb) * CC_t_s(kkk,aaa)

        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    i_tmp = n_kl
    j_tmp = CC_n_occ * CC_n_occ
    l_tmp = n_a * n_b

    Call Dgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,a_tmp,i_tmp, &
               t_tmp,i_tmp,beta,CC_ten_con,j_tmp)


    Deallocate(t_tmp,a_tmp)

    if (i_recv.ne.CC_mpi_did+1) then

      Call MPI_WAIT(req1,stat1,errnum)
      Call MPI_WAIT(req2,stat2,errnum)

      Allocate(a_tmp(CC_mem_kl_D(i_recv),CC_n_occ,CC_n_occ,1),stat=errnum)
      Call check_allocation(errnum,'a_tmp in CC')

      a_tmp = a_recv

      Deallocate(a_recv)

    end if

  end do

  End Subroutine CC_cl_add_a2w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_add_b2w()

  Use timing
  Use CC_cl

  Implicit None

  Integer (kind = 8) :: s_tmp

  Integer :: n_a,a_start,a_end,a_run,n_ij,ij_run,n_kl,ij_start,ij_end
  Integer :: n_cd,n_d,d_start,d_end,d_run,n_b,b_start,b_end,b_run
  Integer :: a_calc,a_finish,i_finish
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,i_tmp,j_tmp,l_tmp
  Integer :: errnum
  Integer :: nth,OMP_GET_NUM_THREADS
  Integer :: n_calc,n_mx,n_iter,i_iter

  Integer :: i_task,i_task_run
  Integer :: req1,req2,i_send,i_recv
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Double precision :: ctime, wtime, c_pre, w_pre, ctime_aux, wtime_aux
  Double precision :: ctime_calc, wtime_calc, ctime_comm, wtime_comm
  Double precision :: ctime_calcb, wtime_calcb, ctime_commb, wtime_commb

  Double precision :: alpha,beta

  Double precision , dimension(:,:,:) , allocatable :: t_tmp
  Double precision , dimension(:,:,:) , allocatable :: RI_tmp,RI_B
  Double precision , dimension(:,:,:,:) , allocatable :: CC_b_aux
  Double precision , dimension(:,:,:) , allocatable :: w_tmp,w_send,w_recv
 
  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  n_ij = CC_n_occ * (CC_n_occ + 1) / 2

  n_d = CC_mem_aa_D(CC_mpi_did+1)
  d_start = CC_index_aa_D(CC_mpi_did+1)
  d_end = d_start - 1 + n_d

  n_cd = n_d * CC_n_vir

  n_mx = Int(CC_tmp_mem/(CC_mem_aa_D(1) * CC_n_vir * CC_n_vir * 2))

  ctime_calc = 0.0D0
  wtime_calc = 0.0D0
  ctime_comm = 0.0D0
  wtime_comm = 0.0D0
  ctime_aux = 0.0D0
  wtime_aux = 0.0D0
  ctime_calcb = 0.0D0
  wtime_calcb = 0.0D0
  ctime_commb = 0.0D0
  wtime_commb = 0.0D0

  !print*,'b2w',n_ab
  Allocate(RI_tmp(CC_n_bas,n_d,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'RI_tmp in CC')

  ! RI_B(bas,ddd,bbb)
  Allocate(RI_B(CC_n_bas,n_d,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'RI_B in CC')

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(bbb,ddd)
  !$OMP DO
  do bbb = 1, CC_n_vir
    do ddd = 1, n_d
      RI_B(:,ddd,bbb) = CC_RI_B(:,ddd,CC_n_occ+bbb)
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(ddd,kkk)
  !$OMP DO
  do kkk = 1, CC_n_occ
    do ddd = 1, n_d
      RI_tmp(:,ddd,kkk) = CC_RI_B(:,ddd,kkk)
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  alpha = -1.0D0
  beta = 1.0D0

  i_tmp = CC_n_bas * n_d

  Call Dgemm('N','N',i_tmp,CC_n_vir,CC_n_occ,alpha,RI_tmp,i_tmp, &
             CC_t_s,CC_n_occ,beta,RI_B,i_tmp)

  Deallocate(RI_tmp)

  Allocate(t_tmp(CC_n_vir,n_d,n_ij),stat=errnum)
  Call check_allocation(errnum,'t_tmp in CC')

  !$OMP PARALLEL Default(shared) &
  !$OMP Private(jjj,iii,ddd,ccc,d_run,ij_run)
  !$OMP DO
  do ij_run = 1, n_ij
    do ccc = 1, CC_n_vir
      do d_run = 1, n_d
        ddd = d_start - 1 + d_run
        Call CC_cl_decode(ij_run,iii,jjj,CC_n_occ,2)
        t_tmp(ccc,d_run,ij_run) = CC_t_d(jjj,iii,d_run,ccc) &
                                + CC_t_s(jjj,ddd) * CC_t_s(iii,ccc)
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  if (n_a.ne.0) then

    if (Mod(n_a,n_mx).eq.0) then
      n_iter = n_a / n_mx
    else
      n_iter = Int(n_a / n_mx) + 1
    end if

  else

    n_iter = 0

  end if

  alpha = 1.0D0
  beta = 0.0D0

  a_finish = a_start - 1
  i_finish = 0

  do i_iter = 1, n_iter

    if (i_iter.ne.n_iter) then
      a_calc = n_mx
    else
      a_calc = n_a - i_finish
    end if

    ! b_aux(c,d,a,b)
    Allocate(CC_b_aux(CC_n_vir,n_d,a_calc,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'CC_b_aux in CC')

    Call get_timestamps(c_pre, w_pre)
    Call CC_cl_aux_b(a_finish,a_calc,CC_b_aux,RI_B,ctime_calcb,wtime_calcb,&
                                                   ctime_commb,wtime_commb)
    Call get_timestamps(ctime, wtime)
    ctime_aux = ctime_aux + ctime - c_pre
    wtime_aux = wtime_aux + wtime - w_pre

    !write(90+myid,*) 'b_aux'
    !do ccc = 1, CC_n_vir
    !  do ddd = 1, CC_n_vir
    !    do aaa = 1, a_calc
    !      do bbb = 1, CC_n_vir
    !        write(90+myid,*) ccc,ddd,aaa,bbb,CC_b_aux(ccc,ddd,aaa,bbb)
    !      end do
    !    end do
    !  end do
    !end do

    i_task = CC_mpi_did + 1
    i_send = CC_mpi_did + 1
    i_recv = CC_mpi_did + 1

    n_kl = CC_mem_ij_D(CC_mpi_did+1)

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

      n_ij = CC_mem_ij_D(i_task)
      ij_start = CC_index_ij_D(i_task)
      ij_end = ij_start - 1 + n_ij

      Allocate(w_tmp(n_ij,a_calc,CC_n_vir),stat=errnum)
      Call check_allocation(errnum,'w_tmp in CC')

      i_tmp = n_d * CC_n_vir
      j_tmp = n_ij
      l_tmp = a_calc * CC_n_vir

      Call get_timestamps(c_pre, w_pre)

      Call Dgemm('T','N',j_tmp,l_tmp,i_tmp,alpha,t_tmp(:,:,ij_start:ij_end),i_tmp, &
                 CC_b_aux,i_tmp,beta,w_tmp,j_tmp)

      Call get_timestamps(ctime, wtime)

      ctime_calc = ctime_calc + ctime - c_pre
      wtime_calc = wtime_calc + wtime - w_pre

      if (i_task_run.ne.1) then

        Call get_timestamps(c_pre, w_pre)

        Call MPI_WAIT(req1,stat1,errnum)
        Call MPI_WAIT(req2,stat2,errnum)

        Call get_timestamps(ctime, wtime)

        ctime_comm = ctime_comm + ctime - c_pre
        wtime_comm = wtime_comm + wtime - w_pre

        !$OMP PARALLEL Default(shared) &
        !$OMP private(ij_run,a_run,aaa,bbb)
        !$OMP DO
        do a_run = 1, a_calc
          do bbb = 1, CC_n_vir
            do ij_run = 1, n_kl
              aaa = i_finish + a_run
              CC_w_d(ij_run,aaa,bbb,1) = CC_w_d(ij_run,aaa,bbb,1) + w_recv(ij_run,a_run,bbb)
            end do
          end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        Deallocate(w_send,w_recv)
      end if

      if (i_task_run.ne.CC_mpi_domain_size) then

        Allocate(w_send(n_ij,a_calc,CC_n_vir),stat=errnum)
        Call check_allocation(errnum,'w_send in CC')

        w_send = w_tmp

        Allocate(w_recv(n_kl,a_calc,CC_n_vir),stat=errnum)
        Call check_allocation(errnum,'w_recv in CC')

        s_tmp = Int(n_ij,8) * Int(a_calc * CC_n_vir,8)

        Call CC_mpi_real_isend(s_tmp,w_send,i_send-1,100,req1,CC_mpi_comm_domain)

        s_tmp = Int(n_kl,8) * Int(a_calc * CC_n_vir,8)
        Call CC_mpi_real_irecv(s_tmp,w_recv,i_recv-1,100,req2,CC_mpi_comm_domain)

      else

        !$OMP PARALLEL Default(shared) &
        !$OMP private(ij_run,a_run,aaa,bbb)
        !$OMP DO
        do a_run = 1, a_calc
          do bbb = 1, CC_n_vir
            do ij_run = 1, n_kl
              aaa = i_finish + a_run
              CC_w_d(ij_run,aaa,bbb,1) = CC_w_d(ij_run,aaa,bbb,1) + w_tmp(ij_run,a_run,bbb)
            end do
          end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL

      end if

      Deallocate(w_tmp)

    end do

    Deallocate(CC_b_aux)

    a_finish = a_finish + a_calc
    i_finish = i_finish + a_calc

  end do

  Deallocate(RI_B,t_tmp)

  !write(90+myid,*) 'w_d'
  !do a_run = 1, CC_mem_aa_G(CC_mpi_gid+1)
  !  do bbb = 1, CC_n_vir
  !    do ij_run = 1, n_kl
  !      aaa = CC_index_aa_G(CC_mpi_gid+1) - 1 + a_run
  !      Call CC_cl_decode(ij_run,iii,jjj,CC_n_occ,2)
  !      write(90+myid,*) iii,jjj,aaa,bbb,CC_w_d(ij_run,aaa,bbb,1)
  !    end do
  !  end do
  !end do
 

  if (myid.eq.0) then
    write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of calc b', &
                                               ctime_calcb, 's (cpu)', &
                                               wtime_calcb, 's (wall clock)'
  end if

  if (myid.eq.0) then
    write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of comm b', &
                                               ctime_commb, 's (cpu)', &
                                               wtime_commb, 's (wall clock)'
  end if

  if (myid.eq.0) then
    write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of evaluating b', &
                                             ctime_aux, 's (cpu)', &
                                             wtime_aux, 's (wall clock)'
  end if

  if (myid.eq.0) then
    write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of adding b', &
                                             ctime_calc, 's (cpu)', &
                                             wtime_calc, 's (wall clock)'
  end if

  if (myid.eq.0) then
    write(use_unit,"(2x,A37,F14.3,A,F14.3,A)") 'Time of comm b', &
                                               ctime_comm, 's (cpu)', &
                                               wtime_comm, 's (wall clock)'
  end if

  End Subroutine CC_cl_add_b2w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_add_tc2w()

  Use CC_cl

  Implicit None

  Integer :: n_ij,ij_start,ij_end,n_a,a_start,a_end,ij_run,code_ij
  Integer :: n_b,b_start,b_end,i_run,j_run,k_run,a_run,b_run,i_tmp
  Integer :: n_ab,ab_start,ab_end,code_ab,i_ab,ab_run
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd
  Integer (kind = 8) :: s_tmp
  Integer :: source_id,errnum

  Integer :: i_task_run,i_task
  Integer :: req1,req2,i_send,i_recv
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2

  Double precision :: alpha,beta
  Double precision , dimension(:,:,:,:) , allocatable :: con_tmp,con_recv
  Double precision , dimension(:,:,:) , allocatable :: w_d_all,w_recv,w_send

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  n_ij = CC_mem_ij_D(CC_mpi_did+1)
  ij_start = CC_index_ij_D(CC_mpi_did+1)
  ij_end = ij_start - 1 + n_ij

  Allocate(w_d_all(n_ij,CC_n_vir,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'w_d_all in CC')

  w_d_all = 0.0D0

  Allocate(con_tmp(CC_n_occ,CC_n_occ,n_a,CC_mem_aa_D(CC_mpi_did+1)),stat=errnum)
  Call check_allocation(errnum,'con_tmp in CC')

  con_tmp = CC_ten_con

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

    n_b = CC_mem_aa_D(i_task)
    b_start = CC_index_aa_D(i_task)
    b_end = b_start - 1 + n_b

    if (i_recv.ne.CC_mpi_did+1) then

      Allocate(con_recv(CC_n_occ,CC_n_occ,n_a,CC_mem_aa_D(i_recv)),stat=errnum)
      Call check_allocation(errnum,'con_recv in CC')

      i_tmp = CC_n_occ * CC_n_occ * n_a
      s_tmp = Int(i_tmp,8) * Int(CC_mem_aa_D(CC_mpi_did+1),8)

      Call CC_mpi_real_isend(s_tmp,CC_ten_con,i_send-1,100,req1,CC_mpi_comm_domain)

      s_tmp = Int(i_tmp,8) * Int(CC_mem_aa_D(i_recv),8)
      Call CC_mpi_real_irecv(s_tmp,con_recv,i_recv-1,100,req2,CC_mpi_comm_domain)

    end if

    do a_run = 1, n_a
      do b_run = 1, n_b
        aaa = a_start - 1 + a_run
        bbb = b_start - 1 + b_run

        !$OMP PARALLEL Default(shared) &
        !$OMP Private(ij_run,code_ij,iii,jjj)
        !$OMP DO
        do ij_run = 1, n_ij
          code_ij = ij_start - 1 + ij_run
          Call CC_cl_decode(code_ij,iii,jjj,CC_n_occ,2)
          w_d_all(ij_run,aaa,bbb) = w_d_all(ij_run,aaa,bbb)      &
                                  + con_tmp(iii,jjj,a_run,b_run)

          w_d_all(ij_run,bbb,aaa) = w_d_all(ij_run,bbb,aaa)      &
                                  + con_tmp(jjj,iii,a_run,b_run) 

        end do
        !$OMP END DO
        !$OMP END PARALLEL

      end do
    end do

    Deallocate(con_tmp)

    if (i_recv.ne.CC_mpi_did+1) then

      Call MPI_WAIT(req1,stat1,errnum)
      Call MPI_WAIT(req2,stat2,errnum)

      Allocate(con_tmp(CC_n_occ,CC_n_occ,n_a,CC_mem_aa_D(i_recv)),stat=errnum)
      Call check_allocation(errnum,'con_tmp in CC')

      con_tmp = con_recv

      Deallocate(con_recv)

    end if

  end do

  CC_w_d(:,:,:,1) = CC_w_d(:,:,:,1) + w_d_all(:,a_start:a_end,:)

  Allocate(w_recv(n_ij,n_a,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'w_recv in CC')

  i_send = CC_mpi_gid + 1
  i_recv = CC_mpi_gid + 1

  do i_task_run = 1, CC_mpi_group_size - 1

    i_recv = i_recv + 1
    if (i_recv.gt.CC_mpi_group_size) then
      i_recv = 1
    end if

    i_send = i_send - 1
    if (i_send.lt.1) then
      i_send = CC_mpi_group_size
    end if

    n_a = CC_mem_aa_G(i_send)
    a_start = CC_index_aa_G(i_send)
    a_end = a_start - 1 + n_a

    Allocate(w_send(n_ij,n_a,CC_n_vir),stat=errnum)
    Call check_allocation(errnum,'w_send in CC')

    w_send = w_d_all(:,a_start:a_end,:)

    s_tmp = Int(n_ij,8) * Int(n_a*CC_n_vir,8)
    Call CC_mpi_real_isend(s_tmp,w_send,i_send-1,100,req1,CC_mpi_comm_group)

    s_tmp = Int(n_ij,8) * Int(CC_mem_aa_G(CC_mpi_gid+1)*CC_n_vir,8)
    Call CC_mpi_real_irecv(s_tmp,w_recv,i_recv-1,100,req2,CC_mpi_comm_group)

    Call MPI_WAIT(req1,stat1,errnum)
    Call MPI_WAIT(req2,stat2,errnum)

    CC_w_d(:,:,:,1) = CC_w_d(:,:,:,1) + w_recv

    Deallocate(w_send)

  end do

  Deallocate(w_recv)

  Deallocate(w_d_all)

  End Subroutine CC_cl_add_tc2w

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

