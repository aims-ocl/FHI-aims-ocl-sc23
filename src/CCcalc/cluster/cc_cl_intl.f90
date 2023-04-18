  Subroutine CC_cl_prepare_intl()

  Use basis
  Use prodbas
  Use hartree_fock
  Use CC_cl
  Use mpi_tasks

  Implicit None

  Integer (kind = 8) :: s_tmp
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,nnn
  Integer :: n_dom,n_grp
  Integer :: i_dom,i_grp
  Integer :: target_id,i_use
  Integer :: errnum,i_basbas,j_basbas,i_task,mpierr
  Integer :: i_col,i_row

  Double precision :: rlt,ddot
  Double precision , dimension(:,:,:,:) , allocatable :: intl_tmp,ovlp_3ks,ovlp_tmp
  Integer , dimension(n_tasks) :: n_loc_tmp
 
  ! Allocate RI coefficients
  Allocate(CC_RI(CC_mem_bas(CC_mpi_did+1), &
                 CC_n_state, CC_n_state,1),stat=errnum)
  Call check_allocation(errnum,'CC_RI in CC')

  ! Get ovlp_3ks matrix
  Allocate(ovlp_3ks(n_loc_prodbas,n_states,n_states,1),stat=errnum)
  Call check_allocation(errnum,'ovlp_3ks(1D) in CC')
  
  if (n_periodic.eq.0) then
    Call transform_ovlp3fn(n_states,ks_eigenvector,ovlp_3fn,ovlp_3ks)
  
    n_loc_tmp = 0
    n_loc_tmp(myid+1) = n_loc_prodbas
    Call CC_mpi_int_allreduce(n_tasks, n_loc_tmp, MPI_COMM_WORLD)
  
    ! Redistribute ovlp_3ks in domains
    CC_RI = 0.0D0
    i_basbas = 0
    do i_task = 1, n_tasks
  
      target_id = i_task - 1
  
      s_tmp = Int(CC_n_state,8) * Int(CC_n_state,8) * Int(n_loc_tmp(i_task),8)

      Allocate(ovlp_tmp(n_loc_tmp(i_task), CC_n_state, CC_n_state,1))
  
      if (myid.eq.target_id) then
        ovlp_tmp(:,:,:,1) = ovlp_3ks(:,CC_valence:n_states,CC_valence:n_states,1)
      else
        ovlp_tmp = 0.0D0
      end if
  
      Call CC_mpi_bcast(s_tmp, ovlp_tmp, target_id, MPI_COMM_WORLD)

      do iii = 1, n_loc_tmp(i_task)
        i_basbas = i_basbas + 1
        jjj = CC_index_bas(CC_mpi_did+1)
        kkk = CC_index_bas(CC_mpi_did+1) + CC_mem_bas(CC_mpi_did+1) - 1
        if ((i_basbas.ge.jjj).and.(i_basbas.le.kkk)) then
          j_basbas = i_basbas - CC_index_bas(CC_mpi_did+1) + 1
          CC_RI(j_basbas,:,:,1) = ovlp_tmp(iii,:,:,1)
        end if
      end do
  
      Deallocate(ovlp_tmp)
    end do
  
    Deallocate(ovlp_3ks)
  
    ! Extract unnecessary memory storage
    Deallocate(ovlp_3fn)
    Deallocate(map_prodbas)
  else
    Call CC_cl_RI_p2c()
  end if

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'RI coefficients are ready.'
  end if

  ! Deallocate the memory about the taburated grids.
  Call deallocate_grid_storage()
  
  ! Deallocate the memory about the auxiliary basis sets
  Call cleanup_basbas()

  ! Allocate RI tensors for b_aux
  ! CC_RI_B(bas,c,n)
  Allocate(CC_RI_B(CC_n_bas,CC_mem_aa_D(CC_mpi_did+1),CC_n_state),stat=errnum)
  Call check_allocation(errnum,'CC_RI_B in CC')

  Call cc_cl_redistr_RI()

  ! Allocate the integral tensors replicated saved in each domain
  ! For (ia|jb)
  ! CC_intl_iajb(i,j,a,b) ---- a
  Allocate(CC_intl_iajb(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1), &
                        CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'CC_intl_iajb in CC')

  Call CC_cl_calc_intl_iajb()

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral (ia|jb) tensor is ready'
  end if

  ! For (ki|ac)
  ! CC_intl_kiac(k,c,i,a)
  Allocate(CC_intl_kiac(CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),CC_n_occ, &
                        CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'CC_intl_kiac in CC')

  Call CC_cl_calc_intl_kiac()

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral (ki|ac) tensor is ready'
  end if

  ! For (ki|lc)
  ! CC_intl_kilc(k,c,i,l)
  Allocate(CC_intl_kilc(CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),CC_n_occ, &
                        CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'CC_intl_kilc in CC')

  Call CC_cl_calc_intl_kilc()

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral (ki|lc) tensor is ready'
  end if

  ! For (ac|kd)
  ! CC_intl_ackd(k,c,d,a)
  Allocate(CC_intl_ackd(CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),CC_n_vir, &
                        CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'CC_intl_ackd in CC')

  Call CC_cl_calc_intl_ackd()

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral (ac|kd) tensor is ready'
  end if

  ! For (ki|lj)
  ! CC_intl_kilj(kl,i,j)
  Allocate(CC_intl_kilj(CC_mem_kl_D(CC_mpi_did+1),CC_n_occ,CC_n_occ,1),stat=errnum)
  Call check_allocation(errnum,'CC_intl_kilj in CC')

  Call CC_cl_calc_intl_kilj()

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Integral (ki|lj) tensor is ready'
  end if

  if (CC_abcd_sv_flag) then

    ! For (ac|bd)
    ! CC_intl_acbd(c,d,a<=b,1)
    Allocate(CC_intl_acbd(CC_mem_aa_D(CC_mpi_did+1),CC_n_vir, &
                          CC_mem_ab_G(CC_mpi_gid+1),1),stat=errnum)
    Call check_allocation(errnum,'CC_intl_acbd in CC')

    Call CC_cl_calc_intl_acbd()

    if (myid.eq.0) then
      write(use_unit,"(2x,A)") 'Integral (ac|bd) tensor is ready'
    end if
  end if

  End Subroutine CC_cl_prepare_intl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine cc_cl_redistr_RI()

  Use CC_cl

  Implicit None

  Integer :: c_start,c_end,n_c,b_start,b_end,nnn

  Integer (kind = 8) :: s_tmp
  Integer :: errnum,source_id

  Integer :: i_task
  
  Double precision , dimension(:,:,:,:) , allocatable :: RI_tmp

  n_c = CC_mem_aa_D(CC_mpi_did+1)
  c_start = CC_index_aa_D(CC_mpi_did+1)
  c_end = c_start - 1 + n_c

  do i_task = 1, CC_mpi_domain_size

    nnn = CC_mem_bas(i_task)
    b_start = CC_index_bas(i_task)
    b_end = b_start + nnn - 1

    Allocate(RI_tmp(nnn,CC_n_vir,CC_n_state,1),stat=errnum)
    Call check_allocation(errnum,'RI_tmp in CC')

    source_id = i_task - 1
    if (CC_mpi_did.eq.source_id) then 
      RI_tmp(:,:,:,1) = CC_RI(:,CC_n_occ+1:CC_n_occ+CC_n_vir,:,1)
    end if

    s_tmp = Int(nnn,8) * Int(CC_n_vir*CC_n_state,8)
    Call CC_mpi_bcast(s_tmp, RI_tmp, source_id, CC_mpi_comm_domain)

    CC_RI_B(b_start:b_end,:,:) = RI_tmp(:,c_start:c_end,:,1)

    Deallocate(RI_tmp)

  end do

  End Subroutine cc_cl_redistr_RI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_calc_intl_iajb()

  Use CC_cl

  Implicit None

  Integer :: n_a,a_start,a_end,n_c,c_start,c_end,c_run,a_run
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,nnn
  Integer :: i_tmp,j_tmp

  Integer (kind = 8) :: s_tmp
  Integer :: i_task,target_id,errnum
  
  Double precision :: alpha,beta
  Double precision , dimension(:,:,:) , allocatable :: mat_A,mat_B
  Double precision , dimension(:,:,:,:) , allocatable :: intl_rlt
 
  ! Calculate integral (ia|jb)
  ! CC_intl_iajb(i,j,a,b) ---- a
  CC_intl_iajb = 0.0D0

  nnn = CC_mem_bas(CC_mpi_did+1)

  alpha = 1.0D0
  beta = 0.0D0

  n_c = CC_mem_aa_G(CC_mpi_gid+1)
  c_start = CC_index_aa_G(CC_mpi_gid+1)
  c_end = c_start - 1 + n_c

  do i_task = 1, CC_mpi_domain_size

    n_a = CC_mem_aa_D(i_task)
    a_start = CC_index_aa_D(i_task) 
    a_end = a_start - 1 + n_a

    i_tmp = CC_n_occ * n_a
    j_tmp = CC_n_occ * n_c

    Allocate(mat_A(nnn,CC_n_occ,n_a), stat=errnum)
    Call check_allocation(errnum,'mat_A in CC')

    Allocate(mat_B(nnn,CC_n_occ,n_c), stat=errnum)
    Call check_allocation(errnum,'mat_B in CC')

    Allocate(intl_rlt(CC_n_occ,n_a,CC_n_occ,n_c), stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')

    mat_A = CC_RI(:,1:CC_n_occ,CC_n_occ+a_start:CC_n_occ+a_end,1)
    mat_B = CC_RI(:,1:CC_n_occ,CC_n_occ+c_start:CC_n_occ+c_end,1)

    Call Dgemm('T','N',i_tmp,j_tmp,nnn,alpha,mat_A,nnn, &
               mat_B,nnn,beta,intl_rlt,i_tmp)

    target_id = i_task - 1
    s_tmp = Int(i_tmp,8) * Int(j_tmp,8)
    Call CC_mpi_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      !(i,j,a,b) <---- (i,a,j,b)
      do jjj = 1, CC_n_occ
        CC_intl_iajb(:,jjj,:,c_start:c_end) = intl_rlt(:,:,jjj,:)
      end do
    end if

    Deallocate(mat_A,mat_B,intl_rlt)

  end do

  i_tmp = CC_n_occ * CC_n_vir
  j_tmp = CC_n_occ * CC_mem_aa_D(CC_mpi_did+1)
 
  s_tmp = Int(i_tmp,8) * Int(j_tmp,8)
  Call CC_mpi_allreduce(s_tmp, CC_intl_iajb, CC_mpi_comm_group)
 
  End Subroutine CC_cl_calc_intl_iajb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_calc_intl_kiac()

  Use CC_cl

  Implicit None

  Integer :: n_c,c_start,c_end,n_a,a_start,a_end,c_run,a_run
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,nnn
  Integer :: i_tmp,j_tmp

  Integer (kind = 8) :: s_tmp
  Integer :: i_task,target_id,errnum
  
  Double precision :: alpha,beta
  Double precision , dimension(:,:,:) , allocatable :: mat_A,mat_B
  Double precision , dimension(:,:,:,:) , allocatable :: intl_rlt
 
  ! Calculate integral (ki|ac)
  ! CC_intl_kiac(k,c,i,a)
  CC_intl_kiac = 0.0D0

  nnn = CC_mem_bas(CC_mpi_did+1)

  alpha = 1.0D0
  beta = 0.0D0

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  do i_task = 1, CC_mpi_domain_size

    n_c = CC_mem_aa_D(i_task)
    c_start = CC_index_aa_D(i_task)
    c_end = c_start - 1 + n_c

    i_tmp = CC_n_occ * CC_n_occ
    j_tmp = n_a * n_c

    Allocate(mat_A(nnn,CC_n_occ,CC_n_occ), stat=errnum)
    Call check_allocation(errnum,'mat_A in CC')

    Allocate(mat_B(nnn,n_c,n_a), stat=errnum)
    Call check_allocation(errnum,'mat_B in CC')

    Allocate(intl_rlt(CC_n_occ,CC_n_occ,n_c,n_a), stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')

    mat_A = CC_RI(:,1:CC_n_occ,1:CC_n_occ,1)
    mat_B = CC_RI(:,CC_n_occ+c_start:CC_n_occ+c_end, &
                    CC_n_occ+a_start:CC_n_occ+a_end,1)

    Call Dgemm('T','N',i_tmp,j_tmp,nnn,alpha,mat_A,nnn, &
               mat_B,nnn,beta,intl_rlt,i_tmp)

    target_id = i_task - 1
    s_tmp = Int(i_tmp,8) * Int(j_tmp,8)
    Call CC_mpi_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

    if (CC_mpi_did.eq.target_id) then
      ! (k,i,c,a) ----> (k,c,i,a)
      do ccc = 1, n_c
        CC_intl_kiac(:,ccc,:,a_start:a_end) = intl_rlt(:,:,ccc,:)
      end do
    end if

    Deallocate(mat_A,mat_B,intl_rlt)

  end do

  i_tmp = CC_mem_aa_D(CC_mpi_did+1) * CC_n_vir
  j_tmp = CC_n_occ * CC_n_occ

  s_tmp = Int(i_tmp,8) * Int(j_tmp,8)
  Call CC_mpi_allreduce(s_tmp, CC_intl_kiac, CC_mpi_comm_group)

  End Subroutine CC_cl_calc_intl_kiac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_calc_intl_kilc()

  Use CC_cl

  Implicit None

  Integer :: n_i,i_start,i_end,n_a,a_start,a_end,i_run,a_run
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,nnn
  Integer :: i_tmp,j_tmp

  Integer (kind = 8) :: s_tmp
  Integer :: i_task,target_id,errnum
  
  Double precision :: alpha,beta
  Double precision , dimension(:,:,:) , allocatable :: mat_A,mat_B
  Double precision , dimension(:,:,:,:) , allocatable :: intl_rlt
 
  ! Calculate integral (ki|lc)
  ! CC_intl_kilc(k,c,i,l)
  CC_intl_kilc = 0.0D0

  nnn = CC_mem_bas(CC_mpi_did+1)

  alpha = 1.0D0
  beta = 0.0D0

  n_i = CC_mem_ii_G(CC_mpi_gid+1)
  i_start = CC_index_ii_G(CC_mpi_gid+1)
  i_end = i_start - 1 + n_i

  if (n_i.ne.0) then

    do i_task = 1, CC_mpi_domain_size

      n_a = CC_mem_aa_D(i_task)
      a_start = CC_index_aa_D(i_task)
      a_end = a_start - 1 + n_a

      i_tmp = CC_n_occ * CC_n_occ
      j_tmp = n_i * n_a

      Allocate(mat_A(nnn,CC_n_occ,CC_n_occ), stat=errnum)
      Call check_allocation(errnum,'mat_A in CC')

      Allocate(mat_B(nnn,n_a,n_i), stat=errnum)
      Call check_allocation(errnum,'mat_B in CC')

      Allocate(intl_rlt(CC_n_occ,CC_n_occ,n_a,n_i), stat=errnum)
      Call check_allocation(errnum,'intl_rlt in CC')

      mat_A = CC_RI(:,1:CC_n_occ,1:CC_n_occ,1)
      mat_B = CC_RI(:,CC_n_occ+a_start:CC_n_occ+a_end,i_start:i_end,1)

      Call Dgemm('T','N',i_tmp,j_tmp,nnn,alpha,mat_A,nnn, &
                 mat_B,nnn,beta,intl_rlt,i_tmp)

      target_id = i_task - 1
      s_tmp = Int(i_tmp,8) * Int(j_tmp,8)
      Call CC_mpi_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

      if (CC_mpi_did.eq.target_id) then
        do ccc = 1, n_a
          !(k,i,c,l) ----> (k,c,i,l)
          CC_intl_kilc(:,ccc,:,i_start:i_end) = intl_rlt(:,:,ccc,:)
        end do
      end if

      Deallocate(mat_A,mat_B,intl_rlt)

    end do

  end if

  i_tmp = CC_n_occ * CC_n_occ
  j_tmp = CC_mem_aa_D(CC_mpi_did+1) * CC_n_occ
 
  s_tmp = Int(i_tmp,8) * Int(j_tmp,8)
  Call CC_mpi_allreduce(s_tmp, CC_intl_kilc, CC_mpi_comm_group)
 
  End Subroutine CC_cl_calc_intl_kilc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_calc_intl_kilj()

  Use CC_cl

  Implicit None

  Integer :: n_i,i_start,i_end,n_k,k_start,k_end,i_run,k_run,kl
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,nnn
  Integer :: i_tmp,j_tmp

  Integer (kind = 8) :: s_tmp
  Integer :: i_task,target_id,errnum
  
  Double precision :: alpha,beta
  Double precision , dimension(:,:,:) , allocatable :: mat_A,mat_B
  Double precision , dimension(:,:,:,:) , allocatable :: intl_rlt
 
  ! Calculate integral (ki|lj)
  ! CC_intl_kilj(kl,i,j)
  CC_intl_kilj = 0.0D0

  nnn = CC_mem_bas(CC_mpi_did+1)

  alpha = 1.0D0
  beta = 0.0D0

  n_i = CC_mem_ii_G(CC_mpi_gid+1)
  i_start = CC_index_ii_G(CC_mpi_gid+1)
  i_end = i_start - 1 + n_i

  if (n_i.ge.1) then

    do i_task = 1, CC_mpi_domain_size

      n_k = CC_mem_kl_D(i_task)
      k_start = CC_index_kl_D(i_task)
      k_end = k_start - 1 + n_k

      Allocate(mat_A(nnn,n_i,1), stat=errnum)
      Call check_allocation(errnum,'mat_A in CC')

      Allocate(mat_B(nnn,CC_n_occ,1), stat=errnum)
      Call check_allocation(errnum,'mat_B in CC')

      Allocate(intl_rlt(n_i,CC_n_occ,n_k,1), stat=errnum)
      Call check_allocation(errnum,'intl_rlt in CC')

      i_tmp = n_i 
      j_tmp = CC_n_occ 

      do k_run = 1, n_k

        kl = k_start - 1 + k_run
      
        Call CC_cl_decode(kl,kkk,lll,CC_n_occ,1)

        mat_A(:,:,1) = CC_RI(:,i_start:i_end,kkk,1)
        mat_B(:,:,1) = CC_RI(:,1:CC_n_occ,lll,1)

        Call Dgemm('T','N',i_tmp,j_tmp,nnn,alpha,mat_A,nnn, &
                   mat_B,nnn,beta,intl_rlt(:,:,k_run,1),i_tmp)
      end do

      target_id = i_task - 1
      s_tmp = Int(i_tmp,8) * Int(j_tmp,8) * Int(n_k)
      Call CC_mpi_reduce(s_tmp,intl_rlt,target_id,CC_mpi_comm_domain)

      if (CC_mpi_did.eq.target_id) then
        ! (i,j,kl) ----> (kl,i,j)
        do kkk = 1, n_k
          CC_intl_kilj(kkk,i_start:i_end,:,1) = intl_rlt(:,:,kkk,1)
        end do
      end if

      Deallocate(mat_A,mat_B,intl_rlt)

    end do
  end if

  End Subroutine CC_cl_calc_intl_kilj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_calc_intl_ackd()

  Use CC_cl

  Implicit None

  Integer :: n_a,a_start,a_end
  Integer :: c_start,c_end,n_c
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,nnn,a_run
  Integer :: i_tmp,j_tmp,l_tmp,m_tmp

  Integer (kind = 8) :: s_tmp
  Integer :: errnum

  Integer :: i_task,i_task_run
  Integer :: req1,req2,i_send,i_recv
  Integer , dimension(MPI_STATUS_SIZE) :: stat1,stat2
  
  Double precision :: alpha,beta
  Double precision , dimension(:,:,:,:) , allocatable :: matA,matB,intl_rlt
  Double precision , dimension(:,:,:,:) , allocatable :: RI_tmp,RI_recv 

  ! Calculate integral (ac|kd)
  ! CC_intl_ackd(k,c,d,a)
  CC_intl_ackd = 0.0D0

  nnn = CC_mem_bas(CC_mpi_did+1)

  alpha = 1.0D0
  beta = 0.0D0

  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  a_start = CC_index_aa_G(CC_mpi_gid+1)
  a_end = a_start - 1 + n_a

  n_c = CC_mem_aa_D(CC_mpi_did+1)
  c_start = CC_index_aa_D(CC_mpi_did+1)
  c_end = c_start - 1 + n_c

  i_task = CC_mpi_did
  i_send = CC_mpi_did + 1
  i_recv = CC_mpi_did + 1

  Allocate(RI_tmp(nnn,CC_n_state,CC_n_state,1),stat=errnum)
  Call check_allocation(errnum,'RI_tmp in CC')

  RI_tmp = CC_RI

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

    nnn = CC_mem_bas(i_task)

    if (i_recv.ne.CC_mpi_did+1) then

      Allocate(RI_recv(CC_mem_bas(i_recv),CC_n_state,CC_n_state,1),stat=errnum)
      Call check_allocation(errnum,'RI_recv in CC')

      i_tmp = CC_n_state * CC_n_state
      s_tmp = Int(i_tmp,8) * Int(CC_mem_bas(CC_mpi_did+1),8)

      Call CC_mpi_real_isend(s_tmp,CC_RI,i_send-1,100,req1,CC_mpi_comm_domain)

      s_tmp = Int(i_tmp,8) * Int(CC_mem_bas(i_recv),8)
      Call CC_mpi_real_irecv(s_tmp,RI_recv,i_recv-1,100,req2,CC_mpi_comm_domain)

    end if

    Allocate(matA(nnn,CC_n_occ,CC_n_vir,1), stat=errnum)
    Call check_allocation(errnum,'mat_A in CC')

    Allocate(matB(nnn,n_c,1,1), stat=errnum)
    Call check_allocation(errnum,'mat_B in CC')

    Allocate(intl_rlt(CC_n_occ,CC_n_vir,n_c,1), stat=errnum)
    Call check_allocation(errnum,'intl_rlt in CC')

    i_tmp = CC_n_occ * CC_n_vir

    do a_run = a_start, a_end

      aaa = a_run + CC_n_occ
      
      matA(:,:,:,1) = RI_tmp(:,1:CC_n_occ,CC_n_occ+1:CC_n_occ+CC_n_vir,1)
      matB(:,:,1,1) = RI_tmp(:,CC_n_occ+c_start:CC_n_occ+c_end,aaa,1)

      Call Dgemm('T','N',i_tmp,n_c,nnn,alpha,matA,nnn, &
                 matB,nnn,beta,intl_rlt,i_tmp)

      do ddd = 1, CC_n_vir
        CC_intl_ackd(:,:,ddd,a_run) = CC_intl_ackd(:,:,ddd,a_run) &
                                    + intl_rlt(:,ddd,:,1)
      end do

    end do

    Deallocate(matA,matB,intl_rlt)

    Deallocate(RI_tmp)

    if (i_recv.ne.CC_mpi_did+1) then

      Call MPI_WAIT(req1,stat1,errnum)
      Call MPI_WAIT(req2,stat2,errnum)

      Allocate(RI_tmp(CC_mem_bas(i_recv),CC_n_state,CC_n_state,1),stat=errnum)
      Call check_allocation(errnum,'RI_tmp in CC')

      RI_tmp = RI_recv

      Deallocate(RI_recv)

    end if

  end do

  i_tmp = CC_mem_aa_D(CC_mpi_did+1) * CC_n_occ
  j_tmp = CC_n_vir * CC_n_vir
 
  s_tmp = Int(i_tmp,8) * Int(j_tmp,8)
  Call CC_mpi_allreduce(s_tmp, CC_intl_ackd, CC_mpi_comm_group)

  End Subroutine CC_cl_calc_intl_ackd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_calc_intl_acbd()

  Use CC_cl

  Implicit None

  Integer :: n_c,c_start,c_end,n_ab,ab_start,ab_end,code_ab,cd_run,ab_run
  Integer :: iii,jjj,kkk,lll,aaa,bbb,ccc,ddd,nnn

  Integer (kind = 8) :: s_tmp
  Integer :: i_task,i_use,i_sv,target_id,errnum,OMP_GET_NUM_THREADS,nth
  
  Double precision :: alpha,beta,ddot
  Double precision , dimension(:,:) , allocatable :: matA,matB
  Double precision , dimension(:,:,:) , allocatable :: intl_rlt
 
  ! Calculate integral (ac|bd)
  ! CC_intl_acbd(c,d,a<=b,1)
  CC_intl_acbd = 0.0D0

  nnn = CC_mem_bas(CC_mpi_did+1)

  alpha = 1.0D0
  beta = 0.0D0

  n_ab = CC_mem_ab_G(CC_mpi_gid+1)
  ab_start = CC_index_ab_G(CC_mpi_gid+1)
  ab_end = ab_start - 1 + n_ab

  n_c = CC_mem_aa_D(CC_mpi_did+1)
  c_start = CC_index_aa_D(CC_mpi_did+1)
  c_end = c_start - 1 + n_c

  s_tmp = Int(CC_n_vir * CC_n_vir,8)

  Allocate(matA(nnn,CC_n_vir), stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  Allocate(matB(nnn,CC_n_vir), stat=errnum)
  Call check_allocation(errnum,'mat_B in CC')

  Allocate(intl_rlt(CC_n_vir,CC_n_vir,2), stat=errnum)
  Call check_allocation(errnum,'mat_A in CC')

  nth = 1

  !$OMP PARALLEL Default(shared)
  !$OMP SINGLE
  !$ nth = OMP_GET_NUM_THREADS()
  !$OMP END SINGLE
  !$OMP END PARALLEL

  !!$ if (nth.gt.1) then
  !!$   Call omp_set_nested(1)
  !!$   Call mkl_set_dynamic(0)
  !!$   Call MKL_SET_NUM_THREADS(nth-1)                                                              
  !!$ end if

  !$OMP PARALLEL Default(Shared) &
  !$OMP Private(i_use,i_sv,code_ab,ab_run,aaa,bbb)
  i_use = 1
  code_ab = ab_start - 1
  do ab_run = 1, n_ab
    code_ab = code_ab + 1

    Call CC_cl_decode(code_ab,aaa,bbb,CC_n_vir,2)

    aaa = aaa + CC_n_occ
    bbb = bbb + CC_n_occ

    !$OMP SINGLE
    matA = CC_RI(:,aaa,CC_n_occ+1:CC_n_occ+CC_n_vir,1)
    matB = CC_RI(:,bbb,CC_n_occ+1:CC_n_occ+CC_n_vir,1)

    Call Dgemm('T','N',CC_n_vir,CC_n_vir,nnn,alpha,matA,nnn, &
               matB,nnn,beta,intl_rlt(:,:,i_use),CC_n_vir)
    !$OMP END SINGLE

    i_sv = i_use
    i_use = 3 - i_use

    !$OMP SINGLE
    Call CC_mpi_allreduce(s_tmp,intl_rlt(:,:,i_sv),CC_mpi_comm_domain)
    CC_intl_acbd(:,:,ab_run,1) = intl_rlt(c_start:c_end,:,i_sv)
    !$OMP END SINGLE NOWAIT
  end do
  !$OMP END PARALLEL

  !!$ if (nth.gt.1) then
  !!$   Call MKL_SET_NUM_THREADS(nth)
  !!$   Call omp_set_nested(0)
  !!$   Call mkl_set_dynamic(1)
  !!$ end if

  Deallocate(matA,matB,intl_rlt)

  End Subroutine CC_cl_calc_intl_acbd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
