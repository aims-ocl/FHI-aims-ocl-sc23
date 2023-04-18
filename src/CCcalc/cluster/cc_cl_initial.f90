  Subroutine CC_cl_initial()

  Use basis
  Use prodbas
  Use hartree_fock
  Use mpi_tasks
  Use CC_cl
  use geometry, only: species


  Implicit None

  Integer , parameter :: ures = 100

  Integer (kind = 8) :: i_tmp,tmp_mem
  Integer :: j_tmp,iii,jjj,mpierr
  Integer :: i_state,errnum,i_task,hdferr
  Double precision , Parameter :: ne_thresh = 1.0D-6
  Double precision :: CC_mem_demand,CC_disk_demand,large_mem
  Double precision , dimension(:) , allocatable :: run_tmp

  Character (len = 20) :: large_mem_dis

  Character (len=80) :: resline,blk

  Logical :: res_ex

! Initialization starts
  CC_res_pt = 0

  Allocate(CC_res_inf(CC_n_domain+1), stat=errnum)
  Call check_allocation(errnum,'CC_res_inf in CC')

  CC_res_inf = 0

! read restart informations if necessary
  if (CC_restart_flag) then
    if (myid.eq.0) then
      Inquire(file='CCrestart',exist=res_ex)
      if (res_ex) then
        open(unit=ures,file='CCrestart')
        read(ures,*) CC_res_pt
        do iii = 1, CC_n_domain
          read(ures,*) CC_res_inf(iii+1)
        end do
        read(ures,*) CC_E_PT
        close(ures)
        CC_res_inf(1) = CC_res_pt
      else
        open(unit=ures,file='CCrestart')
        write(ures,*) 0
        do iii = 1, CC_n_domain
          write(ures,*) 0
        end do
        write(ures,*) 0.0D0
        close(ures)
      end if
    end if

    Call CC_mpi_int_bcast(CC_n_domain+1,CC_res_inf,0,MPI_COMM_WORLD)
    CC_res_pt = CC_res_inf(1)
  end if

! output mission information
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") '--------------------------------------------'
    write(use_unit,"(2x,A)") 'Coupled Cluster(CC) calculation starts'
    write(use_unit,"(2x,A)") '--------------------------------------------'
    write(use_unit,"(A)") ' '
    if (CC_calc_method.eq.1) then
      write(use_unit,"(2x,A)") 'CCSD calculation for closed shell molecules'
    else if (CC_calc_method.eq.2) then
      write(use_unit,"(2x,A)") 'CCSD(T) calculation for closed shell molecules'
    end if
    write(use_unit,"(A)") ' '
    write(use_unit,"(2x,A)") 'Initialization of CCSD calculation'
    if (CC_res_pt.eq.0) then
      write(use_unit,"(2x,A)") 'Restart information will be saved'
    else if ((CC_res_pt.eq.-1).and.(CC_calc_method.eq.2)) then
      write(use_unit,"(2x,A)") 'Restarting perturbative triples calculation'
    else
      write(use_unit,"(2x,A)") 'Restarting unfinished CCSD calculation'
    end if
  end if

! Determine the reference (HF) configuration and number of electrons

  CC_n_elec = 0
  
  do i_state = 1, n_states
    if (dabs(occ_numbers(i_state,1,1)-2.0D0).le.ne_thresh) then
      CC_n_elec = CC_n_elec + 1
    end if
  end do
  
  ! Calculate the number of valence state
  if (flag_frozen_core_postSCF) then
    call count_frozen_core_states(CC_valence)
    CC_n_fc = CC_valence - 1
    if(myid.eq.0) then
      write(use_unit,'(2X,A,I8)') "| Frozen core number      ::", CC_n_fc
    endif
  else if (flag_frozen_core) then
    CC_valence = i_start_mp2
    CC_n_fc = CC_valence - 1
    if(myid.eq.0) then
      write(use_unit,'(2X,A,I8)') "| Frozen core number      ::", CC_n_fc
    endif
  else
    CC_valence = 1
    CC_n_fc = 0
  endif
  
  CC_n_occ = CC_n_elec - CC_n_fc
  CC_n_vir = n_states - CC_n_elec
  CC_n_state = n_states - CC_n_fc
  CC_n_bas = n_basbas

  CC_n_s = CC_n_occ * CC_n_vir

  i_tmp = Int(CC_n_occ,8) * Int(CC_n_vir,8)
  CC_n_d = (i_tmp + 1) * i_tmp / 2

  CC_n_config = CC_n_d + Int(CC_n_s,8)

  ! Create MPI domains and groups
  Call CC_mpi_init()

  if (myid.eq.0) then
    write(use_unit,"(A)") ' '
    write(use_unit,"(2x,A)") 'MPI grouping...'
    write(use_unit,"(2x,A,I3)") 'Number of domains: ',CC_n_domain
    write(use_unit,"(2x,A,I3)") 'Number of tasks in each domain: ',CC_mpi_domain_size
    write(use_unit,"(2x,A)") 'Acceleration method: DIIS'
    write(use_unit,"(2x,A,I7)") 'Number of valence orbitals: ',CC_n_occ
    write(use_unit,"(2x,A,I7)") 'Number of unoccupied orbitals: ',CC_n_vir
    write(use_unit,"(2x,A)") 'Number of configurations'
    write(use_unit,"(10x,A16,8x,A16)") 's-excitation','d-excitation'
    write(use_unit,"(10x,I16,8x,I16)") CC_n_s,CC_n_d
    write(use_unit,"(2x,A,I16,2x,A)") 'Totally ',CC_n_config,'configurations are considered'
  end if

  ! Determine saving strategy
  CC_sv_strategy = 1
  if ((CC_sv_control.eq.'global').or.(CC_sv_control.eq.'GLOBAL')) then
    CC_sv_strategy = 2
  end if

  Call CC_cl_mem_distribution(CC_mem_demand,large_mem,large_mem_dis,tmp_mem)

  if (myid.eq.0) then
    iii = len_trim(large_mem_dis)
    write(use_unit,"(2x,A,F14.6,A)") 'At least ',CC_mem_demand, &
                                     'GB memory space is needed for each task'
    write(use_unit,"(2x,A,A,A,F14.6,A)") 'The largest tensor is ',large_mem_dis(1:iii), &
                                         ', cost',large_mem, 'GB memory for each task'
    write(use_unit,"(2x,A)") 'Preparing four-index integrals...'
  end if

  ! Prepare integrals
  Call CC_cl_prepare_intl()

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Four-index integrals are ready.'
    write(use_unit,"(2x,A)") 'Allocate tensors used in CCSD calculation'
  end if

  ! Allocate all vectors
  ! Allocate run-time temporary vector
  i_tmp = 20 * Int(CC_n_vir * CC_n_vir,8) + tmp_mem

  Allocate(run_tmp(i_tmp), stat=errnum)
  Call check_allocation(errnum,'run-time tmp in CC')

  ! For t_s (i,a)
  Allocate(CC_t_s(CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'CC_t_s in CC')

  ! For t_d (i,j,a,b)
  Allocate(CC_t_d(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'CC_t_s in CC')

  ! For CC_h_ik(i,k)
  Allocate(CC_h_ik(CC_n_occ,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'CC_h_ik in CC')

  ! For CC_g_ik(i,k)
  Allocate(CC_g_ik(CC_n_occ,CC_n_occ),stat=errnum)
  Call check_allocation(errnum,'CC_g_ik in CC')

  ! For CC_h_ca(c,a)
  Allocate(CC_h_ca(CC_n_vir,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'CC_h_ca in CC')

  ! For CC_g_ca(c,a)
  Allocate(CC_g_ca(CC_n_vir,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'CC_g_ca in CC')

  ! For CC_h_ck(k,c)
  Allocate(CC_h_ck(CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'CC_h_ck in CC')

  ! DIIS B matrix
  Allocate(CC_DIIS_Bmat(CC_DIIS_n_sv,CC_DIIS_n_sv),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_Bmat in CC')

  CC_DIIS_Bmat = 0.0D0  

  ! DIIS tau vector
  Allocate(CC_DIIS_tau(CC_DIIS_n_sv),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_tau in CC')

  ! For aux_a (kl,i,j)
  Allocate(CC_a_aux(CC_mem_kl_D(CC_mpi_did+1),CC_n_occ, &
                    CC_n_occ,1),stat=errnum)
  Call check_allocation(errnum,'CC_a_aux in CC')

  ! For aux_j and aux_k (k,c,i,a)
!  Allocate(CC_j_aux(CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),CC_n_occ, &
!                    CC_mem_aa_G(CC_mpi_gid+1)),stat=errnum)
!  Call check_allocation(errnum,'CC_j_aux in CC')
!
!  Allocate(CC_k_aux(CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),CC_n_occ, &
!                    CC_mem_aa_G(CC_mpi_gid+1)),stat=errnum)
!  Call check_allocation(errnum,'CC_k_aux in CC')

  ! For tensor contraction result
  ! (i,j,a,b)
  Allocate(CC_ten_con(CC_n_occ,CC_n_occ,CC_mem_aa_G(CC_mpi_gid+1), &
                      CC_mem_aa_D(CC_mpi_did+1)),stat=errnum)
  Call check_allocation(errnum,'CC_ten_con in CC')

  ! For w_s
  Allocate(CC_w_s(CC_n_occ,CC_n_vir),stat=errnum)
  Call check_allocation(errnum,'CC_w_s in CC')

  ! For w_d (i<=j,a,b)
  Allocate(CC_w_d(CC_mem_ij_D(CC_mpi_did+1),CC_mem_aa_G(CC_mpi_gid+1), &
                  CC_n_vir,1),stat=errnum)
  Call check_allocation(errnum,'CC_w_d in CC')

  Allocate(CC_DIIS_t_s(CC_n_occ,CC_n_vir,CC_DIIS_n_sv,1),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_t_s in CC')
  
  CC_DIIS_t_s = 0.0D0

  Allocate(CC_DIIS_r_s(CC_n_occ,CC_n_vir,CC_DIIS_n_sv,1),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_r_s in CC')
  
  CC_DIIS_r_s = 0.0D0

  Allocate(CC_DIIS_t_d(CC_mem_ij_D(CC_mpi_did+1), &
                       CC_mem_aa_G(CC_mpi_gid+1),CC_n_vir,CC_DIIS_n_sv),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_t_d in CC')
  
  CC_DIIS_t_d = 0.0D0

  Allocate(CC_DIIS_r_d(CC_mem_ij_D(CC_mpi_did+1), &
                       CC_mem_aa_G(CC_mpi_gid+1),CC_n_vir,CC_DIIS_n_sv),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_r_d in CC')
  
  CC_DIIS_r_d = 0.0D0

  Deallocate(run_tmp)

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Tensors allocated successfully.'
  end if

  ! Set up reference (HF) energy
  E_HF = total_energy

  CC_RLE_flag_singular = .false.
  CC_RLE_flag_sv_overflow = .false.

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'End of initialization'
  end if

  End Subroutine CC_cl_initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_mem_distribution(CC_mem_demand,large_mem,la_dis,mt)

  Use CC_cl

  Implicit None

  Double precision :: CC_mem_demand,large_mem
  Integer :: n_i,n_a,n_b,n_ij,n_ab,nnn
  Integer :: i_tmp,i_tmp2,j_tmp,errnum,nth,OMP_GET_NUM_THREADS
  Integer (kind = 8) :: s_tmp,la_mem,la_tmp,mt,mt1
  Character (len = 20) :: la_dis

  ! Determine the distributiuon of vectors on each task
  ! For occupied orbitals
  Allocate(CC_mem_ii_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_ii_D in CC')

  Allocate(CC_index_ii_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_ii_D in CC')

  Call CC_cl_vec_distr(1,CC_n_occ,CC_mpi_domain_size,CC_mem_ii_D,CC_index_ii_D)

  Allocate(CC_mem_ii_G(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_ii_G in CC')

  Allocate(CC_index_ii_G(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_ii_G in CC')

  Call CC_cl_vec_distr(1,CC_n_occ,CC_mpi_group_size,CC_mem_ii_G,CC_index_ii_G)

  Allocate(CC_mem_ij_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_ij_D in CC')

  Allocate(CC_index_ij_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_ij_D in CC')

  j_tmp = CC_n_occ * (CC_n_occ + 1) / 2
  Call CC_cl_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_ij_D,CC_index_ij_D)

  Allocate(CC_mem_ij_G(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_ii_G in CC')

  Allocate(CC_index_ij_G(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_ij_G in CC')

  j_tmp = CC_n_occ * (CC_n_occ - 1) / 2
  Call CC_cl_vec_distr(1,j_tmp,CC_mpi_group_size,CC_mem_ij_G,CC_index_ij_G)

  Allocate(CC_mem_kl_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_kl_D in CC')

  Allocate(CC_index_kl_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_kl_D in CC')

  j_tmp = CC_n_occ * CC_n_occ
  Call CC_cl_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_kl_D,CC_index_kl_D)

  ! For unoccupied orbitals
  Allocate(CC_mem_aa_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_aa_D in CC')

  Allocate(CC_index_aa_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_aa_D in CC')

  Call CC_cl_vec_distr(1,CC_n_vir,CC_mpi_domain_size,CC_mem_aa_D,CC_index_aa_D)

  Allocate(CC_mem_aa_G(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_aa_D in CC')

  Allocate(CC_index_aa_G(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_aa_G in CC')

  Call CC_cl_vec_distr(1,CC_n_vir,CC_mpi_group_size,CC_mem_aa_G,CC_index_aa_G)

  Allocate(CC_mem_ab_G(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_ab_G in CC')

  Allocate(CC_index_ab_G(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_ab_G in CC')

  j_tmp = CC_n_vir * (CC_n_vir + 1) / 2
  Call CC_cl_vec_distr(1,j_tmp,CC_mpi_group_size,CC_mem_ab_G,CC_index_ab_G)

  ! For cd
  Allocate(CC_mem_cd_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_cd_D in CC')

  Allocate(CC_index_cd_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_cd_D in CC')

  j_tmp = CC_n_vir * (CC_n_vir + 1) / 2
  Call CC_cl_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_cd_D,CC_index_cd_D)

  ! For CC_n_bas
  Allocate(CC_mem_bas(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_bas in CC')

  Allocate(CC_index_bas(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_bas in CC')

  Call CC_cl_vec_distr(1,CC_n_bas,CC_mpi_domain_size,CC_mem_bas,CC_index_bas)

  ! Calculate memory needed
  s_tmp = CC_mem_aa_D(1) * CC_n_vir * CC_n_vir * 2
  if (s_tmp < CC_work_tmp_size) then
    CC_tmp_mem = CC_work_tmp_size
  else
    CC_tmp_mem = s_tmp * 4
  end if

  n_i = CC_mem_ii_G(CC_mpi_gid+1)
  n_b = CC_mem_aa_D(CC_mpi_did+1)
  n_ij = CC_mem_ij_D(CC_mpi_did+1)
  n_a = CC_mem_aa_G(CC_mpi_gid+1)
  n_ab = CC_mem_ab_G(CC_mpi_gid+1)
  nnn = CC_mem_bas(CC_mpi_did+1)

  nth = 1
 
  !$OMP PARALLEL Default(shared)
  !$OMP SINGLE
  !$ nth = OMP_GET_NUM_THREADS()
  !$OMP END SINGLE
  !$OMP END PARALLEL

  s_tmp = Int(CC_tmp_mem,8) + Int(CC_MPI_max_len,8)
  la_mem = 0
  la_dis = 'temporary vectors'

  ! RI coefficients
  j_tmp = CC_n_state * CC_n_state
  la_tmp = Int(j_tmp,8) * Int(nnn,8)
  j_tmp = n_b * CC_n_state
  la_tmp = la_tmp + Int(j_tmp,8) * Int(CC_n_bas,8)
  s_tmp = s_tmp + la_tmp
  if (la_tmp.gt.la_mem) then
    la_mem = la_tmp
    la_dis = 'RI coefficients'
  end if

  ! (ia|jb) and (ki|ac)
  j_tmp = CC_n_occ * CC_n_vir * CC_n_occ
  la_tmp = Int(j_tmp,8) * Int(n_b,8)
  s_tmp = s_tmp + la_tmp * 2
  if (la_tmp.gt.la_mem) then
    la_mem = la_tmp
    la_dis = 'integral (ia|jb)'
  end if

  ! (ki|lc)
  j_tmp = CC_n_occ * CC_n_occ * CC_n_occ
  la_tmp = Int(j_tmp,8) * Int(n_b,8)
  s_tmp = s_tmp + la_tmp
  if (la_tmp.gt.la_mem) then
    la_mem = la_tmp
    la_dis = 'integral (ki|lc)'
  end if

  ! (ac|kd)
  j_tmp = CC_n_vir * CC_n_vir
  la_tmp = Int(j_tmp,8) * Int(n_b * CC_n_occ,8)
  s_tmp = s_tmp + la_tmp
  if (la_tmp.gt.la_mem) then
    la_mem = la_tmp
    la_dis = 'integral (ac|kd)'
  end if

  ! (ki|lj)
  j_tmp = CC_n_occ * CC_n_occ
  la_tmp = Int(j_tmp,8) * Int(n_ij,8)
  s_tmp = s_tmp + la_tmp
  if (la_tmp.gt.la_mem) then
    la_mem = la_tmp
    la_dis = 'integral (ki|lj)'
  end if

  ! (ac|bd) if necessary
  !if (CC_abcd_sv_flag) then

  !  j_tmp = n_b * CC_n_vir
  !  la_tmp = Int(j_tmp,8) * Int(n_ab,8)
  !  s_tmp = s_tmp + la_tmp
  !  if (la_tmp.gt.la_mem) then
  !    la_mem = la_tmp
  !    la_dis = 'integral (ac|bd)'
  !  end if

  !end if

  ! CC_t_d
  j_tmp = CC_n_occ * CC_n_vir * CC_n_occ
  la_tmp = Int(j_tmp,8) * Int(n_b,8)
  s_tmp = s_tmp + la_tmp
  if (la_tmp.gt.la_mem) then
    la_mem = la_tmp
    la_dis = 'amplitudes'
  end if

  ! CC_w_d
  la_tmp = Int(n_ij,8) * Int(n_a*CC_n_vir,8)
  s_tmp = s_tmp + la_tmp
  if (la_tmp.gt.la_mem) then
    la_mem = la_tmp
    la_dis = 'amplitudes'
  end if

  ! j_aux and k_aux
  j_tmp = n_a * CC_n_occ * CC_n_occ
  la_tmp = Int(j_tmp,8) * Int(n_b,8)
  s_tmp = s_tmp + la_tmp
  if (la_tmp.gt.la_mem) then
    la_mem = la_tmp
    la_dis = 'intermediate J'
  end if

  ! a_aux
  j_tmp = CC_n_occ * CC_n_occ
  la_tmp = Int(j_tmp,8) * Int(n_ij,8)
  s_tmp = s_tmp + la_tmp
  if (la_tmp.gt.la_mem) then
    la_mem = la_tmp
    la_dis = 'intermediate A'
  end if

  ! tensor contraction intermediate
  j_tmp = n_a * CC_n_occ * CC_n_occ
  la_tmp = Int(j_tmp,8) * Int(n_b,8)
  s_tmp = s_tmp + la_tmp

  ! DIIS vectors
  la_tmp = Int(n_a*CC_n_vir,8) * Int(n_ij*2*CC_DIIS_n_sv,8)
  s_tmp = s_tmp + la_tmp
  if (la_tmp.gt.la_mem) then
    la_mem = la_tmp
    la_dis = 'DIIS vectors'
  end if

  ! temporary tensors
  ! intl(ia|jb)
  mt1 = Int(CC_mem_bas(CC_mpi_did+1),8) &
      * Int(CC_n_occ * CC_n_vir * 2,8)
  mt1 = mt1 + Int(CC_n_occ * CC_n_vir,8) &
      * Int(CC_mem_ii_D(CC_mpi_did+1) * CC_mem_aa_G(CC_mpi_gid+1),8)
  mt = mt1

  ! intl(ac|kd)
  mt1 = Int(CC_mem_bas(CC_mpi_did+1) * (CC_n_occ + CC_mem_ii_D(CC_mpi_did+1) * CC_n_vir),8)
  mt1 = mt1 + Int(CC_n_occ * CC_n_vir,8) * Int(CC_mem_ii_D(CC_mpi_did+1) * 2,8)
  if (mt1.gt.mt) then
    mt = mt1
  end if

  ! intl(ac|bd)
  mt1 = Int(CC_mem_bas(CC_mpi_did+1) * (CC_n_vir + CC_n_vir))
  mt1 = mt1 + Int(CC_n_occ * CC_n_vir,8) * 2
  if (mt1.gt.mt) then
    mt = mt1
  end if

  ! intermediate h & g
  mt1 = Int(n_b * n_i,8) * Int(CC_n_occ * CC_n_vir * 4,8)
  if (mt1.gt.mt) then
    mt = mt1
  end if

  ! intermediate j
  mt1 = Int(n_a * n_b,8) * Int(CC_n_occ * CC_n_occ * 4,8)
  mt1 = mt1 + Int(n_b * CC_n_vir,8) * Int(CC_n_occ * CC_n_occ,8)
  mt1 = mt1 + Int(CC_n_occ,8) * Int(n_b**2,8)
  if (mt1.gt.mt) then
    mt = mt1
  end if

  ! intermediate b
  mt1 = Int(n_b * CC_n_vir,8) * Int(CC_n_bas,8)
  mt1 = mt1 + Int(n_b * CC_n_vir,8) * Int(CC_n_occ * (CC_n_occ+1) /2,8)
  mt1 = mt1 + Int(n_b * n_b,8) * Int(nth,8)
  if (mt1.gt.mt) then
    mt = mt1
  end if

  s_tmp = s_tmp + mt
  if (mt.gt.la_mem) then
    la_mem = mt
    la_dis = 'temporary vectors'
  end if

  ! Calculate memory space demand (MB)
  CC_mem_demand = dble(s_tmp) * 8.0D0 / 1024.0D0 / 1024.0D0 / 1024.0D0
  large_mem = dble(la_mem) * 8.0D0 / 1024.0D0 / 1024.0D0 / 1024.0D0

  End Subroutine CC_cl_mem_distribution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_ini_guess()

! This subroutine evaluate MP1 wave function as the initial guess of CC calculations,
! and calculates MP2 energy

  Use CC_cl

  Implicit None
  Integer :: iii,jjj,aaa,bbb,code_ab,ab_run,ab_start,n_ab,i_start,n_i,i_run
  Integer :: n_a,a_start,a_end,a_run,n_c,c_start,c_end
  Integer :: errnum
  Double precision :: CC_E_diff_0_s,et1,et2,rlt1,rlt2,rlt

  ! Clear w_vector
  Call CC_cl_clean_w()
  CC_t_s = 0.0D0

  ! Initialize w_d
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

  Call CC_cl_add_tc2w()
  Call CC_cl_Jacob()

  !n_a = CC_mem_aa_D(CC_mpi_did+1)
  !a_start = CC_index_aa_D(CC_mpi_did+1)

  !!$OMP PARALLEL Default(Shared) &
  !!$OMP Private(a_run,iii,jjj,aaa,bbb, &
  !!$OMP et1,et2,CC_E_diff_0_s)
  !!$OMP DO
  !do bbb = 1, CC_n_vir
  !  do a_run = 1, n_a
  !    do jjj = 1, CC_n_occ
  !      do iii = 1, CC_n_occ

  !        aaa = a_start - 1 + a_run

  !        ! Calculate the energy difference between reference(HF) state and i_state state.
  !        ! CC_E_diff_0_s = E_0 - E_s = e_i + e_j - e_a - e_b + CC_lv_shift
  !        Call CC_cl_ev_diff(iii,aaa,et1)
  !        Call CC_cl_ev_diff(jjj,bbb,et2)
  !        CC_E_diff_0_s = - et1 - et2 + CC_lv_shift

  !        ! Calculate coefficients
  !        CC_t_d(iii,jjj,a_run,bbb) = CC_intl_iajb(iii,jjj,a_run,bbb) / CC_E_diff_0_s

  !      end do
  !    end do
  !  end do
  !end do
  !!$OMP END DO
  !!$OMP END PARALLEL

  End Subroutine CC_cl_ini_guess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_cl_RI_p2c()

  Use basis
  Use prodbas
  Use hartree_fock
  Use CC_cl
  Use mpi_tasks
  use geometry, only: species
  use pbc_lists, only: Cbasis_to_atom
  Implicit None

  Integer :: i_k_point,id_root,i_k_point_local,k1,k2,k3,k4,k_12,k_34,km1
  Integer :: i_k_pattern,i_k_start,i_run,i_go,i_RI,n_RI,i_k1,i_k2
  Integer :: errnum

  Double complex , dimension(:,:,:,:) , allocatable :: lvl_recip1_tmp, &
                                                       coulomb_matr_recip_tmp, &
                                                       KS_egnvec_tmp, &
                                                       lvl_comm_tmp, coulomb_comm_tmp, &
                                                       KS_egnvec_comm_tmp, &
                                                       RI_tmp_R, RI_tmp_L

  Double complex , dimension(:,:) , allocatable :: lvl_tricoeff_sum_all,lvl_nn_tmp, &
                                                   lvl_nn_tmp_2,lvl_nn_tmp_3,mul_tmp


  Integer (kind = 8) :: s_tmp
  Integer :: i_tmp
  Integer :: n_row,n_col
  Integer :: i_row,i_col
  Integer :: i_basis_1,i_basis_2,i_atom_1,i_atom_2,i_species_1,i_species_2
  Integer :: bboff_1,bboff_2,n_spbb_1,n_spbb_2

  Allocate(lvl_recip1_tmp(max_n_basbas_sp,n_basis,n_basis,n_k_points),stat=errnum)
  Call check_allocation(errnum,'lvl_recip1_tmp in CC')

  lvl_recip1_tmp = 0.0D0

  Allocate(lvl_comm_tmp(max_n_basbas_sp,n_basis,n_basis,1),stat=errnum)
  Call check_allocation(errnum,'lvl_comm_tmp in CC')

  lvl_comm_tmp = 0.0D0

  if (myid .eq. mod(1, n_tasks)) then
    call power_auxmat_lapack_complex(coulomb_matr_recip(:,:,1),0.5d0, '')
  endif 

  Allocate(coulomb_matr_recip_tmp(CC_mem_bas(CC_mpi_did+1),n_basbas, &
                                  n_k_points,1),stat=errnum)
  Call check_allocation(errnum,'coulomb_matr_recip_tmp in CC')

  coulomb_matr_recip_tmp = 0.0D0

  Allocate(coulomb_comm_tmp(n_basbas,n_basbas,1,1),stat=errnum)
  Call check_allocation(errnum,'coulomb_comm_tmp in CC')

  coulomb_comm_tmp = 0.0D0

  Allocate(KS_egnvec_tmp(n_basis,n_states,n_k_points,1),stat=errnum)
  Call check_allocation(errnum,'KS_egnvec_tmp in CC')

  KS_egnvec_tmp = 0.0D0

  Allocate(KS_egnvec_comm_tmp(n_basis,n_states,1,1),stat=errnum)
  Call check_allocation(errnum,'KS_egnvec_comm_tmp in CC')
  
  KS_egnvec_comm_tmp = 0.0D0

  if (flag_KS_eigenfunc_conjg) then
    if (allocated(KS_eigenvector_complex)) then
      KS_eigenvector_complex = conjg(KS_eigenvector_complex)
    end if
  end if

  do i_k_point = 1, n_k_points

    id_root = Mod(i_k_point,n_tasks)
    i_k_point_local = Int((i_k_point - 1) / n_tasks) + 1

    if (myid.eq.id_root) then

      lvl_comm_tmp(:,:,:,1) = lvl_tricoeff_recip1(:,:,:,i_k_point_local)
      coulomb_comm_tmp(:,:,1,1) = coulomb_matr_recip(:,:,i_k_point_local)

      if (real_eigenvectors) then
        KS_egnvec_comm_tmp(:,:,1,1) = KS_eigenvector(:,:,1,i_k_point_local)
      else
        KS_egnvec_comm_tmp(:,:,1,1) = KS_eigenvector_complex(:,:,1,i_k_point_local)
      end if
    end if

    s_tmp = Int(size(lvl_comm_tmp),8)
    Call CC_mpi_complex_bcast(s_tmp, lvl_comm_tmp, id_root, MPI_COMM_WORLD)

    s_tmp = Int(size(coulomb_comm_tmp),8)
    Call CC_mpi_complex_bcast(s_tmp, coulomb_comm_tmp, id_root, MPI_COMM_WORLD)

    s_tmp = Int(size(KS_egnvec_comm_tmp),8)
    Call CC_mpi_complex_bcast(s_tmp, KS_egnvec_comm_tmp, id_root, MPI_COMM_WORLD)
        
    lvl_recip1_tmp(:,:,:,i_k_point) = lvl_comm_tmp(:,:,:,1)

    i_RI = CC_index_bas(CC_mpi_did+1) - 1
    n_RI = CC_mem_bas(CC_mpi_did+1)
    coulomb_matr_recip_tmp(:,:,i_k_point,1) = coulomb_comm_tmp(i_RI+1:i_RI+n_RI,:,1,1)

    KS_egnvec_tmp(:,:,i_k_point,1) = KS_egnvec_comm_tmp(:,:,1,1)

  end do
  
  Deallocate(lvl_comm_tmp)
  Deallocate(coulomb_comm_tmp)

  Deallocate(KS_egnvec_comm_tmp)

  Allocate(RI_tmp_R(CC_mem_bas(CC_mpi_did+1),n_basis,n_basis,1),stat=errnum)
  Call check_allocation(errnum,'RI_tmp_R in CC')

  Allocate(RI_tmp_L(CC_mem_bas(CC_mpi_did+1),n_basis,n_basis,1),stat=errnum)
  Call check_allocation(errnum,'RI_tmp_L in CC')

  Allocate(mul_tmp(CC_mem_bas(CC_mpi_did+1),1),stat=errnum)
  Call check_allocation(errnum,'mul_tmp in CC')

  Allocate(lvl_nn_tmp(n_basis,n_basis),stat=errnum)
  Call check_allocation(errnum,'lvl_nn_tmp in CC')

  Allocate(lvl_nn_tmp_2(n_basis,n_states),stat=errnum)
  Call check_allocation(errnum,'lvl_nn_tmp in CC')

  Allocate(lvl_nn_tmp_3(n_states,n_states),stat=errnum)
  Call check_allocation(errnum,'lvl_nn_tmp in CC')

  Allocate(lvl_tricoeff_sum_all(n_basbas,1),stat=errnum)
  Call check_allocation(errnum,'lvl_tricoeff_sum_all in CC')

 
  k1 = 1
  k2 = 1  
  k_12 = kq_point_list(k2,k1)
  km1 = kq_point_list(1,k1)
  
  !s_tmp = Int(n_basis,8) * Int(n_basis,8)
  
  do i_tmp = 1, n_basis * n_basis
  
    Call CC_cl_decode(i_tmp,i_basis_1,i_basis_2,n_basis,1)
  
    i_atom_1 = Cbasis_to_atom(i_basis_1)
    i_species_1 = species(i_atom_1)
    bboff_1 = atom2basbas_off(i_atom_1)
    n_spbb_1 = sp2n_basbas_sp(i_species_1)
  
    i_atom_2 = Cbasis_to_atom(i_basis_2)
    i_species_2 = species(i_atom_2)
    bboff_2 = atom2basbas_off(i_atom_2)
    n_spbb_2 = sp2n_basbas_sp(i_species_2)
  
    ! Calculate CC_RI_L
    lvl_tricoeff_sum_all = 0.0D0
  
    lvl_tricoeff_sum_all(bboff_1+1:bboff_1+n_spbb_1,1) = &
      lvl_tricoeff_sum_all(bboff_1+1:bboff_1+n_spbb_1,1) + &
      conjg(lvl_recip1_tmp(1:n_spbb_1, i_basis_1, i_basis_2,k1)) + &
      lvl_tricoeff_recip2(1:n_spbb_1, i_basis_2, i_basis_1)
  
    lvl_tricoeff_sum_all(bboff_2+1:bboff_2+n_spbb_2,1) = &
      lvl_tricoeff_sum_all(bboff_2+1:bboff_2+n_spbb_2,1) + &
      lvl_recip1_tmp(1:n_spbb_2, i_basis_2, i_basis_1,k2) + &
      lvl_tricoeff_recip2(1:n_spbb_2, i_basis_1, i_basis_2)
 
 
    i_RI = CC_index_bas(CC_mpi_did+1) - 1
  
    !do i_go = 1, CC_mem_bas(CC_mpi_did+1)
    !  i_RI = i_RI + 1
    !  RI_tmp_L(i_go,i_basis_2,i_basis_1,1) = &
    !           lvl_tricoeff_sum_all(i_RI,1)
    !end do
 
    ! Calculate CC_RI_R
    mul_tmp = MatMul(coulomb_matr_recip_tmp(:,:,k_12,1),lvl_tricoeff_sum_all)
    RI_tmp_R(:,i_basis_2,i_basis_1,1) = mul_tmp(:,1)
  
  end do

  do i_go = 1, CC_mem_bas(CC_mpi_did+1)
  
    !lvl_nn_tmp = RI_tmp_L(i_go,:,:,1)
    !Call zgemm('N', 'N', n_basis, n_states, n_basis, (1.d0,0.d0), &
    !           lvl_nn_tmp, n_basis, KS_egnvec_tmp(:,:,k2,1), n_basis, &
    !           (0.d0,0.d0), lvl_nn_tmp_2, n_basis)
  
    !Call zgemm('C', 'N', n_states, n_states, n_basis, (1.d0,0.d0), &
    !           KS_egnvec_tmp(:,:,k1,1), n_basis, lvl_nn_tmp_2, n_basis, &
    !           (0.d0,0.d0), lvl_nn_tmp_3, n_states)
 
    !CC_RI_L(i_go,:,:,k1,k2) = lvl_nn_tmp_3(CC_n_fc+1:n_states,CC_n_fc+1:n_states)
  
    lvl_nn_tmp = RI_tmp_R(i_go,:,:,1)
    Call zgemm('N', 'N', n_basis, n_states, n_basis, (1.d0,0.d0), &
               lvl_nn_tmp, n_basis, KS_egnvec_tmp(:,:,k2,1), n_basis, &
               (0.d0,0.d0), lvl_nn_tmp_2, n_basis)
  
    Call zgemm('C', 'N', n_states, n_states, n_basis, (1.d0,0.d0), &
               KS_egnvec_tmp(:,:,k1,1), n_basis, lvl_nn_tmp_2, n_basis, &
               (0.d0,0.d0), lvl_nn_tmp_3, n_states)
  
    CC_RI(i_go,:,:,1) = dble(lvl_nn_tmp_3(CC_n_fc+1:n_states,CC_n_fc+1:n_states))
  
  end do

  Deallocate(lvl_tricoeff_sum_all)
  Deallocate(lvl_nn_tmp,lvl_nn_tmp_2,lvl_nn_tmp_3)
  Deallocate(lvl_recip1_tmp,coulomb_matr_recip_tmp,KS_egnvec_tmp)

!  print*,'RI_L'
!  print*,CC_RI_L
!  print*,'RI_R'
!  print*,CC_RI_R



  End Subroutine CC_cl_RI_p2c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

