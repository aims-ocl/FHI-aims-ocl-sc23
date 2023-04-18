  Subroutine CC_3d_initial()

  Use basis
  Use prodbas
  Use hartree_fock
  Use mpi_tasks
  Use CC_3d
  use localorb_io, only: use_unit

  Implicit None

  Integer :: i_k_state,k1,k2,k3,k4,k_12,k_34

  Integer (kind = 8) :: i_tmp,s_tmp,largest_ten
  Integer :: j_tmp,k_tmp,iii,jjj,mpierr
  Integer :: i_state,errnum,i_task,hdferr
  Double precision , Parameter :: ne_thresh = 1.0D-6
  Double precision :: CC_mem_demand,CC_disk_demand
  Double complex , dimension(:) , allocatable :: run_tmp
  Double complex :: rlt

  Character(len=20) :: lt_outline

! Initialization starts
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") '--------------------------------------------'
    write(use_unit,"(2x,A)") 'Coupled Cluster(CC) calculation starts'
    write(use_unit,"(2x,A)") '--------------------------------------------'
    write(use_unit,"(A)") ' '
    write(use_unit,"(2x,A)") 'For periodic closed shell system'
    write(use_unit,"(2x,A)") 'Initialization of CC calculation'
    if (CC_use_disk) then
      if (CC_restart_point.eq.0) then
        write(use_unit,"(2x,A)") 'Restart information will be saved'
      else
        write(use_unit,"(2x,A)") 'Restarting unfinished CC calculation'
      end if
    else
      write(use_unit,"(2x,A)") 'Restart information will not be saved'
    end if
  end if

! Determine the reference (HF) configuration and number of electrons

!  if (CC_restart_point.ne.0) then
!    Call CC_3d_read_restart_info()   
!  else
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
    else
      CC_valence = 1
      CC_n_fc = 0
    endif
  
    CC_n_occ = CC_n_elec - CC_n_fc
    CC_n_vir = n_states - CC_n_elec
    CC_n_state = n_states - CC_n_fc
    CC_n_bas = n_basbas
  
!  end if

  CC_n_s = CC_n_occ * CC_n_vir * n_k_points

  j_tmp = 0
  k_tmp = 0
  do k1 = 1, n_k_points
    do k2 = 1, n_k_points
      do k3 = 1, n_k_points

        Call CC_3d_determine_k_pattern(i_k_state,k1,k2,k3,k4,4)
        Call CC_3d_code(k_12,k1,k2,n_k_points,1)
        Call CC_3d_code(k_34,k3,k4,n_k_points,1)

        if (k_12.lt.k_34) then
          j_tmp = j_tmp + 1
        else if (k_12.eq.k_34) then
          k_tmp = k_tmp + 1
        end if

      end do
    end do
  end do

  CC_n_d = Int(j_tmp,8) * Int(CC_n_occ**2,8) * Int(CC_n_vir**2,8)
  j_tmp = CC_n_occ * CC_n_vir
  CC_n_d = CC_n_d + Int(k_tmp,8) * Int(j_tmp,8) * (Int(j_tmp,8) - 1) / 2

  CC_n_config = Int(CC_n_s,8) + CC_n_d

  ! Determine saving strategy
  CC_sv_strategy = 1
  if ((CC_sv_control.eq.'GLOBAL').or.(CC_sv_control.eq.'global')) then
    CC_sv_strategy = 2
  end if

  ! Create MPI domains and groups
  Call CC_3d_mpi_init()

  if (myid.eq.0) then
    write(use_unit,"(A)") ' '
    write(use_unit,"(2x,A)") 'MPI grouping...'
    write(use_unit,"(2x,A,I3)") 'Number of domains: ',CC_n_domain
    write(use_unit,"(2x,A,I3)") 'Number of tasks in each domain: ',CC_mpi_domain_size
    write(use_unit,"(2x,A)") 'Acceleration method: DIIS'
    write(use_unit,"(2x,A,I6)") 'Number DIIS vectors saved:',CC_DIIS_n_sv
    write(use_unit,"(2x,A,A)") 'Storage strategy: ',CC_sv_control
    if (CC_abcd_sv_flag) then
      write(use_unit,"(2x,A,A)") 'Integral tensor (ac|bd) will be saved in memory.'
    end if
    write(use_unit,"(2x,A,A)") ''
    write(use_unit,"(2x,A,I7)") 'Nuumber of valence orbitals: ',CC_n_occ
    write(use_unit,"(2x,A,I7)") 'Nuumber of unoccupied orbitals: ',CC_n_vir
    write(use_unit,"(2x,A)") 'Number of configurations'
    write(use_unit,"(10x,A16,8x,A16)") 's-excitation','d-excitation'
    write(use_unit,"(10x,I16,8x,I16)") CC_n_s,CC_n_d
    write(use_unit,"(2x,A,I16,2x,A)") 'Totally ',CC_n_config,'configurations are considered'
  end if

  ! Determine the distributiuon of vectors on each task
  ! For two-index distribution
  ! For ij
  Allocate(CC_mem_ij_D(CC_mpi_domain_size,2), stat=errnum)
  Call check_allocation(errnum,'CC_mem_ij_D in CC')

  Allocate(CC_index_ij_D(CC_mpi_domain_size,2), stat=errnum)
  Call check_allocation(errnum,'CC_index_ij_D in CC')

  j_tmp = CC_n_occ * CC_n_occ
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_ij_D(:,1),CC_index_ij_D(:,1))

  j_tmp = CC_n_occ * (CC_n_occ + 1) / 2
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_ij_D(:,2),CC_index_ij_D(:,2))

  ! For ab
  Allocate(CC_mem_ab_D(CC_mpi_domain_size,4), stat=errnum)
  Call check_allocation(errnum,'CC_mem_ab_D in CC')

  Allocate(CC_index_ab_D(CC_mpi_domain_size,4), stat=errnum)
  Call check_allocation(errnum,'CC_index_ab_D in CC')

  j_tmp = CC_n_vir * CC_n_vir
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_ab_D(:,1),CC_index_ab_D(:,1))

  j_tmp = CC_n_vir * (CC_n_vir + 1) / 2
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_ab_D(:,2),CC_index_ab_D(:,2))

  j_tmp = CC_n_vir * (CC_n_vir - 1) / 2
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_ab_D(:,3),CC_index_ab_D(:,3))

  j_tmp = CC_n_vir
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_ab_D(:,4),CC_index_ab_D(:,4))


  ! For ia
  Allocate(CC_mem_ia_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_ia_D in CC')

  Allocate(CC_index_ia_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_ia_D in CC')

  j_tmp = CC_n_occ * CC_n_vir
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_ia_D,CC_index_ia_D)

  ! For cd
  Allocate(CC_mem_cd_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_cd_D in CC')

  Allocate(CC_index_cd_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_cd_D in CC')

  j_tmp = CC_n_vir * CC_n_vir
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_cd_D,CC_index_cd_D)

  ! For kl
  Allocate(CC_mem_kl_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_kl_D in CC')

  Allocate(CC_index_kl_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_kl_D in CC')

  j_tmp = CC_n_occ * CC_n_occ
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_kl_D,CC_index_kl_D)

  ! For i
  Allocate(CC_mem_ii_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_ii_D in CC')

  Allocate(CC_index_ii_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_ii_D in CC')

  j_tmp = CC_n_occ 
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_ii_D,CC_index_ii_D)

  ! For a
  Allocate(CC_mem_aa_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_aa_D in CC')

  Allocate(CC_index_aa_D(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_aa_D in CC')

  j_tmp = CC_n_vir
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_domain_size,CC_mem_aa_D,CC_index_aa_D)

  ! For CC_n_bas
  Allocate(CC_mem_bas(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_bas in CC')

  Allocate(CC_index_bas(CC_mpi_domain_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_bas in CC')

  Call CC_3d_vec_distr(1,CC_n_bas,CC_mpi_domain_size,CC_mem_bas,CC_index_bas)

  ! Determine the distributiuon of k_points on each task
  Allocate(CC_mem_k1(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_k1 in CC')

  Allocate(CC_index_k1(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_k1 in CC')

  j_tmp = n_k_points
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_group_size,CC_mem_k1,CC_index_k1)

  ! Determine the distributiuon of k_points on each task
  Allocate(CC_mem_k2(CC_mpi_group_size,3), stat=errnum)
  Call check_allocation(errnum,'CC_mem_k2 in CC')

  Allocate(CC_index_k2(CC_mpi_group_size,3), stat=errnum)
  Call check_allocation(errnum,'CC_index_k2 in CC')

  j_tmp = n_k_points * n_k_points
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_group_size,CC_mem_k2(:,1),CC_index_k2(:,1))

  j_tmp = n_k_points * (n_k_points + 1) / 2
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_group_size,CC_mem_k2(:,2),CC_index_k2(:,2))

  j_tmp = n_k_points * (n_k_points - 1) / 2
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_group_size,CC_mem_k2(:,3),CC_index_k2(:,3))

  ! Determine the distributiuon of k_patterns on each task
  Allocate(CC_mem_k3(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_mem_k3 in CC')

  Allocate(CC_index_k3(CC_mpi_group_size), stat=errnum)
  Call check_allocation(errnum,'CC_index_k3 in CC')

  j_tmp = n_k_points ** 3
  Call CC_3d_vec_distr(1,j_tmp,CC_mpi_group_size,CC_mem_k3,CC_index_k3)

  ! Determine storage space
  s_tmp = Int(CC_MPI_max_len,8)
  largest_ten = Int(CC_MPI_max_len,8)
  lt_outline = 'MPI comm tmp'

  if (.not.(CC_abcd_sv_flag)) then
    ! Space for RI coefficients
    i_tmp = Int(CC_mem_bas(CC_mpi_did+1)*CC_n_state**2,8)
    j_tmp = n_k_points ** 2
    s_tmp = s_tmp + i_tmp * Int(j_tmp * 2,8)
    if (largest_ten.lt.i_tmp * Int(j_tmp * 2,8)) then
      largest_ten = i_tmp * Int(j_tmp * 2,8)
      lt_outline = 'RI coefficients'
    end if
  end if

  ! Space for run_time temporary vectors
  j_tmp = CC_mem_k1(CC_mpi_gid+1) * n_k_points ** 2
  i_tmp = Int(CC_mem_aa_D(CC_mpi_did+1) * CC_n_vir,8) &
        * Int(CC_n_occ ** 2,8) * Int(j_tmp*6,8) &
        + Int(4 * n_k_points * CC_mem_ab_D(CC_mpi_did+1,3)**2,8)
  s_tmp = s_tmp + i_tmp

  if (largest_ten.lt.i_tmp) then
    largest_ten = i_tmp
    lt_outline = 'Run-time tmp'
  end if

  ! Space for CC_w_d_d, CC_intl_iajb_d, CC_t_d_d, CC_w_d_nd, CC_intl_iajb_nd,
  ! CC_t_d_nd, DIIS_t and DIIS_r
  i_tmp = Int(CC_mem_aa_D(CC_mpi_did+1),8) * Int(CC_n_occ**2,8)
  s_tmp = s_tmp + i_tmp * Int(j_tmp * (3+CC_DIIS_n_sv*2),8)

  if (largest_ten.lt.i_tmp * Int(j_tmp * (3+CC_DIIS_n_sv*2),8)) then
    largest_ten = i_tmp * Int(j_tmp * (3+CC_DIIS_n_sv*2),8)
    lt_outline = 'DIIS storage tensor'
  end if
 
  i_tmp = Int(CC_mem_ab_D(CC_mpi_did+1,3),8) * Int(CC_n_occ**2,8)
  s_tmp = s_tmp + i_tmp * Int(j_tmp * (3+CC_DIIS_n_sv*2),8)

  ! Space for (ai|kc) and (ki|ac)
  i_tmp = Int(CC_n_occ**2,8) * Int(CC_n_vir,8) * Int(CC_mem_aa_D(CC_mpi_did+1),8)
  s_tmp = s_tmp + i_tmp * Int(j_tmp * 2,8)

  ! Space for (ac|kd) and (ai|bc)
  i_tmp = Int(CC_n_occ,8) * Int(CC_n_vir**2,8) * Int(CC_mem_aa_D(CC_mpi_did+1),8)
  s_tmp = s_tmp + i_tmp * Int(j_tmp * 2,8)

  if (largest_ten.lt.i_tmp * Int(j_tmp * 2,8)) then
    largest_ten = i_tmp * Int(j_tmp * 2,8)
    lt_outline = '(ac|kd) integral'
  end if

  ! Space for (ac|kd)_nd
  i_tmp = Int(CC_n_occ,8) * Int(CC_n_vir,8) * Int(CC_mem_ab_D(CC_mpi_did+1,3),8)
  s_tmp = s_tmp + i_tmp * Int(j_tmp,8)

  ! Space for (ac|kd)_d
  i_tmp = Int(CC_n_occ,8) * Int(CC_n_vir,8) * Int(CC_mem_aa_D(CC_mpi_did+1),8)
  s_tmp = s_tmp + i_tmp * Int(j_tmp,8)

  ! Space for (li|kc) and (ki|lc)
  i_tmp = Int(CC_n_occ**2,8) * Int(CC_n_vir,8) * Int(CC_mem_ii_D(CC_mpi_did+1),8)
  s_tmp = s_tmp + i_tmp * Int(j_tmp * 2,8)

  ! Space for (ki|lc)_T
  i_tmp = Int(CC_n_occ,8) * Int(CC_n_vir,8) * Int(CC_mem_ij_D(CC_mpi_did+1,2),8)
  s_tmp = s_tmp + i_tmp * Int(j_tmp * 2,8)

  ! Space for (ki|lj)
  i_tmp = Int(CC_n_occ**2,8) * Int(CC_mem_ij_D(CC_mpi_did+1,2),8)
  s_tmp = s_tmp + i_tmp * Int(j_tmp,8)

  ! Space for (ai|kj)
  i_tmp = Int(CC_n_occ**3,8) * Int(CC_mem_aa_D(CC_mpi_did+1),8)
  s_tmp = s_tmp + i_tmp * Int(j_tmp,8)

  if (CC_abcd_sv_flag) then
    ! Space for (ac|bd)
    i_tmp = Int(CC_mem_ab_D(CC_mpi_did+1,3),8) * Int(CC_n_vir**2,8)
    s_tmp = s_tmp + i_tmp * Int(j_tmp,8)

    if (largest_ten.lt.i_tmp * Int(j_tmp,8)) then
      largest_ten = i_tmp * Int(j_tmp,8)
      lt_outline = '(ac|bd) integral'
    end if

    i_tmp = Int(CC_mem_aa_D(CC_mpi_did+1),8) * Int(CC_n_vir**2,8)
    s_tmp = s_tmp + i_tmp * Int(j_tmp,8)

  end if

  ! Space for CC_aux_a
  i_tmp = Int(CC_n_occ**2,8) * Int(CC_mem_ij_D(CC_mpi_did+1,2),8)
  s_tmp = s_tmp + i_tmp * Int(j_tmp,8)

  ! Space for CC_aux_j , CC_aux_k and contraction intermediate
  i_tmp = Int(CC_n_occ**2,8) * Int(CC_n_vir,8) * Int(CC_mem_aa_D(CC_mpi_did+1),8)
  s_tmp = s_tmp + i_tmp * Int(j_tmp * 3,8)

  ! Space for CC_t_d and (kc|ld)
  if (CC_sv_strategy.eq.1) then
    j_tmp = n_k_points**3
  end if

  s_tmp = s_tmp + i_tmp * Int(j_tmp * 2,8)

  ! Calculate memory space demand (MB)
  CC_mem_demand = dble(s_tmp) * 16.0D0 / 1024.0D0 / 1024.0D0

  if (myid.eq.0) then
    write(use_unit,"(2x,A,F14.3,A)") 'At least ',CC_mem_demand, &
                                     'MB memory space is needed for each task'
    j_tmp = len_trim(lt_outline)
    CC_mem_demand = dble(largest_ten) * 16.0D0 / 1024.0D0 / 1024.0D0
    write(use_unit,"(2x,A,A,A,F14.3,A)") 'The largest tensor ',lt_outline(1:j_tmp),' needs ', &
                                         CC_mem_demand,'MB memory'
    if (CC_use_disk) then
      write(use_unit,"(2x,A,F14.6,A)") 'Totally ',CC_disk_demand,'MB disk space is needed' 
    end if
  end if

  if (CC_mem_check) then
    Call aims_stop('This is just a memory check for CC calculation!')
  end if

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Preparing four-index integrals...'
  end if
  ! Allocate all vectors
  ! Initialize integral vectors
  Call CC_3d_init_intl()

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Four-index integrals are ready.'
    write(use_unit,"(2x,A)") 'Allocating CC tensors...'
  end if

  if (CC_abcd_sv_flag) then
    Deallocate(CC_RI_L,CC_RI_R)
  end if

  ! Allocate run-time temporary vector
  i_tmp = Int(CC_mem_aa_D(CC_mpi_did+1) * CC_n_vir,8) &
        * Int(CC_n_occ ** 2,8) * Int(n_k_points ** 2,8) &
        * Int(CC_mem_k1(CC_mpi_gid+1) * 6,8) &
        + Int(4 * n_k_points * CC_mem_ab_D(CC_mpi_did+1,3)**2,8)

  Allocate(run_tmp(i_tmp), stat=errnum)
  Call check_allocation(errnum,'run-time tmp in CC')

  run_tmp = 0.0D0

  ! DIIS R matrix
  Allocate(CC_DIIS_Bmat(CC_DIIS_n_sv,CC_DIIS_n_sv),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_Bmat in CC')

  ! DIIS tau vector
  Allocate(CC_DIIS_tau(CC_DIIS_n_sv),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_tau in CC')

  ! DIIS storage
  Allocate(CC_DIIS_t_s(CC_n_occ,CC_n_vir,n_k_points,CC_DIIS_n_sv),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_t_s in CC')

  Allocate(CC_DIIS_r_s(CC_n_occ,CC_n_vir,n_k_points,CC_DIIS_n_sv),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_r_s in CC')

  Allocate(CC_DIIS_t_d(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                       n_k_points,CC_mem_k1(CC_mpi_gid+1),CC_DIIS_n_sv),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_t_d in CC')

  Allocate(CC_DIIS_r_d(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                       n_k_points,CC_mem_k1(CC_mpi_gid+1),CC_DIIS_n_sv),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_r_d in CC')

  Allocate(CC_DIIS_t_nd(CC_n_occ,CC_n_occ,CC_mem_ab_D(CC_mpi_did+1,3),n_k_points, &
                        n_k_points,CC_mem_k1(CC_mpi_gid+1),CC_DIIS_n_sv),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_t_nd in CC')

  Allocate(CC_DIIS_r_nd(CC_n_occ,CC_n_occ,CC_mem_ab_D(CC_mpi_did+1,3),n_k_points, &
                        n_k_points,CC_mem_k1(CC_mpi_gid+1),CC_DIIS_n_sv),stat=errnum)
  Call check_allocation(errnum,'CC_DIIS_r_nd in CC')

! Allocate intermediate vectors
  ! For CC_h_ik(i,k,ki)
  Allocate(CC_h_ik(CC_n_occ,CC_n_occ,n_k_points),stat=errnum)
  Call check_allocation(errnum,'CC_h_ik in CC')

  CC_h_ik = 0.0D0

  ! For CC_g_ik(i,k,ki)
  Allocate(CC_g_ik(CC_n_occ,CC_n_occ,n_k_points),stat=errnum)
  Call check_allocation(errnum,'CC_g_ik in CC')

  CC_g_ik = 0.0D0

  ! For CC_h_ca(c,a,ka)
  Allocate(CC_h_ca(CC_n_vir,CC_n_vir,n_k_points),stat=errnum)
  Call check_allocation(errnum,'CC_h_ca in CC')

  CC_h_ca = 0.0D0

  ! For CC_g_ca(c,a,ka)
  Allocate(CC_g_ca(CC_n_vir,CC_n_vir,n_k_points),stat=errnum)
  Call check_allocation(errnum,'CC_g_ca in CC')

  CC_g_ca = 0.0D0

  ! For CC_h_ck(k,c,kk)
  Allocate(CC_h_ck(CC_n_occ,CC_n_vir,n_k_points),stat=errnum)
  Call check_allocation(errnum,'CC_h_ck in CC')

  CC_h_ck = 0.0D0

  ! For aux_a (k<=l,kk,i,j,ki,kk+kl)
  Allocate(CC_a_aux(CC_mem_ij_D(CC_mpi_did+1,2),n_k_points, &
                    CC_n_occ,CC_n_occ,n_k_points, &
                    CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_a_aux in CC')

  CC_a_aux = 0.0D0

  ! For aux_j (k,c,kk,i,a,ki,ka)
  Allocate(CC_j_aux(CC_n_occ,CC_n_vir,n_k_points, &
                    CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                    CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_j_aux in CC')

  CC_j_aux = 0.0D0

  ! For aux_k (k,c,kk,i,a,ki,ka)
  Allocate(CC_k_aux(CC_n_occ,CC_n_vir,n_k_points, &
                    CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                    CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_k_aux in CC')

  CC_k_aux = 0.0D0

  ! For contraction intermediate (j,b,kj,i,a,ki,ka)
  Allocate(CC_ten_con(CC_n_occ,CC_n_vir,n_k_points, &
                      CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                      CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_ten_con in CC')

  CC_ten_con = 0.0D0

  ! Allocate working tensors
  ! For w_s (i,a,ki=ka)
  Allocate(CC_w_s(CC_n_occ,CC_n_vir,n_k_points),stat=errnum)
  Call check_allocation(errnum,'CC_w_s in CC')

  CC_w_s = 0.0D0

  ! For w_d_nd (i,j,a<b,ki,ka,ka+kb)
  Allocate(CC_w_d_nd(CC_n_occ,CC_n_occ,CC_mem_ab_D(CC_mpi_did+1,3), &
                     n_k_points,n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_w_d_nd in CC')

  CC_w_d_nd = 0.0D0

  ! For w_d_d (i,j,a=b,ki,ka,ka+kb)
  Allocate(CC_w_d_d(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1), &
                    n_k_points,n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_w_d_d in CC')

  CC_w_d_d = 0.0D0

  ! Allocate coefficient tensor
  Allocate(CC_t_s(CC_n_occ,CC_n_vir,n_k_points),stat=errnum)
  Call check_allocation(errnum,'CC_t_s in CC')

  CC_t_s = 0.0D0

  ! For CC_t_d (i,j,a<b,ki,ka,ka+kb)
  Allocate(CC_t_d_nd(CC_n_occ,CC_n_occ,CC_mem_ab_D(CC_mpi_did+1,3), &
                     n_k_points,n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_t_d_nd in CC')

  CC_t_d_nd = 0.0D0

  ! For CC_t_d (i,j,a=b,ki,ka,ka+kb)
  Allocate(CC_t_d_d(CC_n_occ,CC_n_occ,CC_mem_aa_D(CC_mpi_did+1), &
                    n_k_points,n_k_points,CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
  Call check_allocation(errnum,'CC_t_d_d in CC')

  CC_t_d_d = 0.0D0

  if (CC_sv_strategy.eq.1) then

    ! CC_t_d_A (j,b,kj,i,a,ki,ka)
    Allocate(CC_t_d_A(CC_n_occ,CC_n_vir,n_k_points, &
                      CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                      n_k_points),stat=errnum)
    Call check_allocation(errnum,'CC_t_d_A in CC')
  
    CC_t_d_A = 0.0D0

  else if (CC_sv_strategy.eq.2) then

    ! CC_t_d (j,b,kj,i,a,ki,ka)
    ! (j,b,kj,i,a,ki,ka)
    Allocate(CC_t_d(CC_n_occ,CC_n_vir,n_k_points, &
                    CC_n_occ,CC_mem_aa_D(CC_mpi_did+1),n_k_points, &
                    CC_mem_k1(CC_mpi_gid+1)),stat=errnum)
    Call check_allocation(errnum,'CC_t_d in CC')
  
    CC_t_d = 0.0D0

  end if

  Deallocate(run_tmp)

  if (CC_restart_point.eq.0) then
  ! Set up reference (HF) energy
    E_HF = total_energy
  end if

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'CC tensors allocated sucessfully.'
    write(use_unit,"(2x,A)") 'End of initialization'
  end if

  End Subroutine CC_3d_initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_3d_ini_guess()

! This subroutine evaluate MP1 wave function as the initial guess of CC calculations,
! and calculates MP2 energy

  Use CC_3d

  Implicit None

  Integer :: k_pattern,k1,k2,k3,k4,k_run,n_k1,k_start,k_end,k_13
  Integer :: errnum

  Integer :: n_dom
  Integer :: iii,jjj,aaa,bbb,a_run,b_run,i_run,code_ab
  Double precision :: CC_E_diff_0_s,et1,et2
  Double precision , dimension(:,:,:,:,:,:,:) , allocatable :: t_tmp

  Integer :: i_tmp,j_tmp,a_start,a_end,b_start,b_end

  ! Clear w_vector
  Call CC_3d_clean_w()
  CC_t_s = 0.0D0

  ! Initialize CC_t_d
  ! for w_d_nd (i,j,a<b,ki,ka,ka+kb) k1 = ka, k2 = ki, k3 = kb, k4 = kj
  n_k1 = CC_mem_k1(CC_mpi_gid+1)
  k_start = CC_index_k1(CC_mpi_gid+1)

  k_13 = k_start - 1
  do k_run = 1, n_k1

    k_13 = k_13 + 1

    do k1 = 1, n_k_points

      Call CC_3d_k_minus(k_13,k1,k3)

      do k2 = 1, n_k_points

        Call CC_3d_k_minus(k_13,k2,k4)

        ! for w_d_nd (i,j,a<b,ki,ka,ka+kb) k1 = ka, k2 = ki, k3 = kb, k4 = kj
        !$OMP PARALLEL Default(Shared) &
        !$OMP Private(code_ab,a_run,aaa,bbb,iii,jjj,et1,et2,CC_E_diff_0_s)
        code_ab = CC_index_ab_D(CC_mpi_did+1,3) - 1

        do a_run = 1, CC_mem_ab_D(CC_mpi_did+1,3)

          code_ab = code_ab + 1
          Call CC_3d_decode(code_ab,aaa,bbb,CC_n_vir,3)

          !$OMP DO COLLAPSE(2)
          do iii = 1, CC_n_occ
            do jjj = 1, CC_n_occ

              ! Calculate the energy difference between reference(HF) state and i_state state.
              ! CC_E_diff_0_s = E_0 - E_s = e_i + e_j - e_a - e_b 
              Call CC_3d_ev_diff(iii,aaa,k2,k1,et1)
              Call CC_3d_ev_diff(jjj,bbb,k4,k3,et2)
              CC_E_diff_0_s = - et1 - et2

              CC_w_d_nd(iii,jjj,a_run,k2,k1,k_run) = &
                    conjg(CC_intl_iajb_nd(iii,jjj,a_run,k2,k1,k_run)) / CC_E_diff_0_s
            end do
          end do
          !$OMP END DO
        end do

        ! for w_d_d (i,j,a=b,ki,ka,ka+kb) k1 = ka, k2 = ki, k3 = kb, k4 = kj
        aaa = CC_index_aa_D(CC_mpi_did+1) - 1
        do a_run = 1, CC_mem_aa_D(CC_mpi_did+1)

          aaa = aaa + 1
          bbb = aaa

          !$OMP DO COLLAPSE(2)
          do iii = 1, CC_n_occ
            do jjj = 1, CC_n_occ

              ! Calculate the energy difference between reference(HF) state and i_state state.
              ! CC_E_diff_0_s = E_0 - E_s = e_i + e_j - e_a - e_b 
              Call CC_3d_ev_diff(iii,aaa,k2,k1,et1)
              Call CC_3d_ev_diff(jjj,bbb,k4,k3,et2)
              CC_E_diff_0_s = - et1 - et2

              CC_w_d_d(iii,jjj,a_run,k2,k1,k_run) = &
                    conjg(CC_intl_iajb_d(iii,jjj,a_run,k2,k1,k_run)) / CC_E_diff_0_s
            end do
          end do
          !$OMP END DO
        end do
        !$OMP END PARALLEL

      end do
    end do
  end do


  End Subroutine CC_3d_ini_guess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

