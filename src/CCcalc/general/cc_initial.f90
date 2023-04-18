  Subroutine CC_initialization()

  Use dimensions
  Use runtime_choices
  Use physics
  Use hartree_fock
  Use basis
  Use prodbas
  Use mpi_tasks
  Use CC_corr

  Implicit None

  Integer :: i_spin,i_state,errnum,i_tmp,j_tmp,k_tmp,i_task
  Double precision , Parameter :: ne_thresh = 1.0D-6
  Double precision :: CC_mem_demand

! Initialization starts
  if (myid.eq.0) then
    write(use_unit,"(2x,A)") '--------------------------------------------'
    write(use_unit,"(2x,A)") 'Coupled Cluster(CC) calculation starts'
    write(use_unit,"(2x,A)") '--------------------------------------------'
    write(use_unit,"(A)") ' '
    write(use_unit,"(2x,A)") 'Initialization of CC calculation'
    write(use_unit,"(2x,A)") 'Preparing four-index integrals...'
  end if

! Determine the reference (HF) configuration and numbers of alpha and beta electrons
  CC_n_elec = 0

  Allocate(CC_config_ref(n_states,2),stat=errnum)
  Call check_allocation(errnum,'CC_config_ref in CC calculation')

  CC_config_ref = 0

  if (n_spin.eq.1) then
    do i_state = 1, n_states
      if (dabs(occ_numbers(i_state,1,1)-1.0D0).le.ne_thresh) then
        CC_config_ref(i_state,1) = 1
        CC_n_elec(1) = CC_n_elec(1) + 1
      else if (dabs(occ_numbers(i_state,1,1)-2.0D0).le.ne_thresh) then
        CC_config_ref(i_state,1) = 1
        CC_config_ref(i_state,2) = 1
        CC_n_elec(1) = CC_n_elec(1) + 1
        CC_n_elec(2) = CC_n_elec(2) + 1
      end if
    end do
  else
    do i_spin = 1, 2
      do i_state = 1, n_states
        if (dabs(occ_numbers(i_state,i_spin,1)-1.0D0).le.ne_thresh) then
          CC_config_ref(i_state,i_spin) = 1
          CC_n_elec(i_spin) = CC_n_elec(i_spin) + 1
        end if
      end do
    end do
  end if

  ! Calculate the number of valence state
  if (flag_frozen_core_postSCF) then
    call count_frozen_core_states(CC_valence)
    if(myid.eq.0) then
       write(use_unit,'(2X,A,I8)') "| Frozen core number      ::", &
           CC_valence - 1
    endif
  else
    CC_valence = 1
  endif

  do i_spin = 1, 2
    CC_n_occ(i_spin) = CC_n_elec(i_spin) - CC_valence + 1
    CC_n_vir(i_spin) = n_states - CC_n_elec(i_spin)
  end do

! get three-center integrals (P|ij)
  Allocate(ovlp_3ks(n_basbas,n_states,n_states,n_spin), stat=errnum)
  Call check_allocation(errnum,'ovlp_3ks (1D) in CC')

  Call CC_get_ovlp3ks()

  ! Extract unnecessary memory storage
  Deallocate(ovlp_3fn)
  Deallocate(map_prodbas)

  ! Deallocate the memory about the taburated grids.
  Call deallocate_grid_storage()

  ! Deallocate the memory about the auxiliary basis sets
  Call cleanup_basbas()

  ! Deallocate the momory related to hartree_fock.mod
  Call cleanup_hartree_fock()

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Four-index integrals are ready'
    write(use_unit,"(2x,A)") 'Calculating memory demanding...'
  end if

! Calculate numbers of configurations
!  CC_n_vir(1) = n_states - CC_n_elec(1)
!  CC_n_vir(2) = n_states - CC_n_elec(2)

  CC_n_1a = CC_n_occ(1) * CC_n_vir(1)
  CC_n_1b = CC_n_occ(2) * CC_n_vir(2)
  CC_n_ab = CC_n_occ(1) * CC_n_vir(1) * CC_n_occ(2) * CC_n_vir(2)
  CC_n_2a = CC_n_occ(1) * (CC_n_occ(1) - 1) * CC_n_vir(1) * (CC_n_vir(1) - 1) / 4
  CC_n_2b = CC_n_occ(2) * (CC_n_occ(2) - 1) * CC_n_vir(2) * (CC_n_vir(2) - 1) / 4

  CC_n_config = CC_n_1a + CC_n_1b + CC_n_ab + CC_n_2a + CC_n_2b

! Initialize vectors
  Allocate(CC_mem_w(n_tasks,5), stat=errnum)
  Call check_allocation(errnum, 'CC_mem_w in CC')

  Allocate(CC_index_w(n_tasks,5), stat=errnum)
  Call check_allocation(errnum, 'CC_w_index in CC')

! Caalculate memory distribution for alpha single-excitaions
  Call CC_mem_distribution(CC_n_1a,n_tasks,CC_mem_w(:,1),CC_index_w(:,1))

! Calculate memory distribution for beta single-excitaions
  Call CC_mem_distribution(CC_n_1b,n_tasks,CC_mem_w(:,2),CC_index_w(:,2))

! Calculate memory distribution for alpha-beta double-excitaions
  Call CC_mem_distribution(CC_n_ab,n_tasks,CC_mem_w(:,3),CC_index_w(:,3))

! Calculate memory distribution for alpha-alpha double-excitaions
  Call CC_mem_distribution(CC_n_2a,n_tasks,CC_mem_w(:,4),CC_index_w(:,4))

! Calculate memory distribution for beta-beta double-excitaions
  Call CC_mem_distribution(CC_n_2b,n_tasks,CC_mem_w(:,5),CC_index_w(:,5))
!  do i_tmp = 1,5
!    write(40+myid,*) 'CC_mem_w(:,:)',CC_mem_w(:,i_tmp)
!    write(40+myid,*) 'CC_index_w(:,:)',CC_index_w(:,i_tmp)
!  end do

! Initialize coefficient vectors
  Allocate(CC_t_1a(CC_n_occ(1),CC_n_vir(1)), stat=errnum)
  Call check_allocation(errnum, 'CC_t_1a in CC')

  Allocate(CC_t_ab(CC_n_ab), stat=errnum)
  Call check_allocation(errnum, 'CC_t_ab in CC')

  Allocate(CC_t_2a(CC_n_2a), stat=errnum)
  Call check_allocation(errnum, 'CC_t_2a in CC')

  if (n_spin.ne.1) then
    Allocate(CC_t_1b(CC_n_occ(2),CC_n_vir(2)), stat=errnum)
    Call check_allocation(errnum, 'CC_t_1b in CC')

    Allocate(CC_t_2b(CC_n_2b), stat=errnum)
    Call check_allocation(errnum, 'CC_t_2b in CC')
  end if

! Initialize delta*t vectors
  Allocate(CC_w_1a(CC_mem_w(myid+1,1)), stat=errnum)
  Call check_allocation(errnum, 'CC_w_1a in CC')

  Allocate(CC_w_ab(CC_mem_w(myid+1,3)), stat=errnum)
  Call check_allocation(errnum, 'CC_w_ab in CC')

  Allocate(CC_w_2a(CC_mem_w(myid+1,4)), stat=errnum)
  Call check_allocation(errnum, 'CC_w_2a in CC')

  if (n_spin.ne.1) then
    Allocate(CC_w_1b(CC_mem_w(myid+1,2)), stat=errnum)
    Call check_allocation(errnum, 'CC_w_1b in CC')

    Allocate(CC_w_2b(CC_mem_w(myid+1,5)), stat=errnum)
    Call check_allocation(errnum, 'CC_w_2b in CC')
  end if

! Initialize vectors for RLE algorithm
! storage of basis functions
  Allocate(CC_t_1a_sv(CC_n_RLE_sv,CC_mem_w(myid+1,1)), stat=errnum)
  Call check_allocation(errnum, 'CC_t_1a_sv in CC')

  Allocate(CC_t_ab_sv(CC_n_RLE_sv,CC_mem_w(myid+1,3)), stat=errnum)
  Call check_allocation(errnum, 'CC_t_ab_sv in CC')

  Allocate(CC_t_2a_sv(CC_n_RLE_sv,CC_mem_w(myid+1,4)), stat=errnum)
  Call check_allocation(errnum, 'CC_t_2a_sv in CC')

  if (n_spin.ne.1) then
    Allocate(CC_t_1b_sv(CC_n_RLE_sv,CC_mem_w(myid+1,2)), stat=errnum)
    Call check_allocation(errnum, 'CC_t_1b_sv in CC')

    Allocate(CC_t_2b_sv(CC_n_RLE_sv,CC_mem_w(myid+1,5)), stat=errnum)
    Call check_allocation(errnum, 'CC_t_2b_sv in CC')
  end if

! Initialize auxiliary vectors

! Calculate memory distribution of auxiliary vectors h_ij, h_ba, h_bj, g_ca, g_ik
  Allocate(CC_mem_h_ij(n_tasks,2), stat=errnum)
  Call check_allocation(errnum, 'CC_mem_h_ij in CC')

  Allocate(CC_mem_h_ba(n_tasks,2), stat=errnum)
  Call check_allocation(errnum, 'CC_mem_h_ba in CC')

  Allocate(CC_mem_h_bj(n_tasks,2), stat=errnum)
  Call check_allocation(errnum, 'CC_mem_h_bj in CC')

  Allocate(CC_mem_g_ca(n_tasks,2), stat=errnum)
  Call check_allocation(errnum, 'CC_mem_g_ca in CC')

  Allocate(CC_mem_g_ik(n_tasks,2), stat=errnum)
  Call check_allocation(errnum, 'CC_mem_g_ik in CC')

  Allocate(CC_index_h_ij(n_tasks,2), stat=errnum)
  Call check_allocation(errnum, 'CC_index_h_ij in CC')

  Allocate(CC_index_h_ba(n_tasks,2), stat=errnum)
  Call check_allocation(errnum, 'CC_index_h_ba in CC')

  Allocate(CC_index_h_bj(n_tasks,2), stat=errnum)
  Call check_allocation(errnum, 'CC_index_h_bj in CC')

  Allocate(CC_index_g_ca(n_tasks,2), stat=errnum)
  Call check_allocation(errnum, 'CC_index_g_ca in CC')

  Allocate(CC_index_g_ik(n_tasks,2), stat=errnum)
  Call check_allocation(errnum, 'CC_index_g_ik in CC')


  do i_tmp = 1, 2

    j_tmp = CC_n_occ(i_tmp) * CC_n_occ(i_tmp)
    Call CC_mem_distribution(j_tmp,n_tasks,CC_mem_h_ij(:,i_tmp),CC_index_h_ij(:,i_tmp))

    j_tmp = CC_n_vir(i_tmp) * CC_n_vir(i_tmp)
    Call CC_mem_distribution(j_tmp,n_tasks,CC_mem_h_ba(:,i_tmp),CC_index_h_ba(:,i_tmp))

    j_tmp = CC_n_vir(i_tmp) * CC_n_occ(i_tmp)
    Call CC_mem_distribution(j_tmp,n_tasks,CC_mem_h_bj(:,i_tmp),CC_index_h_bj(:,i_tmp))

    j_tmp = CC_n_vir(i_tmp) * CC_n_vir(i_tmp)
    Call CC_mem_distribution(j_tmp,n_tasks,CC_mem_g_ca(:,i_tmp),CC_index_g_ca(:,i_tmp))

    j_tmp = CC_n_occ(i_tmp) * CC_n_occ(i_tmp)
    Call CC_mem_distribution(j_tmp,n_tasks,CC_mem_g_ik(:,i_tmp),CC_index_g_ik(:,i_tmp))

  end do

! in alpha space
  Allocate(CC_h_ij_a(CC_n_occ(1),CC_n_occ(1)), stat=errnum)
  Call check_allocation(errnum, 'CC_h_ij_a in CC')

  Allocate(CC_h_ba_a(CC_n_vir(1),CC_n_vir(1)), stat=errnum)
  Call check_allocation(errnum, 'CC_h_ba_a in CC')

  Allocate(CC_h_bj_a(CC_n_vir(1),CC_n_occ(1)), stat=errnum)
  Call check_allocation(errnum, 'CC_h_bj_a in CC')

  Allocate(CC_g_ca_a(CC_n_vir(1),CC_n_vir(1)), stat=errnum)
  Call check_allocation(errnum, 'CC_g_ca_a in CC')

  Allocate(CC_g_ik_a(CC_n_occ(1),CC_n_occ(1)), stat=errnum)
  Call check_allocation(errnum, 'CC_g_ik_a in CC')

! in beta space
  if (n_spin.ne.1) then
    Allocate(CC_h_ij_b(CC_n_occ(2),CC_n_occ(2)), stat=errnum)
    Call check_allocation(errnum, 'CC_h_ij_b in CC')

    Allocate(CC_h_ba_b(CC_n_vir(2),CC_n_vir(2)), stat=errnum)
    Call check_allocation(errnum, 'CC_h_ba_b in CC')
  
    Allocate(CC_h_bj_b(CC_n_vir(2),CC_n_occ(2)), stat=errnum)
    Call check_allocation(errnum, 'CC_h_bj_b in CC')

    Allocate(CC_g_ca_b(CC_n_vir(2),CC_n_vir(2)), stat=errnum)
    Call check_allocation(errnum, 'CC_g_ca_b in CC')

    Allocate(CC_g_ik_b(CC_n_occ(2),CC_n_occ(2)), stat=errnum)
    Call check_allocation(errnum, 'CC_g_ik_b in CC')
  end if

! For auxiliary vectors a and b
  Allocate(CC_mem_a(n_tasks,3), stat=errnum)
  Call check_allocation(errnum, 'CC_mem_a in CC')

  Allocate(CC_index_a(n_tasks,3), stat=errnum)
  Call check_allocation(errnum, 'CC_index_a in CC')

  Allocate(CC_mem_b(n_tasks,3), stat=errnum)
  Call check_allocation(errnum, 'CC_mem_b in CC')

  Allocate(CC_index_b(n_tasks,3), stat=errnum)
  Call check_allocation(errnum, 'CC_index_b in CC')

! in alpha-beta space
  i_tmp = (CC_n_occ(1) * CC_n_occ(2)) ** 2
  Call CC_mem_distribution(i_tmp,n_tasks,CC_mem_a(:,1),CC_index_a(:,1))

  Allocate(CC_aux_a_ab(CC_mem_a(myid+1,1)), stat=errnum)
  Call check_allocation(errnum, 'CC_aux_a_ab in CC')

  i_tmp = (CC_n_vir(1) * CC_n_vir(2)) ** 2
  Call CC_mem_distribution(i_tmp,n_tasks,CC_mem_b(:,1),CC_index_b(:,1))

!  write(40+myid,*) 'CC_mem_b(:,1)',CC_mem_b(:,1)
!  write(40+myid,*) 'CC_mem_b(:,1)',CC_index_b(:,1)

  Allocate(CC_aux_b_ab(CC_mem_b(myid+1,1)), stat=errnum)
  Call check_allocation(errnum, 'CC_aux_b_ab in CC')

! in alpha-alpha space
  i_tmp = ((CC_n_occ(1) - 1) * CC_n_occ(1) / 2) ** 2
  Call CC_mem_distribution(i_tmp,n_tasks,CC_mem_a(:,2),CC_index_a(:,2))

  Allocate(CC_aux_a_aa(CC_mem_a(myid+1,2)), stat=errnum)
  Call check_allocation(errnum, 'CC_aux_a_aa in CC')

  i_tmp = ((CC_n_vir(1) - 1) * CC_n_vir(1) / 2) ** 2
  Call CC_mem_distribution(i_tmp,n_tasks,CC_mem_b(:,2),CC_index_b(:,2))

  Allocate(CC_aux_b_aa(CC_mem_b(myid+1,2)), stat=errnum)
  Call check_allocation(errnum, 'CC_aux_b_aa in CC')

! in beta-beta space if necessary
  if (n_spin.ne.1) then
    i_tmp = ((CC_n_occ(2) - 1) * CC_n_occ(2) / 2) ** 2
    Call CC_mem_distribution(i_tmp,n_tasks,CC_mem_a(:,3),CC_index_a(:,3))

    Allocate(CC_aux_a_bb(CC_mem_a(myid+1,3)), stat=errnum)
    Call check_allocation(errnum, 'CC_aux_a_bb in CC')

    i_tmp = ((CC_n_vir(2) - 1) * CC_n_vir(2) / 2) ** 2
    Call CC_mem_distribution(i_tmp,n_tasks,CC_mem_b(:,3),CC_index_b(:,3))

    Allocate(CC_aux_b_bb(CC_mem_b(myid+1,3)), stat=errnum)
    Call check_allocation(errnum, 'CC_aux_b_bb in CC')
  end if

! For auxiliary vector h
  Allocate(CC_mem_h(n_tasks,6), stat=errnum)
  Call check_allocation(errnum, 'CC_mem_h in CC')

  Allocate(CC_index_h(n_tasks,6), stat=errnum)
  Call check_allocation(errnum, 'CC_index_h in CC')

! icak in alpha-alpha-alpha-alpha space
  i_tmp = (CC_n_occ(1) * CC_n_vir(1)) ** 2
  Call CC_mem_distribution(i_tmp,n_tasks,CC_mem_h(:,1),CC_index_h(:,1))

  Allocate(CC_aux_h_aaaa(CC_mem_h(myid+1,1)), stat=errnum)
  Call check_allocation(errnum, 'CC_aux_h_aaaa in CC')

! icak in alpha-beta-alpha-beta space
  i_tmp = CC_n_occ(2) * CC_n_vir(2) * CC_n_occ(1) * CC_n_vir(1)
  Call CC_mem_distribution(i_tmp,n_tasks,CC_mem_h(:,2),CC_index_h(:,2))

  Allocate(CC_aux_h_abab(CC_mem_h(myid+1,2)), stat=errnum)
  Call check_allocation(errnum, 'CC_aux_h_abab in CC')

! icak in alpha-beta-beta-alpha space
  i_tmp = CC_n_occ(1) * CC_n_vir(2) * CC_n_occ(1) * CC_n_vir(2)
  Call CC_mem_distribution(i_tmp,n_tasks,CC_mem_h(:,3),CC_index_h(:,3))

  Allocate(CC_aux_h_abba(CC_mem_h(myid+1,3)), stat=errnum)
  Call check_allocation(errnum, 'CC_aux_h_abba in CC')

! The counter part of h_icak
  if (n_spin.ne.1) then

! icak in beta-alpha-alpha-beta space
    i_tmp = (CC_n_occ(2) * CC_n_vir(1)) ** 2
    Call CC_mem_distribution(i_tmp,n_tasks,CC_mem_h(:,4),CC_index_h(:,4))

    Allocate(CC_aux_h_baab(CC_mem_h(myid+1,4)), stat=errnum)
    Call check_allocation(errnum, 'CC_aux_h_baab in CC')

! icak in beta-alpha-beta-alpha space
    i_tmp = CC_n_occ(1) * CC_n_vir(1) * CC_n_occ(2) * CC_n_vir(2)
    Call CC_mem_distribution(i_tmp,n_tasks,CC_mem_h(:,5),CC_index_h(:,5))

    Allocate(CC_aux_h_baba(CC_mem_h(myid+1,5)), stat=errnum)
    Call check_allocation(errnum, 'CC_aux_h_baba in CC')

! icak in beta-beta-beta-beta space
    i_tmp = (CC_n_occ(2) * CC_n_vir(2)) ** 2
    Call CC_mem_distribution(i_tmp,n_tasks,CC_mem_h(:,6),CC_index_h(:,6))

    Allocate(CC_aux_h_bbbb(CC_mem_h(myid+1,6)), stat=errnum)
    Call check_allocation(errnum, 'CC_aux_h_bbbb in CC')

  end if 

  Allocate(CC_RLE_t(CC_n_RLE_sv,CC_n_RLE_sv), stat=errnum)
  Call check_allocation(errnum, 'CC_REL_t in CC')

  CC_RLE_t = 0.0D0

! Calculate memory demand for each core
  CC_mem_demand = 0.0D0

  ! memory demand of saving RI integrals
  CC_mem_demand = CC_mem_demand + Real(n_basbas*n_states*n_states*n_spin)

  ! memory demand of saving coefficients
  CC_mem_demand = CC_mem_demand + Real(CC_n_1a)
  CC_mem_demand = CC_mem_demand + Real(CC_n_ab)
  CC_mem_demand = CC_mem_demand + Real(CC_n_2a)
  if (n_spin.ne.1) then
    CC_mem_demand = CC_mem_demand + Real(CC_n_1b)
    CC_mem_demand = CC_mem_demand + Real(CC_n_2b)
  end if

  ! memory demand of saving w vector and RLE vectors
  j_tmp = CC_n_RLE_sv + 1
  CC_mem_demand = CC_mem_demand + Real(CC_mem_w(myid+1,1) * j_tmp)
  CC_mem_demand = CC_mem_demand + Real(CC_mem_w(myid+1,3) * j_tmp)
  CC_mem_demand = CC_mem_demand + Real(CC_mem_w(myid+1,4) * j_tmp)
  if (n_spin.ne.1) then
    CC_mem_demand = CC_mem_demand + Real(CC_mem_w(myid+1,2) * j_tmp)
    CC_mem_demand = CC_mem_demand + Real(CC_mem_w(myid+1,5) * j_tmp)
  end if

  ! memory demand of CC_h_ij_a
  CC_mem_demand = CC_mem_demand + Real(CC_n_occ(1) * CC_n_occ(1))

  ! memory demand of CC_h_ba_a
  CC_mem_demand = CC_mem_demand + Real(CC_n_vir(1) * CC_n_vir(1))

  ! memory demand of CC_h_bj_a
  CC_mem_demand = CC_mem_demand + Real(CC_n_vir(1) * CC_n_occ(1))

  ! memory demand of CC_g_ca_a
  CC_mem_demand = CC_mem_demand + Real(CC_n_vir(1) * CC_n_vir(1))

  ! memory demand of CC_g_ik_a
  CC_mem_demand = CC_mem_demand + Real(CC_n_occ(1) * CC_n_occ(1))

  if (n_spin.ne.1) then
    ! memory demand of CC_h_ij_b
    CC_mem_demand = CC_mem_demand + Real(CC_n_occ(2) * CC_n_occ(2))

    ! memory demand of CC_h_ba_b
    CC_mem_demand = CC_mem_demand + Real(CC_n_vir(2) * CC_n_vir(2))

    ! memory demand of CC_h_bj_b
    CC_mem_demand = CC_mem_demand + Real(CC_n_vir(2) * CC_n_occ(2))

    ! memory demand of CC_g_ca_b
    CC_mem_demand = CC_mem_demand + Real(CC_n_vir(2) * CC_n_vir(2))

    ! memory demand of CC_g_ik_b
    CC_mem_demand = CC_mem_demand + Real(CC_n_occ(2) * CC_n_occ(2))
  end if

  ! memory demand of CC_aux_a_aa
  CC_mem_demand = CC_mem_demand + Real(CC_mem_a(myid+1,2))

  ! memory demand of CC_aux_b_aa
  CC_mem_demand = CC_mem_demand + Real(CC_mem_b(myid+1,2))

  ! memory demand of CC_aux_a_ab
  CC_mem_demand = CC_mem_demand + Real(CC_mem_a(myid+1,1))

  ! memory demand of CC_aux_b_ab
  CC_mem_demand = CC_mem_demand + Real(CC_mem_b(myid+1,1))

  if (n_spin.ne.1) then
    ! memory demand of CC_aux_a_bb
    CC_mem_demand = CC_mem_demand + Real(CC_mem_a(myid+1,3))

    ! memory demand of CC_aux_b_bb
    CC_mem_demand = CC_mem_demand + Real(CC_mem_b(myid+1,3))
  end if

  ! memory demand of CC_aux_h_aaaa
  CC_mem_demand = CC_mem_demand + Real(CC_mem_h(myid+1,1))

  ! memory demand of CC_aux_h_abab
  CC_mem_demand = CC_mem_demand + Real(CC_mem_h(myid+1,2))

  ! memory demand of CC_aux_h_abba
  CC_mem_demand = CC_mem_demand + Real(CC_mem_h(myid+1,3))

  if (n_spin.ne.1) then
    ! memory demand of CC_aux_h_baab
    CC_mem_demand = CC_mem_demand + Real(CC_mem_h(myid+1,4))

    ! memory demand of CC_aux_h_baba
    CC_mem_demand = CC_mem_demand + Real(CC_mem_h(myid+1,5))

    ! memory demand of CC_aux_h_bbbb
    CC_mem_demand = CC_mem_demand + Real(CC_mem_h(myid+1,6))
  end if

  CC_mem_demand = CC_mem_demand * 8.0D0 / 1024.0D0 / 1024.0D0

! Set up reference (HF) energy
  E_HF = total_energy

  Call CC_clear_w()

  CC_RLE_flag_singular = .false.
  CC_RLE_flag_sv_overflow = .false.

  if (myid.eq.0) then
    write(use_unit,"(2x,A)") 'Number of configurations'
    write(use_unit,"(10x,A2,8x,A2,8x,A2,8x,A2,8x,A2,8x)") '1a','1b','ab','2a','2b'
    write(use_unit,"(2x,5I10)") CC_n_1a,CC_n_1b,CC_n_ab,CC_n_2a,CC_n_2b
    write(use_unit,"(2x,A6,I10,2x,A)") 'Totally ',CC_n_config,'configurations are considered'
    write(use_unit,"(2x,F14.6,A)") CC_mem_demand,'MB memory is needed for each core'
    write(use_unit,"(2x,A)") 'End of initialization'
  end if

  End Subroutine CC_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_get_ovlp3ks()

  Use dimensions
  Use runtime_choices
  Use physics
  Use hartree_fock
  Use basis
  Use prodbas
  Use mpi_tasks
  Use CC_corr

  Implicit None

  Integer :: i_basbas, i_spin, i_state, j_state, errnum

  Double precision , dimension(:,:,:,:) , allocatable :: ovlp_3ks_local


  Allocate(ovlp_3ks_local(n_loc_prodbas,n_states,n_states,n_spin),stat=errnum)
  Call check_allocation(errnum,'ovlp_3ks_local(1D) in CC')

  Call transform_ovlp3fn(n_states,ks_eigenvector,ovlp_3fn,ovlp_3ks_local)

  ovlp_3ks = 0.0D0

  do i_spin = 1, n_spin
    do i_state = 1, n_states
      do j_state = 1, n_states
        do i_basbas = 1, n_loc_prodbas
          ovlp_3ks(map_prodbas(i_basbas,myid+1),j_state,i_state,i_spin) = &
                          ovlp_3ks_local(i_basbas,j_state,i_state,i_spin)
        end do
      end do
    end do
  end do

  Deallocate(ovlp_3ks_local)

  Call sync_vector(ovlp_3ks, n_basbas*n_states*n_states*n_spin)

  End Subroutine CC_get_ovlp3ks


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_ini_guess()

! This subroutine evaluate MP1 wave function as the initial guess of CC calculations,
! and calculates MP2 energy

  use dimensions
  use runtime_choices
  use physics
  use mpi_tasks
  Use CC_corr

  Implicit None
  Integer :: i_tmp,j_tmp,k_tmp,i_spin
  Integer (kind = CC_ip) :: i_state,i_state_local
  Integer :: iii,jjj,aaa,bbb
  Integer , dimension(4) :: CC_ab_4_index,CC_xx_4_index
  Double precision :: CC_V_0t,CC_E_diff_0_s,et1,et2

  !Clear the w vector
  Call CC_clear_w()

  !Clear the t vector
  Call CC_clear_t()

  !Clear CC_E_corr
  CC_E_corr = 0.0D0

! Do calculations for alpha-alpha and beta-beta excitations
  do i_spin = 1, n_spin
    i_state = CC_index_w(myid+1,3+i_spin) - 1
  
    do i_state_local = 1, CC_mem_w(myid+1,3+i_spin)
    
      i_state = i_state + 1
    
      ! Calculate the index
      Call CC_decode_index_xx(i_state,i_spin,iii,jjj,aaa,bbb)
  
  ! MP2 calculation
    ! Calculate MP1 wave function
      ! Calculate the integral <ij||ab>
      Call CC_Calc_integral(CC_V_0t,iii,jjj,aaa,bbb,i_spin,i_spin,i_spin,i_spin)
      
      ! Calculate the energy difference between reference state and i_state state.
      ! CC_E_diff_0_s = E_0 - E_s = e_i + e_j - e_a - e_b 
      Call CC_Calc_eigenvalue_diff(aaa,iii,i_spin,i_spin,et1)
      Call CC_Calc_eigenvalue_diff(bbb,jjj,i_spin,i_spin,et2)
      CC_E_diff_0_s =  et1 + et2
  
      ! Calculate coefficients
      if (i_spin.eq.1) then
        CC_t_2a(i_state) = CC_V_0t / CC_E_diff_0_s
        CC_w_2a(i_state_local) = CC_t_2a(i_state) 
      else
        CC_t_2b(i_state) = CC_V_0t / CC_E_diff_0_s
        CC_w_2b(i_state_local) = CC_t_2b(i_state)
      end if
       
    ! Calculate the energy contribution to MP2 calculation for alpha-alpha excitations
      if (i_spin.eq.1) then
        CC_E_corr = CC_E_corr + CC_V_0t * CC_t_2a(i_state)
      else
        CC_E_corr = CC_E_corr + CC_V_0t * CC_t_2b(i_state)
      end if
  
    end do
  end do

  if (n_spin.eq.1) then
    CC_E_corr = CC_E_corr * 2.0D0
  end if

! Do calculations for alpha-beta excitations
  i_state = CC_index_w(myid+1,3) - 1
  Call CC_decode_index_ab(i_state,iii,jjj,aaa,bbb)

  do i_state_local = 1, CC_mem_w(myid+1,3)
    i_state = i_state + 1
! increace the index in alpha-beta space
    Call CC_inc_index_ab(iii,jjj,aaa,bbb)

! MP2 calculation
  ! Calculate MP1 wave-function in alpha-beta space
    ! Calculate the integral <ij||ab>
    Call CC_Calc_integral(CC_V_0t,iii,jjj,aaa,bbb,1,2,1,2)

    ! Calculate the energy difference between reference state and i_state state.
    ! CC_E_diff_0_s = E_0 - E_s = e_i + e_j - e_a - e_b
    Call CC_Calc_eigenvalue_diff(aaa,iii,1,1,et1)
    Call CC_Calc_eigenvalue_diff(bbb,jjj,2,2,et2)
    CC_E_diff_0_s = et1 + et2

    ! Calculate coefficients
    CC_t_ab(i_state) = CC_V_0t / CC_E_diff_0_s

    CC_w_ab(i_state_local) = CC_t_ab(i_state)

  ! Calculate the energy contribution to MP2 calculation for alpha-beta excitations
    CC_E_corr = CC_E_corr + CC_V_0t * CC_t_ab(i_state)

  end do

! Synchronize correlation energy
  Call sync_real_number(CC_E_corr)

! Synchronize coefficients
  Call sync_vector(CC_t_ab,CC_n_ab)
  Call sync_vector(CC_t_2a,CC_n_2a)

  If (n_spin.ne.1) then
    Call sync_vector(CC_t_2b,CC_n_2b)
  end if

! Save the trial function as a basis function
  CC_RLE_ndc = 1

  CC_t_1a_sv(1,:) = CC_w_1a
  CC_t_ab_sv(1,:) = CC_w_ab
  CC_t_2a_sv(1,:) = CC_w_2a

  if (n_spin.ne.1) then
    CC_t_1b_sv(1,:) = CC_w_1b
    CC_t_2b_sv(1,:) = CC_w_2b
  end if

  CC_RLE_t(1,1) = 1.0D0

  End Subroutine CC_ini_guess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Jacob()

! This subroutine does the first iteration of CC SCF procedure

  Use dimensions
  Use mpi_tasks
  Use CC_corr
  
  Implicit None

  ! Generate new t vector
  Call CC_RLE_new_t()

  ! Save new t vector to RLE basis set
  
  Call CC_RLE_adding_t()

  CC_RLE_t(2,2) = 1.0D0

  Call CC_Jacob_sv_t()

  End Subroutine CC_Jacob

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine CC_Jacob_sv_t()

! This subroutine saves the newest trial function(w vector) to t vector

  Use dimensions
  Use mpi_tasks
  Use CC_corr

  Implicit None

  Integer (kind = CC_ip) :: i_run,i_state
  Integer :: i_spin,iii,jjj,aaa,bbb,ainx,iinx,jinx

  ! Synchronize t vector
  Call CC_clear_t()

!  print*,'sv_t'
  do i_spin = 1, n_spin
    i_state = CC_index_w(myid+1,i_spin) - 1
    Call CC_decode_index_s(i_state,i_spin,iii,aaa)

    do i_run = 1, CC_mem_w(myid+1,i_spin)
      Call CC_inc_index_s(i_spin,iii,aaa)
      ainx = aaa - CC_n_elec(i_spin)
      iinx = iii - CC_valence + 1

      if (i_spin.eq.1) then
        CC_t_1a(iinx,ainx) = CC_w_1a(i_run)
      else
        CC_t_1b(iinx,ainx) = CC_w_1b(i_run)
      end if
    end do

    i_state = CC_index_w(myid+1,3+i_spin) - 1

    do i_run = 1, CC_mem_w(myid+1,3+i_spin)

      i_state = i_state + 1

      if (i_spin.eq.1) then
        CC_t_2a(i_state) = CC_w_2a(i_run)
      else
        CC_t_2b(i_state) = CC_w_2b(i_run)
      end if

    end do

  end do

  i_state = CC_index_w(myid+1,3) - 1

  do i_run = 1, CC_mem_w(myid+1,3)
    i_state = i_state + 1
    CC_t_ab(i_state) = CC_w_ab(i_run)
  end do

  Call sync_vector(CC_t_1a,CC_n_1a)
  Call sync_vector(CC_t_ab,CC_n_ab)
  Call sync_vector(CC_t_2a,CC_n_2a)

  if (n_spin.ne.1) then
    Call sync_vector(CC_t_1b,CC_n_1b)
    Call sync_vector(CC_t_2b,CC_n_2b)
  end if


  End Subroutine CC_Jacob_sv_t
