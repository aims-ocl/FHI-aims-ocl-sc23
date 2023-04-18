!****s* FHI-aims/post_scf_correlation_treatment_p0
!  NAME
!   post_scf_correlation_treatment_p0
!  SYNOPSIS

      subroutine post_cpt2_blacs()

!  PURPOSE
!  This subroutine evaluate the necessary matrix elements (3-center overlap,
!  coulomb matrix) used for accurate correlation energy calculations beyond DFT
!  (e.g. MP2, RPA, RPA+, etc).
!
!  USES

      use dimensions
      use species_data
      use runtime_choices
      use pbc_lists, only : k_point_list
      use prodbas, only : OVLP_TYPE_COULOMB
      use hartree_fock
      use hartree_fock_p0
      use physics
      use gw_para
      use my_lvl_triples, only : my_cleanup_lvl_triples, my_initialize_lvl_triples
      use lvl_triples, only : cleanup_lvl_triples, initialize_lvl_triples
      use tight_binding_auxmat
      use timing
      use mpi_tasks
      use sbt_overlap_aims
      use calculate_fock_matrix_p0, only : cleanup_fock_matrix_calculations, &
          init_fock_matrix_calculations, evaluate_exchange_matr_realspace_p0
      use cpt2_blacs
      use lvl_tricoeff_cpt2
      use redistribute_cpt2_kspace
      use pbc_lists, only: k_phase, n_cells
      implicit none

!  ARGUMENTS

!  INPUT
!    none
!  OUTPUT
!    none
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!  local variables
      character*20 filename
      logical :: need_o3fn
      real(kind=8) :: time_prodbas_add, clock_time_prodbas_add
      real(kind=8) :: rpa_tot_en
      real(kind=8) :: pt2_tot_en
      real(kind=8) :: en_se, en_rse
      real(kind=8) :: exx_ene, d_exx_ene(3) !SVL dummy variables, not actually used anywhere but needed for hf_realspace call 

      integer :: n_k_points_min
      integer :: bb_row_s, bb_col_s

      complex(kind=8), allocatable :: coulomb_gamma_blacs(:,:)

! Counter
      integer i_spin
      integer i_state
      integer i_k_point, i_k_point_local
      integer i_band, i_cell
      integer i_k_point_index
      integer i_prodbas_1, i_prodbas_2
      integer i_bb_1, i_bb_2, i_bb_s_1, i_bb_s_2

      integer :: info, mpierr
      character*150 :: info_str
      character(*), parameter :: &
          func='post_scf_correlation_treatment_p0'

      real*8 :: energy_xc

      logical :: gamma_center_cpt2

!      nrThreads   = 1
!      myThreadId = 1
!!$    nrThreads = omp_get_max_threads()
!!$    myThreadId = omp_get_thread_num() + 1

      call perfinit('evpost')
      call perfon('pscf')
      call get_timestamps(time_prodbas_add, clock_time_prodbas_add)
!      call perfon('ini1')
      ! IGOR :: always use the full symmetry, which can dramatically reduce the time cost.
      use_full_symmetry = .true.
      ! IGOR :: ===========================================
!  n_k_points, k_point_list, irk_point_mapping, inv_irk_point_mapping,
!  irk_point_included (all defined in pbc_lists.f90) will be determined here
      call determine_irreducible_k_grid ()

      ! for non-hartree-fock self-consistent calculations, we need to
      ! construct the product (auxiliary) basis functions, and evaluate
      ! the 3-center overlap and coulomb matrix elements here.
!call perfoff

      if (.not.use_hf_kspace) then
!call perfon('ini2')
        if (.not.allocated(n_homo)) then
         allocate (n_homo(n_spin))
          n_homo(:)=0
        endif

        call initialize_prodbas()


!call perfoff
!call perfon('ini3')
        !if(.not.use_hf_realspace)then
           call allocate_hartree_fock_p0()
           if(.not.use_hf_kspace_with_rpa)then
              call cleanup_fock_matrix_calculations
              call init_fock_matrix_calculations(.false., .false.)
              call evaluate_exchange_matr_realspace_p0(KS_eigenvector,KS_eigenvector_complex,occ_numbers, &
                   exx_ene, d_exx_ene, .false., .false.)
              call cleanup_fock_matrix_calculations
           endif
        !endif
!call perfoff
call perfon('ini4')

    !energy_xc=0.0d0
    !call get_exchange_energy_p0 &
    !   (hybrid_coeff, KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue,&
    !    k_weights, occ_numbers,chemical_potential,eex_energy_cpt2, energy_xc)
    !write(use_unit,*) "igor debug, eex_energy_cpt2:", eex_energy_cpt2
    eex_energy_cpt2 = -1.0d0*exx_ene

   call init_cpt2_blacs_distribution(occ_numbers)

   ! check if it is a gamma-only cpt2 calculation
   gamma_center_cpt2 = k_point_list(1,1).lt.1.0d-5 .and. &
                       k_point_list(1,2).lt.1.0d-5 .and. &
                       k_point_list(1,3).lt.1.0d-5
   if (n_k_points.eq.1 .and. gamma_center_cpt2) then
       gamma_only_cpt2 = .true.
   else
       gamma_only_cpt2 = .false.
   end if

!if (myid.eq.0) then
!   print*,'reso: bb,b,st',n_basbas,n_basis,n_states
!   print*,'reso: mnbb,nc,nf',max_n_basis_sp,n_cells,n_full_freq
!end if
        if (myid_col_cpt2.eq.0) then
           allocate(lvl_tricoeff_mod_r(lbb_row:ubb_row,max_n_basis_sp,n_states,n_spin,n_ks_points_task),stat=info)
        else
           allocate(lvl_tricoeff_mod_r(lbb_row:lbb_row,1,1,1,1),stat=info)
        end if
        call check_allocation(info, 'lvl_tricoeff_mod_r', func)

        allocate(coulomb_matr_blacs(lbb_row:ubb_row,lbb_col:ubb_col,n_kq_points_task),stat=info)
        call check_allocation(info, 'coulomb_matr_blacs                 ')

!        allocate(coulomb_coeff_blacs(lbb_row:ubb_row,lbb_col:ubb_col, n_irkq_points_task),stat=info)
!        call check_allocation(info, 'coulomb_matr_blacs                 ')

        allocate(coulomb_gamma_blacs(lbb_row:ubb_row,lbb_col:ubb_col),stat=info)
        call check_allocation(info, 'coulomb_gamma_blacs                 ')

        n_k_points_min=max(n_k_points_task,1)

        allocate(kq_point_list(n_k_points,n_k_points),stat=info) 
        call check_allocation(info, 'kq_point_list', func)
        allocate(kpq_point_list(n_k_points,n_k_points),stat=info) 
        call check_allocation(info, 'kpq_point_list', func)

        if (use_threadsafe_gwinit) then
          call my_initialize_lvl_triples(OVLP_TYPE_COULOMB)
        else
          call initialize_lvl_triples(OVLP_TYPE_COULOMB)
        endif


        call initialize_tb_auxmat(1, OVLP_TYPE_COULOMB)


!        if(n_k_points.gt.1) then
           !ACHTUNG: this section can be removed by removing lvl_tricoef_recip* from evaluate_exchange_matr_kspace_p0!!
!           allocate(lvl_tricoeff_recip1(max_n_basbas_sp,n_basis,n_basis,n_k_points_min),stat=info) 
!           call check_allocation(info, 'lvl_tricoeff_recip1', func)
!           allocate(lvl_tricoeff_recip2(max_n_basbas_sp,n_basis,n_basis),stat=info) 
!           call check_allocation(info, 'lvl_tricoeff_recip2', func)
!           call get_lvl_tricoeff_recip(n_cells_task,lvl_tricoeff_recip1,lvl_tricoeff_recip2)
!        end if


        ! compute the number of real-space unit cells locally. parallelization in comm_col
        ! the entry n_cells+1 contains the contribution from the second atom
        n_cells_task = 0
        do i_cell = 1, n_cells+1
           if(myid_col_cpt2.eq.mod(i_cell,n_tasks_col_cpt2) .and.(myid_col_cpt2 .le. n_cells+1)) then
              n_cells_task = n_cells_task + 1
           endif
        enddo
        
        !compute the lvl_tricoeff_cell array only once
        !if (myid.eq.0) write(use_unit,*) 'igor debug 0a'
        !call get_memory_cpt2('igor debug 0a')
        call cpt2_init_lvl_tricoeff_recip(n_cells,n_cells_task,n_k_points, n_k_points_task, n_ks_points_task,&
             k_phase)
        !call get_memory_cpt2('igor debug 0b')

        if (use_threadsafe_gwinit) then
           call my_cleanup_lvl_triples()
        else
           call cleanup_lvl_triples()
        endif


        call cpt2_get_lvl_tricoeff_recip(n_cells, n_k_points, n_k_points_task, n_ks_points_task,&
             k_phase, KS_eigenvector, KS_eigenvector_complex, lvl_tricoeff_mod_r)
        
        call cpt2_cleanup_lvl_tricoeff_recip

        call deallocate_tb_auxmat()
        if(use_hse .and. hse_omega_hf /= 0.d0 .and. (.not. use_gw_and_hse) &
            .and. (.not. use_dftpt2_and_hse)) then
           call initialize_tb_auxmat(1, OVLP_TYPE_HSE)
        else
           call initialize_tb_auxmat(1, OVLP_TYPE_CUT)
           !call initialize_tb_auxmat(1, OVLP_TYPE_CUT_ANALYTIC)
!          call initialize_periodic_tb_auxmat(1, 1.d0)
        endif

        call determine_k_minus_q_list(kq_point_list,kpq_point_list)

!        call get_coulomb_matr_recip(coulomb_matr_recip,1)
!call perfon('gcmb')
        call get_coulomb_matr_blacs_cpt2(coulomb_matr_blacs,1)
!call perfoff
        if((.not.use_hf_realspace).and.use_hf_kspace_with_rpa)then

            print*,'ACHTUNG: noch nicht an coulomb_matr_blacs angepasst'
            stop
           call evaluate_exchange_matr_kspace_p0 &
                (KS_eigenvalue,KS_eigenvector,KS_eigenvector_complex,occ_numbers)
        endif
!<ADDED CODE>
        call deallocate_tb_auxmat()
!</ADDED CODE>


       if(.not. gamma_cut_coulomb) then
         do i_k_point_local = 1, n_kq_points_task
           i_k_point=n_tasks_kq_cpt2*(i_k_point_local-1) + myid_kq_cpt2 + 1
           !i_k_point = inv_irk_point_mapping(i_irk_point)
           if(i_k_point .eq. 1 .and. all(abs(k_point_list(i_k_point,:)).lt.1.e-2)) then
              coulomb_gamma_blacs(:,:)= coulomb_matr_blacs(:,:,i_k_point_local)        
           endif
         enddo

         call deallocate_tb_auxmat()
         call initialize_periodic_tb_auxmat(1, 1.d0)

         call get_coulomb_matr_blacs_cpt2(coulomb_matr_blacs,1)


         do i_k_point_local = 1, n_kq_points_task
           i_k_point=n_tasks_kq_cpt2*(i_k_point_local-1) + myid_kq_cpt2 + 1
           !i_k_point = inv_irk_point_mapping(i_irk_point)
           if(i_k_point .eq. 1 .and. all(abs(k_point_list(i_k_point,:)).lt.1.e-2)) then
              coulomb_matr_blacs(:,:,i_k_point_local) = coulomb_gamma_blacs(:,:)        
           endif
         enddo

! end of if(.not. gamma_cut_coulomb) 
       endif
       call perfoff
! end of if(.not.use_hf_realspace)
   endif

   deallocate(coulomb_gamma_blacs)

   call deallocate_tb_auxmat()

   if (real_eigenvectors) then
     call redistribute_cpt2_blacs_real()
   else
     call redistribute_cpt2_blacs_complex()
   end if

call perfon('pscpt2')
       call pt2_calculation(pt2_tot_en)
call perfoff

   call deallocate_cpt2_blacs()

   !deallocate(coulomb_matr_blacs)
   call perfoff
   call perfout('pscf')

end subroutine post_cpt2_blacs
!***************
